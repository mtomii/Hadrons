/*
 * A2AMesonFieldProduct.cpp, part of Hadrons
 *
 * Contractor-style named <a2aMatrix> inputs and <product> list. Inputs are
 * read ONCE per trajectory and kept in RAM; there is no EigenDiskVector,
 * Grid_init, lattice Grid, Environment, module registration, or MPI call.
 * XML/trajectory and HDF5 I/O idioms follow Contractor.cpp/A2AMatrix.hpp and
 * A2AMesonFieldSelectHits.cpp. Original Hadrons code: Copyright (C) 2015-2023,
 * Antonin Portelli and contributors. Changes preserve the GPL-2.0-or-later
 * licence. This program is distributed WITHOUT ANY WARRANTY. See LICENSE
 * in the top-level Hadrons distribution.
 */
#include <Hadrons/Global.hpp>
#include <Hadrons/A2AMatrix.hpp>
#include <algorithm>
#include <cassert>
#include <cerrno>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <stdexcept>
#include <utility>
#include <vector>
#include <sys/stat.h>
#include <unistd.h>

#ifndef HAVE_HDF5
#error "A2AMesonFieldProduct requires HDF5 support in Grid/Hadrons"
#endif

using namespace Grid;
using namespace Hadrons;

namespace A2AMesonFieldProduct
{
class TrajRange: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(TrajRange,
                                    long long, start,
                                    long long, end,
                                    long long, step);
    TrajRange(): start(-1), end(-1), step(0) {}
};
class GlobalPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(GlobalPar, TrajRange, trajCounter);
};
class MatrixPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(MatrixPar,
                                    std::string, name,
                                    std::string, file,
                                    std::string, group);
};
class ProductPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(ProductPar,
                                    std::string, mfL,
                                    std::string, mfR,
                                    long long,   deltaL,
                                    long long,   deltaR,
                                    bool,        ifDenseL,
                                    bool,        ifDenseR,
                                    long long,   initDense,
                                    std::string, output,
                                    std::string, outputGroup);
    ProductPar(): deltaL(-1), deltaR(-1), ifDenseL(false), ifDenseR(false), initDense(-1) {}
};

using Matrix = A2AMatrix<ComplexD>;
static_assert(Matrix::IsRowMajor, "Hadrons A2AMatrix is expected to be row-major");
static constexpr hsize_t nsc = Ns*Nc;
static constexpr hsize_t ioTile = 512;  // <= 4 MiB conversion buffer
static constexpr hsize_t chunkSize = 256;
static const char *metadataName = "_A2AProductMetadata";
static const char *sourcesName = "_A2AProductInputs";

// No new matrix/I/O framework: just the named resident inputs and product plan.
struct ResidentMatrix
{
    MatrixPar par;
    std::unique_ptr<H5NS::H5File> file;
    hsize_t nt{0}, rows{0}, cols{0};
    std::size_t realBytes{0};
    std::vector<Matrix> value;
};
struct ProductPlan
{
    ProductPar par;
    std::size_t left, right;
    hsize_t active;
    std::size_t realBytes;
};
struct ComplexValue { double re, im; };

void require(bool ok, const std::string &message)
{
    if (!ok) throw std::runtime_error(message);
}
hsize_t checkedMul(hsize_t a, hsize_t b)
{
    require(a == 0 || b <= std::numeric_limits<hsize_t>::max()/a,
            "Integer overflow in element/byte-count multiplication");
    return a*b;
}
hsize_t checkedAdd(hsize_t a, hsize_t b)
{
    require(b <= std::numeric_limits<hsize_t>::max()-a, "Integer overflow in size addition");
    return a+b;
}
std::size_t asSize(hsize_t n)
{
    require(n <= std::numeric_limits<std::size_t>::max(), "Size exceeds std::size_t");
    return static_cast<std::size_t>(n);
}
Eigen::Index asIndex(hsize_t n)
{
    require(n <= hsize_t(std::numeric_limits<Eigen::Index>::max()), "Dimension exceeds Eigen::Index");
    return static_cast<Eigen::Index>(n);
}
hsize_t matrixBytes(hsize_t rows, hsize_t cols)
{
    const auto n = checkedMul(rows, cols);
    asIndex(n); // Eigen also uses Index for the coefficient count.
    const auto bytes = checkedMul(n, sizeof(ComplexD));
    asSize(bytes);
    return bytes;
}
bool exists(const std::string &path)
{
    struct stat st;
    if (::lstat(path.c_str(), &st) == 0) return true;
    require(errno == ENOENT, "Cannot inspect '" + path + "': " + std::strerror(errno));
    return false;
}
std::string canonicalFile(const std::string &path)
{
    char *resolved = ::realpath(path.c_str(), nullptr);
    require(resolved != nullptr, "Cannot resolve input '" + path + "': " + std::strerror(errno));
    const std::string result(resolved);
    std::free(resolved);
    return result;
}
std::string trajectoryName(std::string name, long long traj)
{
    while (name.find("@traj@") != std::string::npos) tokenReplace(name, "traj", traj);
    return name;
}
std::string groupPath(const std::string &name)
{
    require(!name.empty(), "Empty HDF5 group; use / for the root group");
    std::stringstream in(name);
    std::string part, path, last;
    while (std::getline(in, part, '/'))
    {
        if (part.empty()) continue;
        require(part != "." && part != "..", "Invalid HDF5 group path: " + name);
        path += "/" + part;
        last = part;
    }
    require(last != HADRONS_A2AM_NAME,
            "Specify the parent GROUP, without /" + std::string(HADRONS_A2AM_NAME) + ": " + name);
    return path.empty() ? "/" : path;
}
std::string matrixPath(const std::string &group)
{
    return (group == "/" ? "" : group) + "/" + HADRONS_A2AM_NAME;
}
H5NS::Group createGroupTree(H5NS::H5File &file, const std::string &path)
{
    std::stringstream in(path);
    std::string part, current;
    while (std::getline(in, part, '/'))
    {
        if (part.empty()) continue;
        current += "/" + part;
        if (H5Lexists(file.getId(), current.c_str(), H5P_DEFAULT) == 0)
            file.createGroup(current);
    }
    return file.openGroup(path);
}
H5NS::CompType memoryType(void)
{
    H5NS::CompType type(sizeof(ComplexValue));
    type.insertMember("re", HOFFSET(ComplexValue, re), H5NS::PredType::NATIVE_DOUBLE);
    type.insertMember("im", HOFFSET(ComplexValue, im), H5NS::PredType::NATIVE_DOUBLE);
    return type;
}
H5NS::CompType outputType(std::size_t bytes)
{
    H5NS::CompType type(std::size_t(2*bytes));
    const auto &real = bytes == 4 ? H5NS::PredType::IEEE_F32LE : H5NS::PredType::IEEE_F64LE;
    type.insertMember("re", 0, real);
    type.insertMember("im", bytes, real);
    return type;
}
std::size_t complexPrecision(const H5NS::DataSet &data)
{
    require(data.getTypeClass() == H5T_COMPOUND, "Expected re/im complex compound datatype");
    const auto type = data.getCompType();
    require(type.getNmembers() == 2, "Expected two complex members: re, im");
    const int re = H5Tget_member_index(type.getId(),"re"), im = H5Tget_member_index(type.getId(),"im");
    require(re >= 0 && im >= 0, "Complex datatype must contain re and im");
    const auto rt = type.getMemberDataType(re), it = type.getMemberDataType(im);
    const auto bytes = rt.getSize();
    require(rt.getClass() == H5T_FLOAT && it.getClass() == H5T_FLOAT &&
            bytes == it.getSize() && (bytes == 4 || bytes == 8),
            "Only ComplexF/ComplexD (float/double re and im) are supported");
    return bytes;
}
void inspect(ResidentMatrix &m)
{
    const std::string path = matrixPath(m.par.group);
    std::cout << "Inspect input '" << m.par.name << "': " << m.par.file << " : " << path << std::endl;
    m.file.reset(new H5NS::H5File(m.par.file, H5F_ACC_RDONLY));
    const auto inputGroup = m.file->openGroup(m.par.group);
    require(H5Lexists(inputGroup.getId(), metadataName, H5P_DEFAULT) == 0,
            "An existing Product is not an elementary meson field: outer-index time provenance needs separate handling");
    const auto data = m.file->openDataSet(path);
    const auto space = data.getSpace();
    require(space.getSimpleExtentNdims() == 3, "Expected (time, row, column) for " + m.par.name);
    hsize_t dim[3];
    space.getSimpleExtentDims(dim);
    require(dim[0] > 0 && dim[1] > 0 && dim[2] > 0, "Empty input matrix: " + m.par.name);
    m.nt = dim[0]; m.rows = dim[1]; m.cols = dim[2];
    asIndex(m.nt);
    matrixBytes(m.rows, m.cols);
    require(m.nt <= std::vector<Matrix>().max_size(), "Too many timeslices: " + m.par.name);
    m.realBytes = complexPrecision(data);
}

// Read each numeric element exactly once into the resident matrices. HDF5 may
// read whole chunks internally. Tiles avoid an additional full-timeslice copy.
void load(ResidentMatrix &m)
{
    const auto begin = std::chrono::steady_clock::now();
    const auto path = matrixPath(m.par.group);
    std::cout << "LOAD input '" << m.par.name << "' " << m.nt << " x " << m.rows << " x " << m.cols
              << ": " << m.par.file << " : " << path << std::endl;
    auto data = m.file->openDataSet(path);
    auto space = data.getSpace();
    const auto type = memoryType();
    std::vector<ComplexValue> buf(asSize(checkedMul(std::min(ioTile,m.rows), std::min(ioTile,m.cols))));
    m.value.resize(asSize(m.nt));
    for (hsize_t t = 0; t < m.nt; ++t)
    {
        auto &mat = m.value[asSize(t)];
        mat.resize(asIndex(m.rows), asIndex(m.cols));
        for (hsize_t i = 0; i < m.rows; )
        {
            const auto nr = std::min(ioTile, m.rows-i);
            for (hsize_t j = 0; j < m.cols; )
            {
                const auto nc = std::min(ioTile, m.cols-j);
                const hsize_t count[3] = {1, nr, nc}, offset[3] = {t, i, j};
                H5NS::DataSpace mem(3, count);
                space.selectHyperslab(H5S_SELECT_SET, count, offset);
                data.read(buf.data(), type, mem, space);
                for (hsize_t r = 0; r < nr; ++r)
                for (hsize_t c = 0; c < nc; ++c)
                {
                    const auto &z = buf[std::size_t(r*nc+c)];
                    if (!std::isfinite(z.re) || !std::isfinite(z.im))
                        throw std::runtime_error("Non-finite input '" + m.par.name + "' at t=" + std::to_string(t));
                    // Bounds/products were checked once during inspection.
                    mat(Eigen::Index(i+r), Eigen::Index(j+c)) = ComplexD(z.re, z.im);
                }
                j += nc;
            }
            i += nr;
        }
    }
    std::cout << "  loaded once: " << m.par.name << " ("
              << std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count()
              << " s)" << std::endl;
}

ProductPlan planProduct(ProductPar par, const std::map<std::string,std::size_t> &names,
                        const std::vector<ResidentMatrix> &matrices)
{
    require(names.count(par.mfL) && names.count(par.mfR), "Unknown mfL/mfR name in " + par.output);
    const auto l = names.at(par.mfL), r = names.at(par.mfR);
    const auto &a = matrices[l]; const auto &b = matrices[r];
    require(a.nt == b.nt, "Left/right Nt differ in " + par.output);
    require(par.deltaL >= 0 && par.deltaR >= 0 &&
            hsize_t(par.deltaL) < a.nt && hsize_t(par.deltaR) < a.nt &&
            (par.deltaL == 0 || par.deltaR == 0),
            "Require 0 <= deltaL,deltaR < Nt and at least one offset zero");
    require(!(par.ifDenseL && par.ifDenseR), "Both contracted indices cannot be dense W");
    assert(!(par.ifDenseL && par.ifDenseR));
    require(par.initDense >= 0, "initDense must be non-negative");
    hsize_t active = a.cols;
    if (par.ifDenseL || par.ifDenseR)
    {
        const hsize_t low = par.initDense;
        const hsize_t compact = par.ifDenseL ? a.cols : b.rows;
        const hsize_t expanded = par.ifDenseL ? b.rows : a.cols;
        require(compact >= low && expanded >= low, "initDense exceeds a contracted dimension");
        require((compact-low)%nsc == 0, "Dense high dimension is not a multiple of Ns*Nc");
        require(expanded == checkedAdd(low, checkedMul(a.nt, compact-low)),
                "Dense/expanded inner dimensions do not match in " + par.output);
        active = compact;
    }
    else require(a.cols == b.rows, "Non-dense inner dimensions do not match in " + par.output);
    const auto bytes = std::max(a.realBytes,b.realBytes);
    matrixBytes(a.rows,b.cols);
    checkedMul(checkedMul(checkedMul(a.nt,a.rows),b.cols),2*bytes);
    if (par.ifDenseL) matrixBytes(active,b.cols);
    if (par.ifDenseR) matrixBytes(a.rows,active);
    return {par, l, r, active, bytes};
}

// Copy only metadata belonging to a selected input GROUP, never its matrix or
// unrelated matrices elsewhere in the file. The two sources stay separately labelled.
void copyAttributes(const H5NS::H5Object &src, H5NS::H5Object &dst)
{
    for (int i = 0; i < src.getNumAttrs(); ++i)
    {
        auto attr = src.openAttribute(unsigned(i));
        auto type = attr.getDataType(); auto space = attr.getSpace();
        require(H5Tdetect_class(type.getId(), H5T_REFERENCE) == 0, "Reference-valued metadata unsupported");
        auto out = dst.createAttribute(attr.getName(),type,space);
        const auto n = space.getSimpleExtentNpoints();
        require(n >= 0, "Invalid attribute dataspace");
        if (n == 0) continue;
        std::vector<unsigned char> buf(asSize(checkedMul(hsize_t(n),type.getSize())));
        attr.read(type,buf.data());
        try { out.write(type,buf.data()); }
        catch (...) { H5Dvlen_reclaim(type.getId(),space.getId(),H5P_DEFAULT,buf.data()); throw; }
        require(H5Dvlen_reclaim(type.getId(),space.getId(),H5P_DEFAULT,buf.data()) >= 0,
                "Could not reclaim variable-length metadata");
    }
}
void copyMetadata(const H5NS::Group &src, H5NS::Group &dst, unsigned depth = 0)
{
    require(depth < 64, "Metadata hierarchy too deep/cyclic");
    copyAttributes(src,dst);
    for (hsize_t i = 0; i < src.getNumObjs(); ++i)
    {
        const std::string name = src.getObjnameByIdx(i);
        if (name == HADRONS_A2AM_NAME) continue;
        H5L_info_t info;
        require(H5Lget_info(src.getId(),name.c_str(),&info,H5P_DEFAULT) >= 0, "Cannot inspect metadata link");
        require(info.type == H5L_TYPE_HARD, "Soft/external metadata links unsupported");
        const auto kind = src.getObjTypeByIdx(i);
        if (kind == H5G_GROUP)
        {
            auto s = src.openGroup(name); auto d = dst.createGroup(name);
            copyMetadata(s,d,depth+1);
        }
        else
        {
            require(kind == H5G_DATASET || kind == H5G_TYPE, "Unsupported metadata object");
            if (kind == H5G_DATASET)
            {
                auto d = src.openDataSet(name);
                require(H5Tdetect_class(d.getDataType().getId(),H5T_REFERENCE) == 0, "Reference-valued metadata unsupported");
                for (int a = 0; a < d.getNumAttrs(); ++a)
                    require(H5Tdetect_class(d.openAttribute(unsigned(a)).getDataType().getId(),H5T_REFERENCE) == 0,
                            "Reference-valued metadata attribute unsupported");
            }
            require(H5Ocopy(src.getId(),name.c_str(),dst.getId(),name.c_str(),H5P_DEFAULT,H5P_DEFAULT) >= 0,
                    "Cannot copy input metadata: " + name);
        }
    }
}
void saveMetadata(H5NS::Group &group, const ProductPlan &plan,
                  const ResidentMatrix &left, const ResidentMatrix &right, long long traj)
{
    XmlWriter xml("");
    write(xml,"format",std::string("HadronsA2AMesonFieldProduct"));
    write(xml,"version",1);
    write(xml,"trajectory",traj);
    write(xml,"product",plan.par);
    write(xml,"leftInput",left.par);
    write(xml,"rightInput",right.par);
    write(xml,"leftDimensions",std::vector<unsigned long long>{left.nt,left.rows,left.cols});
    write(xml,"rightDimensions",std::vector<unsigned long long>{right.nt,right.rows,right.cols});
    write(xml,"rowTimeOffset",plan.par.deltaL);
    write(xml,"columnTimeOffset",plan.par.deltaR);
    write(xml,"normalization",std::string("No extra normalization, conjugation or transpose."));
    write(xml,"timeConvention",std::string("C[t]=L[(t+deltaL)%Nt]*R[(t+deltaR)%Nt]; index time follows its own input."));
    const std::string text = xml.docString();
    H5NS::StrType type(H5NS::PredType::C_S1,H5T_VARIABLE); type.setCset(H5T_CSET_UTF8);
    H5NS::DataSpace scalar(H5S_SCALAR);
    auto data = group.createDataSet(metadataName,type,scalar);
    data.write(text,type);
    auto inputs = group.createGroup(sourcesName);
    auto l = inputs.createGroup("left"); auto r = inputs.createGroup("right");
    auto ls = left.file->openGroup(left.par.group); auto rs = right.file->openGroup(right.par.group);
    copyMetadata(ls,l); copyMetadata(rs,r);
}

void product(const ProductPlan &plan, const std::vector<ResidentMatrix> &matrices, long long traj)
{
    const auto &par = plan.par;
    const auto &left = matrices[plan.left]; const auto &right = matrices[plan.right];
    const hsize_t nt = left.nt, low = par.initDense, active = plan.active;
    const hsize_t dim[3] = {nt,left.rows,right.cols};
    const std::string temporary = par.output + ".partial." + std::to_string(::getpid());
    bool ownTemporary = false;
    std::cout << "PRODUCT " << par.mfL << "[t+" << par.deltaL << "] * " << par.mfR
              << "[t+" << par.deltaR << "] -> " << par.output << " : " << matrixPath(par.outputGroup)
              << "; shape " << dim[0] << " x " << dim[1] << " x " << dim[2]
              << "; active K=" << active << "; output Complex" << (plan.realBytes == 4 ? "F" : "D") << std::endl;
    try
    {
        {
            makeFileDir(par.output);
            // Reserve without overwriting another process/stale temporary. Then let
            // Grid initialize its own HDF5 file attributes for standard readers.
            { H5NS::H5File reserve(temporary,H5F_ACC_EXCL); ownTemporary = true; }
            { Hdf5Writer initialize(temporary); }
            H5NS::H5File out(temporary,H5F_ACC_RDWR);
            auto group = createGroupTree(out,par.outputGroup);
            H5NS::DataSpace outSpace(3,dim);
            const hsize_t chunk[3] = {1,std::min(chunkSize,dim[1]),std::min(chunkSize,dim[2])};
            H5NS::DSetCreatPropList plist; plist.setChunk(3,chunk); plist.setFletcher32();
            auto outData = group.createDataSet(HADRONS_A2AM_NAME,outputType(plan.realBytes),outSpace,plist);
            const auto memType = memoryType();
            std::vector<ComplexValue> buf(asSize(checkedMul(std::min(ioTile,dim[1]),std::min(ioTile,dim[2]))));
            Matrix selected, result(asIndex(left.rows),asIndex(right.cols));
            if (par.ifDenseL) selected.resize(asIndex(active),asIndex(right.cols));
            if (par.ifDenseR) selected.resize(asIndex(left.rows),asIndex(active));
            double multiplySeconds = 0., writeSeconds = 0.;
            const double largest = plan.realBytes == 4 ? double(std::numeric_limits<float>::max())
                                                       : std::numeric_limits<double>::max();
            for (hsize_t t = 0; t < nt; ++t)
            {
                // nt <= Eigen::Index::max(); do not add in signed int arithmetic.
                const hsize_t tL = (t+hsize_t(par.deltaL))%nt, tR = (t+hsize_t(par.deltaR))%nt;
                const Matrix &a = left.value[asSize(tL)], &b = right.value[asSize(tR)];
                const Matrix *ap = &a, *bp = &b;
                const auto computeBegin = std::chrono::steady_clock::now();
                if (par.ifDenseL)
                {
                    selected.topRows(asIndex(low)) = b.topRows(asIndex(low));
                    // Dense r=low+h*nsc+a pairs with expanded low+(h*nt+tL)*nsc+a.
                    for (hsize_t h = 0; h < (active-low)/nsc; ++h)
                        selected.middleRows(asIndex(low+h*nsc),asIndex(nsc)) =
                            b.middleRows(asIndex(low+(h*nt+tL)*nsc),asIndex(nsc));
                    bp = &selected;
                }
                else if (par.ifDenseR)
                {
                    selected.leftCols(asIndex(low)) = a.leftCols(asIndex(low));
                    // Right dense W resides at tR, not at the output's reference t.
                    for (hsize_t h = 0; h < (active-low)/nsc; ++h)
                        selected.middleCols(asIndex(low+h*nsc),asIndex(nsc)) =
                            a.middleCols(asIndex(low+(h*nt+tR)*nsc),asIndex(nsc));
                    ap = &selected;
                }
                // The same ordinary Eigen product as the non-MKL mul() path.
                // No adjoint, gamma5, 1/Nhit, global reduction or source-time sum.
                result.noalias() = (*ap)*(*bp);
                multiplySeconds += std::chrono::duration<double>(std::chrono::steady_clock::now()-computeBegin).count();
                const auto writeBegin = std::chrono::steady_clock::now();
                for (hsize_t i = 0; i < dim[1]; )
                {
                    const hsize_t nr = std::min(ioTile,dim[1]-i);
                    for (hsize_t j = 0; j < dim[2]; )
                    {
                        const hsize_t nc = std::min(ioTile,dim[2]-j);
                        for (hsize_t r = 0; r < nr; ++r)
                        for (hsize_t c = 0; c < nc; ++c)
                        {
                            const auto z = result(Eigen::Index(i+r),Eigen::Index(j+c));
                            if (!std::isfinite(z.real()) || !std::isfinite(z.imag()) ||
                                std::abs(z.real()) > largest || std::abs(z.imag()) > largest)
                                throw std::runtime_error("Non-finite result/output overflow at t=" + std::to_string(t));
                            buf[std::size_t(r*nc+c)] = {z.real(),z.imag()};
                        }
                        const hsize_t count[3] = {1,nr,nc}, offset[3] = {t,i,j};
                        H5NS::DataSpace mem(3,count);
                        outSpace.selectHyperslab(H5S_SELECT_SET,count,offset);
                        outData.write(buf.data(),memType,mem,outSpace);
                        j += nc;
                    }
                    i += nr;
                }
                writeSeconds += std::chrono::duration<double>(std::chrono::steady_clock::now()-writeBegin).count();
            }
            saveMetadata(group,plan,left,right,traj);
            out.flush(H5F_SCOPE_GLOBAL);
            std::cout << "  select/GEMM " << multiplySeconds << " s; numeric output " << writeSeconds << " s" << std::endl;
        } // All HDF5 handles are closed before publication.
        if (::link(temporary.c_str(),par.output.c_str()) != 0)
        {
            const int error = errno; ownTemporary = false;
            throw std::runtime_error("Cannot publish output: " + std::string(std::strerror(error)) +
                                     "; completed temporary retained: " + temporary);
        }
        if (::unlink(temporary.c_str()) != 0)
            std::cerr << "Warning: completed output exists, but could not remove " << temporary << std::endl;
        ownTemporary = false;
        std::cout << "Completed: " << par.output << std::endl;
    }
    catch (...) { if (ownTemporary) ::unlink(temporary.c_str()); throw; }
}

void runTrajectory(const std::vector<MatrixPar> &input, const std::vector<ProductPar> &products, long long traj)
{
    std::vector<ResidentMatrix> matrices;
    matrices.reserve(input.size());
    std::map<std::string,std::size_t> names;
    std::map<std::pair<std::string,std::string>,std::size_t> physicalInputs;
    hsize_t inputBytes = 0, maxWorkspace = 0;
    for (auto par: input)
    {
        require(!par.name.empty() && !par.file.empty(), "Empty matrix name/file");
        require(!names.count(par.name), "Duplicate matrix name: " + par.name);
        par.file = canonicalFile(trajectoryName(par.file,traj));
        par.group = groupPath(trajectoryName(par.group,traj));
        const auto key = std::make_pair(par.file,par.group);
        auto existing = physicalInputs.find(key);
        if (existing != physicalInputs.end())
        {
            names.emplace(par.name,existing->second);
            std::cout << "Alias '" << par.name << "' reuses resident input '" << matrices[existing->second].par.name << "'" << std::endl;
            continue;
        }
        ResidentMatrix m; m.par = par; inspect(m);
        inputBytes = checkedAdd(inputBytes,checkedMul(m.nt,matrixBytes(m.rows,m.cols)));
        names.emplace(par.name,matrices.size()); physicalInputs.emplace(key,matrices.size());
        matrices.push_back(std::move(m));
    }
    std::vector<ProductPlan> plans;
    std::set<std::string> outputs;
    for (auto par: products)
    {
        par.output = trajectoryName(par.output,traj);
        par.outputGroup = groupPath(trajectoryName(par.outputGroup,traj));
        require(!par.output.empty() && !exists(par.output), "Empty/existing output (no overwrite): " + par.output);
        require(outputs.insert(par.output).second, "Duplicate product output: " + par.output);
        const auto plan = planProduct(par,names,matrices);
        const auto &l = matrices[plan.left]; const auto &r = matrices[plan.right];
        hsize_t bytes = matrixBytes(l.rows,r.cols);
        if (par.ifDenseL) bytes = checkedAdd(bytes,matrixBytes(plan.active,r.cols));
        if (par.ifDenseR) bytes = checkedAdd(bytes,matrixBytes(l.rows,plan.active));
        maxWorkspace = std::max(maxWorkspace,bytes);
        plans.push_back(plan);
    }
    asSize(checkedAdd(inputBytes,maxWorkspace));
    std::cout << "Resident input payload: " << inputBytes << " bytes (" << double(inputBytes)/(1024.*1024.*1024.)
              << " GiB); " << matrices.size() << " distinct file/group inputs\n"
              << "Largest result/selection workspace: " << maxWorkspace << " bytes; I/O tile <= "
              << ioTile*ioTile*sizeof(ComplexValue) << " bytes (Eigen/HDF5 internal buffers additional)\n"
              << "No automatic disk spill or input eviction." << std::endl;
    for (auto &m: matrices) load(m);          // OUTSIDE the product loop.
    for (const auto &plan: plans) product(plan,matrices,traj);
} // Resident inputs are freed here, before the next trajectory.

void checkXmlFields(const std::string &filename)
{
    XmlReader reader(filename);
    auto fields = [&](const std::vector<std::string> &names)
    {
        for (const auto &name: names)
        {
            require(reader.push(name), "Missing XML element: " + name);
            require(!reader.nextElement(name), "Duplicate XML element: " + name);
            reader.pop();
        }
    };
    fields({"global","a2aMatrix","product"});
    reader.push("global"); fields({"trajCounter"}); reader.push("trajCounter");
    fields({"start","end","step"}); reader.pop(); reader.pop();
    reader.push("a2aMatrix"); require(reader.push("elem"),"No <a2aMatrix><elem> entries");
    do { fields({"name","file","group"}); } while (reader.nextElement("elem"));
    reader.pop(); reader.pop();
    reader.push("product"); require(reader.push("elem"),"No <product><elem> entries");
    do { fields({"mfL","mfR","deltaL","deltaR","ifDenseL","ifDenseR","initDense","output","outputGroup"}); }
    while (reader.nextElement("elem"));
}
void checkSingleProcess(void)
{
    const char *size = nullptr;
    if (std::getenv("OMPI_COMM_WORLD_RANK")) size = std::getenv("OMPI_COMM_WORLD_SIZE");
    else if (std::getenv("PMI_RANK")) size = std::getenv("PMI_SIZE");
    else if (std::getenv("PMIX_RANK")) size = std::getenv("PMIX_UNIV_SIZE");
    else if (std::getenv("SLURM_PROCID")) size = std::getenv("SLURM_NTASKS");
    if (size)
    {
        char *end = nullptr; errno = 0;
        const auto n = std::strtoull(size,&end,10);
        require(errno == 0 && end != size && *end == '\0' && n == 1,
                "Single-process utility: use srun -n 1, not multiple ranks on the same XML");
    }
}
} // namespace A2AMesonFieldProduct

int main(int argc, char *argv[])
{
    using namespace A2AMesonFieldProduct;
    if (argc != 2)
    {
        std::cerr << "usage: " << argv[0] << " <parameter file>" << std::endl;
        return EXIT_FAILURE;
    }
    H5NS::Exception::dontPrint();
    try
    {
        checkSingleProcess(); checkXmlFields(argv[1]);
        GlobalPar global;
        std::vector<MatrixPar> a2aMatrix;
        std::vector<ProductPar> products;
        XmlReader reader(argv[1]);
        read(reader,"global",global); read(reader,"a2aMatrix",a2aMatrix); read(reader,"product",products);
        const auto &range = global.trajCounter;
        require(range.start >= 0 && range.end > range.start && range.step > 0,
                "trajCounter requires 0 <= start < end, step > 0 (end is exclusive)");
        require(!a2aMatrix.empty() && !products.empty(), "a2aMatrix and product lists must be nonempty");
        Eigen::initParallel();
        std::cout << "A2AMesonFieldProduct: one process; Eigen reports " << Eigen::nbThreads()
                  << " thread(s); no lattice Grid or DiskVector" << std::endl;
        for (long long traj = range.start; traj < range.end; )
        {
            std::cout << ":::::::: Trajectory " << traj << std::endl;
            runTrajectory(a2aMatrix,products,traj);
            if (range.step >= range.end-traj) break;
            traj += range.step;
        }
    }
    catch (const H5NS::Exception &e)
    {
        std::cerr << "HDF5 error: " << e.getDetailMsg() << std::endl;
        return EXIT_FAILURE;
    }
    catch (const std::exception &e)
    {
        std::cerr << "Error: " << e.what() << std::endl;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
