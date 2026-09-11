/*
 * A2AMesonFieldSelectHits.cpp, part of Hadrons
 *
 * The XML/trajectory-loop structure is derived from utilities/Contractor.cpp.
 * The HDF5 layout and hyperslab I/O follow Hadrons/A2AMatrix.hpp.
 * Original Hadrons code: Copyright (C) 2015-2023, Antonin Portelli and contributors.
 *
 * Hadrons is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 2 of the License, or (at your option)
 * any later version. This program is distributed WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
 * PARTICULAR PURPOSE. See the GNU General Public License for more details.
 * See LICENSE in the top-level Hadrons distribution.
 */
#include <Hadrons/Global.hpp>
#include <Hadrons/A2AMatrix.hpp>
#include <algorithm>
#include <cerrno>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <limits>
#include <set>
#include <sstream>
#include <stdexcept>
#include <sys/stat.h>
#include <unistd.h>

#ifndef HAVE_HDF5
#error "A2AMesonFieldSelectHits requires HDF5 support in Grid/Hadrons"
#endif

using namespace Grid;
using namespace Hadrons;

namespace A2AMesonFieldSelectHits
{
// Parameters follow Contractor.cpp. No nt, cacheSize, lattice or MPI layout.
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
class AxisPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(AxisPar,
                                    long long,              lowSize,
                                    long long,              modesPerHit,
                                    long long,              nHitIn,
                                    bool,                   isV,
                                    std::vector<long long>, hits);
    AxisPar(): lowSize(-1), modesPerHit(-1), nHitIn(-1), isV(false) {}
};
class SelectionPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(SelectionPar,
                                    std::string, file,
                                    std::string, dataset,
                                    std::string, output,
                                    AxisPar,     left,
                                    AxisPar,     right);
};

// Contiguous input/output runs, not a new I/O abstraction.
struct ModeRun
{
    hsize_t input, output, size;
    double  scale;
};
struct AxisPlan
{
    std::vector<ModeRun> runs;
    hsize_t size;
    double  highScale;
};
struct ComplexValue
{
    double re, im;
};

// These bound our data buffer, not GEMM blocks or physics parameters.
static constexpr hsize_t ioTile     = 512; // <= 4 MiB of ComplexValue
static constexpr hsize_t chunkSize  = 256; // output chunk = (1, <=256, <=256)
static const std::string historyName = "_A2AHitSelectionHistory";

void require(bool ok, const std::string &message)
{
    if (!ok) throw std::runtime_error(message);
}

hsize_t checkedAdd(hsize_t a, hsize_t b)
{
    require(b <= std::numeric_limits<hsize_t>::max() - a,
            "Integer overflow in dimension/offset addition");
    return a + b;
}

hsize_t checkedMul(hsize_t a, hsize_t b)
{
    require(a == 0 || b <= std::numeric_limits<hsize_t>::max()/a,
            "Integer overflow in dimension/byte-count multiplication");
    return a*b;
}

std::size_t asSize(hsize_t n)
{
    require(n <= std::numeric_limits<std::size_t>::max(),
            "Size does not fit std::size_t");
    return static_cast<std::size_t>(n);
}

bool exists(const std::string &path)
{
    struct stat st;
    if (::lstat(path.c_str(), &st) == 0) return true;
    require(errno == ENOENT, "Cannot inspect '" + path + "': " + std::strerror(errno));
    return false;
}

// The token is @traj@ (Contractor convention), not bare 'traj'.
std::string trajectoryName(std::string name, long long traj)
{
    while (name.find("@traj@") != std::string::npos)
        tokenReplace(name, "traj", traj);
    return name;
}

std::string groupPath(const std::string &name)
{
    std::stringstream in(name);
    std::string part, path;
    while (std::getline(in, part, '/'))
    {
        if (part.empty()) continue;
        require(part != "." && part != "..", "dataset must be a normal HDF5 group path");
        path += "/" + part;
    }
    return path.empty() ? "/" : path;
}

void addRun(AxisPlan &plan, hsize_t input, hsize_t size, double scale)
{
    if (size == 0) return;
    if (!plan.runs.empty() && plan.runs.back().scale == scale &&
        checkedAdd(plan.runs.back().input, plan.runs.back().size) == input)
    {
        plan.runs.back().size = checkedAdd(plan.runs.back().size, size);
    }
    else plan.runs.push_back({input, plan.size, size, scale});
    plan.size = checkedAdd(plan.size, size);
}

AxisPlan makePlan(const AxisPar &par, hsize_t inputSize, const std::string &side)
{
    require(par.lowSize >= 0 && par.nHitIn >= 0 && par.modesPerHit >= 0,
            side + ": lowSize, modesPerHit and nHitIn must be non-negative");
    require(par.nHitIn == 0 || par.modesPerHit > 0,
            side + ": modesPerHit must be positive when nHitIn > 0");
    const hsize_t low = static_cast<hsize_t>(par.lowSize);
    const hsize_t q   = static_cast<hsize_t>(par.modesPerHit);
    const hsize_t n   = static_cast<hsize_t>(par.nHitIn);
    require(checkedAdd(low, checkedMul(n, q)) == inputSize,
            side + ": input dimension " + std::to_string(inputSize)
            + " != lowSize + nHitIn*modesPerHit");

    std::set<long long> seen;
    for (const auto h: par.hits)
    {
        require(h >= 0 && h < par.nHitIn, side + ": hit index out of range: " + std::to_string(h));
        require(seen.insert(h).second, side + ": duplicate hit index: " + std::to_string(h));
    }

    AxisPlan plan{{}, 0, 1.};
    if (par.isV && !par.hits.empty())
        plan.highScale = static_cast<double>(par.nHitIn)/static_cast<double>(par.hits.size());
    addRun(plan, 0, low, 1.);
    for (const auto h: par.hits)
        addRun(plan, checkedAdd(low, checkedMul(static_cast<hsize_t>(h), q)), q, plan.highScale);
    require(plan.size > 0, side + ": selection would give a zero-sized output axis");
    return plan;
}

// Copy attributes without assuming a particular meson-field metadata class.
void copyAttributes(const H5NS::H5Object &source, H5NS::H5Object &target)
{
    for (int i = 0; i < source.getNumAttrs(); ++i)
    {
        auto attr = source.openAttribute(static_cast<unsigned int>(i));
        auto type = attr.getDataType();
        auto space = attr.getSpace();
        require(H5Tdetect_class(type.getId(), H5T_REFERENCE) == 0,
                "HDF5 object/region-reference attributes are not supported");
        auto out = target.createAttribute(attr.getName(), type, space);
        const hssize_t points = space.getSimpleExtentNpoints();
        require(points >= 0, "Invalid attribute dataspace");
        if (points == 0) continue;
        const auto bytes = checkedMul(static_cast<hsize_t>(points), type.getSize());
        std::vector<unsigned char> buf(asSize(bytes));
        attr.read(type, buf.data());
        try { out.write(type, buf.data()); }
        catch (...)
        {
            H5Dvlen_reclaim(type.getId(), space.getId(), H5P_DEFAULT, buf.data());
            throw;
        }
        require(H5Dvlen_reclaim(type.getId(), space.getId(), H5P_DEFAULT, buf.data()) >= 0,
                "Could not release variable-length attribute memory");
    }
}

// Standard Hadrons files have one matrix and a metadata tree. No full-matrix copy.
// Reject soft/external links and reference-valued metadata rather than silently
// producing dangling references after changing the mode dimensions.
void copyMetadata(const H5NS::Group &source, H5NS::Group &target,
                  const std::string &path, const std::string &matrixPath,
                  unsigned int depth = 0)
{
    require(depth < 64, "Metadata hierarchy too deep (possibly cyclic hard links)");
    copyAttributes(source, target);
    for (hsize_t i = 0; i < source.getNumObjs(); ++i)
    {
        const std::string name = source.getObjnameByIdx(i);
        const std::string childPath = path + "/" + name;
        H5L_info_t linkInfo;
        require(H5Lget_info(source.getId(), name.c_str(), &linkInfo, H5P_DEFAULT) >= 0,
                "Cannot inspect metadata link: " + childPath);
        require(linkInfo.type == H5L_TYPE_HARD, "Soft/external links are not supported: " + childPath);
        if (childPath == matrixPath) continue;
        require(name != HADRONS_A2AM_NAME,
                "Multiple meson-field datasets in one file are not supported: " + childPath);
        const auto objectType = source.getObjTypeByIdx(i);
        if (objectType == H5G_GROUP)
        {
            auto src = source.openGroup(name);
            auto dst = target.createGroup(name);
            copyMetadata(src, dst, childPath, matrixPath, depth + 1);
        }
        else if (objectType == H5G_DATASET)
        {
            auto data = source.openDataSet(name);
            require(H5Tdetect_class(data.getDataType().getId(), H5T_REFERENCE) == 0,
                    "Reference-valued metadata dataset is not supported: " + childPath);
            for (int a = 0; a < data.getNumAttrs(); ++a)
                require(H5Tdetect_class(data.openAttribute(static_cast<unsigned int>(a))
                                      .getDataType().getId(), H5T_REFERENCE) == 0,
                        "Reference-valued metadata attribute is not supported: " + childPath);
            require(H5Ocopy(source.getId(), name.c_str(), target.getId(), name.c_str(),
                            H5P_DEFAULT, H5P_DEFAULT) >= 0, "Cannot copy metadata: " + childPath);
        }
        else
        {
            require(objectType == H5G_TYPE, "Unsupported HDF5 object: " + childPath);
            require(H5Ocopy(source.getId(), name.c_str(), target.getId(), name.c_str(),
                            H5P_DEFAULT, H5P_DEFAULT) >= 0, "Cannot copy datatype: " + childPath);
        }
    }
}

// Input disk type is retained; memory I/O uses named native-double re/im members.
std::size_t checkComplexType(const H5NS::CompType &type)
{
    require(type.getNmembers() == 2, "Expected a complex compound with two members: re, im");
    const int re = type.getMemberIndex("re"), im = type.getMemberIndex("im");
    require(re >= 0 && im >= 0, "Complex datatype must contain re and im members");
    auto rtype = type.getMemberDataType(static_cast<unsigned int>(re));
    auto itype = type.getMemberDataType(static_cast<unsigned int>(im));
    const std::size_t bytes = rtype.getSize();
    require(rtype.getClass() == H5T_FLOAT && itype.getClass() == H5T_FLOAT &&
            bytes == itype.getSize() && (bytes == 4 || bytes == 8),
            "Only ComplexF/ComplexD files with equal float/double re and im members are supported");
    return bytes;
}

void saveHistory(H5NS::Group &group, const SelectionPar &par, long long traj,
                 const AxisPlan &left, const AxisPlan &right)
{
    H5NS::Group history = H5Lexists(group.getId(), historyName.c_str(), H5P_DEFAULT) > 0
                       ? group.openGroup(historyName) : group.createGroup(historyName);
    std::string name;
    for (unsigned long long i = 0; ; ++i)
    {
        std::ostringstream s;
        s << "step" << std::setw(4) << std::setfill('0') << i;
        name = s.str();
        if (H5Lexists(history.getId(), name.c_str(), H5P_DEFAULT) == 0) break;
    }
    XmlWriter writer("");
    write(writer, "formatVersion", 1);
    write(writer, "trajectory", traj);
    write(writer, "selection", par);
    write(writer, "nHitOutLeft", static_cast<unsigned long long>(par.left.hits.size()));
    write(writer, "nHitOutRight", static_cast<unsigned long long>(par.right.hits.size()));
    write(writer, "highScaleLeft", left.highScale);
    write(writer, "highScaleRight", right.highScale);
    write(writer, "normalization", std::string("Input high V is normalized by 1/nHitIn; high W is unscaled."));
    write(writer, "scope", std::string("External mode axes only; internal loops and time axis are unchanged."));
    const std::string xml = writer.docString();
    H5NS::StrType type(H5NS::PredType::C_S1, H5T_VARIABLE);
    type.setCset(H5T_CSET_UTF8);
    H5NS::DataSpace scalar(H5S_SCALAR);
    auto data = history.createDataSet(name, type, scalar);
    data.write(xml, type);
}

void selectHits(const SelectionPar &par, long long traj)
{
    require(!par.file.empty() && !par.output.empty() && !par.dataset.empty(),
            "file, output and dataset are required");
    require(!exists(par.output), "Refusing to overwrite existing output: " + par.output);
    const std::string groupName = groupPath(par.dataset);
    const std::string matrixPath = (groupName == "/" ? "" : groupName) + "/" + HADRONS_A2AM_NAME;
    const std::string temporary = par.output + ".partial." + std::to_string(::getpid());
    bool ownTemporary = false;

    try
    {
        {
            H5NS::H5File input(par.file, H5F_ACC_RDONLY);
            auto inputData = input.openDataSet(matrixPath);
            auto inputSpace = inputData.getSpace();
            require(inputSpace.getSimpleExtentNdims() == 3, "Expected dimensions (t, left mode, right mode)");
            hsize_t inputDim[3];
            inputSpace.getSimpleExtentDims(inputDim);
            require(inputDim[0] > 0 && inputDim[1] > 0 && inputDim[2] > 0, "Empty input matrix");
            auto diskType = inputData.getCompType();
            const std::size_t realBytes = checkComplexType(diskType);
            const auto left  = makePlan(par.left,  inputDim[1], "left");
            const auto right = makePlan(par.right, inputDim[2], "right");
            const hsize_t outputDim[3] = {inputDim[0], left.size, right.size};
            const hsize_t outputBytes = checkedMul(checkedMul(checkedMul(outputDim[0], outputDim[1]),
                                                             outputDim[2]), diskType.getSize());

            std::cout << "======== Selecting '" << par.file << "' -> '" << par.output << "'\n"
                      << "shape: " << inputDim[0] << " x " << inputDim[1] << " x " << inputDim[2]
                      << " -> " << outputDim[0] << " x " << outputDim[1] << " x " << outputDim[2] << '\n'
                      << "high-mode scales: left=" << left.highScale << ", right=" << right.highScale << '\n'
                      << "output payload: " << outputBytes << " bytes; precision: Complex"
                      << (realBytes == 4 ? "F" : "D") << std::endl;

            makeFileDir(par.output);
            H5NS::H5File output(temporary, H5F_ACC_EXCL);
            ownTemporary = true;
            {
                auto inRoot = input.openGroup("/");
                auto outRoot = output.openGroup("/");
                copyMetadata(inRoot, outRoot, "", matrixPath);
                auto outGroup = output.openGroup(groupName);
                H5NS::DataSpace outputSpace(3, outputDim);
                const hsize_t chunk[3] = {1, std::min(chunkSize, left.size), std::min(chunkSize, right.size)};
                H5NS::DSetCreatPropList plist;
                plist.setChunk(3, chunk);
                plist.setFletcher32();
                auto outputData = outGroup.createDataSet(HADRONS_A2AM_NAME, diskType, outputSpace, plist);
                copyAttributes(inputData, outputData);

                H5NS::CompType memoryType(sizeof(ComplexValue));
                memoryType.insertMember("re", HOFFSET(ComplexValue, re), H5NS::PredType::NATIVE_DOUBLE);
                memoryType.insertMember("im", HOFFSET(ComplexValue, im), H5NS::PredType::NATIVE_DOUBLE);
                const double largest = realBytes == 4 ? static_cast<double>(std::numeric_limits<float>::max())
                                                      : std::numeric_limits<double>::max();
                std::vector<ComplexValue> buf(asSize(checkedMul(std::min(ioTile, left.size),
                                                                std::min(ioTile, right.size))));

                // Each run is entirely low or high (unless the coefficient is identical).
                // No transpose, conjugation, temporal selection, or internal-loop change.
                const auto begin = std::chrono::steady_clock::now();
                for (hsize_t t = 0; t < inputDim[0]; ++t)
                {
                    for (const auto &l: left.runs)
                    for (const auto &r: right.runs)
                    {
                        const double scale = l.scale*r.scale;
                        for (hsize_t i = 0; i < l.size; )
                        {
                            const hsize_t nr = std::min(ioTile, l.size - i);
                            for (hsize_t j = 0; j < r.size; )
                            {
                                const hsize_t nc = std::min(ioTile, r.size - j);
                                const hsize_t count[3] = {1, nr, nc};
                                const hsize_t inOffset[3] = {t, l.input + i, r.input + j};
                                const hsize_t outOffset[3] = {t, l.output + i, r.output + j};
                                H5NS::DataSpace memSpace(3, count);
                                inputSpace.selectHyperslab(H5S_SELECT_SET, count, inOffset);
                                outputSpace.selectHyperslab(H5S_SELECT_SET, count, outOffset);
                                inputData.read(buf.data(), memoryType, memSpace, inputSpace);
                                const std::size_t n = asSize(checkedMul(nr, nc));
                                for (std::size_t k = 0; k < n; ++k)
                                {
                                    if (scale != 1.) { buf[k].re *= scale; buf[k].im *= scale; }
                                    if (!std::isfinite(buf[k].re) || !std::isfinite(buf[k].im) ||
                                        std::abs(buf[k].re) > largest || std::abs(buf[k].im) > largest)
                                        throw std::runtime_error("Non-finite input or normalization overflow at t="
                                                                 + std::to_string(t));
                                }
                                outputData.write(buf.data(), memoryType, memSpace, outputSpace);
                                j += nc;
                            }
                            i += nr;
                        }
                    }
                    const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - begin).count();
                    std::cout << "  timeslice " << t << "/" << inputDim[0] - 1 << " done ("
                              << seconds << " s)" << std::endl;
                }
                saveHistory(outGroup, par, traj, left, right);
            }
            output.flush(H5F_SCOPE_GLOBAL);
            output.close();
        } // Close every HDF5 handle before publishing the completed file.

        // Same-directory hard link: atomic publication with no overwrite race.
        // Unlike rename(), link() fails if another process already made output.
        if (::link(temporary.c_str(), par.output.c_str()) != 0)
        {
            const int error = errno;
            ownTemporary = false;
            throw std::runtime_error("Cannot publish completed output: " + std::string(std::strerror(error))
                                     + ". Completed temporary file retained: " + temporary);
        }
        if (::unlink(temporary.c_str()) != 0)
            std::cerr << "Warning: output is complete, but could not remove " << temporary << std::endl;
        ownTemporary = false;
        std::cout << "Completed: " << par.output << std::endl;
    }
    catch (...)
    {
        // On a normal exception discard our partial result, never an input or
        // someone else's output. SIGKILL can leave a *.partial.PID file behind.
        if (ownTemporary) ::unlink(temporary.c_str());
        throw;
    }
}

// Grid's generic reader warns rather than failing for some missing nodes.
// In particular, an absent <hits> must not silently mean "keep low modes only".
void checkXmlFields(const std::string &filename)
{
    XmlReader reader(filename);
    auto fields = [&](const std::vector<std::string> &names)
    {
        for (const auto &name: names)
        {
            require(reader.push(name), "Missing XML element: " + name);
            reader.pop();
        }
    };
    fields({"global", "selection"});
    reader.push("global");
    require(reader.push("trajCounter"), "Missing XML element: trajCounter");
    fields({"start", "end", "step"});
    reader.pop(); reader.pop();
    reader.push("selection");
    require(reader.push("elem"), "No <selection><elem> entries");
    do
    {
        fields({"file", "dataset", "output", "left", "right"});
        for (const auto &side: {"left", "right"})
        {
            reader.push(side);
            fields({"lowSize", "modesPerHit", "nHitIn", "isV", "hits"});
            reader.pop();
        }
    } while (reader.nextElement("elem"));
}

void checkSingleProcess(void)
{
    // Do not reject a batch allocation merely because SLURM_NTASKS is large.
    // Only check launcher variables when they identify a task in a step.
    const char *size = nullptr;
    if (std::getenv("OMPI_COMM_WORLD_RANK")) size = std::getenv("OMPI_COMM_WORLD_SIZE");
    else if (std::getenv("PMI_RANK")) size = std::getenv("PMI_SIZE");
    else if (std::getenv("SLURM_PROCID")) size = std::getenv("SLURM_NTASKS");
    if (size) require(std::strtoull(size, nullptr, 10) <= 1,
                      "This is a single-process utility. Do not launch the same XML on multiple ranks.");
}
} // namespace A2AMesonFieldSelectHits

int main(int argc, char *argv[])
{
    using namespace A2AMesonFieldSelectHits;
    if (argc != 2)
    {
        std::cerr << "usage: " << argv[0] << " <parameter file>" << std::endl;
        return EXIT_FAILURE;
    }
    H5NS::Exception::dontPrint();
    try
    {
        checkSingleProcess();
        checkXmlFields(argv[1]);
        GlobalPar global;
        std::vector<SelectionPar> selection;
        XmlReader reader(argv[1]);
        read(reader, "global", global);
        read(reader, "selection", selection);
        const auto &range = global.trajCounter;
        require(range.start >= 0 && range.end > range.start && range.step > 0,
                "trajCounter requires 0 <= start < end and step > 0 (end is exclusive)");
        require(!selection.empty(), "No <selection><elem> entries");
        for (long long traj = range.start; traj < range.end; )
        {
            std::cout << ":::::::: Trajectory " << traj << std::endl;
            for (auto par: selection)
            {
                par.file    = trajectoryName(par.file, traj);
                par.dataset = trajectoryName(par.dataset, traj);
                par.output  = trajectoryName(par.output, traj);
                selectHits(par, traj);
            }
            if (range.step >= range.end - traj) break;
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
