/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: utilities/BubbleContractorDWC.cpp

Copyright (C) 2015-2019

Author: Antonin Portelli <antonin.portelli@me.com>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
See the full license in the file "LICENSE" in the top level distribution directory
*************************************************************************************/
/*  END LEGAL */
#include <Hadrons/Global.hpp>
#include <Hadrons/A2AMatrix.hpp>
#include <Hadrons/DiskVector.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/TimerArray.hpp>
#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#ifndef HAVE_HDF5
#error "BubbleContractorDWC requires HDF5 support"
#endif

// Based on BubbleContractor.cpp. Keep one bubble per t0 (no time averaging).
// Local dense-W connection/read helpers follow FlexibleContractorDWC.
// No Grid_init(), lattice Grid, extra conjugation, or hit normalization.

using namespace Grid;
using namespace Hadrons;
#define TIME_MOD(t) timeModulo((t), par.global.nt)
namespace BubbleContractorDWC
{
    class TrajRange: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(TrajRange,
                                        unsigned int, start,
                                        unsigned int, end,
                                        unsigned int, step);
        TrajRange(): start(0), end(0), step(0) {}
    };

    class GlobalPar: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(GlobalPar,
                                        TrajRange, trajCounter,
                                        unsigned int, nt,
                                        std::string, diskVectorDir,
                                        std::string, output);
        GlobalPar(): nt(0) {}
    };
    class A2AMatrixPar: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(A2AMatrixPar,
                                        std::string, file,
                                        std::string, dataset,
                                        unsigned int, cacheSize,
                                        std::string, name,
                                        int, initDenseRow,
                                        int, initDenseCol,
                                        int, timeOffsetRow,
                                        int, timeOffsetCol);
        A2AMatrixPar(): cacheSize(0), initDenseRow(-1), initDenseCol(-1),
                        timeOffsetRow(0), timeOffsetCol(0) {}
    };
    class ProductPar: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(ProductPar,
                                        std::string, terms,
                                        std::vector<unsigned int>, times);
    };
    class CorrelatorResult: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(CorrelatorResult,
                                        std::vector<unsigned int>, times,
                                        std::vector<ComplexD>, correlator);
    };
}
struct BubbleContractorDWCPar
{
    BubbleContractorDWC::GlobalPar                  global;
    std::vector<BubbleContractorDWC::A2AMatrixPar>  a2aMatrix;
    std::vector<BubbleContractorDWC::ProductPar>    product;
};
void saveCorrelator(const BubbleContractorDWC::CorrelatorResult &result,
		    const std::string dir,
		    const std::string fileStem,
                    const unsigned int traj)
{
    std::string filename;
    filename = dir + "/" + ModuleBase::resultFilename(fileStem, traj);
    std::cout << "Saving correlator to '" << filename << "'" << std::endl;
    makeFileDir(filename);
    ResultWriter writer(filename);
    write(writer, fileStem, result);
}
struct Sec
{
    Sec(const double usec)
    {
        seconds = usec/1.0e6;
    }

    double seconds;
};
inline std::ostream & operator<< (std::ostream& s, const Sec &&sec)
{
    s << std::setw(10) << sec.seconds << " sec";

    return s;
}

struct Flops
{
    Flops(const double flops, const double fusec)
    {
        gFlopsPerSec = fusec > 0 ? flops/fusec/1.0e3 : 0.;
    }

    double gFlopsPerSec;
};

inline std::ostream & operator<< (std::ostream& s, const Flops &&f)
{
    s << std::setw(10) << f.gFlopsPerSec << " GFlop/s";

    return s;
}
struct Bytes
{
    Bytes(const double bytes, const double busec)
    {
        gBytesPerSec = busec > 0 ? bytes/busec*1.0e6/1024/1024/1024 : 0.;
    }

    double gBytesPerSec;
};

inline std::ostream & operator<< (std::ostream& s, const Bytes &&b)
{
    s << std::setw(10) << b.gBytesPerSec << " GB/s";

    return s;
}

// Small local helpers; no new shared matrix class or changes to A2AMatrix.hpp.
void require(bool ok, const std::string &message)
{
  if (!ok) { HADRONS_ERROR(Argument, "BubbleContractorDWC: " + message); }
}

unsigned int timeModulo(std::int64_t t, unsigned int nt)
{
  const auto n = static_cast<std::int64_t>(nt);
  return static_cast<unsigned int>((t % n + n) % n);
}

std::size_t checkedMul(std::size_t a, std::size_t b)
{
  require(!b || a <= std::numeric_limits<std::size_t>::max()/b,
          "dimension/byte-count multiplication overflow");
  return a*b;
}

void checkMatrixSize(Eigen::Index rows, Eigen::Index cols)
{
  require(rows > 0 && cols > 0, "empty matrix dimension");
  const auto n = checkedMul(std::size_t(rows), std::size_t(cols));
  require(n <= std::size_t(std::numeric_limits<Eigen::Index>::max()),
          "matrix exceeds Eigen::Index");
  checkedMul(n, sizeof(ComplexD));
}

struct AxisInfo
{
  int initDense;                  // -1 = ordinary; >=0 = dense, including zero lows.
  unsigned int time;             // Already normalized source time, not file time.
};
struct MatrixAxes { AxisInfo row, col; };
struct MatrixDimensions { Eigen::Index rows, cols; };
using Matrix = A2AMatrix<ComplexD>;
using LastMatrix = A2AMatrixTr<ComplexD>;
using IndexList = std::vector<Eigen::Index>;

MatrixAxes matrixAxes(const BubbleContractorDWC::A2AMatrixPar &p,
                      unsigned int t, unsigned int nt)
{
  // Promote BEFORE modulo/addition: a negative offset must not become unsigned.
  return {{p.initDenseRow, timeModulo(std::int64_t(t) + std::int64_t(p.timeOffsetRow) % std::int64_t(nt), nt)},
          {p.initDenseCol, timeModulo(std::int64_t(t) + std::int64_t(p.timeOffsetCol) % std::int64_t(nt), nt)}};
}

// Empty selection means identity (no allocation/copy for ordinary connections).
struct Connection
{
  Eigen::Index active;
  IndexList left, right;
};

Connection connection(Eigen::Index nLeft, AxisInfo left,
                      Eigen::Index nRight, AxisInfo right, unsigned int nt)
{
  require(!(left.initDense >= 0 && right.initDense >= 0),
          "both ends of a contracted connection are dense W");
  if (left.initDense < 0 && right.initDense < 0)
  {
    require(nLeft == nRight, "non-dense contracted dimensions differ");
    return {nLeft, {}, {}};
  }
  const bool denseLeft = left.initDense >= 0;
  const AxisInfo dense = denseLeft ? left : right;
  const auto compact = denseLeft ? nLeft : nRight;
  const auto expanded = denseLeft ? nRight : nLeft;
  const std::size_t low = static_cast<std::size_t>(dense.initDense);
  const std::size_t nsc = std::size_t(Ns)*Nc;
  require(compact >= 0 && expanded >= 0 && std::size_t(compact) >= low &&
          std::size_t(expanded) >= low, "initDense exceeds an axis dimension");
  const auto high = std::size_t(compact) - low;
  require(high % nsc == 0, "dense high dimension is not a multiple of Ns*Nc");
  const auto fullHigh = checkedMul(high, nt);
  require(fullHigh <= std::numeric_limits<std::size_t>::max()-low &&
          std::size_t(expanded) == low + fullHigh,
          "dense/expanded contracted dimensions do not match");
  Connection c{compact, {}, {}};
  if (expanded == compact) return c; // Low-only or Nt=1: identity, still checked.
  IndexList &index = denseLeft ? c.right : c.left;
  index.resize(std::size_t(compact));
  for (std::size_t r = 0; r < low; ++r) index[r] = Eigen::Index(r);
  // Integer divisions/mapping are not in any matrix-element multiplication loop.
  for (std::size_t h = 0; h < high/nsc; ++h)
  for (std::size_t a = 0; a < nsc; ++a)
    index[low+h*nsc+a] = Eigen::Index(low+(h*nt+dense.time)*nsc+a);
  return c;
}

// Copy selected axes only. Branches are outside component loops. Output storage
// order is respected; last terms remain column-major without a mathematical transpose.
template <typename Out, typename In>
void selectMatrix(Out &out, const In &in, const IndexList &rows, const IndexList &cols)
{
  const Eigen::Index nr = rows.empty() ? in.rows() : Eigen::Index(rows.size());
  const Eigen::Index nc = cols.empty() ? in.cols() : Eigen::Index(cols.size());
  checkMatrixSize(nr, nc);
  out.resize(nr, nc);
  if (rows.empty() && cols.empty()) out = in;
  else if (rows.empty())
    for (Eigen::Index c = 0; c < nc; ++c) out.col(c) = in.col(cols[std::size_t(c)]);
  else if (cols.empty())
    for (Eigen::Index r = 0; r < nr; ++r) out.row(r) = in.row(rows[std::size_t(r)]);
  else if constexpr (Out::IsRowMajor)
    for (Eigen::Index r = 0; r < nr; ++r)
    for (Eigen::Index c = 0; c < nc; ++c)
      out(r,c) = in(rows[std::size_t(r)], cols[std::size_t(c)]);
  else
    for (Eigen::Index c = 0; c < nc; ++c)
    for (Eigen::Index r = 0; r < nr; ++r)
      out(r,c) = in(rows[std::size_t(r)], cols[std::size_t(c)]);
}

double multiplyDWC(Matrix &out, const Matrix &a, const Matrix &b,
                   const MatrixAxes &axesA, const MatrixAxes &axesB,
                   unsigned int nt, Matrix &selection)
{
  const auto c = connection(a.cols(), axesA.col, b.rows(), axesB.row, nt);
  checkMatrixSize(a.rows(), b.cols());
  if (!c.left.empty())
  {
    selectMatrix(selection, a, {}, c.left);
    A2AContraction::mul(out, selection, b);
  }
  else if (!c.right.empty())
  {
    selectMatrix(selection, b, c.right, {});
    A2AContraction::mul(out, a, selection);
  }
  else A2AContraction::mul(out, a, b);
  return double(a.rows())*double(b.cols())*(8.*double(c.active)-2.);
}

double traceDWC(ComplexD &acc, const Matrix &a, const LastMatrix &b,
               const MatrixAxes &axesA, const MatrixAxes &axesB,
               unsigned int nt, Matrix &selectedA, LastMatrix &selectedB)
{
  const auto inner = connection(a.cols(), axesA.col, b.rows(), axesB.row, nt);
  const auto closing = connection(b.cols(), axesB.col, a.rows(), axesA.row, nt);
  const Matrix *ap = &a;
  const LastMatrix *bp = &b;
  if (!closing.right.empty() || !inner.left.empty())
  {
    selectMatrix(selectedA, a, closing.right, inner.left);
    ap = &selectedA;
  }
  if (!inner.right.empty() || !closing.left.empty())
  {
    selectMatrix(selectedB, b, inner.right, closing.left);
    bp = &selectedB;
  }
  require(ap->cols() == bp->rows() && ap->rows() == bp->cols(),
          "internal error in selected trace dimensions");
  // This is tr(a*b), NOT an inner product with conjugation. Grid's dot()
  // path conjugates its first argument twice intentionally.
  A2AContraction::accTrMul(acc, *ap, *bp);
  return 8.*double(ap->rows())*double(ap->cols());
}

// Close a single compact matrix without constructing an identity or a full
// expanded matrix. The left endpoint is col, the right endpoint is row.
double traceOneDWC(ComplexD &acc, const Matrix &m, const MatrixAxes &axes,
                   unsigned int nt)
{
  const auto c = connection(m.cols(), axes.col, m.rows(), axes.row, nt);
  if (!c.left.empty())
    for (Eigen::Index r = 0; r < c.active; ++r)
      acc += m(r, c.left[std::size_t(r)]);
  else if (!c.right.empty())
    for (Eigen::Index r = 0; r < c.active; ++r)
      acc += m(c.right[std::size_t(r)], r);
  else
    for (Eigen::Index r = 0; r < c.active; ++r) acc += m(r,r);
  return 2.*double(c.active);
}

// Read F or D independently of HADRONS_A2AM_IO_TYPE. The memory compound is
// explicitly the layout of ComplexValue, never the datatype read from the file.
struct ComplexValue { double re, im; };
MatrixDimensions loadMatrix(const std::string &filename, const std::string &group,
                            unsigned int nt, EigenDiskVector<ComplexD> &v,
                            double &tRead, std::size_t &fileBytes)
{
  H5NS::H5File file(filename, H5F_ACC_RDONLY);
  auto g = file.openGroup(group);
  auto data = g.openDataSet(HADRONS_A2AM_NAME);
  auto space = data.getSpace();
  require(space.getSimpleExtentNdims() == 3, "expected (time,row,col) in " + filename);
  hsize_t d[3]; space.getSimpleExtentDims(d);
  require(d[0] == nt && d[1] > 0 && d[2] > 0, "invalid matrix/time dimensions in " + filename);
  const auto limit = std::uintmax_t(std::numeric_limits<Eigen::Index>::max());
  require(std::uintmax_t(d[1]) <= limit && std::uintmax_t(d[2]) <= limit,
          "matrix dimension exceeds Eigen::Index");
  checkMatrixSize(Eigen::Index(d[1]), Eigen::Index(d[2]));
  require(data.getTypeClass() == H5T_COMPOUND, "expected complex compound in " + filename);
  const auto diskType = data.getCompType();
  const int re = H5Tget_member_index(diskType.getId(), "re");
  const int im = H5Tget_member_index(diskType.getId(), "im");
  require(diskType.getNmembers() == 2 && re >= 0 && im >= 0,
          "expected re/im complex members in " + filename);
  const auto rt = diskType.getMemberDataType(unsigned(re));
  const auto it = diskType.getMemberDataType(unsigned(im));
  require(rt.getClass() == H5T_FLOAT && it.getClass() == H5T_FLOAT &&
          rt.getSize() == it.getSize() && (rt.getSize() == 4 || rt.getSize() == 8),
          "expected float/double re and im in " + filename);
  H5NS::CompType memType(sizeof(ComplexValue));
  memType.insertMember("re", HOFFSET(ComplexValue,re), H5NS::PredType::NATIVE_DOUBLE);
  memType.insertMember("im", HOFFSET(ComplexValue,im), H5NS::PredType::NATIVE_DOUBLE);
  fileBytes = checkedMul(checkedMul(checkedMul(nt, std::size_t(d[1])), std::size_t(d[2])), diskType.getSize());
  const hsize_t tile = 512; // Bounded conversion buffer; not a physics/cache parameter.
  std::vector<ComplexValue> values(checkedMul(std::min(tile,d[1]), std::min(tile,d[2])));
  Matrix buf;
  buf.resize(Eigen::Index(d[1]), Eigen::Index(d[2]));
  std::cout << "  disk Complex" << (rt.getSize()==4 ? "F" : "D")
            << " -> native ComplexD; shape " << nt << " x " << d[1] << " x " << d[2] << std::endl;
  tRead = 0.;
  for (unsigned int tp1 = nt; tp1 > 0; --tp1)
  {
    const unsigned int t = tp1-1;
    tRead -= usecond();
    for (hsize_t i = 0; i < d[1]; )
    {
      const hsize_t nr = std::min(tile,d[1]-i);
      for (hsize_t j = 0; j < d[2]; )
      {
        const hsize_t nc = std::min(tile,d[2]-j);
        const hsize_t count[3] = {1,nr,nc}, offset[3] = {t,i,j};
        H5NS::DataSpace mem(3,count);
        space.selectHyperslab(H5S_SELECT_SET,count,offset);
        data.read(values.data(),memType,mem,space);
        for (hsize_t r = 0; r < nr; ++r)
        for (hsize_t c = 0; c < nc; ++c)
        {
          const auto z = values[std::size_t(r*nc+c)];
          if (!std::isfinite(z.re) || !std::isfinite(z.im))
	  {
	    HADRONS_ERROR(Argument, "non-finite input in " + filename);
	  }
          buf(Eigen::Index(i+r),Eigen::Index(j+c)) = ComplexD(z.re,z.im);
        }
        j += nc;
      }
      i += nr;
    }
    // An assignment copies one time slice into the existing DiskVector cache.
    v[t] = buf;
    tRead += usecond();
  }
  return {Eigen::Index(d[1]), Eigen::Index(d[2])};
}

std::string productStem(const BubbleContractorDWC::ProductPar &p)
{
  const auto term = strToVec<std::string>(p.terms);
  require(!term.empty() && p.times.size() == term.size()-1,
          "each product requires N terms and N-1 times");
  std::string stem;
  for (std::size_t j = 0; j+1 < term.size(); ++j)
    stem += term[j] + "_" + std::to_string(p.times[j]) + "_";
  return stem + term.back();
}

// Check raw text BEFORE Grid's numeric reader can accept only the prefix of
// "4 8" or "4..8". This is scalar validation, not a time-expression parser.
// Use a separate reader because Grid's vector reader marks consumed elem nodes.
void checkXmlFields(const std::string &filename)
{
  XmlReader reader(filename);
  const auto fields = [&](const std::vector<std::string> &names)
  {
    for (const auto &n: names)
    {
      require(reader.push(n), "missing XML element: " + n);
      require(!reader.nextElement(n), "duplicate XML element: " + n);
      reader.pop();
    }
  };
  const auto integer = [&](const std::string &text, const std::string &name, bool isSigned)
  {
    const auto tokens = strToVec<std::string>(text);
    require(tokens.size() == 1, "expected one integer for " + name);
    const auto &s = tokens[0];
    const std::size_t begin = isSigned && (s[0] == '-' || s[0] == '+') ? 1 : 0;
    require(begin < s.size() && std::all_of(s.begin()+begin, s.end(),
                                          [](char c){return c >= '0' && c <= '9';}),
            "expected " + std::string(isSigned ? "a signed" : "a nonnegative")
            + " integer for " + name);
    if (isSigned)
    {
      const auto v = std::stoll(s);
      require(v >= std::numeric_limits<int>::min() && v <= std::numeric_limits<int>::max(),
              "signed integer overflow for " + name);
    }
    else require(std::stoull(s) <= std::numeric_limits<unsigned int>::max(),
                 "unsigned integer overflow for " + name);
  };
  const auto integers = [&](const std::vector<std::string> &names, bool isSigned)
  {
    for (const auto &n: names)
    {
      std::string text; read(reader,n,text); integer(text,n,isSigned);
    }
  };
  fields({"global","a2aMatrix","product"});
  reader.push("global"); fields({"trajCounter","nt","diskVectorDir","output"});
  integers({"nt"},false);
  reader.push("trajCounter"); fields({"start","end","step"});
  integers({"start","end","step"},false); reader.pop(); reader.pop();
  reader.push("a2aMatrix"); require(reader.push("elem"), "empty a2aMatrix list");
  do
  {
    fields({"file","dataset","cacheSize","name","initDenseRow","initDenseCol","timeOffsetRow","timeOffsetCol"});
    integers({"cacheSize"},false);
    integers({"initDenseRow","initDenseCol","timeOffsetRow","timeOffsetCol"},true);
  } while (reader.nextElement("elem"));
  reader.pop(); reader.pop();
  reader.push("product"); require(reader.push("elem"), "empty product list");
  do
  {
    fields({"terms","times"});
    std::string scalar; read(reader,"times",scalar);
    require(strToVec<std::string>(scalar).empty(), "times requires <elem> entries, not a scalar");
    std::vector<std::string> rawTimes;
    read(reader,"times",rawTimes);
    for (const auto &text: rawTimes) integer(text,"times/elem",false);
  } while (reader.nextElement("elem"));
  reader.pop(); reader.pop();
}

void checkSingleProcess()
{
  const char *size = nullptr;
  if (std::getenv("OMPI_COMM_WORLD_RANK")) size = std::getenv("OMPI_COMM_WORLD_SIZE");
  else if (std::getenv("PMI_RANK")) size = std::getenv("PMI_SIZE");
  else if (std::getenv("PMIX_RANK")) size = std::getenv("PMIX_UNIV_SIZE");
  else if (std::getenv("SLURM_PROCID"))
  {
    size = std::getenv("SLURM_STEP_NUM_TASKS");
    if (!size) size = std::getenv("SLURM_NTASKS");
  }
  if (size)
  {
    char *end = nullptr; errno = 0;
    const auto n = std::strtoull(size,&end,10);
    require(errno == 0 && end != size && *end == '\0' && n == 1,
            "single-process utility: launch one task per XML");
  }
}

void validateParameters(const BubbleContractorDWCPar &par)
{
  require(par.global.nt > 0 && par.global.nt <= unsigned(std::numeric_limits<int>::max()), "invalid nt");
  require(par.global.trajCounter.step > 0 && par.global.trajCounter.end > par.global.trajCounter.start,
          "trajCounter requires start < end and step > 0");
  require(!par.global.diskVectorDir.empty() && !par.global.output.empty(), "empty diskVectorDir/output");
  require(!par.a2aMatrix.empty() && !par.product.empty(), "empty matrix/product list");
  std::set<std::string> names, outputs;
  for (const auto &m: par.a2aMatrix)
  {
    require(!m.name.empty() && !m.file.empty() && !m.dataset.empty(), "empty matrix name/file/dataset");
    require(names.insert(m.name).second, "duplicate matrix name: " + m.name);
    require(m.cacheSize >= 1 && m.cacheSize <= par.global.nt, "cacheSize must lie in 1..nt");
    require(m.initDenseRow >= -1 && m.initDenseCol >= -1, "initDense must be -1 or nonnegative");
  }
  for (const auto &p: par.product)
  {
    const auto term = strToVec<std::string>(p.terms);
    require(!term.empty() && p.times.size() == term.size()-1,
            "each product requires N >= 1 terms and N-1 times");
    for (const auto &name: term) require(names.count(name), "unknown term: " + name);
    for (auto t: p.times) require(t < par.global.nt, "times entry must lie in 0..nt-1");
    const auto stem = productStem(p);
    require(outputs.insert(stem).second, "duplicate product output filename: " + stem);
  }
}

int main(int argc, char* argv[])
{
  // parse command line
  std::string   parFilename;
  if (argc != 2)
  {
    std::cerr << "usage: " << argv[0] << " <parameter file>";
    std::cerr << std::endl;

    return EXIT_FAILURE;
  }
  parFilename = argv[1];

  try
  {
  checkSingleProcess();
  H5NS::Exception::dontPrint();
  checkXmlFields(parFilename);
  // parse parameter file
  BubbleContractorDWCPar par;
  XmlReader     reader(parFilename);

  read(reader, "global",    par.global);
  read(reader, "a2aMatrix", par.a2aMatrix);
  read(reader, "product",   par.product);
  validateParameters(par);
  Eigen::initParallel();
  std::cout << "BubbleContractorDWC: one process; Eigen reports "
            << Eigen::nbThreads() << " thread(s); no lattice Grid" << std::endl;
  std::map<std::string, BubbleContractorDWC::A2AMatrixPar> matrixPar;
  for (const auto &m: par.a2aMatrix) matrixPar.emplace(m.name,m);

  // create diskvectors
  std::map<std::string, EigenDiskVector<ComplexD>> a2aMat;

  for (auto &p: par.a2aMatrix)
  {
    std::string dirName = par.global.diskVectorDir + "/" + p.name;

    makeFileDir(dirName); // Parent only; DiskVector creates its own directory.
    a2aMat.emplace(std::piecewise_construct, std::forward_as_tuple(p.name),
                   std::forward_as_tuple(dirName, par.global.nt, p.cacheSize));
  }
  // trajectory loop
  for (unsigned int traj = par.global.trajCounter.start;
       traj < par.global.trajCounter.end; )
  {
    std::cout << ":::::::: Trajectory " << traj << std::endl;

    std::map<std::string, MatrixDimensions> dimensions;
    // load data
    for (auto &p: par.a2aMatrix)
    {
      std::string filename = p.file;
      double      t;

      tokenReplace(filename, "traj", traj);
      std::cout << "======== Loading '" << filename << "'" << std::endl;
      std::size_t fileBytes = 0;
      dimensions.emplace(p.name, loadMatrix(filename, p.dataset, par.global.nt,
                                            a2aMat.at(p.name), t, fileBytes));
      std::cout << "Read " << fileBytes << " disk bytes in " << t/1.0e6
                << " sec (includes conversion/cache assignment), "
                << (t > 0 ? fileBytes/t*1.0e6/1024/1024 : 0.) << " MB/s" << std::endl;
    }
    // contract
    std::cout << "START PRODUCT" << std::endl;
    for (auto &p: par.product)
    {
      std::vector<std::string>               term = strToVec<std::string>(p.terms);
      const auto                            &t = p.times;
      std::vector<A2AMatrixTr<ComplexD>>     lastTerm;
      A2AMatrix<ComplexD>                    prod, buf, tmp;
      LastMatrix                             selectedLast;
      MatrixAxes                             axes;
      const auto                            &matrices = a2aMat; // const []: no write proxy
      TimerArray                             tAr;
      double                                 fusec, busec, flops, bytes;
      BubbleContractorDWC::CorrelatorResult   result;
      tAr.startTimer("Total");
      std::cout << "======== Contraction tr(";
      for (std::size_t g = 0; g < term.size(); ++g)
      {
        std::cout << term[g] << ((g == term.size() - 1) ? ')' : '*');
      }
      std::cout << std::endl;
      const std::string fileStem = productStem(p);
      result.times = term.size() == 1 ? std::vector<unsigned int>{0} : p.times;
      result.correlator.resize(par.global.nt, 0.);

      // Check every edge, including the closing edge, before large temporaries.
      // Dimensions are time independent; the index selections below use real times.
      for (std::size_t j = 0; j < term.size(); ++j)
      {
        const auto k = (j+1)%term.size();
        const auto a = matrixAxes(matrixPar.at(term[j]), j ? t[j-1] : 0, par.global.nt);
        const auto b = matrixAxes(matrixPar.at(term[k]), k ? t[k-1] : 0, par.global.nt);
        connection(dimensions.at(term[j]).cols,a.col,
                   dimensions.at(term[k]).rows,b.row,par.global.nt);
      }

      if (term.size() > 1)
      {
        std::cout << "* Caching column-major last term (not a transpose)" << std::endl;
        const auto d = dimensions.at(term.back());
        const auto cacheBytes = checkedMul(checkedMul(checkedMul(par.global.nt,
                                            std::size_t(d.rows)),std::size_t(d.cols)),sizeof(ComplexD));
        std::cout << "  Additional last-term cache: " << cacheBytes << " bytes" << std::endl;
        lastTerm.resize(par.global.nt);
        for (unsigned int tp = 0; tp < par.global.nt; ++tp)
        {
          tAr.startTimer("Disk vector overhead");
          const A2AMatrix<ComplexD> &ref = matrices.at(term.back())[tp];
          tAr.stopTimer("Disk vector overhead");
          tAr.startTimer("Transpose caching"); // keep the legacy timing label
          lastTerm[tp].resize(ref.rows(), ref.cols());
          thread_for(j,ref.cols(),{
            for (Eigen::Index i = 0; i < ref.rows(); ++i)
              lastTerm[tp](i,j) = ref(i,j);
          });
          tAr.stopTimer("Transpose caching");
        }
        std::cout << Sec(tAr.getDTimer("Transpose caching")) << " "
                  << Bytes(double(cacheBytes), tAr.getDTimer("Transpose caching")) << std::endl;
      }

      for (unsigned int dt = 0; dt < par.global.nt; ++dt)
      {
        std::cout << " -- positions= " << t << ", dt= " << dt << std::endl;
        if (term.size() > 2) std::cout << std::setw(8) << "products";
        flops = 0.;
        bytes = 0.;
        tAr.startTimer("Disk vector overhead");
        prod = matrices.at(term[0])[dt]; // own the copy before another cache access
        tAr.stopTimer("Disk vector overhead");
        axes = matrixAxes(matrixPar.at(term[0]),dt,par.global.nt);
        if (term.size() == 1)
        {
          tAr.startTimer("Linear algebra");
          traceOneDWC(result.correlator[dt],prod,axes,par.global.nt);
          tAr.stopTimer("Linear algebra");
          continue;
        }
        fusec = tAr.getDTimer("A*B algebra");
        busec = tAr.getDTimer("A*B total");
        tAr.startTimer("Linear algebra");
        for (std::size_t j = 1; j+1 < term.size(); ++j)
        {
          const unsigned int tj = TIME_MOD(std::int64_t(dt) + t[j-1]);
          const auto rightAxes = matrixAxes(matrixPar.at(term[j]),tj,par.global.nt);
          tAr.startTimer("Disk vector overhead");
          const A2AMatrix<ComplexD> &ref = matrices.at(term[j])[tj];
          tAr.stopTimer("Disk vector overhead");

          tAr.startTimer("A*B total");
          tAr.startTimer("A*B algebra");
          flops += multiplyDWC(tmp,prod,ref,axes,rightAxes,par.global.nt,buf);
          tAr.stopTimer("A*B algebra");
          prod = tmp;
          axes.col = rightAxes.col; // row keeps the FIRST term's source time
          tAr.stopTimer("A*B total");
          bytes += 3.*double(tmp.rows())*double(tmp.cols())*sizeof(ComplexD);
        }
        if (term.size() > 2)
        {
          std::cout << Sec(tAr.getDTimer("A*B total") - busec) << " "
                    << Flops(flops,tAr.getDTimer("A*B algebra") - fusec) << " "
                    << Bytes(bytes,tAr.getDTimer("A*B total") - busec) << std::endl;
        }
        std::cout << std::setw(8) << "traces";
        fusec = tAr.getDTimer("tr(A*B)");
        busec = fusec;
        const unsigned int tLast = TIME_MOD(std::int64_t(dt) + t.back());
        const auto lastAxes = matrixAxes(matrixPar.at(term.back()),tLast,par.global.nt);
        tAr.startTimer("tr(A*B)");
        flops = traceDWC(result.correlator[dt],prod,lastTerm[tLast],
                         axes,lastAxes,par.global.nt,buf,selectedLast);
        tAr.stopTimer("tr(A*B)");
        bytes = 2.*double(prod.rows())*double(prod.cols())*sizeof(ComplexD);
        tAr.stopTimer("Linear algebra");
        std::cout << Sec(tAr.getDTimer("tr(A*B)") - busec) << " "
                  << Flops(flops,tAr.getDTimer("tr(A*B)") - fusec) << " "
                  << Bytes(bytes,tAr.getDTimer("tr(A*B)") - busec) << std::endl;
        std::cout << t << ":  t0 = " << dt << ", tLast = " << tLast << std::endl;
      }// dt
      saveCorrelator(result,par.global.output,fileStem,traj);
      tAr.stopTimer("Total");
      printTimeProfile(tAr.getTimings(),tAr.getTimer("Total"));
    }
    // Avoid unsigned wrap in the trajectory increment.
    if (par.global.trajCounter.end - traj <= par.global.trajCounter.step) break;
    traj += par.global.trajCounter.step;
  }
  }
  catch (const H5NS::Exception &e)
  {
    std::cerr << "BubbleContractorDWC: HDF5: " << e.getDetailMsg() << std::endl;
    return EXIT_FAILURE;
  }
  catch (const std::exception &e)
  {
    std::cerr << "BubbleContractorDWC: " << e.what() << std::endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
