/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/Utilities/FlexibleContractorDWC.cc

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
#include <cassert>
#include <cerrno>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <limits>
#include <map>
#include <regex>
#include <set>
#include <sstream>
#include <stdexcept>
#include <tuple>
#include <type_traits>
#include <vector>

#ifndef HAVE_HDF5
#error "FlexibleContractorDWC requires HDF5 support"
#endif

// DWC changes: retain the Contractor loops and DiskVector storage, add
// per-axis dense-W pairing, and read HDF5 using an explicit native memory type.
// No Grid_init(), lattice Grid, conjugation, or extra hit normalization.

using namespace Grid;
using namespace Hadrons;

#define TIME_MOD(t) timeModulo((t), par.global.nt)

namespace FlexibleContractorDWC
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

  class TermsPar: Serializable
  {
  public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(TermsPar,
				    std::string, term,
				    bool, tdps,
				    std::string, time);
    TermsPar(): tdps(false) {}
  };

  class ProductPar: Serializable
  {
  public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(ProductPar,
				    std::vector<FlexibleContractorDWC::TermsPar>, terms,
				    std::string, translations,
				    unsigned int, min_t,
				    unsigned int, max_t,
				    bool, translationAverage);
    ProductPar(): min_t(0), max_t(0), translationAverage(false) {}
  };

  class CorrelatorResult: Serializable
  {
  public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(CorrelatorResult,
				    std::vector<FlexibleContractorDWC::A2AMatrixPar>,  a2aMatrix,
				    ProductPar, contraction,
				    std::vector<unsigned int>, times,
				    std::vector<ComplexD>, correlator);
  };
}

struct FlexibleContractorDWCPar
{
    FlexibleContractorDWC::GlobalPar                  global;
    std::vector<FlexibleContractorDWC::A2AMatrixPar>  a2aMatrix;
    std::vector<FlexibleContractorDWC::ProductPar>    product;
};

void makeTimeSeq(std::vector<std::vector<unsigned int>> &timeSeq, 
                 const std::vector<std::set<unsigned int>> &times,
                 std::vector<unsigned int> &current,
                 const unsigned int depth)
{
    if (depth > 0)
    {
        for (auto t: times[times.size() - depth])
        {
            current[times.size() - depth] = t;
            makeTimeSeq(timeSeq, times, current, depth - 1);
        }
    }
    else
    {
        timeSeq.push_back(current);
    }
}

void makeTimeSeq(std::vector<std::vector<unsigned int>> &timeSeq, 
                 const std::vector<std::set<unsigned int>> &times)
{
    std::vector<unsigned int> current(times.size());

    makeTimeSeq(timeSeq, times, current, times.size());
}

std::set<unsigned int> parseTimeRange(const std::string str, const unsigned int nt)
{
  std::regex               rex("([0-9]+)|(([0-9]+)\\.\\.([0-9]+))");
  std::smatch              sm;
  std::vector<std::string> rstr = strToVec<std::string>(str);
  std::set<unsigned int>   tSet;

  for (auto &s: rstr)
  {
    if (!std::regex_match(s, sm, rex))
    {
      HADRONS_ERROR(Argument, "invalid time expression '" + s + "'");
    }
    if (sm[1].matched)
    {
      unsigned int t;
            
      const auto value = std::stoull(sm[1].str());
      if (value >= nt) { HADRONS_ERROR(Range, "time out of range: " + s); }
      t = static_cast<unsigned int>(value);
      if (t >= nt)
      {
	HADRONS_ERROR(Range, "time out of range (from expression '" + str + "')");
      }
      tSet.insert(t);
    }
    else if (sm[2].matched)
    {
      unsigned int ta, tb;

      const auto a = std::stoull(sm[3].str()), b = std::stoull(sm[4].str());
      if (a > b || b >= nt) { HADRONS_ERROR(Range, "invalid time range: " + s); }
      ta = static_cast<unsigned int>(a);
      tb = static_cast<unsigned int>(b);
      if ((ta >= nt) or (tb >= nt))
      {
	HADRONS_ERROR(Range, "time out of range (from expression '" + str + "')");
      }
      for (unsigned int ti = ta; ti <= tb; ++ti)
      {
	tSet.insert(ti);
      }
    }
  }

  return tSet;
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
  if (!ok) { HADRONS_ERROR(Argument, "FlexibleContractorDWC: " + message); }
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

MatrixAxes matrixAxes(const FlexibleContractorDWC::A2AMatrixPar &p,
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

std::string productStem(const FlexibleContractorDWC::ProductPar &p)
{
  std::string stem;
  for (std::size_t j = 0; j < p.terms.size(); ++j)
  {
    const auto &term = p.terms[j];
    if (j) stem += "_";
    stem += term.term + "_" + (term.tdps ? "x" : "y") + term.time;
  }
  return stem;
}

// Do not let missing XML silently turn a dense axis into an ordinary one.
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
  const auto unsignedFields = [&](const std::vector<std::string> &names)
  {
    for (const auto &n: names)
    {
      std::string text; read(reader,n,text);
      const auto token = strToVec<std::string>(text);
      require(token.size() == 1 && std::regex_match(token[0],std::regex("[0-9]+")),
              "expected a nonnegative integer for " + n);
      require(std::stoull(token[0]) <= std::numeric_limits<unsigned int>::max(),
              "unsigned integer overflow for " + n);
    }
  };
  const auto signedFields = [&](const std::vector<std::string> &names)
  {
    for (const auto &n: names)
    {
      std::string text; read(reader,n,text);
      const auto token = strToVec<std::string>(text);
      require(token.size() == 1 && std::regex_match(token[0],std::regex("[-+]?[0-9]+")),
              "expected a signed integer for " + n);
      const auto value = std::stoll(token[0]);
      require(value >= std::numeric_limits<int>::min() && value <= std::numeric_limits<int>::max(),
              "signed integer overflow for " + n);
    }
  };
  fields({"global","a2aMatrix","product"});
  reader.push("global"); fields({"trajCounter","nt","diskVectorDir","output"});
  unsignedFields({"nt"});
  reader.push("trajCounter"); fields({"start","end","step"});
  unsignedFields({"start","end","step"}); reader.pop(); reader.pop();
  reader.push("a2aMatrix"); require(reader.push("elem"), "empty a2aMatrix list");
  do
  {
    fields({"file","dataset","cacheSize","name","initDenseRow","initDenseCol","timeOffsetRow","timeOffsetCol"});
    unsignedFields({"cacheSize"});
    signedFields({"initDenseRow","initDenseCol","timeOffsetRow","timeOffsetCol"});
  }
  while (reader.nextElement("elem"));
  reader.pop(); reader.pop();
  reader.push("product"); require(reader.push("elem"), "empty product list");
  do
  {
    fields({"terms","translations","min_t","max_t","translationAverage"});
    unsignedFields({"min_t","max_t"});
    reader.push("terms"); require(reader.push("elem"), "empty terms list");
    do { fields({"term","tdps","time"}); } while (reader.nextElement("elem"));
    reader.pop(); reader.pop();
  } while (reader.nextElement("elem"));
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

void validateParameters(const FlexibleContractorDWCPar &par)
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
    require(p.terms.size() >= 2 && !p.terms.front().tdps,
            "each product needs at least two terms and first tdps=false");
    require(std::any_of(p.terms.begin(),p.terms.end(),[](const auto &t){return t.tdps;}),
            "each product needs at least one tdps=true term");
    require(p.min_t <= p.max_t && p.max_t < par.global.nt, "invalid min_t/max_t");
    require(!parseTimeRange(p.translations,par.global.nt).empty(), "empty translations");
    for (const auto &term: p.terms)
    {
      require(names.count(term.term), "unknown term: " + term.term);
      require(parseTimeRange(term.time,par.global.nt).size() == 1,
              "one time value per term is required: " + term.term);
    }
    require(outputs.insert(productStem(p)).second, "duplicate product output filename: " + productStem(p));
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
  FlexibleContractorDWCPar par;
  //    unsigned int  nMat,nCont;
  XmlReader     reader(parFilename);

  read(reader, "global",    par.global);
  read(reader, "a2aMatrix", par.a2aMatrix);
  read(reader, "product",   par.product);
  //    nMat  = par.a2aMatrix.size();
  //    nCont = par.product.size();

  validateParameters(par);
  Eigen::initParallel();
  std::cout << "FlexibleContractorDWC: one process; Eigen reports "
            << Eigen::nbThreads() << " thread(s); no lattice Grid" << std::endl;
  std::map<std::string, FlexibleContractorDWC::A2AMatrixPar> matrixPar;
  for (const auto &m: par.a2aMatrix) matrixPar.emplace(m.name,m);

  // create diskvectors
  std::map<std::string, EigenDiskVector<ComplexD>> a2aMat;
  //    unsigned int                                     cacheSize;

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
      //	    double  size;

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
    //EigenDiskVector<ComplexD>::Matrix buf;

    for (auto &p: par.product)
    {
      std::vector<FlexibleContractorDWC::TermsPar>                  terms = p.terms;
      std::vector<std::set<unsigned int>>    times;
      std::vector<std::vector<unsigned int>> timeSeq;
      std::set<unsigned int>                 translations;
      A2AMatrix<ComplexD>                    prod0, prod, buf, tmp;
      LastMatrix                             selectedLast;
      MatrixAxes                             axes0, axes;
      const auto &matrices = a2aMat; // const DiskVector operator[], never a write proxy.
      TimerArray                             tAr;
      double                                 fusec, busec, flops, bytes;
      double                                 fusecp, busecp, flopsp, bytesp;
      std::vector<A2AMatrixTr<ComplexD>>     lastTerm(par.global.nt);
      //	    double  tusec;
      FlexibleContractorDWC::CorrelatorResult   result;

      if ( terms[0].tdps ) {
	std::cout << "You have to specify tdps=False for the 1st meson field of each product" << std::endl;
	assert(!terms[0].tdps);
      }

      tAr.startTimer("Total");
      std::cout << "======== Contraction tr(";
      for (unsigned int g = 0; g < terms.size(); ++g)
      {
	std::cout << terms[g].term << ((g == terms.size() - 1) ? ')' : '*');
	times.push_back(parseTimeRange(terms[g].time, par.global.nt));
      }
      std::cout << std::endl;
      for (auto &m: par.a2aMatrix)
      {
        result.a2aMatrix.push_back(m);
        tokenReplace(result.a2aMatrix.back().file, "traj", traj);
      }
      result.contraction = p;
      result.correlator.resize(par.global.nt, 0.);

      translations = parseTimeRange(p.translations, par.global.nt);
      makeTimeSeq(timeSeq, times);
      require(timeSeq.size() == 1, "exactly one timeSeq is supported");
      // Check every edge, including last -> first, before output creation.
      for (std::size_t j = 0; j < terms.size(); ++j)
      {
        const auto k = (j+1)%terms.size();
        const auto &a = dimensions.at(terms[j].term), &b = dimensions.at(terms[k].term);
        const auto aa = matrixAxes(matrixPar.at(terms[j].term),0,par.global.nt);
        const auto bb = matrixAxes(matrixPar.at(terms[k].term),0,par.global.nt);
        connection(a.cols,aa.col,b.rows,bb.row,par.global.nt);
      }

      std::cout << "* Caching transposed last term " << std::endl;
      for (unsigned int t = 0; t < par.global.nt; ++t)
      {
	tAr.startTimer("Disk vector overhead");
	const A2AMatrix<ComplexD> &ref = matrices.at(terms.back().term)[t];
	tAr.stopTimer("Disk vector overhead");

	tAr.startTimer("Transpose caching");
	lastTerm[t].resize(ref.rows(), ref.cols());
	thread_for( j,ref.cols(),{
	  for (unsigned int i = 0; i < ref.rows(); ++i)
	  {
	    lastTerm[t](i, j) = ref(i, j);
	  }
	});
	tAr.stopTimer("Transpose caching");
      }
      bytes = double(par.global.nt)*double(lastTerm[0].rows())*double(lastTerm[0].cols())*sizeof(ComplexD);
      std::cout << Sec(tAr.getDTimer("Transpose caching")) << " " 
		<< Bytes(bytes, tAr.getDTimer("Transpose caching")) << std::endl;

      std::string              fileStem = "", filename, dataSet;
      std::vector<FlexibleContractorDWC::TermsPar> termsv = result.contraction.terms;
      for (unsigned int i = 0; i < termsv.size() ; i++)
      {
	fileStem += termsv[i].term + "_";
	if ( termsv[i].tdps ) {
	  fileStem += "x";
	} else {
	  fileStem += "y";
	}
	fileStem += termsv[i].time;
	if ( i < termsv.size() - 1 ) fileStem += "_";
      }
      const std::string dir = par.global.output;
      filename = dir + "/" + ModuleBase::resultFilename(fileStem, traj);
      //filename = dir + "/" + RESULT_FILE_NAME(fileStem, traj);
      std::cout << "Saving correlator to '" << filename << "'" << std::endl;
      makeFileDir(filename);
      ResultWriter writer(filename);

      for (unsigned int i = 0; i < timeSeq.size(); ++i)
      {
	unsigned int dti = 0;
	auto         &t = timeSeq[i];

	result.times = t;
	for (unsigned int tLast = 0; tLast < par.global.nt; ++tLast)
	{
	  result.correlator[tLast] = 0.;
	}
	for (auto &dt: translations)
	{
	  std::cout << "* Step " << i*translations.size() + dti + 1
		    << "/" << timeSeq.size()*translations.size()
		    << " -- positions= " << t << ", dt= " << dt << std::endl;
	  flops  = 0.;
	  bytes  = 0.;
	  fusec  = tAr.getDTimer("A*B algebra prec");
	  busec  = tAr.getDTimer("A*B total prec");
	  tAr.startTimer("Linear algebra");
	  unsigned int nlast = terms.size()-1;
	  tAr.startTimer("Disk vector overhead");
          const unsigned int firstTime = TIME_MOD(std::int64_t(t[0]) + dt);
          prod = matrices.at(terms[0].term)[firstTime];
          axes = matrixAxes(matrixPar.at(terms[0].term),firstTime,par.global.nt);
	  tAr.stopTimer("Disk vector overhead");
	  for (unsigned int j = 1; j < terms.size(); ++j)
	  {
	    if ( terms[j].tdps ) break;
	    --nlast;
	    tAr.startTimer("Disk vector overhead");
            const unsigned int refTime = TIME_MOD(std::int64_t(t[j]) + dt);
            const A2AMatrix<ComplexD> &ref = matrices.at(terms[j].term)[refTime];
            const auto refAxes = matrixAxes(matrixPar.at(terms[j].term),refTime,par.global.nt);
	    tAr.stopTimer("Disk vector overhead");

	    tAr.startTimer("A*B total prec");
	    tAr.startTimer("A*B algebra prec");
            flops += multiplyDWC(tmp, prod, ref, axes, refAxes, par.global.nt, buf);
            tAr.stopTimer("A*B algebra prec");
            prod   = tmp;
            axes.col = refAxes.col;
	    tAr.stopTimer("A*B total prec");
	    bytes += 3.*tmp.rows()*tmp.cols()*sizeof(ComplexD);
	  }
	  if (terms.size() - nlast > 1 )
	  {
	    std::cout << std::setw(8) << "t-independent products"
		      << Sec(tAr.getDTimer("A*B total prec") - busec) << " "
		      << Flops(flops, tAr.getDTimer("A*B algebra prec") - fusec) << " "
		      << Bytes(bytes, tAr.getDTimer("A*B total prec") - busec) << std::endl;
	  }
	  std::cout << std::setw(8) << "traces";
	  flops  = 0.;
	  bytes  = 0.;
	  fusec  = tAr.getDTimer("tr(A*B)");
	  busec  = tAr.getDTimer("tr(A*B)");
	  flopsp  = 0.;
	  bytesp  = 0.;
	  fusecp  = tAr.getDTimer("A*B algebra rest");
	  busecp  = tAr.getDTimer("A*B total rest");
	  if ( nlast > 1 ) {
	    tAr.startTimer("Copying A2AMatrix");
            prod0 = prod;
            axes0 = axes;
            tAr.stopTimer("Copying A2AMatrix");
	    for (std::int64_t tLast = std::int64_t(p.min_t)+dt; tLast <= std::int64_t(p.max_t)+dt; ++tLast)
	    {
	      tAr.startTimer("Copying A2AMatrix");
              prod = prod0;
              axes = axes0;
              tAr.stopTimer("Copying A2AMatrix");
	      for (unsigned int j = terms.size() - nlast ; j < terms.size() - 1 ; ++j) {
		tAr.startTimer("Disk vector overhead");
		unsigned int tidx = TIME_MOD(std::int64_t(t[j])+dt);
		if ( terms[j].tdps ) tidx = TIME_MOD(t[j] + tLast);
		const A2AMatrix<ComplexD> &ref = matrices.at(terms[j].term)[tidx];
                const auto refAxes = matrixAxes(matrixPar.at(terms[j].term),tidx,par.global.nt);
		tAr.stopTimer("Disk vector overhead");

		tAr.startTimer("A*B total rest");
		tAr.startTimer("A*B algebra rest");
                flopsp += multiplyDWC(tmp, prod, ref, axes, refAxes, par.global.nt, buf);
                tAr.stopTimer("A*B algebra rest");
                prod   = tmp;
                axes.col = refAxes.col;
		tAr.stopTimer("A*B total rest");
		bytesp += 3.*tmp.rows()*tmp.cols()*sizeof(ComplexD);
	      }
	      tAr.startTimer("tr(A*B)");
	      unsigned int tidx = TIME_MOD(std::int64_t(t[terms.size()-1])+dt);
	      if ( terms[terms.size()-1].tdps )
		tidx = TIME_MOD(t[terms.size()-1] + tLast);
              const auto lastAxes = matrixAxes(matrixPar.at(terms.back().term),tidx,par.global.nt);
              flops += traceDWC(result.correlator[TIME_MOD(tLast-dt)], prod, lastTerm[tidx],
                                axes, lastAxes, par.global.nt, buf, selectedLast);
              tAr.stopTimer("tr(A*B)");
	      bytes += 2.*prod.rows()*prod.cols()*sizeof(ComplexD);
	    }
	    std::cout << Sec(tAr.getDTimer("A*B total rest") - busecp) << " "
		      << Flops(flopsp, tAr.getDTimer("A*B algebra rest") - fusecp) << " "
		      << Bytes(bytesp, tAr.getDTimer("A*B total rest") - busecp) << std::endl;
	  } else {
	    for (std::int64_t tLast = std::int64_t(p.min_t)+dt; tLast <= std::int64_t(p.max_t)+dt; ++tLast)
	    {
	      tAr.startTimer("tr(A*B)");
	      unsigned int tidx = TIME_MOD(t[terms.size()-1]+tLast);
	      assert( terms[terms.size()-1].tdps );
              const auto lastAxes = matrixAxes(matrixPar.at(terms.back().term),tidx,par.global.nt);
              flops += traceDWC(result.correlator[TIME_MOD(tLast-dt)], prod, lastTerm[tidx],
                                axes, lastAxes, par.global.nt, buf, selectedLast);
              tAr.stopTimer("tr(A*B)");
	      bytes += 2.*prod.rows()*prod.cols()*sizeof(ComplexD);
	    }
	  }
	  tAr.stopTimer("Linear algebra");
	  std::cout << Sec(tAr.getDTimer("tr(A*B)") - busec) << " "
		    << Flops(flops, tAr.getDTimer("tr(A*B)") - fusec) << " " 
		    << Bytes(bytes, tAr.getDTimer("tr(A*B)") - busec) << std::endl;
	  if (!p.translationAverage)
	  {
	    dataSet = fileStem + "_dt_" + std::to_string(dt);
	    tAr.startTimer("Write");
            for (const auto z: result.correlator)
              require(std::isfinite(z.real()) && std::isfinite(z.imag()), "non-finite correlator");
            write(writer, dataSet, result);
	    tAr.stopTimer("Write");
	    //saveCorrelator(result, par.global.output, dt, traj);
	    for (unsigned int tLast = 0; tLast < par.global.nt; ++tLast)
	    {
	      result.correlator[tLast] = 0.;
	    }
	  }
	  dti++;
	} // dt
	if (p.translationAverage)
	{
	  for (unsigned int tLast = 0; tLast < par.global.nt; ++tLast)
	  {
	    result.correlator[tLast] /= translations.size();
	  }
	  tAr.startTimer("Write");
          for (const auto z: result.correlator)
            require(std::isfinite(z.real()) && std::isfinite(z.imag()), "non-finite correlator");
          write(writer, fileStem, result);
	  tAr.stopTimer("Write");
	  //saveCorrelator(result, par.global.output, 0, traj);
	}
      }
      tAr.stopTimer("Total");
      printTimeProfile(tAr.getTimings(), tAr.getTimer("Total"));
    }
    if (par.global.trajCounter.step >= par.global.trajCounter.end - traj) break;
    traj += par.global.trajCounter.step;
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
