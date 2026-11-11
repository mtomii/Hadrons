/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: utilities/K2pipiContractorDWC.cpp

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
#include <regex>
#include <set>
#include <string>
#include <tuple>
#include <vector>

#ifndef HAVE_HDF5
#error "K2pipiContractorDWC requires HDF5 support"
#endif

// Based on K2pipiContractor.cpp. Preserve its fixed product shared by all last
// matrices and moving times. DWC helpers follow Flexible/BubbleContractorDWC.
// Each fixed term has one integer base time. Translation averages are separate
// for each last matrix. No lattice Grid, extra conjugation or hit normalization.

using namespace Grid;
using namespace Hadrons;

namespace K2pipiContractorDWC
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
                                        std::vector<unsigned int>, times,
                                        std::vector<std::string>, lastMats,
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
                                        std::vector<unsigned int>, times,
                                        std::vector<ComplexD>, correlator);
    };
}

struct K2pipiContractorDWCPar
{
    K2pipiContractorDWC::GlobalPar                  global;
    std::vector<K2pipiContractorDWC::A2AMatrixPar>   a2aMatrix;
    std::vector<std::string>                       a2aMatrixLast;
    std::vector<K2pipiContractorDWC::ProductPar>     product;
};

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
  if (!ok)
  {
    HADRONS_ERROR(Argument, "K2pipiContractorDWC: " + message);
  }
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

MatrixAxes matrixAxes(const K2pipiContractorDWC::A2AMatrixPar &p,
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

// Only translations retain list/range syntax; fixed times are integers.
std::set<unsigned int> parseTimeRange(const std::string &str, unsigned int nt)
{
    const std::regex rex("([0-9]+)|(([0-9]+)\\.\\.([0-9]+))");
    std::set<unsigned int> values;
    for (const auto &token: strToVec<std::string>(str))
    {
        std::smatch match;
        require(std::regex_match(token, match, rex),
                "invalid translation expression: " + token);
        if (match[1].matched)
        {
            const auto value = std::stoull(match[1].str());
            require(value < nt, "translation out of range: " + token);
            values.insert(static_cast<unsigned int>(value));
        }
        else
        {
            const auto first = std::stoull(match[3].str());
            const auto last  = std::stoull(match[4].str());
            require(first <= last && last < nt,
                    "invalid translation range: " + token);
            for (auto value = first; value <= last; ++value)
            {
                values.insert(static_cast<unsigned int>(value));
            }
        }
    }
    return values;
}

std::string productStem(const K2pipiContractorDWC::ProductPar &p)
{
    const auto terms = strToVec<std::string>(p.terms);
    require(!terms.empty() && terms.size() == p.times.size(),
            "each product requires one time per fixed term");
    std::string stem;
    for (std::size_t j = 0; j < terms.size(); ++j)
    {
        if (j != 0)
        {
            stem += "_";
        }
        stem += terms[j] + "_" + std::to_string(p.times[j]);
    }
    return stem;
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
  fields({"global","a2aMatrix","a2aMatrixLast","product"});
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
    fields({"terms","times","lastMats","translations","min_t","max_t","translationAverage"});
    integers({"min_t","max_t"},false);
    std::string averageText; read(reader,"translationAverage",averageText);
    const auto average = strToVec<std::string>(averageText);
    require(average.size() == 1 && (average[0] == "true" || average[0] == "false"),
            "translationAverage must be true or false");
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

void validateParameters(const K2pipiContractorDWCPar &par)
{
    require(par.global.nt > 0 &&
            par.global.nt <= unsigned(std::numeric_limits<int>::max()), "invalid nt");
    require(par.global.trajCounter.start < par.global.trajCounter.end &&
            par.global.trajCounter.step > 0,
            "trajCounter requires start < end and step > 0");
    require(!par.global.diskVectorDir.empty() && !par.global.output.empty(),
            "empty diskVectorDir/output");
    require(!par.a2aMatrix.empty() && !par.product.empty(), "empty matrix/product list");

    std::set<std::string> names, lastNames, outputs;
    for (const auto &m: par.a2aMatrix)
    {
        require(!m.name.empty() && !m.file.empty() && !m.dataset.empty(),
                "empty matrix name/file/dataset");
        require(names.insert(m.name).second, "duplicate matrix name: " + m.name);
        require(m.cacheSize >= 1 && m.cacheSize <= par.global.nt,
                "cacheSize must lie in 1..nt");
        require(m.initDenseRow >= -1 && m.initDenseCol >= -1,
                "initDense must be -1 or nonnegative");
    }
    require(!par.a2aMatrixLast.empty(), "empty a2aMatrixLast list");
    for (const auto &name: par.a2aMatrixLast)
    {
        require(names.count(name), "unknown a2aMatrixLast name: " + name);
        require(lastNames.insert(name).second, "duplicate a2aMatrixLast name: " + name);
    }
    for (const auto &p: par.product)
    {
        const auto terms = strToVec<std::string>(p.terms);
        require(!terms.empty() && terms.size() == p.times.size(),
                "number of fixed terms (" + std::to_string(terms.size()) +
                ") must equal number of times (" + std::to_string(p.times.size()) + ")");
        for (const auto &name: terms)
        {
            require(names.count(name), "unknown fixed term: " + name);
        }
        for (auto time: p.times)
        {
            require(time < par.global.nt, "times entry must lie in 0..nt-1");
        }
        require(!p.lastMats.empty(), "empty product lastMats list");
        std::set<std::string> productLastNames;
        for (const auto &name: p.lastMats)
        {
            require(lastNames.count(name), "lastMats name is not in a2aMatrixLast: " + name);
            require(productLastNames.insert(name).second, "duplicate lastMats name: " + name);
        }
        require(p.min_t <= p.max_t && p.max_t < par.global.nt, "invalid min_t/max_t");
        require(!parseTimeRange(p.translations, par.global.nt).empty(), "empty translations");
        const auto stem = productStem(p);
        require(outputs.insert(stem).second, "duplicate product output filename: " + stem);
        checkedMul(checkedMul(p.lastMats.size(), std::size_t(par.global.nt)), sizeof(ComplexD));
    }
}

// Check every fixed edge and both closing edges for every last matrix before
// constructing large last caches or opening output files. Dimensions are time
// independent; the actual source times are used later by multiplyDWC/traceDWC.
void checkConnections(const K2pipiContractorDWCPar &par,
                      const std::map<std::string, MatrixDimensions> &dimensions,
                      const std::map<std::string, K2pipiContractorDWC::A2AMatrixPar> &matrixPar)
{
    for (const auto &p: par.product)
    {
        const auto terms = strToVec<std::string>(p.terms);
        for (std::size_t j = 1; j < terms.size(); ++j)
        {
            const auto &left = dimensions.at(terms[j-1]);
            const auto &right = dimensions.at(terms[j]);
            const auto a = matrixAxes(matrixPar.at(terms[j-1]), 0, par.global.nt);
            const auto b = matrixAxes(matrixPar.at(terms[j]), 0, par.global.nt);
            connection(left.cols, a.col, right.rows, b.row, par.global.nt);
        }
        const auto &first = dimensions.at(terms.front());
        const auto &last = dimensions.at(terms.back());
        const auto a = matrixAxes(matrixPar.at(terms.front()), 0, par.global.nt);
        const auto b = matrixAxes(matrixPar.at(terms.back()), 0, par.global.nt);
        for (const auto &name: p.lastMats)
        {
            const auto &d = dimensions.at(name);
            const auto c = matrixAxes(matrixPar.at(name), 0, par.global.nt);
            connection(last.cols, b.col, d.rows, c.row, par.global.nt);
            connection(d.cols, c.col, first.rows, a.row, par.global.nt);
        }
    }
}

void checkCorrelator(const K2pipiContractorDWC::CorrelatorResult &result)
{
    for (const auto z: result.correlator)
    {
        require(std::isfinite(z.real()) && std::isfinite(z.imag()), "non-finite correlator");
    }
}

int main(int argc, char* argv[])
{
    if (argc != 2)
    {
        std::cerr << "usage: " << argv[0] << " <parameter file>" << std::endl;
        return EXIT_FAILURE;
    }
    const std::string parFilename = argv[1];

    try
    {
        checkSingleProcess();
        H5NS::Exception::dontPrint();
        checkXmlFields(parFilename);

        K2pipiContractorDWCPar par;
        XmlReader reader(parFilename);
        read(reader, "global", par.global);
        read(reader, "a2aMatrix", par.a2aMatrix);
        read(reader, "a2aMatrixLast", par.a2aMatrixLast);
        read(reader, "product", par.product);
        validateParameters(par);

        Eigen::initParallel();
        std::cout << "K2pipiContractorDWC: one process; Eigen reports "
                  << Eigen::nbThreads() << " thread(s); no lattice Grid" << std::endl;

        std::map<std::string, K2pipiContractorDWC::A2AMatrixPar> matrixPar;
        for (const auto &m: par.a2aMatrix)
        {
            matrixPar.emplace(m.name, m);
        }
        std::set<std::string> fixedNames;
        for (const auto &p: par.product)
        {
            const auto terms = strToVec<std::string>(p.terms);
            fixedNames.insert(terms.begin(), terms.end());
        }

        // Trajectory loop. Recreate all load destinations on every trajectory:
        // last-only row-major entries are deliberately erased after caching.
        for (unsigned int traj = par.global.trajCounter.start;
             traj < par.global.trajCounter.end; )
        {
            std::cout << ":::::::: Trajectory " << traj << std::endl;
            std::map<std::string, EigenDiskVector<ComplexD>> a2aMat;
            std::map<std::string, std::vector<LastMatrix>> a2aMatLast;
            std::map<std::string, MatrixDimensions> dimensions;

            // Construct DiskVectors directly inside the map, without temporary
            // DiskVectors whose destructor also owns scratch-directory cleanup.
            for (const auto &p: par.a2aMatrix)
            {
                const auto dirName = par.global.diskVectorDir + "/" + p.name;
                makeFileDir(dirName);
                a2aMat.emplace(std::piecewise_construct, std::forward_as_tuple(p.name),
                               std::forward_as_tuple(dirName, par.global.nt, p.cacheSize));
            }

            // Load each input once for this trajectory; all products share it.
            for (const auto &p: par.a2aMatrix)
            {
                std::string filename = p.file;
                tokenReplace(filename, "traj", traj);
                std::cout << "======== Loading '" << filename << "'" << std::endl;
                double tRead;
                std::size_t fileBytes;
                dimensions.emplace(p.name, loadMatrix(filename, p.dataset, par.global.nt,
                                                       a2aMat.at(p.name), tRead, fileBytes));
                std::cout << "Read " << fileBytes << " disk bytes in " << tRead/1.0e6
                          << " sec (includes conversion/cache assignment), "
                          << (tRead > 0 ? fileBytes/tRead*1.0e6/1024/1024 : 0.)
                          << " MB/s" << std::endl;
            }
            checkConnections(par, dimensions, matrixPar);

            // Cache all last matrices once per trajectory, shared by products.
            // Column-major changes storage only; indices and values are unchanged.
            const auto &matrices = a2aMat;
            for (const auto &name: par.a2aMatrixLast)
            {
                const auto &d = dimensions.at(name);
                const auto cacheBytes = checkedMul(checkedMul(checkedMul(par.global.nt,
                                          std::size_t(d.rows)), std::size_t(d.cols)), sizeof(ComplexD));
                std::cout << "* Caching column-major last term " << name
                          << " (" << cacheBytes << " bytes)" << std::endl;
                TimerArray cacheTimer;
                cacheTimer.startTimer("Total");
                auto inserted = a2aMatLast.emplace(name, std::vector<LastMatrix>(par.global.nt));
                auto &lastTerm = inserted.first->second;
                for (unsigned int t = 0; t < par.global.nt; ++t)
                {
                    cacheTimer.startTimer("Disk vector overhead");
                    const Matrix &ref = matrices.at(name)[t];
                    cacheTimer.stopTimer("Disk vector overhead");
                    cacheTimer.startTimer("Column-major caching");
                    lastTerm[t].resize(ref.rows(), ref.cols());
                    thread_for(j, ref.cols(),
                    {
                        for (Eigen::Index i = 0; i < ref.rows(); ++i)
                        {
                            lastTerm[t](i,j) = ref(i,j);
                        }
                    });
                    cacheTimer.stopTimer("Column-major caching");
                }
                cacheTimer.stopTimer("Total");
                std::cout << Sec(cacheTimer.getDTimer("Column-major caching")) << " "
                          << Bytes(double(cacheBytes), cacheTimer.getDTimer("Column-major caching"))
                          << std::endl;
                // A fixed term may also be a last matrix: keep its row-major data.
                if (fixedNames.count(name) == 0)
                {
                    a2aMat.erase(name);
                }
            }

            // Contract. The fixed product is computed ONCE per translation,
            // then reused for all last matrices and all moving times.
            for (const auto &p: par.product)
            {
                const auto terms = strToVec<std::string>(p.terms);
                const auto &times = p.times;
                const auto translations = parseTimeRange(p.translations, par.global.nt);
                const auto stem = productStem(p);
                const auto filename = par.global.output + "/" + ModuleBase::resultFilename(stem, traj);
                Matrix prod, tmp, selection;
                LastMatrix selectedLast;
                MatrixAxes axes;
                TimerArray timers;
                timers.startTimer("Total");

                std::vector<K2pipiContractorDWC::CorrelatorResult> results(p.lastMats.size());
                for (auto &result: results)
                {
                    result.times = times;
                    result.correlator.assign(par.global.nt, ComplexD(0.,0.));
                }

                const auto nProducts = checkedMul(translations.size(), terms.size()-1);
                const auto nTraces = checkedMul(checkedMul(translations.size(), p.lastMats.size()),
                                                std::size_t(p.max_t)-p.min_t+1);
                std::cout << "======== Fixed product " << stem << "; " << nProducts
                          << " A*B, " << nTraces << " tr(A*B); "
                          << p.lastMats.size() << " independent last-matrix results" << std::endl;
                std::cout << "Saving correlators to '" << filename << "'" << std::endl;
                makeFileDir(filename);
                ResultWriter writer(filename);

                std::size_t translationIndex = 0;
                for (const auto dt: translations)
                {
                    std::cout << "* Step " << ++translationIndex << "/" << translations.size()
                              << " -- positions= " << times << ", dt= " << dt << std::endl;
                    // Without averaging, each translation starts from fresh vectors.
                    // With averaging, each last matrix has its own accumulator.
                    if (!p.translationAverage)
                    {
                        for (auto &result: results)
                        {
                            std::fill(result.correlator.begin(), result.correlator.end(), ComplexD(0.,0.));
                        }
                    }

                    const auto productBegin = timers.getDTimer("A*B algebra");
                    double productFlops = 0.;
                    timers.startTimer("Disk vector overhead");
                    const auto firstTime = timeModulo(std::int64_t(times.front())+dt, par.global.nt);
                    // Own this copy before accessing another slice of a DiskVector.
                    prod = matrices.at(terms.front())[firstTime];
                    axes = matrixAxes(matrixPar.at(terms.front()), firstTime, par.global.nt);
                    timers.stopTimer("Disk vector overhead");
                    for (std::size_t j = 1; j < terms.size(); ++j)
                    {
                        const auto refTime = timeModulo(std::int64_t(times[j])+dt, par.global.nt);
                        timers.startTimer("Disk vector overhead");
                        const Matrix &ref = matrices.at(terms[j])[refTime];
                        const auto refAxes = matrixAxes(matrixPar.at(terms[j]), refTime, par.global.nt);
                        timers.stopTimer("Disk vector overhead");
                        timers.startTimer("A*B algebra");
                        productFlops += multiplyDWC(tmp, prod, ref, axes, refAxes, par.global.nt, selection);
                        timers.stopTimer("A*B algebra");
                        prod = tmp;
                        axes.col = refAxes.col;
                    }
                    if (terms.size() > 1)
                    {
                        const auto elapsed = timers.getDTimer("A*B algebra") - productBegin;
                        std::cout << "  fixed products: " << Sec(elapsed) << " "
                                  << Flops(productFlops, elapsed) << std::endl;
                    }

                    const auto traceBegin = timers.getDTimer("tr(A*B)");
                    double traceFlops = 0.;
                    for (std::size_t ell = 0; ell < p.lastMats.size(); ++ell)
                    {
                        const auto &name = p.lastMats[ell];
                        const auto &lastTerm = a2aMatLast.at(name);
                        auto &result = results[ell];
                        for (unsigned int tau = p.min_t; tau <= p.max_t; ++tau)
                        {
                            const auto lastTime = timeModulo(std::int64_t(tau)+dt, par.global.nt);
                            const auto lastAxes = matrixAxes(matrixPar.at(name), lastTime, par.global.nt);
                            timers.startTimer("tr(A*B)");
                            traceFlops += traceDWC(result.correlator[tau], prod, lastTerm[lastTime],
                                                   axes, lastAxes, par.global.nt, selection, selectedLast);
                            timers.stopTimer("tr(A*B)");
                        }
                        if (!p.translationAverage)
                        {
                            const auto dataSet = stem + "_" + name + "_dt_" + std::to_string(dt);
                            timers.startTimer("Write");
                            checkCorrelator(result);
                            write(writer, dataSet, result);
                            timers.stopTimer("Write");
                        }
                    }
                    const auto elapsed = timers.getDTimer("tr(A*B)") - traceBegin;
                    std::cout << "  traces: " << Sec(elapsed) << " "
                              << Flops(traceFlops, elapsed) << std::endl;
                }

                if (p.translationAverage)
                {
                    for (std::size_t ell = 0; ell < p.lastMats.size(); ++ell)
                    {
                        auto &result = results[ell];
                        for (auto &z: result.correlator)
                        {
                            z /= static_cast<double>(translations.size());
                        }
                        timers.startTimer("Write");
                        checkCorrelator(result);
                        write(writer, stem + "_" + p.lastMats[ell], result);
                        timers.stopTimer("Write");
                    }
                }
                timers.stopTimer("Total");
                printTimeProfile(timers.getTimings(), timers.getTimer("Total"));
            }

            // Avoid overflow on the final unsigned trajectory increment.
            if (par.global.trajCounter.step >= par.global.trajCounter.end-traj)
            {
                break;
            }
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
