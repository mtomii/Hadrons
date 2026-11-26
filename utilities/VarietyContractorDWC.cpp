/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: utilities/VarietyContractorDWC.cpp

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
#include <Hadrons/Module.hpp>
#include <algorithm>
#include <atomic>
#include <cerrno>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <exception>
#include <limits>
#include <memory>
#include <new>
#include <stdexcept>
#include <utility>
#include <map>
#include <regex>
#include <set>
#include <string>
#include <tuple>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifndef HAVE_HDF5
#error "VarietyContractorDWC requires HDF5 support"
#endif

// Derived from utilities/K2pipiContractorDWC.cpp at Hadrons commit
// 28d0561d32b3c5bc7368f0812aa367a284a8df98 (local-2025-edits-grid-develop).
// DWC axis mapping and the left-to-right matrix-product order are retained.
// Changes: no DiskVector; retain all fixed products; stream last matrices one
// complete time slice at a time; write flat, cropped, two-dimensional HDF5 data.
// Source translations are rows (sorted, unique); no average or metadata output.
// XML removes diskVectorDir/cacheSize and requires translationAverage=false.
// One process; Eigen threads; no lattice Grid, added conjugation, or hit factors.
// Performance revision of VarietyContractorDWC at Hadrons commit
// 6a0e260b3a629c2397aa18d7d1f72a35b1055365: borrow the first fixed input,
// reuse workspaces, parallelize storage-order selections for fixed products,
// and parallelize independent products when tracing. Each trace uses a serial,
// nonconjugating SIMD inner product; floating-point summation order may differ.
// XML, input IO, time matching, precision, normalization and output are retained.

using namespace Grid;
using namespace Hadrons;

namespace VarietyContractorDWC
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
                                        std::string, output);
        GlobalPar(): nt(0) {}
    };

    class A2AMatrixPar: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(A2AMatrixPar,
                                        std::string, file,
                                        std::string, dataset,
                                        std::string, name,
                                        int, initDenseRow,
                                        int, initDenseCol,
                                        int, timeOffsetRow,
                                        int, timeOffsetCol);
        A2AMatrixPar(): initDenseRow(-1), initDenseCol(-1),
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


}

struct VarietyContractorDWCPar
{
    VarietyContractorDWC::GlobalPar                  global;
    std::vector<VarietyContractorDWC::A2AMatrixPar>   a2aMatrix;
    std::vector<std::string>                       a2aMatrixLast;
    std::vector<VarietyContractorDWC::ProductPar>     product;
};

// Small local helpers; no new shared matrix class or changes to A2AMatrix.hpp.
void require(bool ok, const std::string &message)
{
  if (!ok)
  {
    HADRONS_ERROR(Argument, "VarietyContractorDWC: " + message);
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

MatrixAxes matrixAxes(const VarietyContractorDWC::A2AMatrixPar &p,
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

// Fixed products copy disjoint rows/columns in parallel. A trace worker copies
// serially, without creating a nested or serialized OpenMP region per trace.
template <typename Copy>
void copySlices(Eigen::Index n, bool parallel, const Copy &copy)
{
  if (parallel)
  {
    thread_for(i, n,
    {
      copy(Eigen::Index(i));
    });
  }
  else
  {
    for (Eigen::Index i = 0; i < n; ++i)
    {
      copy(i);
    }
  }
}

// Copy selected axes only, following output storage order. Validation and
// allocation precede any parallel loop; each worker writes disjoint elements.
// Column-major last terms are not mathematically transposed.
template <typename Out, typename In>
void selectMatrix(Out &out, const In &in, const IndexList &rows,
                  const IndexList &cols, bool parallel = false)
{
  const Eigen::Index nr = rows.empty() ? in.rows() : Eigen::Index(rows.size());
  const Eigen::Index nc = cols.empty() ? in.cols() : Eigen::Index(cols.size());
  checkMatrixSize(nr, nc);
  out.resize(nr, nc);
  if (rows.empty() && cols.empty())
  {
    out = in;
  }
  else if constexpr (Out::IsRowMajor)
  {
    if (cols.empty())
    {
      copySlices(nr, parallel, [&](Eigen::Index r)
      {
        out.row(r) = in.row(rows[std::size_t(r)]);
      });
    }
    else if (rows.empty())
    {
      copySlices(nr, parallel, [&](Eigen::Index r)
      {
        for (Eigen::Index c = 0; c < nc; ++c)
        {
          out(r,c) = in(r,cols[std::size_t(c)]);
        }
      });
    }
    else
    {
      copySlices(nr, parallel, [&](Eigen::Index r)
      {
        for (Eigen::Index c = 0; c < nc; ++c)
        {
          out(r,c) = in(rows[std::size_t(r)],cols[std::size_t(c)]);
        }
      });
    }
  }
  else
  {
    if (rows.empty())
    {
      copySlices(nc, parallel, [&](Eigen::Index c)
      {
        out.col(c) = in.col(cols[std::size_t(c)]);
      });
    }
    else if (cols.empty())
    {
      copySlices(nc, parallel, [&](Eigen::Index c)
      {
        for (Eigen::Index r = 0; r < nr; ++r)
        {
          out(r,c) = in(rows[std::size_t(r)],c);
        }
      });
    }
    else
    {
      copySlices(nc, parallel, [&](Eigen::Index c)
      {
        for (Eigen::Index r = 0; r < nr; ++r)
        {
          out(r,c) = in(rows[std::size_t(r)],cols[std::size_t(c)]);
        }
      });
    }
  }
}

double multiplyDWC(Matrix &out, const Matrix &a, const Matrix &b,
                   const MatrixAxes &axesA, const MatrixAxes &axesB,
                   unsigned int nt, Matrix &selection,
                   double &selectionUsec, double &gemmUsec)
{
  const auto c = connection(a.cols(), axesA.col, b.rows(), axesB.row, nt);
  checkMatrixSize(a.rows(), b.cols());
  if (!c.left.empty())
  {
    const auto begin = usecond();
    selectMatrix(selection, a, {}, c.left, true);
    selectionUsec += usecond()-begin;
    const auto gemmBegin = usecond();
    A2AContraction::mul(out, selection, b);
    gemmUsec += usecond()-gemmBegin;
  }
  else if (!c.right.empty())
  {
    const auto begin = usecond();
    selectMatrix(selection, b, c.right, {}, true);
    selectionUsec += usecond()-begin;
    const auto gemmBegin = usecond();
    A2AContraction::mul(out, a, selection);
    gemmUsec += usecond()-gemmBegin;
  }
  else
  {
    const auto gemmBegin = usecond();
    A2AContraction::mul(out, a, b);
    gemmUsec += usecond()-gemmBegin;
  }
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
  // ap is row-major (m,n), bp column-major (n,m), so matching flat positions
  // contain ap(r,c) and bp(c,r). This is tr(ap*bp), with no conjugation.
  // The caller parallelizes independent products; no inner OpenMP or critical.
  using Vector = Eigen::Matrix<ComplexD,Eigen::Dynamic,1>;
  const Eigen::Map<const Vector> av(ap->data(),ap->size());
  const Eigen::Map<const Vector> bv(bp->data(),bp->size());
  acc += (av.array()*bv.array()).sum();
  return 8.*double(ap->rows())*double(ap->cols());
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

std::string productStem(const VarietyContractorDWC::ProductPar &p)
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
  reader.push("global"); fields({"trajCounter","nt","output"});
  require(!reader.push("diskVectorDir"), "diskVectorDir is not used; remove it from VarietyContractorDWC XML");
  integers({"nt"},false);
  reader.push("trajCounter"); fields({"start","end","step"});
  integers({"start","end","step"},false); reader.pop(); reader.pop();
  reader.push("a2aMatrix"); require(reader.push("elem"), "empty a2aMatrix list");
  do
  {
    fields({"file","dataset","name","initDenseRow","initDenseCol","timeOffsetRow","timeOffsetCol"});
    require(!reader.push("cacheSize"), "cacheSize is not used; remove it from VarietyContractorDWC XML");
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

void validateParameters(const VarietyContractorDWCPar &par)
{
    require(par.global.nt > 0 &&
            par.global.nt <= unsigned(std::numeric_limits<int>::max()), "invalid nt");
    require(par.global.trajCounter.start < par.global.trajCounter.end &&
            par.global.trajCounter.step > 0,
            "trajCounter requires start < end and step > 0");
    require(!par.global.output.empty(), "empty output directory");
    require(!par.a2aMatrix.empty() && !par.product.empty(), "empty matrix/product list");

    std::set<std::string> names, lastNames, outputs;
    for (const auto &m: par.a2aMatrix)
    {
        require(!m.name.empty() && !m.file.empty() && !m.dataset.empty(),
                "empty matrix name/file/dataset");
        require(names.insert(m.name).second, "duplicate matrix name: " + m.name);
        require(m.name.find('/') == std::string::npos,
                "matrix name must not contain / (flat output): " + m.name);
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
        require(!p.translationAverage, "translationAverage must be false");
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
void checkConnections(const VarietyContractorDWCPar &par,
                      const std::map<std::string, MatrixDimensions> &dimensions,
                      const std::map<std::string, VarietyContractorDWC::A2AMatrixPar> &matrixPar)
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


std::size_t checkedAdd(std::size_t a, std::size_t b)
{
    require(a <= std::numeric_limits<std::size_t>::max()-b,
            "dimension/byte-count addition overflow");
    return a+b;
}

std::size_t matrixBytes(const MatrixDimensions &d)
{
    checkMatrixSize(d.rows, d.cols);
    return checkedMul(checkedMul(std::size_t(d.rows), std::size_t(d.cols)), sizeof(ComplexD));
}

void printMemory(const std::string &name, std::size_t bytes)
{
    std::cout << "  " << name << ": " << bytes << " bytes ("
              << double(bytes)/1.e9 << " GB, "
              << double(bytes)/(1024.*1024.*1024.) << " GiB)" << std::endl;
}

struct IOStats
{
    std::size_t bytes = 0, slices = 0;
    double readUsec = 0., convertUsec = 0.;
};

void printIO(const std::string &name, const IOStats &s)
{
    const double total = s.readUsec+s.convertUsec;
    std::cout << "  " << name << ": " << s.slices << " full slice(s), "
              << s.bytes << " disk bytes; HDF5 " << s.readUsec/1.e6
              << " sec, native assignment " << s.convertUsec/1.e6
              << " sec, combined " << (total > 0. ? double(s.bytes)/total*1.e6/1048576. : 0.)
              << " MiB/s" << std::endl;
}

// At most one of these buffers owns storage. Readers do not own numeric buffers.
struct ComplexValue { double re, im; };
struct ComplexValueF { float re, im; };
struct ReadBuffer
{
    std::vector<ComplexValueF> f;
    std::vector<ComplexValue> d;

    void *prepare(bool useFloat, std::size_t n)
    {
        if (useFloat)
        {
            if (!d.empty())
            {
                std::vector<ComplexValue>().swap(d);
            }
            f.resize(n);
            return f.data();
        }
        if (!f.empty())
        {
            std::vector<ComplexValueF>().swap(f);
        }
        d.resize(n);
        return d.data();
    }
};

// Header/type checks and F/D conversion follow the baseline loader. The HDF5
// hyperslab is always (1, rows, cols): there is no 512-by-512 IO subdivision.
class MatrixReader
{
public:
    MatrixReader(const std::string &filename, const std::string &group, unsigned int nt)
    : filename_(filename), file_(filename, H5F_ACC_RDONLY),
      data_(file_.openGroup(group).openDataSet(HADRONS_A2AM_NAME)),
      space_(data_.getSpace())
    {
        require(space_.getSimpleExtentNdims() == 3, "expected (time,row,col) in " + filename_);
        space_.getSimpleExtentDims(d_);
        require(d_[0] == nt && d_[1] > 0 && d_[2] > 0,
                "invalid matrix/time dimensions in " + filename_);
        const auto limit = std::uintmax_t(std::numeric_limits<Eigen::Index>::max());
        require(std::uintmax_t(d_[1]) <= limit && std::uintmax_t(d_[2]) <= limit,
                "matrix dimension exceeds Eigen::Index");
        dims_ = {Eigen::Index(d_[1]), Eigen::Index(d_[2])};
        checkMatrixSize(dims_.rows, dims_.cols);
        require(data_.getTypeClass() == H5T_COMPOUND, "expected complex compound in " + filename_);
        const auto diskType = data_.getCompType();
        const int re = H5Tget_member_index(diskType.getId(), "re");
        const int im = H5Tget_member_index(diskType.getId(), "im");
        require(diskType.getNmembers() == 2 && re >= 0 && im >= 0,
                "expected re/im complex members in " + filename_);
        const auto rt = diskType.getMemberDataType(unsigned(re));
        const auto it = diskType.getMemberDataType(unsigned(im));
        require(rt.getClass() == H5T_FLOAT && it.getClass() == H5T_FLOAT &&
                rt.getSize() == it.getSize() && (rt.getSize() == 4 || rt.getSize() == 8),
                "expected float/double re and im in " + filename_);
        readFloat_ =
            (H5Tequal(rt.getId(), H5T_IEEE_F32LE) > 0 ||
             H5Tequal(rt.getId(), H5T_IEEE_F32BE) > 0) &&
            (H5Tequal(it.getId(), H5T_IEEE_F32LE) > 0 ||
             H5Tequal(it.getId(), H5T_IEEE_F32BE) > 0);
        const auto &native = readFloat_ ? H5NS::PredType::NATIVE_FLOAT
                                        : H5NS::PredType::NATIVE_DOUBLE;
        memType_.reset(new H5NS::CompType(readFloat_ ? sizeof(ComplexValueF) : sizeof(ComplexValue)));
        memType_->insertMember("re", readFloat_ ? HOFFSET(ComplexValueF,re)
                                              : HOFFSET(ComplexValue,re), native);
        memType_->insertMember("im", readFloat_ ? HOFFSET(ComplexValueF,im)
                                              : HOFFSET(ComplexValue,im), native);
        nValues_ = checkedMul(std::size_t(d_[1]), std::size_t(d_[2]));
        diskBytes_ = checkedMul(nValues_, diskType.getSize());
        bufferBytes_ = checkedMul(nValues_, readFloat_ ? sizeof(ComplexValueF) : sizeof(ComplexValue));
    }

    const MatrixDimensions &dimensions() const { return dims_; }
    std::size_t bufferBytes() const { return bufferBytes_; }

    template <typename Out>
    void readSlice(unsigned int t, Out &out, ReadBuffer &buffer, IOStats &stats)
    {
        require(t < d_[0], "input time out of range in " + filename_);
        double begin = usecond();
        void *raw = buffer.prepare(readFloat_, nValues_);
        out.resize(dims_.rows, dims_.cols);
        stats.convertUsec += usecond()-begin;
        const hsize_t count[3] = {1,d_[1],d_[2]}, offset[3] = {t,0,0};
        H5NS::DataSpace mem(3,count);
        begin = usecond();
        space_.selectHyperslab(H5S_SELECT_SET,count,offset);
        try
        {
            data_.read(raw,*memType_,mem,space_);
        }
        catch (const H5NS::Exception &e)
        {
            throw std::runtime_error("reading '" + filename_ + "' at time "
                                     + std::to_string(t) + ": " + e.getDetailMsg());
        }
        stats.readUsec += usecond()-begin;
        begin = usecond();
        if (readFloat_)
        {
            assignNative(out, buffer.f.data());
        }
        else
        {
            assignNative(out, buffer.d.data());
        }
        stats.convertUsec += usecond()-begin;
        stats.bytes = checkedAdd(stats.bytes,diskBytes_);
        ++stats.slices;
    }

private:
    template <typename Out, typename Value>
    void assignNative(Out &out, const Value *raw) const
    {
        bool bad = false;
        if constexpr (Out::IsRowMajor)
        {
            // Retain the baseline row-major traversal for fixed inputs.
            for (Eigen::Index r = 0; r < dims_.rows; ++r)
            for (Eigen::Index c = 0; c < dims_.cols; ++c)
            {
                const auto z = raw[std::size_t(r)*std::size_t(dims_.cols)+std::size_t(c)];
                bad = bad || !std::isfinite(z.re) || !std::isfinite(z.im);
                out(r,c) = ComplexD(double(z.re),double(z.im));
            }
        }
        else
        {
            // As in the old last-cache assignment, parallelize over columns.
            // Convert straight into column-major storage without a second D slice.
            // Do not throw from an OpenMP worker; report after its barrier.
            std::atomic<bool> badParallel{false};
            thread_for(c, dims_.cols,
            {
                bool badColumn = false;
                for (Eigen::Index r = 0; r < dims_.rows; ++r)
                {
                    const auto z = raw[std::size_t(r)*std::size_t(dims_.cols)+std::size_t(c)];
                    badColumn = badColumn || !std::isfinite(z.re) || !std::isfinite(z.im);
                    out(r,c) = ComplexD(double(z.re),double(z.im));
                }
                if (badColumn)
                {
                    badParallel.store(true,std::memory_order_relaxed);
                }
            });
            bad = badParallel.load(std::memory_order_relaxed);
        }
        if (bad)
        {
            HADRONS_ERROR(Argument, "VarietyContractorDWC: non-finite input in " + filename_);
        }
    }

    std::string filename_;
    H5NS::H5File file_;
    H5NS::DataSet data_;
    H5NS::DataSpace space_;
    std::unique_ptr<H5NS::CompType> memType_;
    hsize_t d_[3];
    MatrixDimensions dims_;
    bool readFloat_;
    std::size_t nValues_, diskBytes_, bufferBytes_;
};

std::unique_ptr<MatrixReader> openMatrix(const VarietyContractorDWC::A2AMatrixPar &p,
                                       unsigned int traj, unsigned int nt)
{
    std::string filename = p.file;
    tokenReplace(filename,"traj",traj);
    try
    {
        return std::unique_ptr<MatrixReader>(new MatrixReader(filename,p.dataset,nt));
    }
    catch (const H5NS::Exception &e)
    {
        throw std::runtime_error("opening '" + filename + "', group '" + p.dataset
                                 + "': " + e.getDetailMsg());
    }
}

struct TimeMatch
{
    std::size_t translation, offset;
};
struct LastUse
{
    std::size_t product, result;
};
struct ProductPlan
{
    std::string stem;
    std::vector<std::size_t> terms, lasts;
    std::vector<unsigned int> translations;
    std::vector<std::vector<TimeMatch>> matches;
    std::size_t nTime, nValues;
};
struct MatrixUse
{
    std::set<unsigned int> fixedTimes;
    std::vector<bool> lastTimes;
    std::vector<LastUse> lastUses;
};
struct FixedGroup
{
    std::vector<std::size_t> products, firstUses, lastUses;
};
struct Plan
{
    std::vector<ProductPlan> products;
    std::vector<MatrixUse> matrices;
    std::vector<FixedGroup> groups;
};

Plan makePlan(const VarietyContractorDWCPar &par)
{
    const auto nt = par.global.nt;
    Plan plan;
    plan.matrices.resize(par.a2aMatrix.size());
    std::map<std::string,std::size_t> ids;
    for (std::size_t k = 0; k < par.a2aMatrix.size(); ++k)
    {
        ids.emplace(par.a2aMatrix[k].name,k);
        plan.matrices[k].lastTimes.assign(nt,false);
    }
    // Lexical name tuples keep equal first factors adjacent, without changing
    // the factor order inside any product.
    std::map<std::vector<std::string>,std::vector<std::size_t>> groups;
    for (std::size_t ip = 0; ip < par.product.size(); ++ip)
    {
        const auto &p = par.product[ip];
        const auto names = strToVec<std::string>(p.terms);
        ProductPlan q;
        q.stem = productStem(p);
        q.nTime = std::size_t(p.max_t)-p.min_t+1;
        const auto translations = parseTimeRange(p.translations,nt);
        q.translations.assign(translations.begin(),translations.end());
        q.nValues = checkedMul(q.translations.size(),q.nTime);
        checkedMul(q.nValues,sizeof(ComplexD));
        q.matches.resize(nt);
        for (const auto &name: names)
        {
            q.terms.push_back(ids.at(name));
        }
        for (std::size_t j = 0; j < p.lastMats.size(); ++j)
        {
            const auto id = ids.at(p.lastMats[j]);
            q.lasts.push_back(id);
            plan.matrices[id].lastUses.push_back({ip,j});
        }
        for (std::size_t i = 0; i < q.translations.size(); ++i)
        {
            const auto dt = q.translations[i];
            for (std::size_t j = 0; j < q.terms.size(); ++j)
            {
                plan.matrices[q.terms[j]].fixedTimes.insert(
                    timeModulo(std::int64_t(p.times[j])+dt,nt));
            }
            const auto rowOffset = checkedMul(i,q.nTime);
            for (std::size_t j = 0; j < q.nTime; ++j)
            {
                const auto u = timeModulo(std::int64_t(dt)+p.min_t+std::int64_t(j),nt);
                q.matches[u].push_back({i,checkedAdd(rowOffset,j)});
            }
        }
        // Merge required times once per product/last pair, not once per output
        // value; the match table is shared by all last matrices.
        for (const auto id: q.lasts)
        for (unsigned int u = 0; u < nt; ++u)
        {
            if (!q.matches[u].empty())
            {
                plan.matrices[id].lastTimes[u] = true;
            }
        }
        groups[names].push_back(ip);
        plan.products.push_back(std::move(q));
    }
    const auto absent = std::numeric_limits<std::size_t>::max();
    std::vector<std::size_t> first(par.a2aMatrix.size(),absent), last(first);
    for (const auto &g: groups)
    {
        const auto ig = plan.groups.size();
        FixedGroup out;
        out.products = g.second;
        for (const auto &name: g.first)
        {
            const auto id = ids.at(name);
            if (first[id] == absent)
            {
                first[id] = ig;
            }
            last[id] = ig;
        }
        plan.groups.push_back(std::move(out));
    }
    for (std::size_t id = 0; id < first.size(); ++id)
    {
        if (first[id] != absent)
        {
            plan.groups[first[id]].firstUses.push_back(id);
            plan.groups[last[id]].lastUses.push_back(id);
        }
    }
    return plan;
}

struct CachedProduct
{
    Matrix value;
    MatrixAxes axes;
};
struct ProductData
{
    std::vector<CachedProduct> fixed;
    std::vector<std::vector<ComplexD>> results;
    std::vector<std::size_t> written;
};

std::size_t traceWorkerCount(const Plan &plan)
{
    std::size_t jobs = 0;
    for (const auto &use: plan.matrices)
    {
        jobs = std::max(jobs,use.lastUses.size());
    }
#ifdef _OPENMP
    const auto threads = std::size_t(std::max(1,omp_get_max_threads()));
#else
    const std::size_t threads = 1;
#endif
    return std::max(std::size_t(1),std::min(jobs,threads));
}

struct TraceWorkspace
{
    Matrix selectedA;
    LastMatrix selectedB;
    double flops = 0.;
    std::size_t traces = 0;
    std::exception_ptr error;
};

// Estimate logical array payloads with conservative group overlap. Retained
// vector capacities, transient reallocations and Eigen/HDF5 workspaces are
// additional. No large numerical arrays are constructed by this estimate.
void reportMemory(const VarietyContractorDWCPar &par, const Plan &plan,
                  const std::vector<MatrixDimensions> &dims,
                  const std::vector<std::size_t> &readBuffers)
{
    std::vector<std::size_t> fixedInput(dims.size()), productBytes(plan.products.size());
    std::vector<std::size_t> work(plan.products.size());
    std::size_t resultBytes = 0, productTotal = 0, resultCount = 0;
    std::size_t nProducts = 0, nDatasets = 0, matchingBytes = 0;
    std::size_t fixedRead = 0, lastRead = 0, lastMatrix = 0, selectedA = 0, selectedB = 0;
    for (std::size_t id = 0; id < dims.size(); ++id)
    {
        fixedInput[id] = checkedMul(plan.matrices[id].fixedTimes.size(),matrixBytes(dims[id]));
        if (!plan.matrices[id].fixedTimes.empty())
        {
            fixedRead = std::max(fixedRead,readBuffers[id]);
        }
        if (!plan.matrices[id].lastUses.empty())
        {
            lastRead = std::max(lastRead,readBuffers[id]);
            lastMatrix = std::max(lastMatrix,matrixBytes(dims[id]));
        }
    }
    for (std::size_t ip = 0; ip < plan.products.size(); ++ip)
    {
        const auto &q = plan.products[ip];
        const auto &p = par.product[ip];
        MatrixDimensions current = dims[q.terms.front()];
        auto axes = matrixAxes(par.a2aMatrix[q.terms.front()],p.times.front(),par.global.nt);
        std::size_t selection = 0, intermediate = 0;
        for (std::size_t j = 1; j < q.terms.size(); ++j)
        {
            const auto &right = dims[q.terms[j]];
            const auto rightAxes = matrixAxes(par.a2aMatrix[q.terms[j]],p.times[j],par.global.nt);
            const auto c = connection(current.cols,axes.col,right.rows,rightAxes.row,par.global.nt);
            if (!c.left.empty())
            {
                selection = std::max(selection,matrixBytes({current.rows,c.active}));
            }
            else if (!c.right.empty())
            {
                selection = std::max(selection,matrixBytes({c.active,right.cols}));
            }
            const MatrixDimensions next{current.rows,right.cols};
            if (j+1 < q.terms.size())
            {
                intermediate = std::max(intermediate,matrixBytes(next));
            }
            current = next;
            axes.col = rightAxes.col;
        }
        // Final results are already counted in productTotal. At most two owned
        // non-final intermediates coexist, including the reused tmp from the
        // preceding translation. The first input remains a borrowed reference.
        work[ip] = checkedAdd(checkedMul(2,intermediate),selection);
        productBytes[ip] = checkedMul(q.translations.size(),matrixBytes(current));
        productTotal = checkedAdd(productTotal,productBytes[ip]);
        nProducts = checkedAdd(nProducts,q.translations.size());
        nDatasets = checkedAdd(nDatasets,q.lasts.size());
        const auto count = checkedMul(q.nValues,q.lasts.size());
        resultCount = checkedAdd(resultCount,count);
        resultBytes = checkedAdd(resultBytes,checkedMul(count,sizeof(ComplexD)));
        matchingBytes = checkedAdd(matchingBytes,checkedMul(q.nValues,sizeof(TimeMatch)));
        for (const auto id: q.lasts)
        {
            const auto &last = dims[id];
            const auto lastAxes = matrixAxes(par.a2aMatrix[id],0,par.global.nt);
            const auto inner = connection(current.cols,axes.col,last.rows,lastAxes.row,par.global.nt);
            const auto closing = connection(last.cols,lastAxes.col,current.rows,axes.row,par.global.nt);
            const auto bytes = matrixBytes({closing.active,inner.active});
            if (!closing.right.empty() || !inner.left.empty())
            {
                selectedA = std::max(selectedA,bytes);
            }
            if (!inner.right.empty() || !closing.left.empty())
            {
                selectedB = std::max(selectedB,bytes);
            }
        }
    }
    std::size_t alive = 0, accumulated = 0, buildPeak = 0;
    for (const auto &g: plan.groups)
    {
        for (const auto id: g.firstUses)
        {
            alive = checkedAdd(alive,fixedInput[id]);
        }
        std::size_t groupWork = 0;
        for (const auto ip: g.products)
        {
            accumulated = checkedAdd(accumulated,productBytes[ip]);
            groupWork = std::max(groupWork,work[ip]);
        }
        buildPeak = std::max(buildPeak,checkedAdd(checkedAdd(checkedAdd(accumulated,alive),groupWork),fixedRead));
        for (const auto id: g.lastUses)
        {
            require(alive >= fixedInput[id],"internal fixed-input lifetime estimate error");
            alive -= fixedInput[id];
        }
    }
    require(alive == 0 && accumulated == productTotal,"internal memory estimate error");
    auto tracePeak = checkedAdd(productTotal,resultBytes);
    tracePeak = checkedAdd(tracePeak,lastRead);
    tracePeak = checkedAdd(tracePeak,lastMatrix);
    const auto workers = traceWorkerCount(plan);
    const auto traceWork = checkedMul(workers,checkedAdd(selectedA,selectedB));
    tracePeak = checkedAdd(tracePeak,traceWork);
    std::cout << "======== Array memory estimate (decimal GB; capacity/transient/library overhead excluded)"
              << std::endl;
    std::cout << "  " << nProducts << " fixed matrices; " << par.product.size()
              << " output files; " << nDatasets << " datasets; " << resultCount
              << " complex output values" << std::endl;
    printMemory("fixed-product cache",productTotal);
    printMemory("cropped results",resultBytes);
    printMemory("fixed-build numeric-array estimate",buildPeak);
    printMemory("trace-phase numeric-array estimate",tracePeak);
    std::cout << "  trace workers: at most " << workers
              << " (independent XML products)" << std::endl;
    printMemory("trace selection workspaces (all workers)",traceWork);
    printMemory("time-match entries (vector capacity/headers additional)",matchingBytes);
}

void buildFixedProducts(const VarietyContractorDWCPar &par, const Plan &plan,
                        unsigned int traj, std::vector<ProductData> &output)
{
    const auto begin = usecond();
    std::map<std::size_t,std::vector<Matrix>> live;
    ReadBuffer buffer;
    double algebraUsec = 0., algebraFlops = 0.;
    double selectionUsec = 0., gemmUsec = 0.;
    std::size_t built = 0;
    for (std::size_t ig = 0; ig < plan.groups.size(); ++ig)
    {
        const auto &g = plan.groups[ig];
        std::cout << "======== Fixed-input group " << ig+1 << "/" << plan.groups.size()
                  << std::endl;
        for (const auto id: g.firstUses)
        {
            const auto &mp = par.a2aMatrix[id];
            std::cout << "Loading fixed input '" << mp.name << "' ("
                      << plan.matrices[id].fixedTimes.size() << " required times)" << std::endl;
            auto input = openMatrix(mp,traj,par.global.nt);
            auto inserted = live.emplace(id,std::vector<Matrix>(par.global.nt));
            require(inserted.second,"internal duplicate fixed-input load");
            auto &slices = inserted.first->second;
            IOStats stats;
            for (const auto t: plan.matrices[id].fixedTimes)
            {
                input->readSlice(t,slices[t],buffer,stats);
            }
            printIO(mp.name,stats);
        }
        // Grouped term names have fixed dimensions. Reuse these work arrays
        // across translations/products, then release them at the group boundary.
        Matrix tmp, selection;
        for (const auto ip: g.products)
        {
            const auto &p = par.product[ip];
            const auto &q = plan.products[ip];
            auto &cached = output[ip].fixed;
            cached.reserve(q.translations.size());
            for (const auto dt: q.translations)
            {
                const auto algebraBegin = usecond();
                const auto firstTime = timeModulo(std::int64_t(p.times.front())+dt,par.global.nt);
                const Matrix &first = live.at(q.terms.front())[firstTime];
                Matrix prod;
                auto axes = matrixAxes(par.a2aMatrix[q.terms.front()],firstTime,par.global.nt);
                if (q.terms.size() == 1)
                {
                    // A one-factor fixed result must survive release of live.
                    prod = first;
                }
                for (std::size_t j = 1; j < q.terms.size(); ++j)
                {
                    const auto t = timeModulo(std::int64_t(p.times[j])+dt,par.global.nt);
                    const Matrix &right = live.at(q.terms[j])[t];
                    const auto rightAxes = matrixAxes(par.a2aMatrix[q.terms[j]],t,par.global.nt);
                    if (j == 1)
                    {
                        algebraFlops += multiplyDWC(prod,first,right,axes,rightAxes,
                                                    par.global.nt,selection,selectionUsec,gemmUsec);
                    }
                    else
                    {
                        algebraFlops += multiplyDWC(tmp,prod,right,axes,rightAxes,
                                                    par.global.nt,selection,selectionUsec,gemmUsec);
                        // Preserve left-to-right multiplication without copying prod.
                        prod.swap(tmp);
                    }
                    axes.col = rightAxes.col;
                }
                cached.push_back({std::move(prod),axes});
                ++built;
                algebraUsec += usecond()-algebraBegin;
            }
        }
        for (const auto id: g.lastUses)
        {
            require(live.erase(id) == 1,"internal fixed-input release error");
            std::cout << "Released fixed input '" << par.a2aMatrix[id].name << "'" << std::endl;
        }
    }
    require(live.empty(),"internal fixed-input lifetime error");
    std::cout << "Fixed products complete: " << built << "; phase " << (usecond()-begin)/1.e6
              << " sec, algebra/owned assignment " << algebraUsec/1.e6 << " sec, "
              << (algebraUsec > 0. ? algebraFlops/algebraUsec/1.e3 : 0.) << " GFlop/s"
              << std::endl;
    std::cout << "  Fixed selection: " << selectionUsec/1.e6
              << " sec; GEMM: " << gemmUsec/1.e6 << " sec, "
              << (gemmUsec > 0. ? algebraFlops/gemmUsec/1.e3 : 0.)
              << " GFlop/s (GEMM only)" << std::endl;
    // live and its conversion buffer die before the result arrays are allocated.
}

void contractLastMatrices(const VarietyContractorDWCPar &par, const Plan &plan,
                          std::vector<std::unique_ptr<MatrixReader>> &readers,
                          std::vector<ProductData> &output)
{
    for (std::size_t ip = 0; ip < plan.products.size(); ++ip)
    {
        const auto &q = plan.products[ip];
        auto &d = output[ip];
        require(d.fixed.size() == q.translations.size(),"internal missing fixed products");
        d.results.resize(q.lasts.size());
        d.written.assign(q.lasts.size(),0);
        for (auto &result: d.results)
        {
            result.assign(q.nValues,ComplexD(0.,0.));
        }
    }
    ReadBuffer buffer;
    LastMatrix last;
    // Buffers survive EMF/time batches. The full last slice is shared read-only;
    // only selected matrices and counters belong to individual trace workers.
    std::vector<TraceWorkspace> workspace(traceWorkerCount(plan));
    std::vector<IOStats> stats(par.a2aMatrix.size());
    double traceUsec = 0., traceFlops = 0.;
    std::size_t traces = 0;
    const auto begin = usecond();
    for (unsigned int u = 0; u < par.global.nt; ++u)
    {
        const auto stepBegin = usecond();
        std::cout << "======== Last-matrix absolute time " << u << "/" << par.global.nt-1
                  << std::endl;
        // Keep the XML last-matrix order; numeric buffers are shared, not per-reader.
        for (const auto &name: par.a2aMatrixLast)
        {
            // Lookup outside all element/trace loops; names are unique.
            const auto it = std::find_if(par.a2aMatrix.begin(),par.a2aMatrix.end(),
                [&](const VarietyContractorDWC::A2AMatrixPar &m){return m.name == name;});
            const auto id = std::size_t(it-par.a2aMatrix.begin());
            const auto &use = plan.matrices.at(id);
            if (!use.lastTimes[u])
            {
                continue;
            }
            require(bool(readers[id]),"internal missing last-matrix reader");
            readers[id]->readSlice(u,last,buffer,stats[id]);
            const auto axes = matrixAxes(par.a2aMatrix[id],u,par.global.nt);
            const auto traceBegin = usecond();
            for (auto &w: workspace)
            {
                w.flops = 0.;
                w.traces = 0;
                w.error = nullptr;
            }
            const auto nUses = use.lastUses.size();
#ifdef _OPENMP
            const int workers = int(std::min(workspace.size(),nUses));
#pragma omp parallel for schedule(dynamic) num_threads(workers)
#endif
            for (std::size_t iu = 0; iu < nUses; ++iu)
            {
#ifdef _OPENMP
                const auto worker = std::size_t(omp_get_thread_num());
#else
                const std::size_t worker = 0;
#endif
                auto &w = workspace[worker];
                if (w.error)
                {
                    continue;
                }
                try
                {
                    const auto &useProduct = use.lastUses[iu];
                    const auto ip = useProduct.product;
                    auto &d = output[ip];
                    auto &result = d.results[useProduct.result];
                    const auto &matches = plan.products[ip].matches[u];
                    double productFlops = 0.;
                    // A product/last pair occurs only once in use.lastUses.
                    // Its results and written counter therefore have one owner.
                    for (const auto &match: matches)
                    {
                        const auto &fixed = d.fixed[match.translation];
                        productFlops += traceDWC(result[match.offset],fixed.value,last,
                                                fixed.axes,axes,par.global.nt,
                                                w.selectedA,w.selectedB);
                    }
                    d.written[useProduct.result] += matches.size();
                    w.flops += productFlops;
                    w.traces += matches.size();
                }
                catch (...)
                {
                    // No exception may escape an OpenMP worker. Each worker
                    // owns its error slot; rethrow on the main thread below.
                    w.error = std::current_exception();
                }
            }
            // Includes the parallel-loop barrier. Use elapsed wall time, not
            // the sum of individual worker durations, for the node's GFlop/s.
            traceUsec += usecond()-traceBegin;
            for (const auto &w: workspace)
            {
                if (w.error)
                {
                    std::rethrow_exception(w.error);
                }
                traceFlops += w.flops;
                traces += w.traces;
            }
        }
        std::cout << "  time step completed in " << (usecond()-stepBegin)/1.e6
                  << " sec; cumulative traces " << traces << std::endl;
    }
    for (std::size_t ip = 0; ip < plan.products.size(); ++ip)
    {
        for (const auto n: output[ip].written)
        {
            require(n == plan.products[ip].nValues,"internal missing/duplicate trace count");
        }
    }
    std::cout << "======== Last-matrix IO totals" << std::endl;
    for (std::size_t id = 0; id < stats.size(); ++id)
    {
        if (stats[id].slices)
        {
            printIO(par.a2aMatrix[id].name,stats[id]);
        }
    }
    std::cout << "Traces complete: " << traces << "; phase " << (usecond()-begin)/1.e6
              << " sec, trace algebra " << traceUsec/1.e6 << " sec, "
              << (traceUsec > 0. ? traceFlops/traceUsec/1.e3 : 0.) << " GFlop/s" << std::endl;
}

void writeResults(const VarietyContractorDWCPar &par, const Plan &plan,
                  unsigned int traj, const std::vector<ProductData> &output)
{
    const auto begin = usecond();
    // Direct HDF5 only: do NOT instantiate ResultWriter/CorrelatorResult.
    // Grid's complex memory type is reused, but no Grid serializer is invoked.
    H5NS::CompType diskType(sizeof(ComplexValue));
    diskType.insertMember("re",HOFFSET(ComplexValue,re),H5NS::PredType::IEEE_F64LE);
    diskType.insertMember("im",HOFFSET(ComplexValue,im),H5NS::PredType::IEEE_F64LE);
    for (std::size_t ip = 0; ip < plan.products.size(); ++ip)
    {
        const auto &q = plan.products[ip];
        const auto filename = par.global.output + "/"
            + ModuleBase::resultFilename(q.stem,traj,"h5");
        makeFileDir(filename);
        std::cout << "Saving " << q.lasts.size() << " flat datasets to '" << filename << "'"
                  << std::endl;
        // Match the existing writer's replacement semantics for an existing file.
        H5NS::H5File file(filename,H5F_ACC_TRUNC);
        const hsize_t shape[2] = {hsize_t(q.translations.size()),hsize_t(q.nTime)};
        H5NS::DataSpace space(2,shape);
        for (std::size_t ell = 0; ell < q.lasts.size(); ++ell)
        {
            const auto &values = output[ip].results[ell];
            for (const auto z: values)
            {
                if (!std::isfinite(z.real()) || !std::isfinite(z.imag()))
                {
                    HADRONS_ERROR(Argument,"VarietyContractorDWC: non-finite correlator in " + q.stem);
                }
            }
            const auto name = q.stem + "_" + par.a2aMatrix[q.lasts[ell]].name;
            // Default creation properties: contiguous, uncompressed, fixed shape.
            auto dataset = file.createDataSet(name,diskType,space);
            dataset.write(values.data(),Hdf5Type<ComplexD>::type());
        }
    }
    std::cout << "Output complete: " << (usecond()-begin)/1.e6 << " sec" << std::endl;
}

int main(int argc, char *argv[])
{
    if (argc != 2)
    {
        std::cerr << "usage: " << argv[0] << " <parameter file>" << std::endl;
        return EXIT_FAILURE;
    }
    try
    {
        checkSingleProcess();
        H5NS::Exception::dontPrint();
        checkXmlFields(argv[1]);
        VarietyContractorDWCPar par;
        XmlReader reader(argv[1]);
        read(reader,"global",par.global);
        read(reader,"a2aMatrix",par.a2aMatrix);
        read(reader,"a2aMatrixLast",par.a2aMatrixLast);
        read(reader,"product",par.product);
        validateParameters(par);
        const auto plan = makePlan(par);
        std::map<std::string,VarietyContractorDWC::A2AMatrixPar> matrixPar;
        for (const auto &m: par.a2aMatrix)
        {
            matrixPar.emplace(m.name,m);
        }
        Eigen::initParallel();
        std::cout << "VarietyContractorDWC: one process; Eigen reports "
                  << Eigen::nbThreads() << " thread(s); no lattice Grid; no DiskVector"
                  << std::endl;
        for (unsigned int traj = par.global.trajCounter.start;
             traj < par.global.trajCounter.end; )
        {
            const auto begin = usecond();
            std::cout << ":::::::: Trajectory " << traj << std::endl;
            // Every reader and numeric destination belongs to this trajectory.
            std::vector<std::unique_ptr<MatrixReader>> lastReaders(par.a2aMatrix.size());
            std::vector<MatrixDimensions> dims(par.a2aMatrix.size());
            std::vector<std::size_t> readBuffers(par.a2aMatrix.size());
            std::map<std::string,MatrixDimensions> dimensions;
            for (std::size_t id = 0; id < par.a2aMatrix.size(); ++id)
            {
                const auto &m = par.a2aMatrix[id];
                auto input = openMatrix(m,traj,par.global.nt);
                dims[id] = input->dimensions();
                readBuffers[id] = input->bufferBytes();
                dimensions.emplace(m.name,dims[id]);
                std::cout << "Header '" << m.name << "': " << par.global.nt << " x "
                          << dims[id].rows << " x " << dims[id].cols << std::endl;
                if (!plan.matrices[id].lastUses.empty())
                {
                    lastReaders[id] = std::move(input);
                }
            }
            checkConnections(par,dimensions,matrixPar);
            reportMemory(par,plan,dims,readBuffers);
            std::vector<ProductData> output(par.product.size());
            buildFixedProducts(par,plan,traj,output);
            contractLastMatrices(par,plan,lastReaders,output);
            // No input handles or fixed-product arrays are needed while writing.
            lastReaders.clear();
            for (auto &p: output)
            {
                std::vector<CachedProduct>().swap(p.fixed);
            }
            writeResults(par,plan,traj,output);
            std::cout << "Trajectory " << traj << " complete in " << (usecond()-begin)/1.e6
                      << " sec" << std::endl;
            if (par.global.trajCounter.step >= par.global.trajCounter.end-traj)
            {
                break;
            }
            traj += par.global.trajCounter.step;
        }
    }
    catch (const std::bad_alloc &)
    {
        std::cerr << "VarietyContractorDWC: memory allocation failed; see the array estimate above"
                  << std::endl;
        return EXIT_FAILURE;
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
