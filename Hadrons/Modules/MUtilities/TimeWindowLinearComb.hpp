/*
 * TimeWindowLinearComb.hpp, part of Hadrons.
 *
 * Module structure follows MUtilities/PropagatorToTimeWindow.hpp;
 * matrix type and ordering follow MContraction/A2ANewVxW.hpp.
 *
 * Hadrons is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * Hadrons is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * See the file LICENSE in the top level Hadrons directory.
 */

#ifndef Hadrons_MUtilities_TimeWindowLinearComb_hpp_
#define Hadrons_MUtilities_TimeWindowLinearComb_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <cstddef>
#include <limits>
#include <string>
#include <vector>

/* out[i] = alpha*input1[i] + beta*input2[i].
 *
 * All three objects have exactly the standard-allocator type
 * std::vector<SpinColourMatrix_v> used by A2ANewVxW and the MT modules.
 * A2ANewVxW orders elements as i = u*nSpace + p, with u the window-time
 * index and p the local spatial SIMD-object index.
 * XML nTime is the full input/output window length, not a prefix to select.
 * setup() uses nTime and the standard full 4D Grid only; it never retrieves
 * input objects. This avoids recursively constructing their producers during
 * memory profiling. Requires MPI_t = SIMD_t = 1, as for A2ANewVxW.
 *
 * The vector carries no Grid or time-window metadata. The caller must
 * ensure that corresponding indices have the intended physical meaning.
 * This module preserves the order and combines corresponding elements;
 * it does not shift times, conjugate, normalize, or perform an MPI sum.
 *
 * Inputs are read-only; output is a separate vector named getName().
 * Identical input names and, with nTime=0, two empty inputs are supported.
 * CPU thread_for operates directly on the inline SIMD matrix values.
 */

BEGIN_HADRONS_NAMESPACE

BEGIN_MODULE_NAMESPACE(MUtilities)

class TimeWindowLinearCombPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(TimeWindowLinearCombPar,
                                    std::string, input1,
                                    std::string, input2,
                                    ComplexD,    alpha,
                                    ComplexD,    beta,
                                    int,         nTime);
};

template <typename FImpl>
class TTimeWindowLinearComb: public Module<TimeWindowLinearCombPar>
{
public:
    typedef typename FImpl::SiteSpinor vobj;
    typedef typename vobj::vector_type vector_type;
    typedef iSpinColourMatrix<vector_type> SpinColourMatrix_v;
    typedef std::vector<SpinColourMatrix_v> MatrixVector;
    typedef typename SpinColourMatrix_v::scalar_type Coefficient;
public:
    TTimeWindowLinearComb(const std::string name);
    virtual ~TTimeWindowLinearComb(void) {};
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    virtual void setup(void);
    virtual void execute(void);
private:
    std::size_t nOutput_{0};
};

MODULE_REGISTER_TMP(TimeWindowLinearComb,
                    TTimeWindowLinearComb<FIMPL>, MUtilities);

// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TTimeWindowLinearComb<FImpl>::TTimeWindowLinearComb(const std::string name)
: Module<TimeWindowLinearCombPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TTimeWindowLinearComb<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().input1};
    if (par().input2 != par().input1)
    {
        in.push_back(par().input2);
    }
    return in;
}

template <typename FImpl>
std::vector<std::string> TTimeWindowLinearComb<FImpl>::getOutput(void)
{
    return {getName()};
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TTimeWindowLinearComb<FImpl>::setup(void)
{
    if (par().input1.empty() || par().input2.empty())
    {
        HADRONS_ERROR(Argument, "TimeWindowLinearComb: empty input name");
    }
    if (getName() == par().input1 || getName() == par().input2)
    {
        HADRONS_ERROR(Argument,
                      "TimeWindowLinearComb: output name must differ from input names");
    }

    // Shape comes from XML and the standard Grid, not from input objects.
    GridCartesian *grid = this->env().template getGrid<vector_type>();
    if (grid->_ndimension != 4 || grid->_isCheckerBoarded)
    {
        HADRONS_ERROR(Argument, "TimeWindowLinearComb requires a full 4D Grid");
    }
    if (grid->_processors[Tp] != 1 || grid->_simd_layout[Tp] != 1)
    {
        HADRONS_ERROR(Argument, "TimeWindowLinearComb requires MPI_t = SIMD_t = 1");
    }
    const int nt = grid->_gdimensions[Tp];
    if (nt <= 0 || par().nTime < 0 || par().nTime > nt)
    {
        HADRONS_ERROR(Argument, "TimeWindowLinearComb requires 0 <= nTime <= Nt");
    }
    const std::size_t nSites = grid->oSites();
    if (nSites == 0 || nSites % std::size_t(nt) != 0)
    {
        HADRONS_ERROR(Size, "TimeWindowLinearComb: inconsistent standard Grid layout");
    }
    const std::size_t nSpace = nSites/std::size_t(nt);
    // nTime <= Nt ensures that this product is no larger than nSites.
    nOutput_ = std::size_t(par().nTime)*nSpace;
    if (nOutput_ > std::numeric_limits<std::size_t>::max()/sizeof(SpinColourMatrix_v) ||
        nOutput_ > MatrixVector().max_size())
    {
        HADRONS_ERROR(Size, "TimeWindowLinearComb: output size exceeds addressable range");
    }

    // Construct the output by size; do not copy either input or a fill value.
    // All elements are overwritten by execute(). Element construction still
    // has its ordinary initialization cost.
    envCreate(MatrixVector, getName(), 1, nOutput_);

    // Standard-vector allocation is not necessarily counted by Grid's
    // automatic MemoryProfiler. This message reports its payload explicitly.
    LOG(Message) << "TimeWindowLinearComb: output '" << getName()
                 << "', " << par().nTime << " time slices, "
                 << nOutput_ << " SIMD matrix elements, "
                 << sizeString(nOutput_*sizeof(SpinColourMatrix_v)) << "/rank"
                 << " (standard-vector allocation; reported explicitly)"
                 << std::endl;
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TTimeWindowLinearComb<FImpl>::execute(void)
{
    const auto &a = envGet(MatrixVector, par().input1);
    const auto &b = envGet(MatrixVector, par().input2);
    auto       &out = envGet(MatrixVector, getName());
    if (a.size() != nOutput_ || b.size() != nOutput_ || out.size() != nOutput_)
    {
        HADRONS_ERROR(Size,
                      "TimeWindowLinearComb: sizes do not match nTime and Grid; expected "
                      + std::to_string(nOutput_) + ", input1="
                      + std::to_string(a.size()) + ", input2="
                      + std::to_string(b.size()) + ", output="
                      + std::to_string(out.size()));
    }

    // Convert once to the scalar precision used by the input/output matrices.
    const Coefficient alpha(par().alpha);
    const Coefficient beta (par().beta);

    LOG(Message) << "TimeWindowLinearComb: '" << getName() << "' = "
                 << alpha << " * '" << par().input1 << "' + "
                 << beta << " * '" << par().input2 << "'" << std::endl;

    startTimer("Time-window linear combination");
    thread_for(i, nOutput_,
    {
        out[i] = alpha*a[i] + beta*b[i];
    });
    stopTimer("Time-window linear combination");
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MUtilities_TimeWindowLinearComb_hpp_
