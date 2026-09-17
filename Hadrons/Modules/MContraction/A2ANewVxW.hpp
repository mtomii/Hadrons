/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: Hadrons/Modules/MContraction/A2ANewVxW.hpp

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
#ifndef Hadrons_MContraction_A2ANewVxW_hpp_
#define Hadrons_MContraction_A2ANewVxW_hpp_
#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <limits>

/* Based on MContraction/A2AVxW.hpp, with full-field inputs as in A2ANewVxPi.
 *
 * mat[u*nSpace+p] = sum_r outerProduct(left[leftIndex[r]](x),
 *                                    right[rightIndex[r]](x)).
 * Grid's outerProduct conjugates the RIGHT spinor. Inputs are ordinary
 * column spinors, NOT the old already-conjugated right-hand intermediates.
 * x has time (tref+delt_min+u) mod Nt and local spatial SIMD index p.
 *
 * Dense means that the mode index is compact, possibly inherited from a
 * meson field. tDense>=0 selects a fixed source time; tDense=-1 uses x's
 * time (for a raw dense W). No hit normalization or gamma5 is inserted.
 * Output: exactly std::vector<SpinColourMatrix_v>, with the old window
 * ordering, not a full PropagatorField. No I/O or global sum is performed.
 * CPU, full non-checkerboard 4D input fields, MPI_t=SIMD_t=1.
 */

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                Close A2A mode indices in a selected time window             *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MContraction)
class A2ANewVxWPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(A2ANewVxWPar,
                                    std::string, left,
                                    std::string, right,
                                    int,         tref,
                                    int,         delt_min,
                                    int,         delt_max,
                                    bool,        ifDenseLeft,
                                    bool,        ifDenseRight,
                                    int,         nLow,
                                    int,         tDense);
};

template <typename FImpl>
class TA2ANewVxW: public Module<A2ANewVxWPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    typedef typename FImpl::SiteSpinor vobj;
    typedef typename vobj::vector_type vector_type;
    typedef iSpinColourMatrix<vector_type> SpinColourMatrix_v;
public:
    // constructor
    TA2ANewVxW(const std::string name);
    // destructor
    virtual ~TA2ANewVxW(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
    int nt_{0};
    std::size_t nTime_{0}, nSpace_{0}, nOutput_{0};
    std::size_t nLeft_{0}, nRight_{0}, nActive_{0};
};
MODULE_REGISTER_TMP(A2ANewVxW, TA2ANewVxW<FIMPL>, MContraction);

/******************************************************************************
 *                        TA2ANewVxW implementation                           *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TA2ANewVxW<FImpl>::TA2ANewVxW(const std::string name)
: Module<A2ANewVxWPar>(name)
{}
// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TA2ANewVxW<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().left, par().right};

    return in;
}

template <typename FImpl>
std::vector<std::string> TA2ANewVxW<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};

    return out;
}
// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TA2ANewVxW<FImpl>::setup(void)
{
    if (par().left.empty() || par().right.empty() || par().nLow < 0 ||
        par().tDense < -1)
    {
        HADRONS_ERROR(Argument, "A2ANewVxW: invalid input name, nLow or tDense");
    }
    assert(!(par().ifDenseLeft && par().ifDenseRight));
    if (par().ifDenseLeft && par().ifDenseRight)
    {
        HADRONS_ERROR(Argument, "A2ANewVxW: both contracted indices cannot be dense W");
    }

    const auto &leftV  = envGet(std::vector<FermionField>, par().left);
    const auto &rightW = envGet(std::vector<FermionField>, par().right);
    if (leftV.empty() || rightW.empty())
    {
        HADRONS_ERROR(Size, "A2ANewVxW: empty input vector set");
    }
    GridBase *grid = leftV[0].Grid();
    if (grid->_ndimension != 4 || grid->_isCheckerBoarded ||
        env().getObjectLs(par().left) != 1 || env().getObjectLs(par().right) != 1)
    {
        HADRONS_ERROR(Argument, "A2ANewVxW requires full, non-checkerboard 4D fields");
    }
    if (grid->_processors[Tp] != 1 || grid->_simd_layout[Tp] != 1)
    {
        HADRONS_ERROR(Argument, "A2ANewVxW requires MPI_t = SIMD_t = 1");
    }
    for (const auto &field: leftV)
    {
        if (field.Grid() != grid)
        {
            HADRONS_ERROR(Size, "A2ANewVxW: left fields have different Grids");
        }
    }
    for (const auto &field: rightW)
    {
        if (field.Grid() != grid)
        {
            HADRONS_ERROR(Size, "A2ANewVxW: left/right fields have different Grids");
        }
    }

    nt_ = grid->_gdimensions[Tp];
    const std::int64_t nTime = std::int64_t(par().delt_max)
                              - std::int64_t(par().delt_min) + 1;
    if (nt_ <= 0 || nTime < 1 || nTime > nt_)
    {
        HADRONS_ERROR(Argument, "A2ANewVxW requires 1 <= time-window length <= Nt");
    }
    nTime_ = static_cast<std::size_t>(nTime);
    nLeft_ = leftV.size();
    nRight_ = rightW.size();
    const std::size_t nLow = par().nLow;
    const std::size_t nsc = std::size_t(Ns)*FImpl::Dimension;
    if (nLeft_ < nLow || nRight_ < nLow)
    {
        HADRONS_ERROR(Size, "A2ANewVxW: nLow exceeds an input vector count");
    }
    nActive_ = nLeft_;
    if (par().ifDenseLeft || par().ifDenseRight)
    {
        const std::size_t compact = par().ifDenseLeft ? nLeft_ : nRight_;
        const std::size_t expanded = par().ifDenseLeft ? nRight_ : nLeft_;
        const std::size_t high = compact - nLow;
        if (high % nsc != 0)
        {
            HADRONS_ERROR(Size, "A2ANewVxW: compact high size is not a multiple of Ns*Nc");
        }
        if (high > (std::numeric_limits<std::size_t>::max() - nLow)/nt_)
        {
            HADRONS_ERROR(Size, "A2ANewVxW: expanded dimension overflow");
        }
        if (expanded != nLow + high*std::size_t(nt_))
        {
            HADRONS_ERROR(Size, "A2ANewVxW: dense/expanded dimensions do not match");
        }
        nActive_ = compact;
    }
    else if (nLeft_ != nRight_)
    {
        HADRONS_ERROR(Size, "A2ANewVxW: non-dense left/right sizes differ");
    }

    if (grid->_slice_nblock[Tp] <= 0 || grid->_slice_block[Tp] <= 0 ||
        grid->_slice_stride[Tp] <= 0 || grid->_ostride[Tp] <= 0)
    {
        HADRONS_ERROR(Size, "A2ANewVxW: invalid spatial slice layout");
    }
    const std::size_t e1 = grid->_slice_nblock[Tp];
    const std::size_t e2 = grid->_slice_block[Tp];
    if (e1 > std::numeric_limits<std::size_t>::max()/e2)
    {
        HADRONS_ERROR(Size, "A2ANewVxW: spatial size overflow");
    }
    nSpace_ = e1*e2;
    if (nSpace_ > std::numeric_limits<std::size_t>::max()/std::size_t(nt_) ||
        nSpace_*std::size_t(nt_) != std::size_t(grid->oSites()))
    {
        HADRONS_ERROR(Size, "A2ANewVxW: spatial slice size disagrees with input Grid");
    }
    // nTime_<=nt_, so the preceding check also protects this product.
    nOutput_ = nTime_*nSpace_;
    if (nOutput_ > std::numeric_limits<std::size_t>::max()/sizeof(SpinColourMatrix_v) ||
        nOutput_ > std::vector<SpinColourMatrix_v>().max_size())
    {
        HADRONS_ERROR(Size, "A2ANewVxW: output size exceeds addressable range");
    }

    // Keep the EXACT old vector type (standard allocator) for the MT modules.
    // Grid's automatic MemoryProfiler need not account for this allocation.
    envCreate(std::vector<SpinColourMatrix_v>, getName(), 1, nOutput_, Zero());
    LOG(Message) << "A2ANewVxW: " << nActive_ << " active modes, "
                 << nTime_ << " time slices, " << nSpace_
                 << " local spatial SIMD objects/time; output "
                 << sizeString(nOutput_*sizeof(SpinColourMatrix_v)) << "/rank"
                 << " (standard-vector allocation; reported explicitly)" << std::endl;
}
// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TA2ANewVxW<FImpl>::execute(void)
{
    const auto &leftV  = envGet(std::vector<FermionField>, par().left);
    const auto &rightW = envGet(std::vector<FermionField>, par().right);
    auto       &mat = envGet(std::vector<SpinColourMatrix_v>, getName());
    assert(leftV.size() == nLeft_ && rightW.size() == nRight_);
    assert(mat.size() == nOutput_);
    if (leftV.size() != nLeft_ || rightW.size() != nRight_ || mat.size() != nOutput_)
    {
        HADRONS_ERROR(Size, "A2ANewVxW: object sizes changed after setup");
    }
    GridBase *grid = leftV[0].Grid();
    const int tref = int((std::int64_t(par().tref) % nt_ + nt_) % nt_);
    const bool denseLeft = par().ifDenseLeft;
    const bool denseRight = par().ifDenseRight;
    const bool localDenseTime = (denseLeft || denseRight) && par().tDense == -1;
    const int fixedDenseTime = par().tDense >= 0 ? par().tDense % nt_ : 0;
    const std::size_t nLow = par().nLow;
    const std::size_t nsc = std::size_t(Ns)*FImpl::Dimension;
    const std::size_t e2 = grid->_slice_block[Tp];
    const std::size_t stride = grid->_slice_stride[Tp];

    // r includes low modes; in a dense case r=nLow+h*nsc+a above nLow.
    // The expanded partner is nLow+(h*nt+sourceTime)*nsc+a.
    std::vector<std::size_t> leftIndex(nActive_), rightIndex(nActive_);
    for (std::size_t r = 0; r < nActive_; ++r)
    {
        leftIndex[r] = r;
        rightIndex[r] = r;
    }
    std::vector<LatticeView<vobj>> leftView, rightView;
    leftView.reserve(nActive_);
    rightView.reserve(nActive_);

    startTimer("V x W");
    try
    {
        for (std::size_t u = 0; u < nTime_; ++u)
        {
            const std::int64_t delta = std::int64_t(par().delt_min) + std::int64_t(u);
            const int tx = int(((std::int64_t(tref) + delta) % nt_ + nt_) % nt_);
            // Fixed pairing: map/open once. Local pairing: map/open per time.
            // Branches and integer divisions stay OUTSIDE spatial/mode sums.
            if (u == 0 || localDenseTime)
            {
                const std::size_t sourceTime = localDenseTime ? tx : fixedDenseTime;
                if (denseLeft)
                {
                    for (std::size_t r = nLow; r < nActive_; ++r)
                    {
                        const std::size_t h = (r - nLow)/nsc;
                        const std::size_t a = (r - nLow)%nsc;
                        rightIndex[r] = nLow + (h*nt_ + sourceTime)*nsc + a;
                    }
                }
                else if (denseRight)
                {
                    for (std::size_t r = nLow; r < nActive_; ++r)
                    {
                        const std::size_t h = (r - nLow)/nsc;
                        const std::size_t a = (r - nLow)%nsc;
                        leftIndex[r] = nLow + (h*nt_ + sourceTime)*nsc + a;
                    }
                }
                // Views are pointers to existing fields, not field copies.
                // Store them in active-mode order; the hot loop uses only r.
                for (std::size_t r = 0; r < nActive_; ++r)
                {
                    leftView.push_back(leftV[leftIndex[r]].View(CpuRead));
                    rightView.push_back(rightW[rightIndex[r]].View(CpuRead));
                }
            }
            LOG(Message) << "A2ANewVxW: tref=" << tref << ", delta=" << delta
                         << ", t_op=" << tx << ", window index=" << u << std::endl;
            const std::size_t so = std::size_t(tx)*grid->_ostride[Tp];

            // Each thread owns the whole SIMD matrix at one spatial position.
            // No packing/GEMM, dense branch, or global reduction in this loop.
            thread_for(p, nSpace_, {
                const std::size_t ss = so + (p/e2)*stride + p%e2;
                SpinColourMatrix_v sum = Zero();
                for (std::size_t r = 0; r < nActive_; ++r)
                {
                    sum += outerProduct(leftView[r][ss], rightView[r][ss]);
                }
                mat[u*nSpace_ + p] = sum;
            });

            if (localDenseTime)
            {
                for (auto &view: rightView) view.ViewClose();
                for (auto &view: leftView)  view.ViewClose();
                rightView.clear();
                leftView.clear();
            }
        }
    }
    catch (...)
    {
        for (auto &view: rightView) view.ViewClose();
        for (auto &view: leftView)  view.ViewClose();
        stopTimer("V x W");
        throw;
    }
    for (auto &view: rightView) view.ViewClose();
    for (auto &view: leftView)  view.ViewClose();
    stopTimer("V x W");
}
END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE
#endif // Hadrons_MContraction_A2ANewVxW_hpp_
