/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: Hadrons/Modules/MContraction/A2ANewPixW.hpp

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
#ifndef Hadrons_MContraction_A2ANewPixW_hpp_
#define Hadrons_MContraction_A2ANewPixW_hpp_
#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/A2AMatrix.hpp>
#include <Hadrons/DiskVector.hpp>
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <limits>

/* Based on MContraction/A2APixW.hpp.
 *
 * out_i(x) = sum_r right[rightIndex[r]](x)
 *                 * conjugate(M(tmes)(i, mesIndex[r])).
 * This is a COLUMN spinor, i.e. the complex conjugate of the old
 * A2APixW output. Close it with sum_i outerProduct(left_i, out_i).
 *
 * Only the contracted (column) meson index is mapped. The row index is
 * preserved, including any compact dense-W index and its time tmes.
 * No additional hit normalization is applied.
 *
 * CPU execution; full 4D fields, MPI_t = SIMD_t = 1. Spatial SIMD lanes
 * are explicitly packed into scalar ComplexD matrices. No Grid changes.
 */

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                   Meson field times an A2A vector set                      *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MContraction)
class A2ANewPixWPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(A2ANewPixWPar,
                                    std::string, right,
                                    std::string, mes,
                                    int,         tmes,
                                    int,         delt_max,
                                    int,         delt_min,
                                    bool,        ifDenseRight,
                                    bool,        ifDenseMes,
                                    int,         nLow,
                                    int,         size,
                                    int,         siteBlock,
                                    int,         modeBlock);
};

template <typename FImpl>
class TA2ANewPixW: public Module<A2ANewPixWPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    typedef typename FImpl::SiteSpinor vobj;
    typedef typename vobj::scalar_object sobj;
    typedef typename vobj::vector_type vector_type;
    typedef A2AMatrix<ComplexD> Matrix;
    typedef A2AMatrixTr<ComplexD> BufferMatrix;
    typedef hostVector<ComplexD> Buffer;
public:
    // constructor
    TA2ANewPixW(const std::string name);
    // destructor
    virtual ~TA2ANewPixW(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
    int nt_{0}, nSimd_{0}, nThread_{1};
    std::size_t nSpace_{0}, siteBlock_{0}, modeBlock_{0}, qMax_{0};
    std::size_t nActive_{0}, nRight_{0}, nMesCols_{0};
    std::size_t rightStride_{0}, resultStride_{0}, mesStride_{0};
};
MODULE_REGISTER_TMP(A2ANewPixW, TA2ANewPixW<FIMPL>, MContraction);

/******************************************************************************
 *                       TA2ANewPixW implementation                           *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TA2ANewPixW<FImpl>::TA2ANewPixW(const std::string name)
: Module<A2ANewPixWPar>(name)
{}
// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TA2ANewPixW<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().right, par().mes};

    return in;
}

template <typename FImpl>
std::vector<std::string> TA2ANewPixW<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};

    return out;
}
// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TA2ANewPixW<FImpl>::setup(void)
{
    if (par().right.empty() || par().mes.empty() || par().nLow < 0 ||
        par().size <= 0 || par().siteBlock <= 0 || par().modeBlock <= 0)
    {
        HADRONS_ERROR(Argument, "A2ANewPixW: invalid input name, size or block size");
    }
    if (par().ifDenseRight && par().ifDenseMes)
    {
        HADRONS_ERROR(Argument, "A2ANewPixW: both contracted indices cannot be dense W");
    }

    const auto &rightW = envGet(std::vector<FermionField>, par().right);
    if (rightW.empty())
    {
        HADRONS_ERROR(Size, "A2ANewPixW: empty right vector set");
    }
    GridBase *grid = rightW[0].Grid();
    if (grid->_ndimension != 4 || grid->_isCheckerBoarded ||
        env().getObjectLs(par().right) != 1)
    {
        HADRONS_ERROR(Argument, "A2ANewPixW requires full, non-checkerboard 4D fields");
    }
    if (grid->_processors[Tp] != 1 || grid->_simd_layout[Tp] != 1)
    {
        HADRONS_ERROR(Argument, "A2ANewPixW requires MPI_t = SIMD_t = 1");
    }
    for (const auto &field: rightW)
    {
        if (field.Grid() != grid)
        {
            HADRONS_ERROR(Size, "A2ANewPixW: right fields have different Grids");
        }
    }

    nt_ = grid->_gdimensions[Tp];
    const std::int64_t nTime = std::int64_t(par().delt_max)
                              - std::int64_t(par().delt_min) + 1;
    if (nt_ <= 0 || nTime < 1 || nTime > nt_)
    {
        HADRONS_ERROR(Argument, "A2ANewPixW requires 1 <= time-window length <= Nt");
    }
    const std::size_t nLow = par().nLow;
    const std::size_t nsc = std::size_t(Ns)*FImpl::Dimension;
    nRight_ = rightW.size();
    if (nRight_ < nLow)
    {
        HADRONS_ERROR(Size, "A2ANewPixW: nLow exceeds right.size()");
    }

    // Arithmetic only for dimensions/allocation, never in the GEMM loop.
    const auto checkedProduct = [](std::size_t a, std::size_t b)
    {
        if (b && a > std::numeric_limits<std::size_t>::max()/b)
        {
            HADRONS_ERROR(Size, "A2ANewPixW: size multiplication overflow");
        }
        return a*b;
    };
    nActive_ = nRight_;
    nMesCols_ = nRight_;
    if (par().ifDenseRight)
    {
        if ((nRight_ - nLow) % nsc != 0)
        {
            HADRONS_ERROR(Size, "A2ANewPixW: dense right high size is not a multiple of Ns*Nc");
        }
        const std::size_t high = checkedProduct(nRight_ - nLow, nt_);
        if (high > std::numeric_limits<std::size_t>::max() - nLow)
        {
            HADRONS_ERROR(Size, "A2ANewPixW: expanded dimension overflow");
        }
        nMesCols_ = nLow + high;
    }
    else if (par().ifDenseMes)
    {
        if ((nRight_ - nLow) % (std::size_t(nt_)*nsc) != 0)
        {
            HADRONS_ERROR(Size, "A2ANewPixW: expanded right high size is not a multiple of Nt*Ns*Nc");
        }
        nActive_ = nLow + (nRight_ - nLow)/nt_;
        nMesCols_ = nActive_;
    }
    if (nMesCols_ > std::size_t(std::numeric_limits<Eigen::Index>::max()))
    {
        HADRONS_ERROR(Size, "A2ANewPixW: matrix dimension exceeds Eigen::Index");
    }

    nSimd_ = vector_type::Nsimd();
    nSpace_ = std::size_t(grid->_slice_nblock[Tp])*grid->_slice_block[Tp];
    if (nSpace_ == 0 || nSimd_ < 1 ||
        nSpace_*std::size_t(nt_) != std::size_t(grid->oSites()))
    {
        HADRONS_ERROR(Size, "A2ANewPixW: invalid spatial slice layout");
    }
    siteBlock_ = std::min<std::size_t>(par().siteBlock, nSpace_);
    modeBlock_ = std::min<std::size_t>(par().modeBlock, par().size);
    const std::size_t nBlock = (nSpace_ - 1)/siteBlock_ + 1;
    nThread_ = int(std::min<std::size_t>(thread_max(), nBlock));
    qMax_ = checkedProduct(checkedProduct(siteBlock_, nSimd_), nsc);
    rightStride_ = checkedProduct(qMax_, nActive_);
    resultStride_ = checkedProduct(qMax_, modeBlock_);
    mesStride_ = par().ifDenseRight ? checkedProduct(modeBlock_, nActive_) : 0;
    const std::size_t nr = checkedProduct(nThread_, rightStride_);
    const std::size_t ny = checkedProduct(nThread_, resultStride_);
    const std::size_t nm = checkedProduct(nThread_, mesStride_);
    const std::size_t maxElements = std::min(
        std::numeric_limits<std::size_t>::max()/sizeof(ComplexD),
        std::size_t(std::numeric_limits<Eigen::Index>::max()));
    if (nr > maxElements || ny > maxElements - nr || nm > maxElements - nr - ny)
    {
        HADRONS_ERROR(Size, "A2ANewPixW: workspace size exceeds addressable range");
    }
    const std::size_t outBytes = checkedProduct(
        checkedProduct(par().size, grid->oSites()), sizeof(vobj));

    // Do not inspect mes[tmes] here: memory profiling may precede its load.
    // hostVector allocations use Grid's profiler; Eigen::Map owns no memory.
    envCreate(std::vector<FermionField>, getName(), 1, par().size, grid);
    envTmp(Buffer, "rightBuf", 1, nr);
    envTmp(Buffer, "resultBuf", 1, ny);
    if (mesStride_ > 0)
    {
        envTmp(Buffer, "mesBuf", 1, nm);
    }

    LOG(Message) << "A2ANewPixW: output " << par().size << " full 4D fields ("
                 << sizeString(outBytes) << "/rank); active contracted modes "
                 << nActive_ << std::endl;
    LOG(Message) << "A2ANewPixW: siteBlock=" << siteBlock_
                 << ", modeBlock=" << modeBlock_ << ", workers=" << nThread_
                 << ", explicit workspace=" << sizeString((nr + ny + nm)*sizeof(ComplexD))
                 << "/rank (Eigen internal workspace is additional)" << std::endl;
}
// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TA2ANewPixW<FImpl>::execute(void)
{
    auto       &vec = envGet(std::vector<FermionField>, getName());
    const auto &rightW = envGet(std::vector<FermionField>, par().right);
    const auto &meson_f = envGet(EigenDiskVector<ComplexD>, par().mes);
    GridBase  *grid = rightW[0].Grid();
    const int tmes = int((std::int64_t(par().tmes) % nt_ + nt_) % nt_);
    const std::size_t nLow = par().nLow;
    const std::size_t nsc = std::size_t(Ns)*FImpl::Dimension;
    assert(vec.size() == std::size_t(par().size));
    assert(rightW.size() == nRight_);
    if (meson_f.getSize() != unsigned(nt_))
    {
        HADRONS_ERROR(Size, "A2ANewPixW: meson time extent differs from Nt");
    }

    // One collective DiskVector access, outside OpenMP; no matrix copy.
    startTimer("Meson field access");
    const auto &meson = meson_f[tmes];
    stopTimer("Meson field access");
    assert(meson.rows() == par().size);
    if (meson.rows() != par().size || meson.cols() < 0 ||
        std::size_t(meson.cols()) != nMesCols_)
    {
        HADRONS_ERROR(Size, "A2ANewPixW: meson dimensions do not match size/right/dense flags");
    }

    // r is the ACTIVE contracted index, including low modes.
    // Low: r=0,...,nLow-1 maps to itself on both sides.
    // Dense high: r=nLow+h*nsc+a, hence h=(r-nLow)/nsc, a=(r-nLow)%nsc.
    // Expanded partner: nLow+(h*nt+t)*nsc+a. Use t=tmes for dense Mes,
    // and t=tx for dense Right. No d()/k() functions are called in site loops.
    std::vector<std::size_t> rightIndex(nActive_), mesIndex(nActive_);
    for (std::size_t r = 0; r < nActive_; ++r)
    {
        rightIndex[r] = r;
        mesIndex[r] = r;
        if (par().ifDenseMes && r >= nLow)
        {
            const std::size_t h = (r - nLow)/nsc;
            const std::size_t a = (r - nLow)%nsc;
            rightIndex[r] = nLow + (h*nt_ + tmes)*nsc + a;
        }
    }

    envGetTmp(Buffer, rightBuf);
    envGetTmp(Buffer, resultBuf);
    ComplexD *mesData = nullptr;
    if (mesStride_ > 0)
    {
        envGetTmp(Buffer, mesBuf);
        mesData = mesBuf.data();
    }
    const std::size_t e2 = grid->_slice_block[Tp];
    const std::size_t stride = grid->_slice_stride[Tp];
    const std::size_t nBlock = (nSpace_ - 1)/siteBlock_ + 1;
    Eigen::initParallel();

    startTimer("Output initialization");
    for (auto &field: vec)
    {
        field = Zero();
    }
    stopTimer("Output initialization");

    // Views contain pointers, not copies of the lattice fields. Open/close
    // outside OpenMP. Only the active right vectors need an open view.
    std::vector<LatticeView<vobj>> rightView, outView;
    rightView.reserve(nActive_);
    outView.reserve(vec.size());
    try
    {
        for (std::size_t r = 0; r < nActive_; ++r)
        {
            rightView.push_back(rightW[rightIndex[r]].View(CpuRead));
        }
        for (auto &field: vec)
        {
            outView.push_back(field.View(CpuWrite));
        }
        startTimer("Pi x W");
        for (std::int64_t delta = par().delt_min; delta <= par().delt_max; ++delta)
        {
            const int tx = int(((std::int64_t(tmes) + delta) % nt_ + nt_) % nt_);
            if (par().ifDenseRight)
            {
                for (std::size_t r = nLow; r < nActive_; ++r)
                {
                    const std::size_t h = (r - nLow)/nsc;
                    const std::size_t a = (r - nLow)%nsc;
                    mesIndex[r] = nLow + (h*nt_ + tx)*nsc + a;
                }
            }
            LOG(Message) << "A2ANewPixW: tmes=" << tmes << ", delta=" << delta
                         << ", t_op=" << tx << std::endl;
            const std::size_t so = std::size_t(tx)*grid->_ostride[Tp];

            // Same spatial site belongs to one worker, including all lanes.
            // The explicit team size cannot exceed the setup workspace count.
            DO_PRAGMA(omp parallel for schedule(static) num_threads(nThread_))
            for (std::uint64_t ib = 0; ib < nBlock; ++ib)
            {
                const std::size_t tid = thread_num();
                const std::size_t site0 = ib*siteBlock_;
                const std::size_t nb = std::min(siteBlock_, nSpace_ - site0);
                const std::size_t nq = nb*nSimd_*nsc;
                Eigen::Map<BufferMatrix> x(rightBuf.data() + tid*rightStride_, qMax_, nActive_);
                Eigen::Map<BufferMatrix> y(resultBuf.data() + tid*resultStride_, qMax_, modeBlock_);

                // Pack each active right vector once per spatial block.
                for (std::size_t r = 0; r < nActive_; ++r)
                for (std::size_t b = 0; b < nb; ++b)
                {
                    const std::size_t sp = site0 + b;
                    const std::size_t ss = so + (sp/e2)*stride + sp%e2;
                    const auto &rv = rightView[r][ss];
                    for (int lane = 0; lane < nSimd_; ++lane)
                    {
                        const sobj s = extractLane(lane, rv);
                        const std::size_t q0 = (b*nSimd_ + lane)*nsc;
                        for (int spin = 0; spin < Ns; ++spin)
                        for (int colour = 0; colour < FImpl::Dimension; ++colour)
                        {
                            x(q0 + spin*FImpl::Dimension + colour, r) = s()(spin)(colour);
                        }
                    }
                }
                for (std::size_t i0 = 0; i0 < vec.size(); i0 += modeBlock_)
                {
                    const std::size_t ni = std::min(modeBlock_, vec.size() - i0);
                    if (par().ifDenseRight)
                    {
                        Eigen::Map<Matrix> coeff(mesData + tid*mesStride_, modeBlock_, nActive_);
                        for (std::size_t i = 0; i < ni; ++i)
                        for (std::size_t r = 0; r < nActive_; ++r)
                        {
                            coeff(i, r) = meson(i0 + i, mesIndex[r]);
                        }
                        y.topLeftCorner(nq, ni).noalias() =
                            x.topRows(nq)*coeff.topRows(ni).adjoint();
                    }
                    else
                    {
                        y.topLeftCorner(nq, ni).noalias() =
                            x.topRows(nq)*meson.block(i0, 0, ni, nActive_).adjoint();
                    }

                    // One complete SIMD object is assigned at a time.
                    for (std::size_t i = 0; i < ni; ++i)
                    for (std::size_t b = 0; b < nb; ++b)
                    {
                        const std::size_t sp = site0 + b;
                        const std::size_t ss = so + (sp/e2)*stride + sp%e2;
                        vobj value;
                        for (int lane = 0; lane < nSimd_; ++lane)
                        {
                            sobj s;
                            const std::size_t q0 = (b*nSimd_ + lane)*nsc;
                            for (int spin = 0; spin < Ns; ++spin)
                            for (int colour = 0; colour < FImpl::Dimension; ++colour)
                            {
                                s()(spin)(colour) = y(q0 + spin*FImpl::Dimension + colour, i);
                            }
                            insertLane(lane, value, s);
                        }
                        outView[i0 + i][ss] = value;
                    }
                }
            }
        }
        stopTimer("Pi x W");
    }
    catch (...)
    {
        for (auto &view: outView)   view.ViewClose();
        for (auto &view: rightView) view.ViewClose();
        throw;
    }
    for (auto &view: outView)   view.ViewClose();
    for (auto &view: rightView) view.ViewClose();
}
END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE
#endif // Hadrons_MContraction_A2ANewPixW_hpp_
