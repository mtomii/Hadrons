/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: Hadrons/Modules/MContraction/A2ANewVxPi.hpp

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
#ifndef Hadrons_MContraction_A2ANewVxPi_hpp_
#define Hadrons_MContraction_A2ANewVxPi_hpp_
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

/* Derived from A2ANewPixW.hpp, using the convention of A2AVxPi.hpp.
 *
 * out_j(x) = sum_r left[leftIndex[r]](x) * M(tmes)(mesIndex[r], j).
 * Neither factor is conjugated. Close with sum_j outerProduct(out_j, right_j)
 * and the appropriate mode/time mapping to obtain L M R^dagger.
 *
 * Only the contracted (row) meson index is mapped. The column index is kept
 * verbatim, including any compact dense-W index. For an elementary meson
 * field that index refers to tmes; a product M1(t1)M2(t2) retains t2 on its
 * right edge, which the downstream closing module must know separately.
 * No hit normalization is added and no index metadata object is created.
 *
 * CPU execution; full 4D fields, MPI_t = SIMD_t = 1. Spatial SIMD lanes
 * are packed into scalar ComplexD matrices. No Grid changes are required.
 */

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                   A2A vector set times a meson field                      *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MContraction)
class A2ANewVxPiPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(A2ANewVxPiPar,
                                    std::string, left,
                                    std::string, mes,
                                    int,         tmes,
                                    int,         delt_max,
                                    int,         delt_min,
                                    bool,        ifDenseLeft,
                                    bool,        ifDenseMes,
                                    int,         nLow,
                                    int,         size,
                                    int,         siteBlock,
                                    int,         modeBlock);
};

template <typename FImpl>
class TA2ANewVxPi: public Module<A2ANewVxPiPar>
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
    TA2ANewVxPi(const std::string name);
    // destructor
    virtual ~TA2ANewVxPi(void) {};
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
    std::size_t nActive_{0}, nLeft_{0}, nMesRows_{0};
    std::size_t leftStride_{0}, resultStride_{0}, mesStride_{0};
};
MODULE_REGISTER_TMP(A2ANewVxPi, TA2ANewVxPi<FIMPL>, MContraction);

/******************************************************************************
 *                       TA2ANewVxPi implementation                           *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TA2ANewVxPi<FImpl>::TA2ANewVxPi(const std::string name)
: Module<A2ANewVxPiPar>(name)
{}
// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TA2ANewVxPi<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().left, par().mes};

    return in;
}

template <typename FImpl>
std::vector<std::string> TA2ANewVxPi<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};

    return out;
}
// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TA2ANewVxPi<FImpl>::setup(void)
{
    if (par().left.empty() || par().mes.empty() || par().nLow < 0 ||
        par().size <= 0 || par().siteBlock <= 0 || par().modeBlock <= 0)
    {
        HADRONS_ERROR(Argument, "A2ANewVxPi: invalid input name, size or block size");
    }
    if (par().ifDenseLeft && par().ifDenseMes)
    {
        HADRONS_ERROR(Argument, "A2ANewVxPi: both contracted indices cannot be dense W");
    }

    const auto &leftV = envGet(std::vector<FermionField>, par().left);
    if (leftV.empty())
    {
        HADRONS_ERROR(Size, "A2ANewVxPi: empty left vector set");
    }
    GridBase *grid = leftV[0].Grid();
    if (grid->_ndimension != 4 || grid->_isCheckerBoarded ||
        env().getObjectLs(par().left) != 1)
    {
        HADRONS_ERROR(Argument, "A2ANewVxPi requires full, non-checkerboard 4D fields");
    }
    if (grid->_processors[Tp] != 1 || grid->_simd_layout[Tp] != 1)
    {
        HADRONS_ERROR(Argument, "A2ANewVxPi requires MPI_t = SIMD_t = 1");
    }
    for (const auto &field: leftV)
    {
        if (field.Grid() != grid)
        {
            HADRONS_ERROR(Size, "A2ANewVxPi: left fields have different Grids");
        }
    }

    nt_ = grid->_gdimensions[Tp];
    const std::int64_t nTime = std::int64_t(par().delt_max)
                              - std::int64_t(par().delt_min) + 1;
    if (nt_ <= 0 || nTime < 1 || nTime > nt_)
    {
        HADRONS_ERROR(Argument, "A2ANewVxPi requires 1 <= time-window length <= Nt");
    }
    const std::size_t nLow = par().nLow;
    const std::size_t nsc = std::size_t(Ns)*FImpl::Dimension;
    nLeft_ = leftV.size();
    if (nLeft_ < nLow)
    {
        HADRONS_ERROR(Size, "A2ANewVxPi: nLow exceeds left.size()");
    }

    // Arithmetic only for dimensions/allocation, never in the GEMM loop.
    const auto checkedProduct = [](std::size_t a, std::size_t b)
    {
        if (b && a > std::numeric_limits<std::size_t>::max()/b)
        {
            HADRONS_ERROR(Size, "A2ANewVxPi: size multiplication overflow");
        }
        return a*b;
    };
    nActive_ = nLeft_;
    nMesRows_ = nLeft_;
    if (par().ifDenseLeft)
    {
        if ((nLeft_ - nLow) % nsc != 0)
        {
            HADRONS_ERROR(Size, "A2ANewVxPi: dense left high size is not a multiple of Ns*Nc");
        }
        const std::size_t high = checkedProduct(nLeft_ - nLow, nt_);
        if (high > std::numeric_limits<std::size_t>::max() - nLow)
        {
            HADRONS_ERROR(Size, "A2ANewVxPi: expanded dimension overflow");
        }
        nMesRows_ = nLow + high;
    }
    else if (par().ifDenseMes)
    {
        if ((nLeft_ - nLow) % (std::size_t(nt_)*nsc) != 0)
        {
            HADRONS_ERROR(Size, "A2ANewVxPi: expanded left high size is not a multiple of Nt*Ns*Nc");
        }
        nActive_ = nLow + (nLeft_ - nLow)/nt_;
        nMesRows_ = nActive_;
    }
    if (nMesRows_ > std::size_t(std::numeric_limits<Eigen::Index>::max()))
    {
        HADRONS_ERROR(Size, "A2ANewVxPi: matrix dimension exceeds Eigen::Index");
    }

    nSimd_ = vector_type::Nsimd();
    nSpace_ = std::size_t(grid->_slice_nblock[Tp])*grid->_slice_block[Tp];
    if (nSpace_ == 0 || nSimd_ < 1 ||
        nSpace_*std::size_t(nt_) != std::size_t(grid->oSites()))
    {
        HADRONS_ERROR(Size, "A2ANewVxPi: invalid spatial slice layout");
    }
    siteBlock_ = std::min<std::size_t>(par().siteBlock, nSpace_);
    modeBlock_ = std::min<std::size_t>(par().modeBlock, par().size);
    const std::size_t nBlock = (nSpace_ - 1)/siteBlock_ + 1;
    nThread_ = int(std::min<std::size_t>(thread_max(), nBlock));
    qMax_ = checkedProduct(checkedProduct(siteBlock_, nSimd_), nsc);
    leftStride_ = checkedProduct(qMax_, nActive_);
    resultStride_ = checkedProduct(qMax_, modeBlock_);
    mesStride_ = par().ifDenseLeft ? checkedProduct(nActive_, modeBlock_) : 0;
    const std::size_t nx = checkedProduct(nThread_, leftStride_);
    const std::size_t ny = checkedProduct(nThread_, resultStride_);
    const std::size_t nm = checkedProduct(nThread_, mesStride_);
    const std::size_t maxElements = std::min(
        std::numeric_limits<std::size_t>::max()/sizeof(ComplexD),
        std::size_t(std::numeric_limits<Eigen::Index>::max()));
    if (nx > maxElements || ny > maxElements - nx || nm > maxElements - nx - ny)
    {
        HADRONS_ERROR(Size, "A2ANewVxPi: workspace size exceeds addressable range");
    }
    const std::size_t outBytes = checkedProduct(
        checkedProduct(par().size, grid->oSites()), sizeof(vobj));

    // Do not inspect mes[tmes] here: memory profiling may precede its load.
    // hostVector allocations use Grid's profiler; Eigen::Map owns no memory.
    envCreate(std::vector<FermionField>, getName(), 1, par().size, grid);
    envTmp(Buffer, "leftBuf", 1, nx);
    envTmp(Buffer, "resultBuf", 1, ny);
    if (mesStride_ > 0)
    {
        envTmp(Buffer, "mesBuf", 1, nm);
    }

    LOG(Message) << "A2ANewVxPi: output " << par().size << " full 4D fields ("
                 << sizeString(outBytes) << "/rank); active contracted modes "
                 << nActive_ << std::endl;
    LOG(Message) << "A2ANewVxPi: siteBlock=" << siteBlock_
                 << ", modeBlock=" << modeBlock_ << ", workers=" << nThread_
                 << ", explicit workspace=" << sizeString((nx + ny + nm)*sizeof(ComplexD))
                 << "/rank (Eigen internal workspace is additional)" << std::endl;
}
// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TA2ANewVxPi<FImpl>::execute(void)
{
    auto       &vec = envGet(std::vector<FermionField>, getName());
    const auto &leftV = envGet(std::vector<FermionField>, par().left);
    const auto &meson_f = envGet(EigenDiskVector<ComplexD>, par().mes);
    GridBase  *grid = leftV[0].Grid();
    const int tmes = int((std::int64_t(par().tmes) % nt_ + nt_) % nt_);
    const std::size_t nLow = par().nLow;
    const std::size_t nsc = std::size_t(Ns)*FImpl::Dimension;
    assert(vec.size() == std::size_t(par().size));
    assert(leftV.size() == nLeft_);
    if (meson_f.getSize() != unsigned(nt_))
    {
        HADRONS_ERROR(Size, "A2ANewVxPi: meson time extent differs from Nt");
    }

    // One collective DiskVector access, outside OpenMP; no matrix copy.
    startTimer("Meson field access");
    const auto &meson = meson_f[tmes];
    stopTimer("Meson field access");
    assert(meson.cols() == par().size);
    if (meson.cols() != par().size || meson.rows() < 0 ||
        std::size_t(meson.rows()) != nMesRows_)
    {
        HADRONS_ERROR(Size, "A2ANewVxPi: meson dimensions do not match size/left/dense flags");
    }

    // r is the ACTIVE contracted index, including low modes.
    // Low: r=0,...,nLow-1 maps to itself on both sides.
    // Dense high: r=nLow+h*nsc+a, hence h=(r-nLow)/nsc, a=(r-nLow)%nsc.
    // Expanded partner: nLow+(h*nt+t)*nsc+a. Use t=tmes for dense Mes,
    // and t=tx for dense Left. Mapping arithmetic is outside spatial loops.
    const bool denseLeft = par().ifDenseLeft;
    const bool denseMes  = par().ifDenseMes;
    std::vector<std::size_t> leftIndex(nActive_), mesIndex(nActive_);
    for (std::size_t r = 0; r < nActive_; ++r)
    {
        leftIndex[r] = r;
        mesIndex[r]  = r;
    }
    if (denseMes)
    {
        for (std::size_t r = nLow; r < nActive_; ++r)
        {
            const std::size_t h = (r - nLow)/nsc;
            const std::size_t a = (r - nLow)%nsc;
            leftIndex[r] = nLow + (h*nt_ + tmes)*nsc + a;
        }
    }

    envGetTmp(Buffer, leftBuf);
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
    // outside OpenMP. Only the active left vectors need an open view.
    std::vector<LatticeView<vobj>> leftView, outView;
    leftView.reserve(nActive_);
    outView.reserve(vec.size());
    try
    {
        for (std::size_t r = 0; r < nActive_; ++r)
        {
            leftView.push_back(leftV[leftIndex[r]].View(CpuRead));
        }
        for (auto &field: vec)
        {
            outView.push_back(field.View(CpuWrite));
        }
        startTimer("V x Pi");
        for (std::int64_t delta = par().delt_min; delta <= par().delt_max; ++delta)
        {
            const int tx = int(((std::int64_t(tmes) + delta) % nt_ + nt_) % nt_);
            if (denseLeft)
            {
                for (std::size_t r = nLow; r < nActive_; ++r)
                {
                    const std::size_t h = (r - nLow)/nsc;
                    const std::size_t a = (r - nLow)%nsc;
                    mesIndex[r] = nLow + (h*nt_ + tx)*nsc + a;
                }
            }
            LOG(Message) << "A2ANewVxPi: tmes=" << tmes << ", delta=" << delta
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
                Eigen::Map<BufferMatrix> x(leftBuf.data() + tid*leftStride_, qMax_, nActive_);
                Eigen::Map<BufferMatrix> y(resultBuf.data() + tid*resultStride_, qMax_, modeBlock_);

                // Pack each active left vector once per spatial block.
                for (std::size_t r = 0; r < nActive_; ++r)
                for (std::size_t b = 0; b < nb; ++b)
                {
                    const std::size_t sp = site0 + b;
                    const std::size_t ss = so + (sp/e2)*stride + sp%e2;
                    const auto &lv = leftView[r][ss];
                    for (int lane = 0; lane < nSimd_; ++lane)
                    {
                        const sobj s = extractLane(lane, lv);
                        const std::size_t q0 = (b*nSimd_ + lane)*nsc;
                        for (int spin = 0; spin < Ns; ++spin)
                        for (int colour = 0; colour < FImpl::Dimension; ++colour)
                        {
                            x(q0 + spin*FImpl::Dimension + colour, r) = s()(spin)(colour);
                        }
                    }
                }
                for (std::size_t j0 = 0; j0 < vec.size(); j0 += modeBlock_)
                {
                    const std::size_t nj = std::min(modeBlock_, vec.size() - j0);
                    // One branch per GEMM, not per site/spin/colour/mode term.
                    // Keeping this choice here avoids duplicating pack/write loops.
                    if (denseLeft)
                    {
                        Eigen::Map<Matrix> coeff(mesData + tid*mesStride_, nActive_, modeBlock_);
                        for (std::size_t r = 0; r < nActive_; ++r)
                        for (std::size_t j = 0; j < nj; ++j)
                        {
                            coeff(r, j) = meson(mesIndex[r], j0 + j);
                        }
                        y.topLeftCorner(nq, nj).noalias() =
                            x.topRows(nq)*coeff.leftCols(nj);
                    }
                    else
                    {
                        y.topLeftCorner(nq, nj).noalias() =
                            x.topRows(nq)*meson.block(0, j0, nActive_, nj);
                    }

                    // One complete SIMD object is assigned at a time.
                    for (std::size_t j = 0; j < nj; ++j)
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
                                s()(spin)(colour) = y(q0 + spin*FImpl::Dimension + colour, j);
                            }
                            insertLane(lane, value, s);
                        }
                        outView[j0 + j][ss] = value;
                    }
                }
            }
        }
        stopTimer("V x Pi");
    }
    catch (...)
    {
        for (auto &view: outView)   view.ViewClose();
        for (auto &view: leftView) view.ViewClose();
        throw;
    }
    for (auto &view: outView)   view.ViewClose();
    for (auto &view: leftView) view.ViewClose();
}
END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE
#endif // Hadrons_MContraction_A2ANewVxPi_hpp_
