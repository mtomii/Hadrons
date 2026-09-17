/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: Hadrons/Modules/MUtilities/PropagatorToTimeWindow.hpp

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
#ifndef Hadrons_MUtilities_PropagatorToTimeWindow_hpp_
#define Hadrons_MUtilities_PropagatorToTimeWindow_hpp_
#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <limits>

/* Derived from MContraction/A2ANewVxW.hpp (window layout) and MIO/WriteProp.hpp
 * (PropagatorField input). Pure copy, without conjugation or normalization.
 *
 * tmin/tmax are INCLUSIVE, UNWRAPPED integer times. Only the input lookup is
 * reduced modulo Nt; endpoints are NOT reduced before computing the length.
 * mat[u*nSpace+p] = prop(x_p, (tmin+u) mod Nt), u=0,...,tmax-tmin.
 * Negative times and times >= Nt are allowed; require 1 <= length <= Nt.
 *
 * Output is exactly std::vector<SpinColourMatrix_v> (standard allocator),
 * matching A2ANewVxW's local spatial SIMD ordering. No Grid/time metadata
 * is embedded in the vector. CPU, full non-checkerboard 4D input, MPI_t=SIMD_t=1.
 * The input is unchanged. No I/O, gamma5, dense-mode mapping, or global sum.
 */

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *              Copy a propagator field into a selected time window            *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MUtilities)
class PropagatorToTimeWindowPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(PropagatorToTimeWindowPar,
                                    std::string, prop,
                                    int,         tmin,
                                    int,         tmax);
};

template <typename FImpl>
class TPropagatorToTimeWindow: public Module<PropagatorToTimeWindowPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    typedef typename FImpl::SiteSpinor::vector_type vector_type;
    typedef iSpinColourMatrix<vector_type> SpinColourMatrix_v;
public:
    // constructor
    TPropagatorToTimeWindow(const std::string name);
    // destructor
    virtual ~TPropagatorToTimeWindow(void) {};
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
};
MODULE_REGISTER_TMP(PropagatorToTimeWindow, TPropagatorToTimeWindow<FIMPL>, MUtilities);

/******************************************************************************
 *                    TPropagatorToTimeWindow implementation                   *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TPropagatorToTimeWindow<FImpl>::TPropagatorToTimeWindow(const std::string name)
: Module<PropagatorToTimeWindowPar>(name)
{}
// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TPropagatorToTimeWindow<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().prop};

    return in;
}

template <typename FImpl>
std::vector<std::string> TPropagatorToTimeWindow<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};

    return out;
}
// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TPropagatorToTimeWindow<FImpl>::setup(void)
{
    if (par().prop.empty())
    {
        HADRONS_ERROR(Argument, "PropagatorToTimeWindow: empty input name");
    }
    const auto &prop = envGet(PropagatorField, par().prop);
    GridBase *grid = prop.Grid();
    if (grid->_ndimension != 4 || grid->_isCheckerBoarded ||
        env().getObjectLs(par().prop) != 1)
    {
        HADRONS_ERROR(Argument, "PropagatorToTimeWindow requires full, non-checkerboard 4D fields");
    }
    if (grid->_processors[Tp] != 1 || grid->_simd_layout[Tp] != 1)
    {
        HADRONS_ERROR(Argument, "PropagatorToTimeWindow requires MPI_t = SIMD_t = 1");
    }

    nt_ = grid->_gdimensions[Tp];
    // Cast BEFORE subtraction: INT_MIN/INT_MAX are valid endpoint values.
    const std::int64_t nTime = std::int64_t(par().tmax)
                              - std::int64_t(par().tmin) + 1;
    if (nt_ <= 0 || nTime < 1 || nTime > nt_)
    {
        HADRONS_ERROR(Argument, "PropagatorToTimeWindow requires 1 <= tmax-tmin+1 <= Nt");
    }
    nTime_ = static_cast<std::size_t>(nTime);

    if (grid->_slice_nblock[Tp] <= 0 || grid->_slice_block[Tp] <= 0 ||
        grid->_slice_stride[Tp] <= 0 || grid->_ostride[Tp] <= 0)
    {
        HADRONS_ERROR(Size, "PropagatorToTimeWindow: invalid spatial slice layout");
    }
    const std::size_t e1 = grid->_slice_nblock[Tp];
    const std::size_t e2 = grid->_slice_block[Tp];
    if (e1 > std::numeric_limits<std::size_t>::max()/e2)
    {
        HADRONS_ERROR(Size, "PropagatorToTimeWindow: spatial size overflow");
    }
    nSpace_ = e1*e2;
    if (nSpace_ > std::numeric_limits<std::size_t>::max()/std::size_t(nt_) ||
        nSpace_*std::size_t(nt_) != std::size_t(grid->oSites()))
    {
        HADRONS_ERROR(Size, "PropagatorToTimeWindow: spatial slice size disagrees with input Grid");
    }
    // nTime_<=nt_, so the preceding check also protects this product.
    nOutput_ = nTime_*nSpace_;
    if (nOutput_ > std::numeric_limits<std::size_t>::max()/sizeof(SpinColourMatrix_v) ||
        nOutput_ > std::vector<SpinColourMatrix_v>().max_size())
    {
        HADRONS_ERROR(Size, "PropagatorToTimeWindow: output size exceeds addressable range");
    }

    // Keep the exact vector type used by A2ANewVxW and the old MT modules.
    // Grid's automatic MemoryProfiler need not account for this allocation.
    envCreate(std::vector<SpinColourMatrix_v>, getName(), 1, nOutput_, Zero());
    LOG(Message) << "PropagatorToTimeWindow: [" << par().tmin << ", " << par().tmax
                 << "] (inclusive, unwrapped), " << nTime_ << " time slices, "
                 << nSpace_ << " local spatial SIMD objects/time; output "
                 << sizeString(nOutput_*sizeof(SpinColourMatrix_v)) << "/rank"
                 << " (standard-vector allocation; reported explicitly)" << std::endl;
}
// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TPropagatorToTimeWindow<FImpl>::execute(void)
{
    const auto &prop = envGet(PropagatorField, par().prop);
    auto       &mat = envGet(std::vector<SpinColourMatrix_v>, getName());
    assert(mat.size() == nOutput_);
    if (mat.size() != nOutput_)
    {
        HADRONS_ERROR(Size, "PropagatorToTimeWindow: output size changed after setup");
    }
    GridBase *grid = prop.Grid();
    const std::size_t e2     = grid->_slice_block[Tp];
    const std::size_t stride = grid->_slice_stride[Tp];

    startTimer("Time window extraction");
    {
        // One CPU view for the entire window; closes automatically on scope exit.
        autoView(prop_v, prop, CpuRead);

        for (std::size_t u = 0; u < nTime_; ++u)
        {
            const std::int64_t t = std::int64_t(par().tmin) + std::int64_t(u);
            const int tx = int((t % nt_ + nt_) % nt_);
            const std::size_t so = std::size_t(tx)*grid->_ostride[Tp];
            LOG(Message) << "PropagatorToTimeWindow: time=" << t << ", t_op=" << tx
                         << ", window index=" << u << std::endl;

            thread_for(p, nSpace_, {
                const std::size_t ss = so + (p/e2)*stride + p%e2;
                mat[u*nSpace_ + p] = prop_v[ss];
            });
        }
    }
    stopTimer("Time window extraction");
}
END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE
#endif // Hadrons_MUtilities_PropagatorToTimeWindow_hpp_
