/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: Hadrons/Modules/MIO/SavePropagatorTimeWindow.hpp

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
#ifndef Hadrons_MIO_SavePropagatorTimeWindow_hpp_
#define Hadrons_MIO_SavePropagatorTimeWindow_hpp_
#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/FieldIo.hpp>
#include <Hadrons/TimeWindowIo.hpp>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <sstream>
#include <type_traits>

/* Based on MIO/SaveSpinColorDiagonalNoise.hpp and MIO/SaveField.hpp.
 * Input: std::vector<SpinColourMatrix_v>, NOT a PropagatorField.
 * Input order: prop[u*nSpace+p], as in A2ANewVxW/PropagatorToTimeWindow.
 * nTime is inferred from prop.size()/nSpace. No time origin is inferred.
 * A short I/O Grid stores those slices at u=0,...,nTime-1, unchanged.
 * Exactly one same-precision SciDAC field is written to filestem.traj.bin.
 * Requires the Environment's standard full 4D layout and MPI_t=SIMD_t=1.
 * All ranks call the writer. No numerical transformation or A2A operation.
 */

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                      Save a propagator time-window vector                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MIO)
class SavePropagatorTimeWindowPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(SavePropagatorTimeWindowPar,
                                    std::string, prop,
                                    std::string, filestem);
};

template <typename FImpl>
class TSavePropagatorTimeWindow: public Module<SavePropagatorTimeWindowPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    typedef typename FImpl::SiteSpinor::vector_type vector_type;
    typedef iSpinColourMatrix<vector_type> SpinColourMatrix_v;
    static_assert(std::is_same<typename PropagatorField::vector_object,
                               SpinColourMatrix_v>::value,
                  "Window matrix and PropagatorField site types must agree");
public:
    // constructor
    TSavePropagatorTimeWindow(const std::string name);
    // destructor
    virtual ~TSavePropagatorTimeWindow(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    virtual std::vector<std::string> getOutputFiles(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
    std::size_t nTime_{0}, nSpace_{0}, nElements_{0};
    std::string gridName_;
};
MODULE_REGISTER_TMP(SavePropagatorTimeWindow,
                    TSavePropagatorTimeWindow<FIMPL>, MIO);

/******************************************************************************
 *                   TSavePropagatorTimeWindow implementation                 *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TSavePropagatorTimeWindow<FImpl>::
TSavePropagatorTimeWindow(const std::string name)
: Module<SavePropagatorTimeWindowPar>(name)
{}
// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TSavePropagatorTimeWindow<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().prop};

    return in;
}

template <typename FImpl>
std::vector<std::string> TSavePropagatorTimeWindow<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {};

    return out;
}

template <typename FImpl>
std::vector<std::string> TSavePropagatorTimeWindow<FImpl>::getOutputFiles(void)
{
    std::vector<std::string> out = {resultFilename(par().filestem, "bin")};

    return out;
}
// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TSavePropagatorTimeWindow<FImpl>::setup(void)
{
    if (par().prop.empty() || par().filestem.empty())
    {
        HADRONS_ERROR(Argument, "SavePropagatorTimeWindow: empty prop or filestem");
    }
    const auto &prop = envGet(std::vector<SpinColourMatrix_v>, par().prop);
    GridCartesian *grid = envGetGrid(PropagatorField);
    if (grid->_ndimension != 4 || grid->_isCheckerBoarded ||
        env().getObjectLs(par().prop) != 1)
    {
        HADRONS_ERROR(Argument, "SavePropagatorTimeWindow requires a full 4D layout");
    }
    if (grid->_processors[Tp] != 1 || grid->_simd_layout[Tp] != 1)
    {
        HADRONS_ERROR(Argument, "SavePropagatorTimeWindow requires MPI_t = SIMD_t = 1");
    }
    const int nt = grid->_gdimensions[Tp];
    if (nt <= 0 || grid->_slice_nblock[Tp] <= 0 || grid->_slice_block[Tp] <= 0)
    {
        HADRONS_ERROR(Size, "SavePropagatorTimeWindow: invalid spatial slice layout");
    }
    const std::size_t e1 = grid->_slice_nblock[Tp];
    const std::size_t e2 = grid->_slice_block[Tp];
    if (e1 > std::numeric_limits<std::size_t>::max()/e2)
    {
        HADRONS_ERROR(Size, "SavePropagatorTimeWindow: spatial size overflow");
    }
    nSpace_ = e1*e2;
    if (nSpace_ > std::numeric_limits<std::size_t>::max()/std::size_t(nt) ||
        nSpace_*std::size_t(nt) != std::size_t(grid->oSites()))
    {
        HADRONS_ERROR(Size, "SavePropagatorTimeWindow: inconsistent standard Grid layout");
    }
    if (prop.empty() || prop.size() % nSpace_ != 0)
    {
        HADRONS_ERROR(Size, "SavePropagatorTimeWindow: prop size is not a positive number of spatial slices");
    }
    nElements_ = prop.size();
    nTime_ = nElements_/nSpace_;
    if (nTime_ > std::size_t(nt) ||
        nElements_ > std::numeric_limits<std::size_t>::max()/sizeof(SpinColourMatrix_v))
    {
        HADRONS_ERROR(Size, "SavePropagatorTimeWindow: window exceeds Nt or byte size overflows");
    }
    // All ranks must construct the same short global Grid.
    std::uint64_t bossTime = nTime_;
    grid->Broadcast(grid->BossRank(), &bossTime, sizeof(bossTime));
    std::uint32_t mismatch = (bossTime != nTime_);
    grid->GlobalSum(mismatch);
    if (mismatch != 0)
    {
        HADRONS_ERROR(Size, "SavePropagatorTimeWindow: window lengths differ between ranks");
    }

    Coordinate dimensions = grid->_gdimensions;
    dimensions[Tp] = static_cast<int>(nTime_);
    // One cached Grid per geometry/window length, NOT per module/trajectory.
    // No cached raw pointer is kept across profiling's env().freeAll().
    std::ostringstream key;
    key << "__PropagatorTimeWindowGrid";
    for (int mu = 0; mu < 4; ++mu)
    {
        key << "_" << grid->_gdimensions[mu]
            << "_" << grid->_simd_layout[mu] << "_" << grid->_processors[mu];
    }
    key << "_nt" << nTime_;
    gridName_ = key.str();
    if (!env().hasCreatedObject(gridName_))
    {
        envCache(GridCartesian, gridName_, 1,
                 dimensions, grid->_simd_layout, grid->_processors, *grid);
    }
    auto &windowGrid = envGet(GridCartesian, gridName_);
    if (!(windowGrid._gdimensions == dimensions) ||
        !(windowGrid._simd_layout == grid->_simd_layout) ||
        !(windowGrid._processors == grid->_processors) ||
        !(windowGrid._processor_coor == grid->_processor_coor) ||
        std::size_t(windowGrid.oSites()) != nElements_)
    {
        HADRONS_ERROR(Size, "SavePropagatorTimeWindow: incompatible cached window Grid");
    }
    // Field has a raw Grid pointer. Keep the Grid alive until this temporary
    // is destroyed, including during memory-profile cleanup.
    envTmp(PropagatorField, "ioBuf", 1, &windowGrid);
    env().addObjectDependency(env().getObjectAddress(gridName_),
                             env().getObjectAddress(getName() + "_tmp_ioBuf"));
    LOG(Message) << "SavePropagatorTimeWindow: inferred " << nTime_
                 << " time slices, " << nSpace_ << " local spatial SIMD objects/time; "
                 << "I/O buffer " << sizeString(nElements_*sizeof(SpinColourMatrix_v))
                 << "/rank (Grid I/O internal buffers are additional)" << std::endl;
}
// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TSavePropagatorTimeWindow<FImpl>::execute(void)
{
    const auto &prop = envGet(std::vector<SpinColourMatrix_v>, par().prop);
    auto &windowGrid = envGet(GridCartesian, gridName_);
    envGetTmp(PropagatorField, ioBuf);
    GridCartesian *grid = envGetGrid(PropagatorField);
    assert(prop.size() == nElements_);
    if (prop.size() != nElements_ || ioBuf.Grid() != &windowGrid ||
        std::size_t(windowGrid.oSites()) != nElements_)
    {
        HADRONS_ERROR(Size, "SavePropagatorTimeWindow: input or buffer changed after setup");
    }
    const std::size_t e2 = windowGrid._slice_block[Tp];
    const std::size_t stride = windowGrid._slice_stride[Tp];
    startTimer("Pack time window");
    {
        autoView(buf_v, ioBuf, CpuWrite);
        for (std::size_t u = 0; u < nTime_; ++u)
        {
            const std::size_t so = u*windowGrid._ostride[Tp];
            thread_for(p, nSpace_, {
                const std::size_t ss = so + (p/e2)*stride + p%e2;
                buf_v[ss] = prop[u*nSpace_ + p];
            });
        }
    }
    stopTimer("Pack time window");

    PropagatorTimeWindowMetadata md;
    md.format         = "HadronsPropagatorTimeWindow";
    md.version        = 1;
    md.fullDimensions = grid->_gdimensions.toVector();
    md.nTime          = static_cast<unsigned int>(nTime_);
    md.trajectory     = vm().getTrajectory();
    md.sourceObject   = par().prop;
    XmlWriter xmlWriter("", "hadronsPropagatorTimeWindow");
    write(xmlWriter, "metadata", md);
    const std::string filename = resultFilename(par().filestem, "bin");
    LOG(Message) << "Saving time-window vector '" << par().prop << "' to '"
                 << filename << "' (" << nTime_ << " slices, original times not stored)"
                 << std::endl;
    // Collective: do NOT put this inside IsBoss(). Existing FieldWriter handles
    // metadata ownership, coordinate-order conversion, precision and checksums.
    // Standard FieldWriter semantics: an existing filename is overwritten.
    startTimer("Write time window");
    {
        FieldWriter<PropagatorField> writer(&windowGrid);
        writer.open(filename);
        writer.writeHeader(xmlWriter, "hadronsPropagatorTimeWindow");
        writer.writeField(ioBuf, md);
        writer.close();
    }
    stopTimer("Write time window");
}

END_MODULE_NAMESPACE
END_HADRONS_NAMESPACE
#endif // Hadrons_MIO_SavePropagatorTimeWindow_hpp_
