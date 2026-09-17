/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: Hadrons/Modules/MIO/LoadPropagatorTimeWindow.hpp

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
#ifndef Hadrons_MIO_LoadPropagatorTimeWindow_hpp_
#define Hadrons_MIO_LoadPropagatorTimeWindow_hpp_
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

/* Based on MIO/SavePropagatorTimeWindow.hpp and MIO/LoadField.hpp.
 * Output: exactly std::vector<SpinColourMatrix_v>, in the original window order.
 * setup() uses XML nTime only: NO file access and NO trajectory lookup.
 * execute() verifies the stored header before reading one same-precision field
 * from filestem.traj.bin, then copies short-Grid slices u=0,...,nTime-1.
 * nTime is the STORED window length, not a prefix/subwindow to select.
 * Original tmin/tmax/meson times are deliberately not stored or reconstructed.
 * CPU, Environment's standard full 4D layout, MPI_t=SIMD_t=1.
 * Uses the saver cache key/lifetime convention; all ranks call the reader.
 */

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                      Load a propagator time-window vector                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MIO)
class LoadPropagatorTimeWindowPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(LoadPropagatorTimeWindowPar,
                                    std::string, filestem,
                                    int,         nTime);
};

template <typename FImpl>
class TLoadPropagatorTimeWindow: public Module<LoadPropagatorTimeWindowPar>
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
    TLoadPropagatorTimeWindow(const std::string name);
    // destructor
    virtual ~TLoadPropagatorTimeWindow(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
    std::size_t nTime_{0}, nSpace_{0}, nElements_{0};
    std::string gridName_;
};
MODULE_REGISTER_TMP(LoadPropagatorTimeWindow,
                    TLoadPropagatorTimeWindow<FIMPL>, MIO);

/******************************************************************************
 *                   TLoadPropagatorTimeWindow implementation                 *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TLoadPropagatorTimeWindow<FImpl>::
TLoadPropagatorTimeWindow(const std::string name)
: Module<LoadPropagatorTimeWindowPar>(name)
{}
// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TLoadPropagatorTimeWindow<FImpl>::getInput(void)
{
    std::vector<std::string> in = {};

    return in;
}

template <typename FImpl>
std::vector<std::string> TLoadPropagatorTimeWindow<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};

    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TLoadPropagatorTimeWindow<FImpl>::setup(void)
{
    if (par().filestem.empty() || par().nTime <= 0)
    {
        HADRONS_ERROR(Argument, "LoadPropagatorTimeWindow: empty filestem or non-positive nTime");
    }
    GridCartesian *grid = envGetGrid(PropagatorField);
    if (grid->_ndimension != 4 || grid->_isCheckerBoarded)
    {
        HADRONS_ERROR(Argument, "LoadPropagatorTimeWindow requires a full 4D layout");
    }
    if (grid->_processors[Tp] != 1 || grid->_simd_layout[Tp] != 1)
    {
        HADRONS_ERROR(Argument, "LoadPropagatorTimeWindow requires MPI_t = SIMD_t = 1");
    }
    const int nt = grid->_gdimensions[Tp];
    if (nt <= 0 || grid->_slice_nblock[Tp] <= 0 || grid->_slice_block[Tp] <= 0)
    {
        HADRONS_ERROR(Size, "LoadPropagatorTimeWindow: invalid spatial slice layout");
    }
    const std::size_t e1 = grid->_slice_nblock[Tp];
    const std::size_t e2 = grid->_slice_block[Tp];
    if (e1 > std::numeric_limits<std::size_t>::max()/e2)
    {
        HADRONS_ERROR(Size, "LoadPropagatorTimeWindow: spatial size overflow");
    }
    nSpace_ = e1*e2;
    if (nSpace_ > std::numeric_limits<std::size_t>::max()/std::size_t(nt) ||
        nSpace_*std::size_t(nt) != std::size_t(grid->oSites()))
    {
        HADRONS_ERROR(Size, "LoadPropagatorTimeWindow: inconsistent standard Grid layout");
    }
    if (par().nTime > nt)
    {
        HADRONS_ERROR(Argument, "LoadPropagatorTimeWindow: nTime exceeds the full Grid Nt");
    }
    nTime_ = static_cast<std::size_t>(par().nTime);
    // nTime_ <= nt, so the preceding check also protects this product.
    nElements_ = nTime_*nSpace_;
    if (nElements_ > std::numeric_limits<std::size_t>::max()/sizeof(SpinColourMatrix_v) ||
        nElements_ > std::vector<SpinColourMatrix_v>().max_size())
    {
        HADRONS_ERROR(Size, "LoadPropagatorTimeWindow: output size exceeds addressable range");
    }
    // All ranks must construct the same short global Grid.
    std::uint64_t bossTime = nTime_;
    grid->Broadcast(grid->BossRank(), &bossTime, sizeof(bossTime));
    std::uint32_t mismatch = (bossTime != nTime_);
    grid->GlobalSum(mismatch);
    if (mismatch != 0)
    {
        HADRONS_ERROR(Size, "LoadPropagatorTimeWindow: window lengths differ between ranks");
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
        HADRONS_ERROR(Size, "LoadPropagatorTimeWindow: incompatible cached window Grid");
    }
    // Standard-vector output matches A2ANewVxW/old MT contractions exactly.
    // Its allocation need not appear in Grid's automatic MemoryProfiler.
    envCreate(std::vector<SpinColourMatrix_v>, getName(), 1, nElements_, Zero());
    // Field has a raw Grid pointer. Keep the Grid alive until this temporary
    // is destroyed, including during memory-profile cleanup.
    envTmp(PropagatorField, "ioBuf", 1, &windowGrid);
    env().addObjectDependency(env().getObjectAddress(gridName_),
                             env().getObjectAddress(getName() + "_tmp_ioBuf"));
    LOG(Message) << "LoadPropagatorTimeWindow: " << nTime_
                 << " time slices, " << nSpace_ << " local spatial SIMD objects/time; "
                 << "output " << sizeString(nElements_*sizeof(SpinColourMatrix_v))
                 << "/rank plus an equally-sized I/O buffer"
                 << " (Grid I/O internal buffers are additional)" << std::endl;
}
// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TLoadPropagatorTimeWindow<FImpl>::execute(void)
{
    auto &prop = envGet(std::vector<SpinColourMatrix_v>, getName());
    auto &windowGrid = envGet(GridCartesian, gridName_);
    envGetTmp(PropagatorField, ioBuf);
    GridCartesian *grid = envGetGrid(PropagatorField);
    assert(prop.size() == nElements_);
    if (prop.size() != nElements_ || ioBuf.Grid() != &windowGrid ||
        std::size_t(windowGrid.oSites()) != nElements_ || par().nTime <= 0 ||
        static_cast<std::size_t>(par().nTime) != nTime_)
    {
        HADRONS_ERROR(Size, "LoadPropagatorTimeWindow: output, buffer or nTime changed after setup");
    }

    // Trajectory is used ONLY at execution, never during memory profiling.
    const std::string filename = resultFilename(par().filestem, "bin");
    LOG(Message) << "Loading time-window vector '" << getName() << "' from '"
                 << filename << "' (" << nTime_ << " slices)" << std::endl;
    PropagatorTimeWindowMetadata header, record;
    startTimer("Read time window");
    try
    {
        // Collective field read: do NOT enclose this in IsBoss().
        FieldReader<PropagatorField> reader(&windowGrid);
        reader.open(filename);
        std::string xmlString;
        reader.readHeader(xmlString);
        if (xmlString.empty())
        {
            HADRONS_ERROR(Io, "LoadPropagatorTimeWindow: empty file header in '" + filename + "'");
        }
        XmlReader xmlReader(xmlString, true, "hadronsPropagatorTimeWindow");
        if (!xmlReader.push("metadata"))
        {
            HADRONS_ERROR(Io, "LoadPropagatorTimeWindow: missing time-window metadata in '" + filename + "'");
        }
        xmlReader.pop();
        read(xmlReader, "metadata", header);

        // Reject a wrong window BEFORE reading the binary field. Never truncate,
        // pad, or resize according to the file after setup() allocated the output.
        if (header.format != "HadronsPropagatorTimeWindow" || header.version != 1)
        {
            HADRONS_ERROR(Io, "LoadPropagatorTimeWindow: unsupported format/version in '" + filename + "'");
        }
        if (header.fullDimensions != grid->_gdimensions.toVector())
        {
            HADRONS_ERROR(Size, "LoadPropagatorTimeWindow: full lattice dimensions differ in '" + filename + "'");
        }
        if (header.nTime != nTime_)
        {
            HADRONS_ERROR(Size, "LoadPropagatorTimeWindow: nTime mismatch in '" + filename +
                          "': XML requests " + std::to_string(nTime_) +
                          " slices, file contains " + std::to_string(header.nTime));
        }
        if (header.trajectory != vm().getTrajectory())
        {
            HADRONS_ERROR(Io, "LoadPropagatorTimeWindow: trajectory mismatch in '" + filename +
                          "': expected " + std::to_string(vm().getTrajectory()) +
                          ", file contains " + std::to_string(header.trajectory));
        }

        // Same-precision SciDAC read, including the existing Grid checksums.
        reader.readField(ioBuf, record);
        if (record.format != header.format || record.version != header.version ||
            record.fullDimensions != header.fullDimensions || record.nTime != header.nTime ||
            record.trajectory != header.trajectory || record.sourceObject != header.sourceObject)
        {
            HADRONS_ERROR(Io, "LoadPropagatorTimeWindow: header/record metadata disagree in '" + filename + "'");
        }
        reader.close();
    }
    catch (...)
    {
        stopTimer("Read time window");
        throw;
    }
    stopTimer("Read time window");
    // A different output module name is allowed; sourceObject is provenance only.
    LOG(Message) << "Saved source object: '" << header.sourceObject
                 << "'; original times are not stored" << std::endl;

    const std::size_t e2 = windowGrid._slice_block[Tp];
    const std::size_t stride = windowGrid._slice_stride[Tp];
    startTimer("Unpack time window");
    {
        autoView(buf_v, ioBuf, CpuRead);
        for (std::size_t u = 0; u < nTime_; ++u)
        {
            const std::size_t so = u*windowGrid._ostride[Tp];
            thread_for(p, nSpace_, {
                const std::size_t ss = so + (p/e2)*stride + p%e2;
                prop[u*nSpace_ + p] = buf_v[ss];
            });
        }
    }
    stopTimer("Unpack time window");
}

END_MODULE_NAMESPACE
END_HADRONS_NAMESPACE
#endif // Hadrons_MIO_LoadPropagatorTimeWindow_hpp_
