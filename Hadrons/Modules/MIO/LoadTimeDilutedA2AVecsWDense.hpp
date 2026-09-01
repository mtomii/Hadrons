/*
 * LoadTimeDilutedA2AVecsWDense.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2026
 *
 * Author: Jonas Hildebrand <jonas.hildebrand@uconn.edu>
 *
 * Hadrons is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * Hadrons is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Hadrons.  If not, see <http://www.gnu.org/licenses/>.
 *
 * See the full license in the file "LICENSE" in the top level distribution
 * directory.
 */

/*  END LEGAL */
#ifndef Hadrons_MIO_LoadTimeDilutedA2AVecsWDense_hpp_
#define Hadrons_MIO_LoadTimeDilutedA2AVecsWDense_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/A2AVectors.hpp>
#include <Hadrons/DilutedNoise.hpp>
#include <Grid/qcd/utils/A2Autils.h>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *  Assemble a full A2A W vector array in the dense/combined high-mode        *
 *  representation: an optional block of low modes (nLow > 0) read binned     *
 *  from disk, followed by Nsc = 12 dense fermion fields per hit expanded     *
 *  from the raw noise (one field per spin-color index, all timeslice         *
 *  sources superimposed), instead of Nsc*nt time-diluted fields per hit.     *
 *  Since the meson-field GEMM is batched over timeslices and never sums      *
 *  across them, the dense fields pass through it unchanged and exactly:      *
 *  slot (hit, sc) restricted to timeslice t reproduces the fully expanded    *
 *  W vector (hit*nt + t)*Nsc + sc bit for bit. Downstream consumers must     *
 *  respect that pairing (and apply nothing with temporal extent to these     *
 *  fields).                                                                  *
 *                                                                            *
 *  Output layout: slot = nLow + hit*Nsc + sc (sc fastest); hit order is      *
 *  the field order of the input noise object, i.e. the fileStems order of    *
 *  MIO::LoadTimeDilutedSpinColorDiagonalNoise. This is the normative         *
 *  convention for all contractor-side index math.                            *
 *                                                                            *
 *  No 1/nHit normalization is applied here: by convention the hit-average    *
 *  factor lives entirely on the V side (see MIO::LoadBinnedA2AVecsV), and    *
 *  W stays the raw unit-modulus noise. The fully expanded fallback for       *
 *  this module is MIO::LoadBinnedA2AVecsW; the file-based V/W loader is      *
 *  MIO::LoadCombinedA2AVecsV.                                                *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MIO)

class LoadTimeDilutedA2AVecsWDensePar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(LoadTimeDilutedA2AVecsWDensePar,
                                    std::string,  lowFilestem,
                                    unsigned int, nLow,
                                    std::string,  noise);
};

template <typename FImpl, int lowBinSize>
class TLoadTimeDilutedA2AVecsWDense: public Module<LoadTimeDilutedA2AVecsWDensePar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    typedef iVector<typename FImpl::SiteSpinor, lowBinSize> LowBinnedSpinor;
public:
    // constructor
    TLoadTimeDilutedA2AVecsWDense(const std::string name);
    // destructor
    virtual ~TLoadTimeDilutedA2AVecsWDense(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(LoadTimeDilutedA2AVecsWDense200, ARG(TLoadTimeDilutedA2AVecsWDense<FIMPL, 200>), MIO);
MODULE_REGISTER_TMP(LoadTimeDilutedA2AVecsWDense100, ARG(TLoadTimeDilutedA2AVecsWDense<FIMPL, 100>), MIO);
MODULE_REGISTER_TMP(LoadZTimeDilutedA2AVecsWDense200, ARG(TLoadTimeDilutedA2AVecsWDense<ZFIMPL, 200>), MIO);
MODULE_REGISTER_TMP(LoadZTimeDilutedA2AVecsWDense100, ARG(TLoadTimeDilutedA2AVecsWDense<ZFIMPL, 100>), MIO);

/******************************************************************************
 *                 TLoadTimeDilutedA2AVecsWDense implementation               *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl, int lowBinSize>
TLoadTimeDilutedA2AVecsWDense<FImpl, lowBinSize>::TLoadTimeDilutedA2AVecsWDense(const std::string name)
: Module<LoadTimeDilutedA2AVecsWDensePar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl, int lowBinSize>
std::vector<std::string> TLoadTimeDilutedA2AVecsWDense<FImpl, lowBinSize>::getInput(void)
{
    std::vector<std::string> in = {par().noise};

    return in;
}

template <typename FImpl, int lowBinSize>
std::vector<std::string> TLoadTimeDilutedA2AVecsWDense<FImpl, lowBinSize>::getOutput(void)
{
    std::vector<std::string> out = {getName()};

    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl, int lowBinSize>
void TLoadTimeDilutedA2AVecsWDense<FImpl, lowBinSize>::setup(void)
{
    auto &noise = envGet(SpinColorDiagonalNoise<FImpl>, par().noise);

    if (par().nLow % lowBinSize != 0)
    {
        HADRONS_ERROR(Size, "nLow (" + std::to_string(par().nLow)
                            + ") is not a multiple of the low-mode bin size ("
                            + std::to_string(lowBinSize) + ")");
    }
    if ((par().nLow > 0) && par().lowFilestem.empty())
    {
        HADRONS_ERROR(Argument, "nLow > 0 but lowFilestem is empty");
    }

    const int    nsc   = Ns*FImpl::Dimension;
    unsigned int total = par().nLow + noise.size()*nsc;
    auto         grid  = envGetGrid(FermionField);

    envCreate(std::vector<FermionField>, getName(), 1, 0, grid);
    auto &out = envGet(std::vector<FermionField>, getName());
    out.reserve(total);
    for (unsigned int i = 0; i < total; ++i)
    {
        out.emplace_back(grid, CpuWrite);
    }
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl, int lowBinSize>
void TLoadTimeDilutedA2AVecsWDense<FImpl, lowBinSize>::execute(void)
{
    auto         &noise = envGet(SpinColorDiagonalNoise<FImpl>, par().noise);
    auto         &out   = envGet(std::vector<FermionField>, getName());
    const int    nsc    = Ns*FImpl::Dimension;
    unsigned int offset = 0;

    if (par().nLow > 0)
    {
        int Nb = par().nLow / lowBinSize;
        std::vector<Lattice<LowBinnedSpinor>> bvec(Nb, envGetGrid(Lattice<LowBinnedSpinor>));

        LOG(Message) << "Loading " << par().nLow << " low-mode A2A vectors from "
                     << Nb << " files of " << lowBinSize << " binned vectors" << std::endl;
        A2AVectorsIo::read(bvec, par().lowFilestem, true, vm().getTrajectory());
        A2Autils<FImpl>::template UnpackBinnedVectors<lowBinSize>(out, offset, bvec);
        offset += par().nLow;
    }

    LOG(Message) << "Expanding noise '" << par().noise << "' into "
                 << noise.size()*nsc << " dense high-mode W vectors ("
                 << noise.size() << " hits, " << nsc
                 << " spin-color indices each)" << std::endl;
    for (int h = 0; h < noise.size(); ++h)
    for (int sc = 0; sc < nsc; ++sc)
    {
        startTimer("Dense W expansion");
        out[offset] = noise.getFullVolumeFerm(h, sc);
        offset++;
        stopTimer("Dense W expansion");
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MIO_LoadTimeDilutedA2AVecsWDense_hpp_
