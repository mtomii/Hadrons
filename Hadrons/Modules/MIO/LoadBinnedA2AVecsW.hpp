/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: Hadrons/Modules/MIO/LoadBinnedA2AVecsW.hpp

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
#ifndef Hadrons_MIO_LoadBinnedA2AVecsW_hpp_
#define Hadrons_MIO_LoadBinnedA2AVecsW_hpp_
#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/A2AVectors.hpp>
#include <Hadrons/DilutedNoise.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                 Module to load all-to-all W vectors                        *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MIO)
class LoadBinnedA2AVecsWPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(LoadBinnedA2AVecsWPar,
                                    std::string,  filestem,
                                    bool,         multiFile,
                                    unsigned int, lowSize,
                                    std::string,  noise);
};
template <typename FImpl, int binSize>
class TLoadBinnedA2AVecsW: public Module<LoadBinnedA2AVecsWPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    typedef typename FImpl::SiteSpinor::vector_type vector_type;
    typedef iVector<iVector<iVector<vector_type, Nc>, Ns>, binSize> SiteSpinorSet;
public:
    // constructor
    TLoadBinnedA2AVecsW(const std::string name);
    // destructor
    virtual ~TLoadBinnedA2AVecsW(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
    void makeHighModeW(FermionField &wout, const FermionField &noise);
};
MODULE_REGISTER_TMP(LoadBinnedA2AVecsW100,
                    ARG(TLoadBinnedA2AVecsW<FIMPL, 100>), MIO);
MODULE_REGISTER_TMP(LoadBinnedA2AVecsW200,
                    ARG(TLoadBinnedA2AVecsW<FIMPL, 200>), MIO);
MODULE_REGISTER_TMP(LoadZBinnedA2AVecsW100,
                    ARG(TLoadBinnedA2AVecsW<ZFIMPL, 100>), MIO);
MODULE_REGISTER_TMP(LoadZBinnedA2AVecsW200,
                    ARG(TLoadBinnedA2AVecsW<ZFIMPL, 200>), MIO);
/******************************************************************************
 *                  TLoadBinnedA2AVecsW implementation                        *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl, int binSize>
TLoadBinnedA2AVecsW<FImpl, binSize>::
TLoadBinnedA2AVecsW(const std::string name)
: Module<LoadBinnedA2AVecsWPar>(name)
{}
// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl, int binSize>
std::vector<std::string> TLoadBinnedA2AVecsW<FImpl, binSize>::getInput(void)
{
    std::vector<std::string> in = {par().noise};

    return in;
}

template <typename FImpl, int binSize>
std::vector<std::string> TLoadBinnedA2AVecsW<FImpl, binSize>::getOutput(void)
{
    std::vector<std::string> out = {getName()};

    return out;
}
// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl, int binSize>
void TLoadBinnedA2AVecsW<FImpl, binSize>::setup(void)
{
    auto &noise = envGet(SpinColorDiagonalNoise<FImpl>, par().noise);

    assert(par().lowSize % binSize == 0);
    assert((par().lowSize == 0) or (!par().filestem.empty()));
    envCreate(std::vector<FermionField>, getName(), 1,
              par().lowSize + noise.fermSize(), envGetGrid(FermionField));
}

// high mode W /////////////////////////////////////////////////////////////////
template <typename FImpl, int binSize>
void TLoadBinnedA2AVecsW<FImpl, binSize>::
makeHighModeW(FermionField &wout, const FermionField &noise)
{
    wout = noise;
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl, int binSize>
void TLoadBinnedA2AVecsW<FImpl, binSize>::execute(void)
{
    auto &noise = envGet(SpinColorDiagonalNoise<FImpl>, par().noise);
    auto &w     = envGet(std::vector<FermionField>, getName());
    int  Nb     = par().lowSize / binSize;

    assert(w.size() == par().lowSize + noise.fermSize());

    // Low modes
    if (par().lowSize > 0)
    {
        Lattice<SiteSpinorSet> bvec(envGetGrid(Lattice<SiteSpinorSet>));

        LOG(Message) << "Loading " << par().lowSize << " low-mode W vectors from "
                     << Nb << " files of " << binSize << " binned vectors"
                     << std::endl;
        if (par().multiFile)
        {
            for (int i = 0; i < Nb; ++i)
            {
                A2AVectorsIo::readElement(par().filestem, bvec, i,
                                          vm().getTrajectory());
                for (int j = 0; j < binSize; ++j)
                {
                    w[i*binSize + j] = peekLorentz(bvec, j);
                }
            }
        }
        else
        {
            ScidacReader binReader;

            A2AVectorsIo::openReader(binReader, par().filestem,
                                     vm().getTrajectory());
            for (int i = 0; i < Nb; ++i)
            {
                A2AVectorsIo::readRecord(binReader, bvec, i);
                for (int j = 0; j < binSize; ++j)
                {
                    w[i*binSize + j] = peekLorentz(bvec, j);
                }
            }
            binReader.close();
        }
    }

    // High modes
    for (unsigned int ih = 0; ih < noise.fermSize(); ih++)
    {
        startTimer("W high mode");
        LOG(Message) << "W vector i = " << par().lowSize + ih
                     << " (" << ((par().lowSize > 0) ? "high " : "")
                     << "stochastic mode)" << std::endl;
        makeHighModeW(w[par().lowSize + ih], noise.getFerm(ih));
        stopTimer("W high mode");
    }
}

END_MODULE_NAMESPACE
END_HADRONS_NAMESPACE
#endif // Hadrons_MIO_LoadBinnedA2AVecsW_hpp_
