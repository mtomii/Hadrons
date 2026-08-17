/*
 * A2ALowModeBinned.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2026
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Fionn O hOgain <fionn.o.hogain@ed.ac.uk>
 * Author: Fionn Ó hÓgáin <fionnoh@gmail.com>
 * Author: fionnoh <fionnoh@gmail.com>
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
#ifndef Hadrons_MUtilities_A2ALowModeBinned_hpp_
#define Hadrons_MUtilities_A2ALowModeBinned_hpp_
#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/EigenPack.hpp>
#include <Hadrons/A2AVectors.hpp>

BEGIN_HADRONS_NAMESPACE
/******************************************************************************
 *                  Create low-mode all-to-all V & W vectors                  *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MUtilities)
class A2ALowModeBinnedPar: Serializable
{
public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(A2ALowModeBinnedPar,
                                  std::string, action,
                                  std::string, eigenPack,
                                  std::string, output,
                                  bool,        multiFile);
};
template <typename FImpl, typename Pack, int binSize>
class TA2ALowModeBinned : public Module<A2ALowModeBinnedPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    typedef HADRONS_DEFAULT_SCHUR_A2A<FImpl> A2A;
    typedef typename FImpl::SiteSpinor::vector_type vector_type;
    typedef iVector<iVector<iVector<vector_type, Nc>, Ns>, binSize> SiteSpinorSet;
public:
    // constructor
    TA2ALowModeBinned(const std::string name);
    // destructor
    virtual ~TA2ALowModeBinned(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
    unsigned int Nl_{0};
};
MODULE_REGISTER_TMP(A2ALowModeBinned100,
    ARG(TA2ALowModeBinned<FIMPL, BaseFermionEigenPack<FIMPL>, 100>), MUtilities);
MODULE_REGISTER_TMP(A2ALowModeBinned200,
    ARG(TA2ALowModeBinned<FIMPL, BaseFermionEigenPack<FIMPL>, 200>), MUtilities);
MODULE_REGISTER_TMP(ZA2ALowModeBinned100,
    ARG(TA2ALowModeBinned<ZFIMPL, BaseFermionEigenPack<ZFIMPL>, 100>), MUtilities);
MODULE_REGISTER_TMP(ZA2ALowModeBinned200,
    ARG(TA2ALowModeBinned<ZFIMPL, BaseFermionEigenPack<ZFIMPL>, 200>), MUtilities);
/******************************************************************************
 *                    TA2ALowModeBinned implementation                        *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl, typename Pack, int binSize>
TA2ALowModeBinned<FImpl, Pack, binSize>::TA2ALowModeBinned(const std::string name)
: Module<A2ALowModeBinnedPar>(name)
{}
// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl, typename Pack, int binSize>
std::vector<std::string> TA2ALowModeBinned<FImpl, Pack, binSize>::getInput(void)
{
    assert(!par().eigenPack.empty());
    std::vector<std::string> in = {par().action, par().eigenPack};

    return in;
}

template <typename FImpl, typename Pack, int binSize>
std::vector<std::string> TA2ALowModeBinned<FImpl, Pack, binSize>::getOutput(void)
{
    std::vector<std::string> out = {};

    return out;
}
// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl, typename Pack, int binSize>
void TA2ALowModeBinned<FImpl, Pack, binSize>::setup(void)
{
    assert(!par().eigenPack.empty());
    auto &action = envGet(FMat, par().action);
    auto &epack  = envGet(Pack, par().eigenPack);
    int  Ls      = env().getObjectLs(par().action);

    Nl_ = epack.evec.size();
    assert(Nl_ % binSize == 0);
    envTmpLat(FermionField, "v");
    envTmpLat(FermionField, "w");
    envTmpLat(Lattice<SiteSpinorSet>, "vBin");
    envTmpLat(Lattice<SiteSpinorSet>, "wBin");
    if (Ls > 1)
    {
        envTmpLat(FermionField, "f5", Ls);
    }
    envTmp(A2A, "a2a", 1, action);
}
// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl, typename Pack, int binSize>
void TA2ALowModeBinned<FImpl, Pack, binSize>::execute(void)
{
    int Ls = env().getObjectLs(par().action);
    envGetTmp(FermionField, v);
    envGetTmp(FermionField, w);
    envGetTmp(Lattice<SiteSpinorSet>, vBin);
    envGetTmp(Lattice<SiteSpinorSet>, wBin);
    envGetTmp(A2A, a2a);
    ScidacWriter vWriter(vBin.Grid()->IsBoss());
    ScidacWriter wWriter(wBin.Grid()->IsBoss());

    LOG(Message) << "Getting low-mode part of all-to-all vectors "
                 << " from '" << par().eigenPack << "' ("
                 << Nl_ << " eigenPack modes)" << std::endl;

    if ((!par().output.empty()) && (!par().multiFile))
    {
        A2AVectorsIo::openWriter(vWriter, par().output + "_v", vBin.Grid(),
                                 vm().getTrajectory());
        A2AVectorsIo::openWriter(wWriter, par().output + "_w", wBin.Grid(),
                                 vm().getTrajectory());
    }

    // Low modes
    for (unsigned int il = 0; il < Nl_; il++)
    {
        auto &epack = envGet(Pack, par().eigenPack);
        startTimer("V low mode");
        LOG(Message) << "V vector i = " << il << " (low mode)" << std::endl;
        if (Ls == 1)
        {
            a2a.makeLowModeV(v, epack.evec[il], epack.eval[il]);
        }
        else
        {
            envGetTmp(FermionField, f5);
            a2a.makeLowModeV5D(v, f5, epack.evec[il], epack.eval[il]);
        }
        stopTimer("V low mode");
        startTimer("W low mode");
        LOG(Message) << "W vector i = " << il << " (low mode)" << std::endl;
        if (Ls == 1)
        {
            a2a.makeLowModeW(w, epack.evec[il], epack.eval[il]);
        }
        else
        {
            envGetTmp(FermionField, f5);
            a2a.makeLowModeW5D(w, f5, epack.evec[il], epack.eval[il]);
        }
        stopTimer("W low mode");

        pokeLorentz(vBin, v, il % binSize);
        pokeLorentz(wBin, w, il % binSize);

        if (((il + 1) % binSize == 0) && (!par().output.empty()))
        {
            unsigned int ib = il / binSize;

            startTimer("V I/O");
            if (par().multiFile)
            {
                A2AVectorsIo::writeElement(par().output + "_v", vBin, ib,
                                           vm().getTrajectory());
            }
            else
            {
                A2AVectorsIo::writeRecord(vWriter, vBin, ib);
            }
            stopTimer("V I/O");
            startTimer("W I/O");
            if (par().multiFile)
            {
                A2AVectorsIo::writeElement(par().output + "_w", wBin, ib,
                                           vm().getTrajectory());
            }
            else
            {
                A2AVectorsIo::writeRecord(wWriter, wBin, ib);
            }
            stopTimer("W I/O");
        }
    }

    if ((!par().output.empty()) && (!par().multiFile))
    {
        vWriter.close();
        wWriter.close();
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE
#endif // Hadrons_MUtilities_A2ALowModeBinned_hpp_
