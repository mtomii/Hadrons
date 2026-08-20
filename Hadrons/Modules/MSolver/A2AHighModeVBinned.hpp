/*
 * A2AHighModeVBinned.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
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
#ifndef Hadrons_MSolver_A2AHighModeVBinned_hpp_
#define Hadrons_MSolver_A2AHighModeVBinned_hpp_
#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Solver.hpp>
#include <Hadrons/A2AVectors.hpp>
#include <Hadrons/DilutedNoise.hpp>

BEGIN_HADRONS_NAMESPACE
/******************************************************************************
 *                    Create high-mode all-to-all V vectors                   *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MSolver)
class A2AHighModeVBinnedPar: Serializable
{
public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(A2AHighModeVBinnedPar,
                                  std::string, noise,
                                  std::string, action,
                                  std::string, solver,
                                  std::string, output,
                                  bool,        multiFile);
};
template <typename FImpl, int binSize>
class TA2AHighModeVBinned : public Module<A2AHighModeVBinnedPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    SOLVER_TYPE_ALIASES(FImpl,);
    typedef HADRONS_DEFAULT_SCHUR_A2A<FImpl> A2A;
    typedef typename FImpl::SiteSpinor::vector_type vector_type;
    typedef iVector<iVector<iVector<vector_type, Nc>, Ns>, binSize> SiteSpinorSet;
public:
    // constructor
    TA2AHighModeVBinned(const std::string name);
    // destructor
    virtual ~TA2AHighModeVBinned(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};
MODULE_REGISTER_TMP(A2AHighModeVBinned96,
    ARG(TA2AHighModeVBinned<FIMPL, 96>), MSolver);
MODULE_REGISTER_TMP(A2AHighModeVBinned192,
    ARG(TA2AHighModeVBinned<FIMPL, 192>), MSolver);
MODULE_REGISTER_TMP(ZA2AHighModeVBinned96,
    ARG(TA2AHighModeVBinned<ZFIMPL, 96>), MSolver);
MODULE_REGISTER_TMP(ZA2AHighModeVBinned192,
    ARG(TA2AHighModeVBinned<ZFIMPL, 192>), MSolver);
/******************************************************************************
 *                    TA2AHighModeVBinned implementation                      *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl, int binSize>
TA2AHighModeVBinned<FImpl, binSize>::TA2AHighModeVBinned(const std::string name)
: Module<A2AHighModeVBinnedPar>(name)
{}
// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl, int binSize>
std::vector<std::string> TA2AHighModeVBinned<FImpl, binSize>::getInput(void)
{
    std::vector<std::string> in = {par().solver, par().noise};
    return in;
}

template <typename FImpl, int binSize>
std::vector<std::string> TA2AHighModeVBinned<FImpl, binSize>::getOutput(void)
{
    std::vector<std::string> out = {};
    return out;
}
// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl, int binSize>
void TA2AHighModeVBinned<FImpl, binSize>::setup(void)
{
    auto &noise  = envGet(SpinColorDiagonalNoise<FImpl>, par().noise);
    auto &action = envGet(FMat, par().action);
    auto &solver = envGet(Solver, par().solver);
    int  Ls      = env().getObjectLs(par().action);
    assert(noise.fermSize() % binSize == 0);
    envTmpLat(FermionField, "v");
    envTmpLat(Lattice<SiteSpinorSet>, "vBin");
    if (Ls > 1)
    {
        envTmpLat(FermionField, "f5", Ls);
    }
    envTmp(A2A, "a2a", 1, action, solver);
}
// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl, int binSize>
void TA2AHighModeVBinned<FImpl, binSize>::execute(void)
{
    auto &noise = envGet(SpinColorDiagonalNoise<FImpl>, par().noise);
    int  Ls     = env().getObjectLs(par().action);
    envGetTmp(FermionField, v);
    envGetTmp(Lattice<SiteSpinorSet>, vBin);
    envGetTmp(A2A, a2a);
    ScidacWriter vWriter(vBin.Grid()->IsBoss());
    LOG(Message) << "Computing high-mode part of all-to-all V vectors "
                 << "using noise '" << par().noise << "' ("
                 << noise.fermSize() << " noise vectors)" << std::endl;
    if ((!par().output.empty()) && (!par().multiFile))
    {
        A2AVectorsIo::openWriter(vWriter, par().output, vBin.Grid(),
                                 vm().getTrajectory());
    }
    // High modes
    for (unsigned int ih = 0; ih < noise.fermSize(); ih++)
    {
        startTimer("V high mode");
        LOG(Message) << "V vector i = " << ih
                     << " (stochastic mode)" << std::endl;
        if (Ls == 1)
        {
            a2a.makeHighModeV(v, noise.getFerm(ih));
        }
        else
        {
            envGetTmp(FermionField, f5);
            a2a.makeHighModeV5D(v, f5, noise.getFerm(ih));
        }
        stopTimer("V high mode");
        pokeLorentz(vBin, v, ih % binSize);

        if (((ih + 1) % binSize == 0) && (!par().output.empty()))
        {
            unsigned int ib = ih / binSize;
            startTimer("V I/O");
            if (par().multiFile)
            {
                A2AVectorsIo::writeElement(par().output, vBin, ib,
                                           vm().getTrajectory());
            }
            else
            {
                A2AVectorsIo::writeRecord(vWriter, vBin, ib);
            }
            stopTimer("V I/O");
        }
    }
    if ((!par().output.empty()) && (!par().multiFile))
    {
        vWriter.close();
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE
#endif // Hadrons_MSolver_A2AHighModeVBinned_hpp_
