/*
 * LoadTimeDilutedSpinColorDiagonalNoise.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2023
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Fionn O hOgain <fionn.o.hogain@ed.ac.uk>
 * Author: Fionn Ó hÓgáin <fionnoh@gmail.com>
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
#ifndef Hadrons_MIO_LoadTimeDilutedSpinColorDiagonalNoise_hpp_
#define Hadrons_MIO_LoadTimeDilutedSpinColorDiagonalNoise_hpp_
#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/DilutedNoise.hpp>
#include <Hadrons/A2AVectors.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *              Load time diluted spin-color diagonal noise                   *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MIO)
class LoadTimeDilutedSpinColorDiagonalNoisePar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(
        LoadTimeDilutedSpinColorDiagonalNoisePar,
        std::vector<std::string>, fileStems,
        unsigned int,             nNoisePerFileStem);
};

template <typename FImpl>
class TLoadTimeDilutedSpinColorDiagonalNoise:
    public Module<LoadTimeDilutedSpinColorDiagonalNoisePar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
public:
    // constructor
    TLoadTimeDilutedSpinColorDiagonalNoise(const std::string name);
    // destructor
    virtual ~TLoadTimeDilutedSpinColorDiagonalNoise(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(
    LoadTimeDilutedSpinColorDiagonalNoise,
    TLoadTimeDilutedSpinColorDiagonalNoise<FIMPL>, MIO);
MODULE_REGISTER_TMP(
    LoadZTimeDilutedSpinColorDiagonalNoise,
    TLoadTimeDilutedSpinColorDiagonalNoise<ZFIMPL>, MIO);

/******************************************************************************
 *          TLoadTimeDilutedSpinColorDiagonalNoise implementation             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TLoadTimeDilutedSpinColorDiagonalNoise<FImpl>::
TLoadTimeDilutedSpinColorDiagonalNoise(const std::string name)
: Module<LoadTimeDilutedSpinColorDiagonalNoisePar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string>
TLoadTimeDilutedSpinColorDiagonalNoise<FImpl>::getInput(void)
{
    std::vector<std::string> in;

    return in;
}

template <typename FImpl>
std::vector<std::string>
TLoadTimeDilutedSpinColorDiagonalNoise<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};

    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TLoadTimeDilutedSpinColorDiagonalNoise<FImpl>::setup(void)
{
    assert(!par().fileStems.empty());
    assert(par().nNoisePerFileStem > 0);

    unsigned int nNoise =
        par().fileStems.size()*par().nNoisePerFileStem;

    LOG(Message) << "setting up time-diluted, spin-color diagonal noise"
                 << std::endl;
    envCreateDerived(SpinColorDiagonalNoise<FImpl>,
                     TimeDilutedNoise<FImpl>,
                     getName(), 1, envGetGrid(FermionField), nNoise);
    LOG(Message) << "time-diluted, spin-color diagonal noise set up!"
                 << std::endl;
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TLoadTimeDilutedSpinColorDiagonalNoise<FImpl>::execute(void)
{
    auto &noise = envGet(SpinColorDiagonalNoise<FImpl>, getName());
    auto &v     = noise.getNoise();

    assert(v.size() ==
           par().fileStems.size()*par().nNoisePerFileStem);

    std::vector<ComplexField> vTmp(par().nNoisePerFileStem,
                                   envGetGrid(ComplexField));

    for (unsigned int i = 0; i < par().fileStems.size(); ++i)
    {
        LOG(Message) << "Loading " << par().nNoisePerFileStem
                     << " time-diluted, spin-color diagonal noise fields from "
                     << par().fileStems[i] << std::endl;
        A2AVectorsIo::read(vTmp, par().fileStems[i], false,
                           vm().getTrajectory());

        for (unsigned int j = 0; j < par().nNoisePerFileStem; ++j)
        {
            v[i*par().nNoisePerFileStem + j] = vTmp[j];
        }
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE
#endif // Hadrons_MIO_LoadTimeDilutedSpinColorDiagonalNoise_hpp_
