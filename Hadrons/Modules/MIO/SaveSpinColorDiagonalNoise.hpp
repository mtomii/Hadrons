/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: Hadrons/Modules/MIO/SaveSpinColorDiagonalNoise.hpp

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
#ifndef Hadrons_MIO_SaveSpinColorDiagonalNoise_hpp_
#define Hadrons_MIO_SaveSpinColorDiagonalNoise_hpp_
#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/DilutedNoise.hpp>
#include <Hadrons/A2AVectors.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                  Module to save spin-color diagonal noise                  *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MIO)
class SaveSpinColorDiagonalNoisePar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(SaveSpinColorDiagonalNoisePar,
                                    std::string, filestem,
                                    std::string, noise);
};
template <typename FImpl>
class TSaveSpinColorDiagonalNoise: public Module<SaveSpinColorDiagonalNoisePar>
{
public:
    // constructor
    TSaveSpinColorDiagonalNoise(const std::string name);
    // destructor
    virtual ~TSaveSpinColorDiagonalNoise(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};
MODULE_REGISTER_TMP(SaveSpinColorDiagonalNoise,
                    TSaveSpinColorDiagonalNoise<FIMPL>, MIO);
MODULE_REGISTER_TMP(SaveZSpinColorDiagonalNoise,
                    TSaveSpinColorDiagonalNoise<ZFIMPL>, MIO);

/******************************************************************************
 *               TSaveSpinColorDiagonalNoise implementation                  *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TSaveSpinColorDiagonalNoise<FImpl>::
TSaveSpinColorDiagonalNoise(const std::string name)
: Module<SaveSpinColorDiagonalNoisePar>(name)
{}
// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TSaveSpinColorDiagonalNoise<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().noise};

    return in;
}

template <typename FImpl>
std::vector<std::string> TSaveSpinColorDiagonalNoise<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {};

    return out;
}
// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TSaveSpinColorDiagonalNoise<FImpl>::setup(void)
{
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TSaveSpinColorDiagonalNoise<FImpl>::execute(void)
{
    auto &noise = envGet(SpinColorDiagonalNoise<FImpl>, par().noise);
    auto &v     = noise.getNoise();

    assert(!v.empty());

    LOG(Message) << "Saving " << v.size()
                 << " spin-color diagonal noise fields to "
                 << par().filestem << std::endl;
    A2AVectorsIo::write(par().filestem, v, false, vm().getTrajectory());
}

END_MODULE_NAMESPACE
END_HADRONS_NAMESPACE
#endif // Hadrons_MIO_SaveSpinColorDiagonalNoise_hpp_
