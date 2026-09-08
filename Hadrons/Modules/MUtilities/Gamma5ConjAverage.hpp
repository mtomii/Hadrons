/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: Hadrons/Modules/MUtilities/Gamma5ConjAverage.hpp

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

/* Derived from Hadrons/Modules/MIO/WriteProp.hpp.
 * out(x) = (prop(x) + gamma5 * prop(x)^dagger * gamma5)/2.
 * The dagger is the adjoint in both spin and colour indices.
 * This is a site-local operation: coordinates are not shifted or reflected.
 * The input is left unchanged; the output has the same Grid and Ls metadata.
 */
#ifndef Hadrons_MUtilities_Gamma5ConjAverage_hpp_
#define Hadrons_MUtilities_Gamma5ConjAverage_hpp_
#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <cassert>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *            Average with the gamma5-conjugated propagator field             *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MUtilities)
class Gamma5ConjAveragePar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(Gamma5ConjAveragePar,
                                    std::string, prop);
};
template <typename FImpl>
class TGamma5ConjAverage: public Module<Gamma5ConjAveragePar>
{
public:
    typedef typename FImpl::PropagatorField Field;
public:
    // constructor
    TGamma5ConjAverage(const std::string name);
    // destructor
    virtual ~TGamma5ConjAverage(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};
MODULE_REGISTER_TMP(Gamma5ConjAverage, TGamma5ConjAverage<FIMPL>, MUtilities);

/******************************************************************************
 *                     TGamma5ConjAverage implementation                      *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TGamma5ConjAverage<FImpl>::TGamma5ConjAverage(const std::string name)
: Module<Gamma5ConjAveragePar>(name)
{}
// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TGamma5ConjAverage<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().prop};

    return in;
}

template <typename FImpl>
std::vector<std::string> TGamma5ConjAverage<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};

    return out;
}
// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TGamma5ConjAverage<FImpl>::setup(void)
{
    assert(!par().prop.empty());

    const auto &prop = envGet(Field, par().prop);

    envCreate(Field, getName(), env().getObjectLs(par().prop), prop.Grid());
}
// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TGamma5ConjAverage<FImpl>::execute(void)
{
    const auto &prop = envGet(Field, par().prop);
    auto       &out = envGet(Field, getName());
    Gamma      g5(Gamma::Algebra::Gamma5);

    LOG(Message) << "Averaging '"
                 << par().prop << "' with its gamma5 conjugate, output '" << getName()
                 << "'" << std::endl;

    startTimer("Gamma5 conjugation average");
    out = Real(0.5)*(prop + g5*adj(prop)*g5);
    stopTimer("Gamma5 conjugation average");
}
END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE
#endif // Hadrons_MUtilities_Gamma5ConjAverage_hpp_
