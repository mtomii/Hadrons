/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/Modules/MGauge/APESmear.hpp

Copyright (C) 2015-2020

Author: Michael Abramczyk

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
#ifndef Hadrons_MGauge_APESmear_hpp_
#define Hadrons_MGauge_APESmear_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         APESmear                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MGauge)

class APESmearPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(APESmearPar,
                                    std::string, gauge,
                                    double, alpha,
                                    unsigned int, N,
                                    unsigned int, orthog_axis
                                    );
};

template <typename GImpl>
class TAPESmear: public Module<APESmearPar>
{
public:
    GAUGE_TYPE_ALIASES(GImpl,);
public:
    // constructor
    TAPESmear(const std::string name);
    // destructor
    virtual ~TAPESmear(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(APESmear, TAPESmear<GIMPL>, MGauge);

/******************************************************************************
 *                 TAPESmear implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename GImpl>
TAPESmear<GImpl>::TAPESmear(const std::string name)
: Module<APESmearPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename GImpl>
std::vector<std::string> TAPESmear<GImpl>::getInput(void)
{
    std::vector<std::string> in = {par().gauge};
    
    return in;
}

template <typename GImpl>
std::vector<std::string> TAPESmear<GImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename GImpl>
void TAPESmear<GImpl>::setup(void)
{
    envCreateLat(GaugeField, getName());
    envTmpLat(GaugeField, "Utmp");
}

// execution ///////////////////////////////////////////////////////////////////
template <typename GImpl>
void TAPESmear<GImpl>::execute(void)
{
    LOG(Message) << "APE smearing " << par().gauge << std::endl;
    
    auto &U = envGet(GaugeField, par().gauge);
    auto &Usmeared = envGet(GaugeField, getName());
    envGetTmp(GaugeField, Utmp);
    
    Usmeared = U;
    Utmp = Zero();
    
    //LOG(Message) << "Printing original gauge field" << std::endl;
    //LOG(Message) << peekLorentz(U, 3) << std::endl;
    
    Smear_APE<GImpl> APE(par().alpha); // alpha is the relative smearing weight [reference needed]
    
    // Perform N APE smearing steps
    for(int i = 0; i < par().N; i++)
    {
		APE.smearAlt(Utmp, Usmeared, par().orthog_axis); // compute sums of staples and store in Utmp
		
		// Conventional factors to ensure alpha is the relative weight
		Usmeared = (1.0 - par().alpha) * Usmeared;
				
		Usmeared = Usmeared + Utmp; // update Usmeared by one step
            Usmeared = ProjectOnGroup(Usmeared); // project onto SU(3)
    }
    
    //LOG(Message) << "Printing smeared gauge field" << std::endl;
    //LOG(Message) << peekLorentz(Usmeared, 3) << std::endl;

	LOG(Message) << "Parameters: alpha = " << par().alpha << ", N = " << par().N << std::endl;
    LOG(Message) << "plaquette= " << WilsonLoops<GImpl>::avgPlaquette(U)
                 << std::endl;
    // copy back time link since it was fouled up by 1-alpha factor
    LatticeColourMatrix Ut(U.Grid());
    Ut = PeekIndex<LorentzIndex>(U, 3);
    PokeIndex<LorentzIndex>(Usmeared, Ut, 3);
    LOG(Message) << "plaquette= " << std::setprecision(16) <<  WilsonLoops<GImpl>::avgPlaquette(Usmeared)
                 << std::endl;
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MGauge_APESmear_hpp_
