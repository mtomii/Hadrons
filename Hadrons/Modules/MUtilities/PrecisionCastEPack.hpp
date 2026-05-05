/*
 * PrecisionCastEPack.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2023
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Fionn O hOgain <fionn.o.hogain@ed.ac.uk>
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
#ifndef Hadrons_MUtilities_PrecisionCastEPack_hpp_
#define Hadrons_MUtilities_PrecisionCastEPack_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/EigenPack.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                          Precision cast module                             *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MUtilities)

class PrecisionCastEPackPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(PrecisionCastEPackPar,
                                    std::string, in,
				    bool, redBlack);
};

template <typename FImplIn, typename FImplOut>
class TPrecisionCastEPack: public Module<PrecisionCastEPackPar>
{
public:
  typedef BaseFermionEigenPack<FImplIn>   PackIn;
  //typedef typename PackIn::Field          FieldIn;
  //typedef typename PackIn::FieldIo        FieldIoIn;
  typedef BaseFermionEigenPack<FImplOut>  BasePackOut;
  typedef FermionEigenPack<FImplOut>      PackOut;
  typedef typename PackOut::Field         FieldOut;

public:
    // constructor
    TPrecisionCastEPack(const std::string name);
    // destructor
    virtual ~TPrecisionCastEPack(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(FermionDoublePrecisionCastEPack, 
                    ARG(TPrecisionCastEPack<FIMPLF, FIMPLD>), 
		    MUtilities);

/******************************************************************************
 *                    TPrecisionCastEPack implementation                      *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImplIn, typename FImplOut>
TPrecisionCastEPack<FImplIn, FImplOut>::TPrecisionCastEPack(const std::string name)
: Module<PrecisionCastEPackPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImplIn, typename FImplOut>
std::vector<std::string> TPrecisionCastEPack<FImplIn, FImplOut>::getInput(void)
{
    std::vector<std::string> in = {par().in};
    
    return in;
}

template <typename FImplIn, typename FImplOut>
std::vector<std::string> TPrecisionCastEPack<FImplIn, FImplOut>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImplIn, typename FImplOut>
void TPrecisionCastEPack<FImplIn, FImplOut>::setup(void)
{
    GridBase *gridIo = nullptr;

    /*
    if (typeHash<Field>() != typeHash<FieldIo>())
    {
        gridIo = envGetRbGrid(FieldIo, env().getObjectLs(par().in));
    }
    // learn from EigenPackLCDecompress.hpp as needed
    */

    auto &in  = envGet(PackIn,  par().in);
    GridBase* grid        = in.evec[0].Grid();
    //GridBase* grid        = getGrid<FieldIn>(true, env().getObjectLs(par().in));

    //auto        &action    = envGet(FermionOperator<FImplOut>, par().action);
    //GridBase* frbGrid = action.FermionRedBlackGrid();
    GridBase* frbGrid = getGrid<FieldOut>(par().redBlack, env().getObjectLs(par().in));

    envCreateDerived(BasePackOut, PackOut, getName(), 
		     env().getObjectLs(par().in),
		     in.evec.size(), frbGrid, gridIo);

    //envCreateLat(BasePackOut, getName());
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImplIn, typename FImplOut>
void TPrecisionCastEPack<FImplIn, FImplOut>::execute(void)
{
    LOG(Message) << "Casting Eigen Pack '" << par().in << "'" << std::endl;
    LOG(Message) << "In  type: " << typeName<PackIn>() << std::endl;
    LOG(Message) << "Out type: " << typeName<PackOut>() << std::endl;

    auto &in  = envGet(PackIn,  par().in);
    auto &out = envGetDerived(BasePackOut, PackOut, getName());

    out.record = in.record;

    for(unsigned int i=0; i<in.evec.size(); ++i)
    {
      //precisionChange(out.eval[i], in.eval[i]);
      precisionChange(out.evec[i], in.evec[i]);
      out.eval[i] = in.eval[i];
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MUtilities_PrecisionCastEPack_hpp_
