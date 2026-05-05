/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/Modules/MIO/SaveBinnedA2AVecs.hpp

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
#ifndef Hadrons_MIO_SaveBinnedA2AVecs_hpp_
#define Hadrons_MIO_SaveBinnedA2AVecs_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/A2AVectors.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                    Module to load all-to-all vectors                       *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MIO)

class SaveBinnedA2AVecsPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(SaveBinnedA2AVecsPar,
                                    std::string,  filestem,
                                    bool,         multiFile,
				    std::string,  a2a);
};

template <typename FImpl, int binSize>
class TSaveBinnedA2AVecs: public Module<SaveBinnedA2AVecsPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
  typedef typename FImpl::SiteSpinor::vector_type vector_type;
  typedef iVector<iVector<iVector<vector_type, Nc>, Ns>, binSize > SiteSpinorSet;
public:
    // constructor
    TSaveBinnedA2AVecs(const std::string name);
    // destructor
    virtual ~TSaveBinnedA2AVecs(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(SaveBinnedA2AVecs48, ARG(TSaveBinnedA2AVecs<FIMPL, 48>), MIO);
MODULE_REGISTER_TMP(SaveBinnedA2AVecs64, ARG(TSaveBinnedA2AVecs<FIMPL, 64>), MIO);
MODULE_REGISTER_TMP(SaveBinnedA2AVecs96, ARG(TSaveBinnedA2AVecs<FIMPL, 96>), MIO);
MODULE_REGISTER_TMP(SaveBinnedA2AVecs144, ARG(TSaveBinnedA2AVecs<FIMPL, 144>), MIO);
MODULE_REGISTER_TMP(SaveBinnedA2AVecs149, ARG(TSaveBinnedA2AVecs<FIMPL, 149>), MIO);
MODULE_REGISTER_TMP(SaveBinnedA2AVecs173, ARG(TSaveBinnedA2AVecs<FIMPL, 173>), MIO);
MODULE_REGISTER_TMP(SaveBinnedA2AVecs192, ARG(TSaveBinnedA2AVecs<FIMPL, 192>), MIO);
MODULE_REGISTER_TMP(SaveBinnedA2AVecs196, ARG(TSaveBinnedA2AVecs<FIMPL, 196>), MIO);
MODULE_REGISTER_TMP(SaveBinnedA2AVecs197, ARG(TSaveBinnedA2AVecs<FIMPL, 197>), MIO);
MODULE_REGISTER_TMP(SaveBinnedA2AVecs216, ARG(TSaveBinnedA2AVecs<FIMPL, 216>), MIO);
MODULE_REGISTER_TMP(SaveBinnedA2AVecs298, ARG(TSaveBinnedA2AVecs<FIMPL, 298>), MIO);


/******************************************************************************
 *                     TSaveBinnedA2AVecs implementation                      *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl, int binSize>
TSaveBinnedA2AVecs<FImpl, binSize>::TSaveBinnedA2AVecs(const std::string name)
: Module<SaveBinnedA2AVecsPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl, int binSize>
std::vector<std::string> TSaveBinnedA2AVecs<FImpl, binSize>::getInput(void)
{
  std::vector<std::string> in = {par().a2a};
    
    return in;
}

template <typename FImpl, int binSize>
std::vector<std::string> TSaveBinnedA2AVecs<FImpl, binSize>::getOutput(void)
{
    std::vector<std::string> out = {};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl, int binSize>
void TSaveBinnedA2AVecs<FImpl, binSize>::setup(void)
{
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl, int binSize>
void TSaveBinnedA2AVecs<FImpl, binSize>::execute(void)
{
  auto &v = envGet(std::vector<FermionField>, par().a2a );
  assert( v.size() % binSize == 0 );
  const int Nb = v.size() / binSize;

  std::vector<Lattice<SiteSpinorSet> > bvec(Nb,envGetGrid(Lattice<SiteSpinorSet>));

  for( int i = 0 ; i < v.size() ; i += binSize ) {
    for( int j = 0 ; j < binSize ; ++j ) {
      pokeLorentz(bvec[i/binSize],v[i+j],j);
    }
  }
  LOG(Message) << "Saving " << v.size() << " A2A vectors in a binned format to "
	       << Nb << " files" << std::endl;
  A2AVectorsIo::write(par().filestem, bvec, par().multiFile, vm().getTrajectory());
}

END_MODULE_NAMESPACE
END_HADRONS_NAMESPACE

#endif // Hadrons_MIO_SaveBinnedA2AVecs_hpp_
