/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/Modules/MIO/LoadBinnedA2AVecs.hpp

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
#ifndef Hadrons_MIO_LoadBinnedA2AVecs_hpp_
#define Hadrons_MIO_LoadBinnedA2AVecs_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/A2AVectors.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                    Module to load all-to-all vectors                       *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MIO)

class LoadBinnedA2AVecsPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(LoadBinnedA2AVecsPar,
                                    std::string,  filestem,
                                    bool,         multiFile,
                                    unsigned int, size);
};

template <typename FImpl, int binSize>
class TLoadBinnedA2AVecs: public Module<LoadBinnedA2AVecsPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
  typedef typename FImpl::SiteSpinor::vector_type vector_type;
  typedef iVector<iVector<iVector<vector_type, Nc>, Ns>, binSize > SiteSpinorSet;
public:
    // constructor
    TLoadBinnedA2AVecs(const std::string name);
    // destructor
    virtual ~TLoadBinnedA2AVecs(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
  // read & peek
  //template <typename BinnedField>
  //void readpeek(std::vector<FermionField> &vec, BinnedField &bvec);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(LoadBinnedA2AVecs48, ARG(TLoadBinnedA2AVecs<FIMPL, 48>), MIO);
MODULE_REGISTER_TMP(LoadBinnedA2AVecs64, ARG(TLoadBinnedA2AVecs<FIMPL, 64>), MIO);
MODULE_REGISTER_TMP(LoadBinnedA2AVecs96, ARG(TLoadBinnedA2AVecs<FIMPL, 96>), MIO);
MODULE_REGISTER_TMP(LoadBinnedA2AVecs144, ARG(TLoadBinnedA2AVecs<FIMPL, 144>), MIO);
MODULE_REGISTER_TMP(LoadBinnedA2AVecs149, ARG(TLoadBinnedA2AVecs<FIMPL, 149>), MIO);
MODULE_REGISTER_TMP(LoadBinnedA2AVecs173, ARG(TLoadBinnedA2AVecs<FIMPL, 173>), MIO);
MODULE_REGISTER_TMP(LoadBinnedA2AVecs192, ARG(TLoadBinnedA2AVecs<FIMPL, 192>), MIO);
MODULE_REGISTER_TMP(LoadBinnedA2AVecs196, ARG(TLoadBinnedA2AVecs<FIMPL, 196>), MIO);
MODULE_REGISTER_TMP(LoadBinnedA2AVecs197, ARG(TLoadBinnedA2AVecs<FIMPL, 197>), MIO);
MODULE_REGISTER_TMP(LoadBinnedA2AVecs216, ARG(TLoadBinnedA2AVecs<FIMPL, 216>), MIO);
MODULE_REGISTER_TMP(LoadBinnedA2AVecs298, ARG(TLoadBinnedA2AVecs<FIMPL, 298>), MIO);

/******************************************************************************
 *                     TLoadBinnedA2AVecs implementation                      *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl, int binSize>
TLoadBinnedA2AVecs<FImpl, binSize>::TLoadBinnedA2AVecs(const std::string name)
: Module<LoadBinnedA2AVecsPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl, int binSize>
std::vector<std::string> TLoadBinnedA2AVecs<FImpl, binSize>::getInput(void)
{
    std::vector<std::string> in;
    
    return in;
}

template <typename FImpl, int binSize>
std::vector<std::string> TLoadBinnedA2AVecs<FImpl, binSize>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl, int binSize>
void TLoadBinnedA2AVecs<FImpl, binSize>::setup(void)
{
  envCreate(std::vector<FermionField>, getName(), 1, par().size, 
	    envGetGrid(FermionField));
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl, int binSize>
void TLoadBinnedA2AVecs<FImpl, binSize>::execute(void)
{
  auto &v = envGet(std::vector<FermionField>, getName());

  assert( v.size() % binSize == 0 );
  int Nb = v.size() / binSize;

  std::vector<Lattice<SiteSpinorSet> > bvec(Nb,envGetGrid(Lattice<SiteSpinorSet>));

  LOG(Message) << "Loading " << v.size() << " A2A vectors from "
	       << Nb << " files of " << binSize << " binned vectors" << std::endl;
  A2AVectorsIo::read(bvec, par().filestem, par().multiFile, vm().getTrajectory());
  for( int i = 0 ; i < v.size() ; i += binSize ) {
    for( int j = 0 ; j < binSize ; ++j ) {
      v[i+j] = peekLorentz(bvec[i/binSize],j);
    }
  }
}

END_MODULE_NAMESPACE
END_HADRONS_NAMESPACE

#endif // Hadrons_MIO_LoadBinnedA2AVecs_hpp_
