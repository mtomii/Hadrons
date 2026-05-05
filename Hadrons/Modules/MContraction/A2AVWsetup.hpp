/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/Modules/MContraction/A2AVWsetup.hpp

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
#ifndef Hadrons_MContraction_A2AVWsetup_hpp_
#define Hadrons_MContraction_A2AVWsetup_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/A2AMatrix.hpp>
#include <Hadrons/DiskVector.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                  From closed loop from all-to-all vectors                  *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MContraction)

class A2AVWsetupPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(A2AVWsetupPar,
                                    std::string, vw,
				    int, t_mes_base,
				    int, delt_max,
				    int, delt_min,
				    bool, ifconj);
};

template <typename FImpl>
class TA2AVWsetup: public Module<A2AVWsetupPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
public:
    // constructor
    TA2AVWsetup(const std::string name);
    // destructor
    virtual ~TA2AVWsetup(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(A2AVWsetup, TA2AVWsetup<FIMPL>, MContraction);

/******************************************************************************
 *                        TA2AVWsetup implementation                          *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TA2AVWsetup<FImpl>::TA2AVWsetup(const std::string name)
: Module<A2AVWsetupPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TA2AVWsetup<FImpl>::getInput(void)
{
  std::vector<std::string> in = {par().vw};

  return in;
}

template <typename FImpl>
std::vector<std::string> TA2AVWsetup<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TA2AVWsetup<FImpl>::setup(void)
{
  typedef typename FImpl::SiteSpinor vobj;
  typedef typename vobj::vector_type vector_type;
  typedef iSpinColourVector<vector_type> SpinColourVector_v;
  envCreate(std::vector<SpinColourVector_v>, getName(), 1, 0, Zero());
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TA2AVWsetup<FImpl>::execute(void)
{
  typedef typename FImpl::SiteSpinor vobj;
  typedef typename vobj::vector_type vector_type;
  typedef iSpinColourMatrix<vector_type> SpinColourMatrix_v;
  typedef iSpinColourVector<vector_type> SpinColourVector_v;
  typedef iSinglet<vector_type> Scalar_v;
  auto &vec   = envGet(std::vector<SpinColourVector_v>, getName());
  auto &VW = envGet(std::vector<FermionField>, par().vw);

  GridBase *grid = VW[0].Grid();
  int nt         = env().getDim().back();
  int stepsize = par().delt_max - par().delt_min + 1;
  int v_os  = par().t_mes_base + par().delt_min;

  int orthogdim = 3;

  assert( grid->_ldimensions[orthogdim] >= stepsize );

  std::cout << "# t_mes_base: " << par().t_mes_base << std::endl;

  int tmin_rep = ( v_os  + nt ) % grid->_ldimensions[orthogdim];

  int N_i = VW.size();

  int rd= grid->_rdimensions[orthogdim];//2
  int e1= grid->_slice_nblock[orthogdim];//1
  int e2= grid->_slice_block [orthogdim];//64 must be 4^3
  int stride=grid->_slice_stride[orthogdim];//128
  int MFrvol = grid->_rdimensions[0]
    *          grid->_rdimensions[1]
    *          grid->_rdimensions[2]
    *          grid->_rdimensions[3]
    /          grid->_rdimensions[orthogdim]
    *          stepsize * N_i;
  vec.assign(MFrvol,Zero());
  std::cout << "Allocated " << vec.size() << " spin-color vectors for V x Pi" << std::endl;
  thread_for(i,N_i,{
    //auto rhs_w = VW[i].View();
    autoView(rhs_w,VW[i],CpuRead);
    for(int it=0;it<stepsize;it++){
      int r = ( tmin_rep + it ) % grid->_ldimensions[orthogdim];
      int so=r*grid->_ostride[orthogdim];
      for(int n=0;n<e1;n++){
      for(int b=0;b<e2;b++){
	int ss = so+n*stride+b;
	int sv = b+e2*(n+e1*(it+stepsize*i));
	if (par().ifconj) {
	  vec[sv] = conjugate(rhs_w[ss]);
	} else {
	  vec[sv] = rhs_w[ss];
	}
      }}
    }
  });
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_A2AVWsetup_hpp_
