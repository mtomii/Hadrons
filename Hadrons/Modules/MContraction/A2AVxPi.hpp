/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/Modules/MContraction/A2AVxPi.hpp

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
#ifndef Hadrons_MContraction_A2AVxPi_hpp_
#define Hadrons_MContraction_A2AVxPi_hpp_

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

class A2AVxPiPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(A2AVxPiPar,
                                    std::string, left,
				    std::string, mes,
				    int, t_mes_base,
				    int, delt_max,
				    int, delt_min);
};

template <typename FImpl>
class TA2AVxPi: public Module<A2AVxPiPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
public:
    // constructor
    TA2AVxPi(const std::string name);
    // destructor
    virtual ~TA2AVxPi(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(A2AVxPi, TA2AVxPi<FIMPL>, MContraction);

/******************************************************************************
 *                         TA2AVxPi implementation                            *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TA2AVxPi<FImpl>::TA2AVxPi(const std::string name)
: Module<A2AVxPiPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TA2AVxPi<FImpl>::getInput(void)
{
  std::vector<std::string> in = {par().left, par().mes};

  return in;
}

template <typename FImpl>
std::vector<std::string> TA2AVxPi<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TA2AVxPi<FImpl>::setup(void)
{
  typedef typename FImpl::SiteSpinor vobj;
  typedef typename vobj::vector_type vector_type;
  typedef iSpinColourVector<vector_type> SpinColourVector_v;
  //envCreateLat(PropagatorField, getName());
  envCreate(std::vector<SpinColourVector_v>, getName(), 1, 0, Zero());
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TA2AVxPi<FImpl>::execute(void)
{
  typedef typename FImpl::SiteSpinor vobj;
  typedef typename vobj::vector_type vector_type;
  typedef iSpinColourMatrix<vector_type> SpinColourMatrix_v;
  typedef iSpinColourVector<vector_type> SpinColourVector_v;
  typedef iSinglet<vector_type> Scalar_v;
  auto &vec   = envGet(std::vector<SpinColourVector_v>, getName());
  auto &leftV = envGet(std::vector<FermionField>, par().left);
  auto &meson_f = envGet(EigenDiskVector<Complex>, par().mes);

  GridBase *grid = leftV[0].Grid();
  int nt         = env().getDim().back();
  int stepsize = par().delt_max - par().delt_min + 1;
  int v_os  = par().t_mes_base + par().delt_min;

  grid->show_decomposition();

  int orthogdim = 3;

  assert( grid->_ldimensions[orthogdim] >= stepsize );

  std::vector<A2AMatrix<Complex> > meson;
  std::cout << "# t_mes_base: " << par().t_mes_base << std::endl;

  int tmin_rep = ( v_os  + nt ) % grid->_ldimensions[orthogdim];
  int tmax_rep = tmin_rep + stepsize - 1;
  for( int tl = tmin_rep ; tl <= tmax_rep ; ++tl ) {
    int t = tl % grid->_ldimensions[orthogdim] + grid->_lstart[orthogdim];
    int tmes;
    int t_tmes;
    for ( int t0 = 0 ; t0 < nt ; t0 += grid->_ldimensions[orthogdim] ) {
      tmes = ( t0 + par().t_mes_base ) % nt;
      t_tmes = t - tmes;
      if ( t_tmes >= par().delt_min && t_tmes <= par().delt_max ) break;
    }
    if( grid->_lstart[0] + grid->_lstart[1] + grid->_lstart[2] == 0 )
      printf("##  t_op = %d, t_mes = %d\n",t,tmes);
    meson.push_back(meson_f[tmes]);
  }

  int N_i = meson[0].rows();
  assert ( N_i == leftV.size() );
  int N_j = meson[0].cols();

  int rd= grid->_rdimensions[orthogdim];//2
  int e1= grid->_slice_nblock[orthogdim];//1
  int e2= grid->_slice_block [orthogdim];//64 must be 4^3
  int stride=grid->_slice_stride[orthogdim];//128
  std::cout << "rdim: " << rd << std::endl;
  std::cout << "slice_nblock: " << e1 << std::endl;
  std::cout << "slice_block: " << e2 << std::endl;
  std::cout << "slice_stride: " << stride << std::endl;
  int MFrvol = grid->_rdimensions[0]
    *          grid->_rdimensions[1]
    *          grid->_rdimensions[2]
    *          grid->_rdimensions[3]
    /          grid->_rdimensions[orthogdim]
    *          stepsize * N_j;
  vec.assign(MFrvol,Zero());
  std::cout << "Allocated " << vec.size() << " spin-color vectors for V x Pi" << std::endl;
  thread_for(j,N_j,{
    for(int it=0;it<stepsize;it++){
      int r = ( tmin_rep + it ) % grid->_ldimensions[orthogdim];
      int so=r*grid->_ostride[orthogdim];
      for(int i=0;i<N_i;++i){
	//auto lhs_v = leftV[i].View();
	autoView(lhs_v,leftV[i],CpuRead);
	for(int n=0;n<e1;n++){
	for(int b=0;b<e2;b++){
	  int ss = so+n*stride+b;
	  int sv = b+e2*(n+e1*(it+stepsize*j));
	  auto left = lhs_v[ss];
	  for(int s1=0;s1<Ns;s1++)
	  for(int c1=0;c1<Nc;c1++){
	    vec[sv]()(s1)(c1) += left()(s1)(c1) * meson[it](i,j);
	  }
	}}
      }
    }
  });
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_A2AVxPi_hpp_
