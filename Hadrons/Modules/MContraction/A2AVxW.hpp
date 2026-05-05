/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/Modules/MContraction/A2AVxW.hpp

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
#ifndef Hadrons_MContraction_A2AVxW_hpp_
#define Hadrons_MContraction_A2AVxW_hpp_

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

class A2AVxWPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(A2AVxWPar,
				    int, Nmode,
				    int, indentL,
				    int, indentR,
                                    std::string, field,
                                    std::string, right,
				    std::string, left);
};

template <typename FImpl>
class TA2AVxW: public Module<A2AVxWPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
public:
    // constructor
    TA2AVxW(const std::string name);
    // destructor
    virtual ~TA2AVxW(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(A2AVxW, TA2AVxW<FIMPL>, MContraction);

/******************************************************************************
 *                         TA2AVxW implementation                           *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TA2AVxW<FImpl>::TA2AVxW(const std::string name)
: Module<A2AVxWPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TA2AVxW<FImpl>::getInput(void)
{
  std::vector<std::string> in = {par().left, par().right, par().field};

  return in;
}

template <typename FImpl>
std::vector<std::string> TA2AVxW<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TA2AVxW<FImpl>::setup(void)
{
  typedef typename FImpl::SiteSpinor vobj;
  typedef typename vobj::vector_type vector_type;
  typedef iSpinColourVector<vector_type> SpinColourVector_v;
  typedef iSpinColourMatrix<vector_type> SpinColourMatrix_v;
  //envCreateLat(PropagatorField, getName());
  envCreate(std::vector<SpinColourMatrix_v>, getName(), 1, 0, Zero());
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TA2AVxW<FImpl>::execute(void)
{
  typedef typename FImpl::SiteSpinor vobj;
  typedef typename vobj::vector_type vector_type;
  typedef iSpinColourMatrix<vector_type> SpinColourMatrix_v;
  typedef iSpinColourVector<vector_type> SpinColourVector_v;
  typedef iSinglet<vector_type> Scalar_v;

  auto &field = envGet(std::vector<FermionField>, par().field);
  GridBase *grid = field[0].Grid();

  auto &mat    = envGet(std::vector<SpinColourMatrix_v>, getName());
  auto &rightW = envGet(std::vector<SpinColourVector_v>, par().right);
  auto &leftV  = envGet(std::vector<SpinColourVector_v>, par().left);
  auto &indentL = par().indentL;
  auto &indentR = par().indentR;

  int orthogdim = 3;
  int e1= grid->_slice_nblock[orthogdim];//1
  int e2= grid->_slice_block [orthogdim];//64 must be 4^3

  int N_i = par().Nmode;

  std::cout << "## Vsize: " << leftV.size() << ", Wsize: " << rightW.size() << std::endl;
  int MFrvolL = int(leftV.size() / N_i);
  int MFrvolR = int(rightW.size() / N_i);
  assert ( leftV.size() % N_i == 0 );
  assert ( rightW.size() % N_i == 0 );
  int ntL = int( MFrvolL / e1 / e2 );
  int ntR = int( MFrvolR / e1 / e2 );
  assert( ntL - indentL == ntR - indentR );
  assert ( MFrvolL % (e1*e2) == 0 );
  assert ( MFrvolR % (e1*e2) == 0 );

  int vol = ( ntL - indentL ) * e1 * e2;

  mat.assign(vol,Zero());
  std::cout << "Allocated " << mat.size() << " spin-color matrices for V x W" << std::endl;
  thread_for(ix,vol,{
    int ixL = ix + indentL * e1 * e2;
    int ixR = ix + indentR * e1 * e2;
    for(int i=0;i<N_i;i++){
      int svL = ixL+MFrvolL*i;
      int svR = ixR+MFrvolR*i;
      for(int s1=0;s1<Ns;s1++)
      for(int s2=0;s2<Ns;s2++)
      for(int c1=0;c1<Nc;c1++)
      for(int c2=0;c2<Nc;c2++){
	mat[ix]()(s1,s2)(c1,c2) += leftV[svL]()(s1)(c1) * rightW[svR]()(s2)(c2);
      }
    }
  });
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_A2AVxW_hpp_
