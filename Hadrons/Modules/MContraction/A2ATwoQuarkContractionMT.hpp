/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/Modules/MContraction/A2ATwoQuarkContractionMT.hpp

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
#ifndef Hadrons_MContraction_A2ATwoQuarkContractionMT_hpp_
#define Hadrons_MContraction_A2ATwoQuarkContractionMT_hpp_

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

class A2ATwoQuarkContractionMTPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(A2ATwoQuarkContractionMTPar,
				    int, ntmat,
                                    std::string, field,
                                    std::string, output,
				    std::string, mat,
				    std::string, gammas);
};

template <typename FImpl>
class TA2ATwoQuarkContractionMT: public Module<A2ATwoQuarkContractionMTPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
public:
    // constructor
    TA2ATwoQuarkContractionMT(const std::string name);
    // destructor
    virtual ~TA2ATwoQuarkContractionMT(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
  std::vector<std::vector<Gamma::Algebra> >       gamma_;
  std::vector<std::string> nameg_;
  std::vector<int> types_;
  std::string mat_;
};

MODULE_REGISTER_TMP(A2ATwoQuarkContractionMT, TA2ATwoQuarkContractionMT<FIMPL>, MContraction);

class CorrelatorResult2: Serializable
{
public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(CorrelatorResult2,
				  std::vector<ComplexD>, correlator);
};


/******************************************************************************
 *                 TA2ATwoQuarkContractionMT implementation                   *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TA2ATwoQuarkContractionMT<FImpl>::TA2ATwoQuarkContractionMT(const std::string name)
: Module<A2ATwoQuarkContractionMTPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TA2ATwoQuarkContractionMT<FImpl>::getInput(void)
{
  std::vector<std::string> in = {par().mat, par().field};

  return in;
}

template <typename FImpl>
std::vector<std::string> TA2ATwoQuarkContractionMT<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TA2ATwoQuarkContractionMT<FImpl>::setup(void)
{
  gamma_.clear();
  std::vector<std::string> tmp = strToVec<std::string>(par().gammas);
  nameg_ = tmp;
  for ( int ig = 0; ig < tmp.size(); ++ig ) {
    std::vector<Gamma::Algebra> vec;
    vec.clear();
    if ( tmp[ig] == "GammaMU" ) {
      vec = {
	Gamma::Algebra::GammaX,
	Gamma::Algebra::GammaY,
	Gamma::Algebra::GammaZ,
	Gamma::Algebra::GammaT
      };
    } else if ( tmp[ig] == "GammaMUGamma5" ) {
      vec = {
	Gamma::Algebra::GammaXGamma5,
	Gamma::Algebra::GammaYGamma5,
	Gamma::Algebra::GammaZGamma5,
	Gamma::Algebra::GammaTGamma5
      };
    } else {
      vec = strToVec<Gamma::Algebra>(tmp[ig]);
    }
    gamma_.push_back(vec);
  }
  mat_ = par().mat;
  typedef typename FImpl::SiteSpinor vobj;
  typedef typename vobj::vector_type vector_type;
  typedef iSpinColourVector<vector_type> SpinColourVector_v;
  typedef iSpinColourMatrix<vector_type> SpinColourMatrix_v;
  //envCreateLat(PropagatorField, getName());
  envCreate(std::vector<SpinColourMatrix_v>, getName(), 1, 0, Zero());
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TA2ATwoQuarkContractionMT<FImpl>::execute(void)
{
  typedef typename FImpl::SiteSpinor vobj;
  typedef typename vobj::vector_type vector_type;
  typedef typename vobj::scalar_type scalar_type;
  typedef iSpinColourMatrix<vector_type> SpinColourMatrix_v;
  typedef iSpinMatrix<vector_type> SpinMatrix_v;
  typedef iColourMatrix<vector_type> ColourMatrix_v;
  typedef iSpinColourVector<vector_type> SpinColourVector_v;
  typedef iSinglet<vector_type> Scalar_v;
  typedef iSinglet<scalar_type> Scalar_s;

  auto &field = envGet(std::vector<FermionField>, par().field);
  GridBase *grid = field[0].Grid();

  auto &mat  = envGet(std::vector<SpinColourMatrix_v>, par().mat);
  auto &nt   = par().ntmat;

  const int Nsimd = grid->Nsimd();

  int vol3d = mat.size() / nt;

  Scalar_v Cv0 = Zero();
  std::vector<Scalar_v> corr_v0(nt,Cv0);
  int num_corr = gamma_.size();
  std::vector<std::vector<Scalar_v> > corr_v(num_corr,corr_v0);

  int thread_vol = nt * gamma_.size();

  thread_for(itg,thread_vol,{
    int ig = itg % gamma_.size();
    int it  = int(itg / gamma_.size());
    //int it1 = it + indent1;
    //int itg = ig + gamma1_.size() * isct;
    std::vector<Gamma::Algebra> gvec = gamma_[ig];
    for(int ix3d=0;ix3d<vol3d;ix3d++){
      int ix  = ix3d + vol3d*it;
      for(int igg=0;igg<gvec.size();igg++){
	SpinColourMatrix_v WM = mat[ix] * Gamma(gvec[igg]);
	Scalar_v val = trace(WM);
	corr_v[ig][it]()()() += val()()();
      }// igg
    }// ix3d
  });

  ComplexD C0(0.,0.);
  std::vector<ComplexD> corr0(nt,C0);
  std::vector<std::vector<ComplexD> > corr(num_corr,corr0);
  thread_for(itg,thread_vol,{
    int ig = itg % gamma_.size();
    int it  = int(itg / gamma_.size());
    //int itg = ig + gamma1_.size() * isct;
    ExtractBuffer<Scalar_s> extracted(Nsimd);
    extract(corr_v[ig][it],extracted);
    for(int isimd=0;isimd<Nsimd;isimd++){
      corr[ig][it]=corr[ig][it]+extracted[isimd];
    }
  });
  for(int i=0;i<thread_vol;i++){
    int it = i % nt;
    int ig = int(i / nt);
    grid->GlobalSum(corr[ig][it]);
  }

  std::string filename = par().output + "/" + ModuleBase::resultFilename(getName(), vm().getTrajectory());
  LOG(Message) << "Saving correlator to '" << filename << "'" << std::endl;
  if( grid->_lstart[0] + grid->_lstart[1] + grid->_lstart[2] + grid->_lstart[3] == 0 ) {

    ResultWriter writer(filename);

    std::vector<std::string> gam = strToVec<std::string>(par().gammas);
    for(int ig=0;ig<num_corr;ig++){
      CorrelatorResult2 out;
      out.correlator = corr[ig];
      std::string dataSet = getName() + "_" + gam[ig];
      write(writer, dataSet, out);
    }
  }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_A2ATwoQuarkContractionMT_hpp_
