/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/Modules/MContraction/A2AFourQuarkContractionMT.hpp

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
#ifndef Hadrons_MContraction_A2AFourQuarkContractionMT_hpp_
#define Hadrons_MContraction_A2AFourQuarkContractionMT_hpp_

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

class A2AFourQuarkContractionMTPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(A2AFourQuarkContractionMTPar,
				    int, ntmat1,
				    int, ntmat2,
				    int, indent1,
				    int, indent2,
                                    std::string, field,
				    std::string, sctypes,
                                    std::string, output,
                                    std::string, mat1,
				    std::string, mat2,
				    std::string, gammas1,
				    std::string, gammas2);
};

template <typename FImpl>
class TA2AFourQuarkContractionMT: public Module<A2AFourQuarkContractionMTPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
public:
    // constructor
    TA2AFourQuarkContractionMT(const std::string name);
    // destructor
    virtual ~TA2AFourQuarkContractionMT(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
  std::vector<std::vector<Gamma::Algebra> >       gamma1_;
  std::vector<std::vector<Gamma::Algebra> >       gamma2_;
  std::vector<std::string> nameg1_;
  std::vector<std::string> nameg2_;
  std::vector<int> types_;
  std::string mat1_;
  std::string mat2_;
};

MODULE_REGISTER_TMP(A2AFourQuarkContractionMT, TA2AFourQuarkContractionMT<FIMPL>, MContraction);

class CorrelatorResult: Serializable
{
public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(CorrelatorResult,
				  std::vector<ComplexD>, correlator);
};

/******************************************************************************
 *                 TA2AFourQuarkContractionMT implementation                  *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TA2AFourQuarkContractionMT<FImpl>::TA2AFourQuarkContractionMT(const std::string name)
: Module<A2AFourQuarkContractionMTPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TA2AFourQuarkContractionMT<FImpl>::getInput(void)
{
  std::vector<std::string> in = {par().mat1, par().mat2, par().field};

  return in;
}

template <typename FImpl>
std::vector<std::string> TA2AFourQuarkContractionMT<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TA2AFourQuarkContractionMT<FImpl>::setup(void)
{
  gamma1_.clear();
  gamma2_.clear();
  std::vector<std::string> tmp1 = strToVec<std::string>(par().gammas1);
  std::vector<std::string> tmp2 = strToVec<std::string>(par().gammas2);
  nameg1_ = tmp1;
  nameg2_ = tmp2;
  assert( tmp1.size() == tmp2.size() );
  for ( int ig = 0; ig < tmp1.size(); ++ig ) {
    std::vector<Gamma::Algebra> vec;
    vec.clear();
    if ( tmp1[ig] == "GammaMU" ) {
      vec = {
	Gamma::Algebra::GammaX,
	Gamma::Algebra::GammaY,
	Gamma::Algebra::GammaZ,
	Gamma::Algebra::GammaT
      };
    } else if ( tmp1[ig] == "GammaMUGamma5" ) {
      vec = {
	Gamma::Algebra::GammaXGamma5,
	Gamma::Algebra::GammaYGamma5,
	Gamma::Algebra::GammaZGamma5,
	Gamma::Algebra::GammaTGamma5
      };
    } else {
      vec = strToVec<Gamma::Algebra>(tmp1[ig]);
    }
    gamma1_.push_back(vec);

    vec.clear();
    if ( tmp2[ig] == "GammaMU" ) {
      vec = {
	Gamma::Algebra::GammaX,
	Gamma::Algebra::GammaY,
	Gamma::Algebra::GammaZ,
	Gamma::Algebra::GammaT
      };
    } else if ( tmp2[ig] == "GammaMUGamma5" ) {
      vec = {
	Gamma::Algebra::GammaXGamma5,
	Gamma::Algebra::GammaYGamma5,
	Gamma::Algebra::GammaZGamma5,
	Gamma::Algebra::GammaTGamma5
      };
    } else {
      vec = strToVec<Gamma::Algebra>(tmp2[ig]);
    }
    gamma2_.push_back(vec);
  }
  types_ = strToVec<int>(par().sctypes);
  mat1_ = par().mat1;
  mat2_ = par().mat2;
  typedef typename FImpl::SiteSpinor vobj;
  typedef typename vobj::vector_type vector_type;
  typedef iSpinColourVector<vector_type> SpinColourVector_v;
  typedef iSpinColourMatrix<vector_type> SpinColourMatrix_v;
  //envCreateLat(PropagatorField, getName());
  envCreate(std::vector<SpinColourMatrix_v>, getName(), 1, 0, Zero());
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TA2AFourQuarkContractionMT<FImpl>::execute(void)
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

  auto &mat1    = envGet(std::vector<SpinColourMatrix_v>, par().mat1);
  auto &mat2    = envGet(std::vector<SpinColourMatrix_v>, par().mat2);
  auto &ntmat1  = par().ntmat1;
  auto &ntmat2  = par().ntmat2;
  auto &indent1 = par().indent1;
  auto &indent2 = par().indent2;

  const int Nsimd = grid->Nsimd();

  int nt = ntmat2 - indent2;
  if ( nt > ntmat1 - indent1 ) nt = ntmat1 - indent1;

  int vol3d = mat1.size() / ntmat1;
  assert ( vol3d == mat2.size() / ntmat2 );

  Scalar_v Cv0 = Zero();
  std::vector<Scalar_v> corr_v0(nt,Cv0);
  int num_corr = gamma1_.size() * types_.size();
  std::vector<std::vector<Scalar_v> > corr_v(num_corr,corr_v0);

  int thread_vol = nt * gamma1_.size() * types_.size();

  thread_for(ittg,thread_vol,{
    int ig = ittg % gamma1_.size();
    int itt  = int(ittg / gamma1_.size());
    int isct = itt % types_.size();
    int it   = int(itt / types_.size());
    int it1 = it + indent1;
    int it2 = it + indent2;
    int itg = ig + gamma1_.size() * isct;
    std::vector<Gamma::Algebra> gvec1 = gamma1_[ig];
    std::vector<Gamma::Algebra> gvec2 = gamma2_[ig];
    for(int ix3d=0;ix3d<vol3d;ix3d++){
      int ix1  = ix3d + vol3d*it1;
      int ix2  = ix3d + vol3d*it2;
      for(int igg=0;igg<gvec1.size();igg++){
	SpinColourMatrix_v WM1 = mat1[ix1] * Gamma(gvec1[igg]);
	SpinColourMatrix_v WM2 = mat2[ix2] * Gamma(gvec2[igg]);
	Scalar_v val = Zero();
	if ( isct == 0 ) {
	  val = trace(WM1) * trace(WM2);
	} else if ( isct == 1 ) {
	  for(int s1=0;s1<Ns;++s1)
	  for(int s2=0;s2<Ns;++s2)
	  for(int c1=0;c1<Nc;++c1)
	  for(int c2=0;c2<Nc;++c2){
	    val()()() += WM1()(s1,s2)(c1,c2) * WM2()(s2,s1)(c2,c1);
	  }
	} else if ( isct == 2 ) {
	  ColourMatrix_v CM1 = Zero();
	  ColourMatrix_v CM2 = Zero();
	  for(int s1=0;s1<Ns;++s1)
	  for(int c1=0;c1<Nc;++c1)
	  for(int c2=0;c2<Nc;++c2){
	    CM1()()(c1,c2) += WM1()(s1,s1)(c1,c2);
	    CM2()()(c1,c2) += WM2()(s1,s1)(c1,c2);
	  }
	  for(int c1=0;c1<Nc;++c1)
	  for(int c2=0;c2<Nc;++c2){
	    val()()() += CM1()()(c1,c2) * CM2()()(c2,c1);
	  }
	} else if ( isct == 3 ) {
	  SpinMatrix_v SM1 = Zero();
	  SpinMatrix_v SM2 = Zero();
	  for(int s1=0;s1<Ns;++s1)
	  for(int s2=0;s2<Ns;++s2)
	  for(int c1=0;c1<Nc;++c1){
	    SM1()(s1,s2)() += WM1()(s1,s2)(c1,c1);
	    SM2()(s1,s2)() += WM2()(s1,s2)(c1,c1);
	  }
	  for(int s1=0;s1<Ns;++s1)
	  for(int s2=0;s2<Ns;++s2){
	    val()()() += SM1()(s1,s2)() * SM2()(s2,s1)();
	  }
	}
	corr_v[itg][it]()()() += val()()();
      }// igg
    }// ix3d
  });

  ComplexD C0(0.,0.);
  std::vector<ComplexD> corr0(nt,C0);
  std::vector<std::vector<ComplexD> > corr(num_corr,corr0);
  thread_for(ittg,thread_vol,{
    int ig = ittg % gamma1_.size();
    int itt  = int(ittg / gamma1_.size());
    int isct = itt % types_.size();
    int it   = int(itt / types_.size());
    int itg = ig + gamma1_.size() * isct;
    ExtractBuffer<Scalar_s> extracted(Nsimd);
    extract(corr_v[itg][it],extracted);
    for(int isimd=0;isimd<Nsimd;isimd++){
      corr[itg][it]=corr[itg][it]+extracted[isimd];
    }
  });
  for(int i=0;i<thread_vol;i++){
    int it = i % nt;
    int itg = int(i / nt);
    grid->GlobalSum(corr[itg][it]);
  }

  std::string filename = par().output + "/" + ModuleBase::resultFilename(getName(), vm().getTrajectory());
  LOG(Message) << "Saving correlator to '" << filename << "'" << std::endl;
  if( grid->_lstart[0] + grid->_lstart[1] + grid->_lstart[2] + grid->_lstart[3] == 0 ) {

    ResultWriter writer(filename);

    std::vector<std::string> gam1 = strToVec<std::string>(par().gammas1);
    std::vector<std::string> gam2 = strToVec<std::string>(par().gammas2);
    for(int itg=0;itg<num_corr;itg++){
      CorrelatorResult out;
      int ig = itg % gamma1_.size();
      int isct = int( itg / gamma1_.size() );
      out.correlator = corr[itg];
      std::string dataSet = getName() + "_" + gam1[ig] + "_" + gam2[ig] + "_sort" + std::to_string(isct);
      write(writer, dataSet, out);
    }
  }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_A2AFourQuarkContractionMT_hpp_
