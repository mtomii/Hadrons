/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/Modules/MContraction/A2AExtendedMesonField.hpp

Copyright (C) 2015-2019

Author: Masaaki Tomii <masaaki.tomii@uconn.edu>
*************************************************************************************/
/*  END LEGAL */
#ifndef Hadrons_MContraction_A2AExtendedMesonField_hpp_
#define Hadrons_MContraction_A2AExtendedMesonField_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/A2AMatrix.hpp>

BEGIN_HADRONS_NAMESPACE
//   _
//  / \
//  \ /
// --x--
//

/******************************************************************************
 *                All-to-all extended meson field creation                    *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MContraction)

class A2AExtendedMesonFieldPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(A2AExtendedMesonFieldPar,
                                    int, cacheBlock,
				    std::string, types,
                                    std::string, left,
                                    std::string, right,
				    std::string, loop_vw1,
				    std::string, loop_vw2,
                                    std::string, output,
                                    std::string, gammas1,
				    std::string, gammas2);
};

class A2AExtendedMesonFieldMetadata: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(A2AExtendedMesonFieldMetadata,
                                    std::string, gamma1,
				    std::string, gamma2);
};

template <typename FImpl>
class TA2AExtendedMesonField : public Module<A2AExtendedMesonFieldPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
  //typedef typename FImpl::FermionField FermionField;
  typedef typename FImpl::SiteSpinor vobj;
  typedef typename vobj::scalar_object sobj;
  typedef typename vobj::scalar_type scalar_type;
  typedef typename vobj::vector_type vector_type;
public:
    // constructor
    TA2AExtendedMesonField(const std::string name);
    // destructor
    virtual ~TA2AExtendedMesonField(void){};
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
};

MODULE_REGISTER(A2AExtendedMesonField, TA2AExtendedMesonField<FIMPL>, MContraction);

/******************************************************************************
*               TA2AExtendedMesonField implementation                         *
******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TA2AExtendedMesonField<FImpl>::TA2AExtendedMesonField(const std::string name)
: Module<A2AExtendedMesonFieldPar>(name)
{
}

// dependencies/products ///////////////////////////////////////////////////////
// Modification needed??
template <typename FImpl>
std::vector<std::string> TA2AExtendedMesonField<FImpl>::getInput(void)
{
  std::vector<std::string> in = {par().left, par().right, par().loop_vw1, par().loop_vw2};

    return in;
}

template <typename FImpl>
std::vector<std::string> TA2AExtendedMesonField<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {};

    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TA2AExtendedMesonField<FImpl>::setup(void)
{
    gamma1_.clear();
    gamma2_.clear();
    std::vector<std::string> tmp1 = strToVec<std::string>(par().gammas1);
    std::vector<std::string> tmp2 = strToVec<std::string>(par().gammas2);
    nameg1_ = tmp1;
    nameg2_ = tmp2;
    types_ = strToVec<int>(par().types);
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
      } else if ( tmp1[ig] == "SigmaMUNU" ) {
	vec = {
	       Gamma::Algebra::SigmaXY,
	       Gamma::Algebra::SigmaXZ,
	       Gamma::Algebra::SigmaXT,
	       Gamma::Algebra::SigmaYZ,
	       Gamma::Algebra::SigmaYT,
	       Gamma::Algebra::SigmaZT
	};
      } else if ( tmp1[ig] == "SigmaMUNUGamma5" ) {
	vec = {
	       Gamma::Algebra::MinusSigmaZT,
	       Gamma::Algebra::SigmaYT,
	       Gamma::Algebra::MinusSigmaYZ,
	       Gamma::Algebra::MinusSigmaXT,
	       Gamma::Algebra::SigmaXZ,
	       Gamma::Algebra::MinusSigmaXY
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
      } else if ( tmp2[ig] == "SigmaMUNU" ) {
	vec = {
	       Gamma::Algebra::SigmaXY,
	       Gamma::Algebra::SigmaXZ,
	       Gamma::Algebra::SigmaXT,
	       Gamma::Algebra::SigmaYZ,
	       Gamma::Algebra::SigmaYT,
	       Gamma::Algebra::SigmaZT
	};
      } else if ( tmp2[ig] == "SigmaMUNUGamma5" ) {
	vec = {
	       Gamma::Algebra::MinusSigmaZT,
	       Gamma::Algebra::SigmaYT,
	       Gamma::Algebra::MinusSigmaYZ,
	       Gamma::Algebra::MinusSigmaXT,
	       Gamma::Algebra::SigmaXZ,
	       Gamma::Algebra::MinusSigmaXY
	};
      } else {
	vec = strToVec<Gamma::Algebra>(tmp2[ig]);
      }
      gamma2_.push_back(vec);
    }
    auto &left  = envGet(std::vector<FermionField>, par().left);
    GridBase *grid = left[0].Grid();
    envCreateLat(PropagatorField, "propagatorLoop");
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TA2AExtendedMesonField<FImpl>::execute(void)
{
  typedef iSpinColourVector<vector_type> SpinColourVector_v;
  typedef iSpinColourMatrix<vector_type> SpinColourMatrix_v;
  typedef iSpinMatrix<vector_type> SpinMatrix_v;
  typedef iSinglet<vector_type> Scalar_v;
  typedef iSinglet<scalar_type> Scalar_s;
    auto &left  = envGet(std::vector<FermionField>, par().left);
    auto &right = envGet(std::vector<FermionField>, par().right);
    auto &loop1 = envGet(std::vector<FermionField>, par().loop_vw1);
    auto &loop2 = envGet(std::vector<FermionField>, par().loop_vw2);

    GridBase *grid = left[0].Grid();

    int orthogdim = 3;
    int rd=grid->_rdimensions[orthogdim];//2
    int ld=grid->_ldimensions[orthogdim];
    int Nd=grid->_ndimension;
    int Nsimd=grid->Nsimd();
    int e1=    grid->_slice_nblock[orthogdim];//1
    int e2=    grid->_slice_block [orthogdim];//64 must be 4^3
    int stride=grid->_slice_stride[orthogdim];//128
    LOG(Message) << "Computing all-to-all EXTENDED meson fields" << std::endl;
    LOG(Message) << "R dimension: " << rd << std::endl;
    LOG(Message) << "Slice nblock: " << e1 << std::endl;
    LOG(Message) << "Slice block: " << e2 << std::endl;
    LOG(Message) << "Slice stride: " << stride << std::endl;

    int nt         = env().getDim().back();
    int N_i        = left.size();
    int N_j        = right.size();
    //int block      = par().block;
    int cacheBlock = par().cacheBlock;
    Vector<HADRONS_A2AM_IO_TYPE> mBuf; mBuf.resize(nt*N_i*N_j);
    Vector<Complex> mCacheBuf; mCacheBuf.resize(nt*cacheBlock*cacheBlock);

  int ngamma = 0;
    for ( int ig = 0 ; ig < gamma1_.size() ; ++ig ) {
      ngamma += gamma1_[ig].size();
    }

    LOG(Message) << "Left: '" << par().left << "' Right: '"
		 << par().right << "'" << std::endl;
    LOG(Message) << "A2AVectors for loop: '" << par().loop_vw1
		 << "' and '" << par().loop_vw2 << "'" << std::endl;
    LOG(Message) << "Spin bilinears1:" << std::endl;
    for (auto &g: gamma1_)
    {
        LOG(Message) << "  " << g << std::endl;
    }
    LOG(Message) << "Spin bilinears2:" << std::endl;
    for (auto &g: gamma2_)
    {
        LOG(Message) << "  " << g << std::endl;
    }
    LOG(Message) << "Meson field size: " << nt << "*" << N_i << "*" << N_j 
                 << " (filesize " << sizeString(nt*N_i*N_j*sizeof(HADRONS_A2AM_IO_TYPE)) 
                 << "/momentum/bilinear)" << std::endl;

    auto &loop = envGet(PropagatorField, "propagatorLoop");
    loop = Zero();
    for (unsigned int i = 0; i < loop1.size(); ++i)
    {
      loop += outerProduct(loop1[i], loop2[i]);
    }
    autoView(loopv,loop,CpuRead);
    LOG(Message) << "Quark loop calculated" << std::endl;

    std::vector<SpinColourVector_v> zeroVec(rd*grid->_ostride[orthogdim],Zero());
    std::vector<std::vector<SpinColourVector_v> > leftv(N_i,zeroVec);
    for ( unsigned int i = 0; i < N_i; i++ ){
      autoView(tmp,left[i],CpuRead);
      thread_for(r,rd,{
	int so=r*grid->_ostride[orthogdim];
	for(int n=0;n<e1;n++)
	for(int b=0;b<e2;b++){
	  int ss = so+n*stride+b;
	  leftv[i][ss] = conjugate(tmp[ss]);
	}
      });
    }
    LOG(Message) << "Memory for left vectors allocated" << std::endl;

    std::function<void(SpinColourVector_v&, const SpinColourMatrix_v&, const SpinColourVector_v&, const std::vector<Gamma::Algebra>, const std::vector<Gamma::Algebra>)> addLoopRight;
    for (int &type: types_){
      std::vector<SpinColourMatrix_v> tloopv(rd*grid->_ostride[orthogdim],Zero());
      LOG(Message) << "Starting calculation with type " << type << std::endl;
      if ( type == 0 ) {
	thread_for(r,rd,{
	  int so=r*grid->_ostride[orthogdim];
	  for(int n=0;n<e1;n++)
	  for(int b=0;b<e2;b++){
	    int ss = so+n*stride+b;
	    for(int s1=0;s1<Ns;++s1)
	    for(int s2=0;s2<Ns;++s2){
	      tloopv[ss]()(s1,s2)(0,0) += loopv[ss]()(s1,s2)(0,0) + loopv[ss]()(s1,s2)(1,1) + loopv[ss]()(s1,s2)(2,2);
	    }
	  }
	});
	LOG(Message) << "Color contraction for quark loop done" << std::endl;

	addLoopRight = [this](SpinColourVector_v& loopRight, const SpinColourMatrix_v& loopm, const SpinColourVector_v rightv,
			      const std::vector<Gamma::Algebra> gamma1, const std::vector<Gamma::Algebra> gamma2){
	  SpinMatrix_v spinLoop = Zero();
	  for(int s1=0;s1<Ns;++s1)
	  for(int s2=0;s2<Ns;++s2){
	    spinLoop()(s1,s2)() = loopm()(s1,s2)(0,0);
	  }
	  for(int mu=0;mu<gamma1.size();++mu){
	    SpinMatrix_v GLoop = Gamma(gamma2[mu]) * spinLoop;
	    Scalar_v trGLoop;
	    trGLoop()()() = GLoop()(0,0)() + GLoop()(1,1)() + GLoop()(2,2)() + GLoop()(3,3)();
	    SpinColourVector_v Grightv = Gamma(gamma1[mu]) * rightv;
	    for(int s=0;s<Ns;++s)
	    for(int c=0;c<Nc;++c){
	      loopRight()(s)(c) += Grightv()(s)(c) * trGLoop()()();
	    }
	  }
	};
      }
      if ( type == 1 ) {
	addLoopRight = [this](SpinColourVector_v& loopRight, const SpinColourMatrix_v& loopm, const SpinColourVector_v rightv,
			      const std::vector<Gamma::Algebra> gamma1, const std::vector<Gamma::Algebra> gamma2){
          for(int s=0;s<Ns;++s)
          for(int c=0;c<Nc;++c){
	    loopRight()(s)(c)
	      += loopm()(s,0)(c,0) * rightv()(0)(0)
	      +  loopm()(s,0)(c,1) * rightv()(0)(1)
	      +  loopm()(s,0)(c,2) * rightv()(0)(2)
	      +  loopm()(s,1)(c,0) * rightv()(1)(0)
	      +  loopm()(s,1)(c,1) * rightv()(1)(1)
	      +  loopm()(s,1)(c,2) * rightv()(1)(2)
	      +  loopm()(s,2)(c,0) * rightv()(2)(0)
	      +  loopm()(s,2)(c,1) * rightv()(2)(1)
	      +  loopm()(s,2)(c,2) * rightv()(2)(2)
	      +  loopm()(s,3)(c,0) * rightv()(3)(0)
	      +  loopm()(s,3)(c,1) * rightv()(3)(1)
	      +  loopm()(s,3)(c,2) * rightv()(3)(2);
	  }
	};
      }
      if ( type == 2 ) {
	addLoopRight = [this](SpinColourVector_v& loopRight, const SpinColourMatrix_v& loopm, const SpinColourVector_v rightv,
			      const std::vector<Gamma::Algebra> gamma1, const std::vector<Gamma::Algebra> gamma2){
	  for(int mu=0;mu<gamma1.size();++mu){
	    int s1 = int(mu/Ns);
	    int s2 = mu%Ns;
	    SpinColourVector_v Grightv = Gamma(gamma1[mu]) * rightv;
	    for(int s=0;s<Ns;++s)
	    for(int c=0;c<Nc;++c){
	      loopRight()(s)(c)
		+= loopm()(s1,s2)(c,0) * Grightv()(s)(0)
		+  loopm()(s1,s2)(c,1) * Grightv()(s)(1)
		+  loopm()(s1,s2)(c,2) * Grightv()(s)(2);
	    }
	  }
	};
      }
      if ( type == 3 ) {
	addLoopRight = [this](SpinColourVector_v& loopRight, const SpinColourMatrix_v& loopm, const SpinColourVector_v rightv,
			      const std::vector<Gamma::Algebra> gamma1, const std::vector<Gamma::Algebra> gamma2){
          for(int s=0;s<Ns;++s)
          for(int c=0;c<Nc;++c){
	    loopRight()(s)(c)
	      += loopm()(s,0)(0,0) * rightv()(0)(c)
	      +  loopm()(s,1)(0,0) * rightv()(1)(c)
	      +  loopm()(s,2)(0,0) * rightv()(2)(c)
	      +  loopm()(s,3)(0,0) * rightv()(3)(c);
	  }
	};
      }

      for (int ig = 0 ; ig < gamma1_.size() ; ++ig ){
	LOG(Message) << "Starting calculation of " << nameg1_[ig] << "_" << nameg2_[ig] << std::endl;
	if ( type == 1 ) {
	  std::vector<Gamma::Algebra> gamma1 = gamma1_[ig];
	  std::vector<Gamma::Algebra> gamma2 = gamma2_[ig];
	  thread_for(r,rd,{
	    int so=r*grid->_ostride[orthogdim];
	    for(int n=0;n<e1;n++)
	    for(int b=0;b<e2;b++){
	      int ss = so+n*stride+b;
	      tloopv[ss] = Zero();
	      for(int mu=0;mu<gamma1.size();++mu)
		tloopv[ss] = tloopv[ss] + Gamma(gamma1[mu]) * loopv[ss] * Gamma(gamma2[mu]);
	    }
	  });
	  LOG(Message) << "Sum over mu done with quark loop" << std::endl;
	}
	if ( type == 2 ) {
	  std::vector<Gamma::Algebra> gamma2 = gamma2_[ig];
	  thread_for(r,rd,{
	    int so=r*grid->_ostride[orthogdim];
	    for(int n=0;n<e1;n++)
	    for(int b=0;b<e2;b++){
	      int ss = so+n*stride+b;
	      tloopv[ss] = Zero();
	      for(int mu=0;mu<gamma2.size();++mu) {
		SpinColourMatrix_v tmp = Gamma(gamma2[mu])*loopv[ss];
		int s1 = int(mu/Ns);
		int s2 = mu%Ns;
		for(int c1=0;c1<Nc;++c1)
		for(int c2=0;c2<Nc;++c2){
		  tloopv[ss]()(s1,s2)(c1,c2)
		    = tmp()(0,0)(c1,c2)
		    + tmp()(1,1)(c1,c2)
		    + tmp()(2,2)(c1,c2)
		    + tmp()(3,3)(c1,c2);
		}
	      }
	    }
	  });
	  LOG(Message) << "Spin contraction with quark loop done for " << nameg2_[ig] << std::endl;
	}
	if ( type == 3 ) {
	  std::vector<Gamma::Algebra> gamma1 = gamma1_[ig];
	  std::vector<Gamma::Algebra> gamma2 = gamma2_[ig];
	  thread_for(r,rd,{
	    int so=r*grid->_ostride[orthogdim];
	    for(int n=0;n<e1;n++)
	    for(int b=0;b<e2;b++){
	      int ss = so+n*stride+b;
	      SpinMatrix_v tmp = Zero();
	      for(int s1=0;s1<Ns;++s1)
	      for(int s2=0;s2<Ns;++s2){
		tmp()(s1,s2)() = loopv[ss]()(s1,s2)(0,0) + loopv[ss]()(s1,s2)(1,1) + loopv[ss]()(s1,s2)(2,2);
	      }
	      tloopv[ss] = Zero();
	      for(int mu=0;mu<gamma1.size();++mu) {
		SpinMatrix_v tmp2 = Gamma(gamma1[mu]) * tmp * Gamma(gamma2[mu]);
		for(int s1=0;s1<Ns;++s1)
		for(int s2=0;s2<Ns;++s2){
		  tloopv[ss]()(s1,s2)(0,0) = tloopv[ss]()(s1,s2)(0,0) + tmp2()(s1,s2)();
		}
	      }
	    }
	  });
	  LOG(Message) << "Color contraction for quark loop and corresponding sum over mu done" << std::endl;
	}
	A2AMatrixSet<HADRONS_A2AM_IO_TYPE> emf(mBuf.data(),1,1,nt,N_i,N_j);
	std::vector<std::vector<SpinColourVector_v> > loopRight(N_j,zeroVec);
	LOG(Message) << "Making loop x right vec" << std::endl;
	for ( unsigned int j = 0; j < N_j; j++ ){
	  // Left x Loop
	  autoView(right_v,right[j],CpuRead);

	  thread_for(r,rd,{
	    int so=r*grid->_ostride[orthogdim];
	    for(int n=0;n<e1;n++)
	    for(int b=0;b<e2;b++){
	      int ss = so+n*stride+b;
	      addLoopRight(loopRight[j][ss],tloopv[ss],right_v[ss],gamma1_[ig],gamma2_[ig]);
	    }
	  });
	}
	LOG(Message) << "Done" << std::endl;
	LOG(Message) << "Making EMF" << std::endl;
	for ( unsigned int i = 0; i < N_i; i += cacheBlock )
	for ( unsigned int j = 0; j < N_j; j += cacheBlock ){
	  int Nii = MIN(N_i-i,cacheBlock);
	  int Njj = MIN(N_j-j,cacheBlock);

	  int MFrvol = rd * Nii * Njj;
	  int MFlvol = ld * Nii * Njj;

	  Vector<Scalar_v > lvSum(MFrvol);
	  thread_for( r, MFrvol,{
	    lvSum[r] = Zero();
	  });
	  // potentially wasting cores here if local time extent too small
	  thread_for(r,rd,{
	    int so=r*grid->_ostride[orthogdim]; // base offset for start of plane
	    int base=Nii*Njj*r;
	    for(int n=0;n<e1;n++)
	    for(int b=0;b<e2;b++){
	      int ss= so+n*stride+b;
	      for(int ii=0;ii<Nii;ii++)
	      for(int jj=0;jj<Njj;jj++){
		int idx = jj+Njj*ii+base;
		SpinMatrix_v tmp = Zero();
		for(int c=0;c<Nc;++c)
		for(int s=0;s<Ns;++s){
		  lvSum[idx]()()() += leftv[i+ii][ss]()(s)(c) * loopRight[j+jj][ss]()(s)(c);
		}
	      }
	    }
	  });

	  Vector<Scalar_s > lsSum(MFlvol);
	  thread_for(r,MFlvol,{
	    lsSum[r]=scalar_type(0.0);
	  });
	  thread_for(rt,rd,{
	    Coordinate icoor(Nd);
	    ExtractBuffer<Scalar_s> extracted(Nsimd);
	    for(int ii=0;ii<Nii;ii++)
	    for(int jj=0;jj<Njj;jj++){
	      int ij_rdx = jj+Njj*(ii+Nii*rt);
	      extract(lvSum[ij_rdx],extracted);
	      for(int idx=0;idx<Nsimd;idx++){
		grid->iCoorFromIindex(icoor,idx);
		int ldx    = rt+icoor[orthogdim]*rd;
		int ij_ldx = jj+Njj*(ii+Nii*ldx);
		lsSum[ij_ldx]=lsSum[ij_ldx]+extracted[idx];
	      }
	    }
	  });

	  int pd = grid->_processors[orthogdim];
	  int pc = grid->_processor_coor[orthogdim];
	  A2AMatrixSet<Complex> emfCache(mCacheBuf.data(),1,1,nt,Nii,Njj);
	  thread_for_collapse(2,lt,ld,{
	    for(int pt=0;pt<pd;pt++){
	      int t = lt + pt*ld;
	      if (pt == pc){
		for(int ii=0;ii<Nii;ii++)
		for(int jj=0;jj<Njj;jj++){
		  int ij_dx = jj+Njj*(ii+Nii*lt);//m+Nmom*i + Nmom*cacheBlock * j + Nmom*lock * Rblock * lt;
		  emfCache(0,0,t,ii,jj) = lsSum[ij_dx]()()();
		}
	      } else {
		const scalar_type zz(0.0);
		for(int ii=0;ii<Nii;ii++)
		for(int jj=0;jj<Njj;jj++){
		  emfCache(0,0,t,ii,jj) = zz;
		}
	      }
	    }
	  });
	  grid->GlobalSumVector(&emfCache(0,0,0,0,0),nt*Nii*Njj);
	  thread_for_collapse( 5,e,1,{
	    for(int s =0;s< 1;s++)
	    for(int t =0;t< nt;t++)
	    for(int ii=0;ii< Nii;ii++)
	    for(int jj=0;jj< Njj;jj++)
	    {
	      emf(0,0,t,i+ii,j+jj) = emfCache(0,0,t,ii,jj);
	    }
	  });
	}// i,j
	LOG(Message) << "Done" << std::endl;

        double       blockSize, ioTime;

	std::string ioname = "type" + std::to_string(type) + "_" + nameg1_[ig] + "_" + nameg2_[ig];
	std::string filename = par().output + "." + std::to_string(vm().getTrajectory()) + "/" + ioname + ".h5";
        LOG(Message) << "Writing block to " << filename << std::endl;
        makeFileDir(filename, grid);

        //ioTime = -GET_TIMER("IO: write block");
        //START_TIMER("IO: total");
#ifdef HADRONS_A2AM_PARALLEL_IO
        grid->Barrier();
        LOG(Message) << "HADRONS_A2AM_PARALLEL_IO" << std::endl;
	if ( grid->ThisRank() == 0 ) {
#endif
	  // make task list for current node
	  A2AMatrixIo<HADRONS_A2AM_IO_TYPE> io(filename,ioname,nt, N_i, N_j);
	  A2AExtendedMesonFieldMetadata md;
	  md.gamma1 = nameg1_[ig];
	  md.gamma2 = nameg2_[ig];

	  // memory consuming
	  io.initFile(md, MAX(N_i,N_j));
	  //START_TIMER("IO: write block");
	  io.saveBlock(emf, 0, 0, 0, 0);
	  //STOP_TIMER("IO: write block");

#ifdef HADRONS_A2AM_PARALLEL_IO
	}
	grid->Barrier();
#endif
        //STOP_TIMER("IO: total");
      }// ig
    }// type

}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_A2AExtendedMesonField_hpp_
