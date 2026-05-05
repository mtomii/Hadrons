/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/Modules/MContraction/A2AFewMesonField.hpp

Copyright (C) 2015-2019

Author: Masaaki Tomii <masaaki.tomii@uconn.edu>
*************************************************************************************/
/*  END LEGAL */
#ifndef Hadrons_MContraction_A2AFewMesonField_hpp_
#define Hadrons_MContraction_A2AFewMesonField_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/A2AMatrix.hpp>

BEGIN_HADRONS_NAMESPACE
/******************************************************************************
 *                All-to-all extended meson field creation                    *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MContraction)

class A2AFewMesonFieldPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(A2AFewMesonFieldPar,
                                    int, cacheBlock,
                                    int, block,
                                    std::string, left,
                                    std::string, right,
                                    std::string, output,
                                    std::string, gammas,
				    std::vector<std::string>, mom);
};

class A2AFewMesonFieldMetadata: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(A2AFewMesonFieldMetadata,
                                    std::vector<RealF>, momentum,
                                    Gamma::Algebra, gamma);
};

template <typename FImpl>
class TA2AFewMesonField : public Module<A2AFewMesonFieldPar>
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
    TA2AFewMesonField(const std::string name);
    // destructor
    virtual ~TA2AFewMesonField(void){};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
    bool                               hasPhase_{false};
    std::string                        momphName_;
    std::vector<Gamma::Algebra> gamma_;
    std::vector<std::vector<Real>>     mom_;
  //const std::vector<LatticeComplex> &mom_;
};

MODULE_REGISTER(A2AFewMesonField, TA2AFewMesonField<FIMPL>, MContraction);

/******************************************************************************
*                     TA2AFewMesonField implementation                        *
******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TA2AFewMesonField<FImpl>::TA2AFewMesonField(const std::string name)
: Module<A2AFewMesonFieldPar>(name)
, momphName_(name + "_momph")
{
}

// dependencies/products ///////////////////////////////////////////////////////
// Modification needed??
template <typename FImpl>
std::vector<std::string> TA2AFewMesonField<FImpl>::getInput(void)
{
  std::vector<std::string> in = {par().left, par().right};
  return in;
}

template <typename FImpl>
std::vector<std::string> TA2AFewMesonField<FImpl>::getOutput(void)
{
  std::vector<std::string> out = {};
  return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TA2AFewMesonField<FImpl>::setup(void)
{
    gamma_.clear();
    mom_.clear();
    if (par().gammas == "all")
    {
        gamma_ = {
            Gamma::Algebra::Gamma5,
            Gamma::Algebra::Identity,    
            Gamma::Algebra::GammaX,
            Gamma::Algebra::GammaY,
            Gamma::Algebra::GammaZ,
            Gamma::Algebra::GammaT,
            Gamma::Algebra::GammaXGamma5,
            Gamma::Algebra::GammaYGamma5,
            Gamma::Algebra::GammaZGamma5,
            Gamma::Algebra::GammaTGamma5,
            Gamma::Algebra::SigmaXY,
            Gamma::Algebra::SigmaXZ,
            Gamma::Algebra::SigmaXT,
            Gamma::Algebra::SigmaYZ,
            Gamma::Algebra::SigmaYT,
            Gamma::Algebra::SigmaZT
        };
    }
    else
    {
        gamma_ = strToVec<Gamma::Algebra>(par().gammas);
    }
    for (auto &pstr: par().mom)
    {
        auto p = strToVec<Real>(pstr);

        if (p.size() != env().getNd() - 1)
        {
            HADRONS_ERROR(Size, "Momentum has " + std::to_string(p.size())
                                + " components instead of " 
                                + std::to_string(env().getNd() - 1));
        }
        mom_.push_back(p);
    }
    envCache(std::vector<ComplexField>, momphName_, 1, 
             par().mom.size(), envGetGrid(ComplexField));
    envTmpLat(ComplexField, "coor");
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TA2AFewMesonField<FImpl>::execute(void)
{
  typedef iSpinColourVector<vector_type> SpinColourVector_v;
  typedef iSpinColourMatrix<vector_type> SpinColourMatrix_v;
  typedef iSpinMatrix<vector_type> SpinMatrix_v;
  typedef iSinglet<vector_type> Scalar_v;
  typedef iSinglet<scalar_type> Scalar_s;
  auto &left  = envGet(std::vector<FermionField>, par().left);
  auto &right = envGet(std::vector<FermionField>, par().right);

  GridBase *grid = left[0].Grid();

  int orthogdim = 3;
  int rd=grid->_rdimensions[orthogdim];//2
  int ld=grid->_ldimensions[orthogdim];
  int Nd=grid->_ndimension;
  int Nsimd=grid->Nsimd();
  int e1=    grid->_slice_nblock[orthogdim];//1
  int e2=    grid->_slice_block [orthogdim];//64 must be 4^3
  int stride=grid->_slice_stride[orthogdim];//128
  LOG(Message) << "Computing all-to-all meson fields" << std::endl;
  LOG(Message) << "R dimension: " << rd << std::endl;
  LOG(Message) << "Slice nblock: " << e1 << std::endl;
  LOG(Message) << "Slice block: " << e2 << std::endl;
  LOG(Message) << "Slice stride: " << stride << std::endl;
  LOG(Message) << "O stride: " << grid->_ostride[orthogdim] << std::endl;
  LOG(Message) << "Nsimd: " << Nsimd << std::endl;
  grid->show_decomposition();

  int nt         = env().getDim().back();
  int N_i        = left.size();
  int N_j        = right.size();
  int ngamma     = gamma_.size();
  int nmom       = mom_.size();
  int block      = par().block;
  int cacheBlock = par().cacheBlock;
  Vector<HADRONS_A2AM_IO_TYPE> mBuf; mBuf.resize(nmom*nt*block*block);
  Vector<Complex> mCacheBuf; mCacheBuf.resize(nmom*nt*cacheBlock*cacheBlock);

  LOG(Message) << "Left: '" << par().left << "' Right: '"
	       << par().right << "'" << std::endl;
    LOG(Message) << "Momenta:" << std::endl;
    for (auto &p: mom_)
    {
        LOG(Message) << "  " << p << std::endl;
    }
    LOG(Message) << "Spin bilinears:" << std::endl;
    for (auto &g: gamma_)
    {
        LOG(Message) << "  " << g << std::endl;
    }
  LOG(Message) << "Meson field size: " << nt << "*" << N_i << "*" << N_j 
	       << " (filesize " << sizeString(nt*N_i*N_j*sizeof(HADRONS_A2AM_IO_TYPE)) << ")" << std::endl;

  auto &ph = envGet(std::vector<ComplexField>, momphName_);
  LOG(Message) << "Phase object gotten" << std::endl;
  if (!hasPhase_)
  {
    startTimer("Momentum phases");
    for (unsigned int j = 0; j < nmom; ++j)
    {
      Complex           i(0.0,1.0);
      std::vector<Real> p;

      envGetTmp(ComplexField, coor);
      ph[j] = Zero();
      for(unsigned int mu = 0; mu < mom_[j].size(); mu++)
      {
	LatticeCoordinate(coor, mu);
	ph[j] = ph[j] + (mom_[j][mu]/env().getDim(mu))*coor;
      }
      ph[j] = exp((Real)(2*M_PI)*i*ph[j]);
    }
    hasPhase_ = true;
    stopTimer("Momentum phases");
    LOG(Message) << "Momentum phases set" << std::endl;
  }

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

  for (auto &gamma: gamma_){
    LOG(Message) << "Starting calculation with: " << gamma << std::endl;
    for ( unsigned int j = 0; j < N_j; j += block ){
      int Njj = MIN(N_j-j,block);
      std::vector<std::vector<SpinColourVector_v> > rightv(Njj,zeroVec);
      for ( unsigned int jj = 0; jj < Njj; jj++ ){
	autoView(tmp,right[j+jj],CpuRead);
	thread_for(r,rd,{
	  int so=r*grid->_ostride[orthogdim];
	  for(int n=0;n<e1;n++)
	  for(int b=0;b<e2;b++){
	    int ss = so+n*stride+b;
	    rightv[jj][ss] = Gamma(gamma) * tmp[ss];
	  }
	});
      }

      for ( unsigned int i = 0; i < N_i; i += block ){
	int Nii = MIN(N_i-i,block);
	A2AMatrixSet<HADRONS_A2AM_IO_TYPE> mf(mBuf.data(),nmom,1,nt,Nii,Njj);

	for ( unsigned int ii = 0; ii < Nii; ii += cacheBlock )
	for ( unsigned int jj = 0; jj < Njj; jj += cacheBlock ){
	  int Niii = MIN(Nii-ii,cacheBlock);
	  int Njjj = MIN(Njj-jj,cacheBlock);
	  int MFrvol = rd * nmom * Niii * Njjj;
	  int MFlvol = ld * nmom * Niii * Njjj;

	  A2AMatrixSet<Complex> mfCache(mCacheBuf.data(),nmom,1,nt,Niii,Njjj);
	  Vector<Scalar_v > lvSum(MFrvol);
	  thread_for( r, MFrvol,{
	    lvSum[r] = Zero();
	  });
	  
	  // potentially wasting cores here if local time extent too small
	  thread_for(r,rd,{
	    int so=r*grid->_ostride[orthogdim]; // base offset for start of plane
	    int base=nmom*Niii*Njjj*r;
	    for(int n=0;n<e1;n++)
	    for(int b=0;b<e2;b++){
	      int ss= so+n*stride+b;

	      for(int iii=0;iii<Niii;iii++)
	      for(int jjj=0;jjj<Njjj;jjj++){
		Scalar_v tmp = Zero();
		tmp()()()
		  = leftv[i+ii+iii][ss]()(0)(0) * rightv[jj+jjj][ss]()(0)(0)
		  + leftv[i+ii+iii][ss]()(0)(1) * rightv[jj+jjj][ss]()(0)(1)
		  + leftv[i+ii+iii][ss]()(0)(2) * rightv[jj+jjj][ss]()(0)(2)
		  + leftv[i+ii+iii][ss]()(1)(0) * rightv[jj+jjj][ss]()(1)(0)
		  + leftv[i+ii+iii][ss]()(1)(1) * rightv[jj+jjj][ss]()(1)(1)
		  + leftv[i+ii+iii][ss]()(1)(2) * rightv[jj+jjj][ss]()(1)(2)
		  + leftv[i+ii+iii][ss]()(2)(0) * rightv[jj+jjj][ss]()(2)(0)
		  + leftv[i+ii+iii][ss]()(2)(1) * rightv[jj+jjj][ss]()(2)(1)
		  + leftv[i+ii+iii][ss]()(2)(2) * rightv[jj+jjj][ss]()(2)(2)
		  + leftv[i+ii+iii][ss]()(3)(0) * rightv[jj+jjj][ss]()(3)(0)
		  + leftv[i+ii+iii][ss]()(3)(1) * rightv[jj+jjj][ss]()(3)(1)
		  + leftv[i+ii+iii][ss]()(3)(2) * rightv[jj+jjj][ss]()(3)(2);
		for (int m=0;m<nmom;++m){
		  int idx = m+nmom*(jjj+Njjj*(iii+Niii*r));
		  autoView(mom_v,ph[m],CpuRead);
		  auto phase = mom_v[ss];
		  mac(&lvSum[idx],&tmp,&phase);
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
	    for(int iii=0;iii<Niii;iii++)
	    for(int jjj=0;jjj<Njjj;jjj++)
	    for(int m=0;m<nmom;++m){
	      int ij_rdx = m+nmom*(jjj+Njjj*(iii+Niii*rt));
	      extract(lvSum[ij_rdx],extracted);
	      for(int idx=0;idx<Nsimd;idx++){
		grid->iCoorFromIindex(icoor,idx);
		int ldx    = rt+icoor[orthogdim]*rd;
		int ij_ldx = m+nmom*(jjj+Njjj*(iii+Niii*ldx));
		lsSum[ij_ldx]=lsSum[ij_ldx]+extracted[idx];
	      }
	    }
	  });

	  int pd = grid->_processors[orthogdim];
	  int pc = grid->_processor_coor[orthogdim];
	  thread_for_collapse(2,lt,ld,{
	    for(int pt=0;pt<pd;pt++){
	      int t = lt + pt*ld;
	      if (pt == pc){
		for(int iii=0;iii<Niii;iii++)
		for(int jjj=0;jjj<Njjj;jjj++)
		for(int m=0;m<nmom;++m){
		  int ij_dx = m+nmom*(jjj+Njjj*(iii+Niii*lt));
		  mfCache(m,0,t,iii,jjj) = lsSum[ij_dx]()()();
		}
	      } else {
		const scalar_type zz(0.0);
		for(int iii=0;iii<Niii;iii++)
		for(int jjj=0;jjj<Njjj;jjj++)
		for(int m=0;m<nmom;++m){
		  mfCache(m,0,t,iii,jjj) = zz;
		}
	      }
	    }
	  });
	  grid->GlobalSumVector(&mfCache(0,0,0,0,0),nmom*nt*Niii*Njjj);
	  //grid->GlobalSumVector(mfCache.data(),nmom*nt*Niii*Njjj);

	  thread_for_collapse( 5,e,nmom,{
            for(int s =0;s< 1;s++)
            for(int t =0;t< nt;t++)
            for(int iii=0;iii< Niii;iii++)
            for(int jjj=0;jjj< Njjj;jjj++)
            {
	      mf(e,s,t,ii+iii,jj+jjj) = mfCache(e,s,t,iii,jjj);
	    }
          });
	}

	{
	  std::string filename = par().output + "." + std::to_string(vm().getTrajectory()) + "/" + "mkdir" + ".h5";
	  makeFileDir(filename, grid);
	}
#ifdef HADRONS_A2AM_PARALLEL_IO
	grid->Barrier();
        unsigned int myRank = grid->ThisRank(), nRank  = grid->RankCount();
        for(int f = myRank; f < nmom; f += nRank){
	  std::stringstream ss;
	  ss << gamma << "_";
	  for(unsigned int mu=0; mu<mom_[f].size(); ++mu){
	    ss << mom_[f][mu] << ((mu == mom_[f].size() - 1) ? "" : "_");
	  }
	  std::string ioname = ss.str();
	  std::string filename = par().output + "." + std::to_string(vm().getTrajectory()) + "/" + ioname + ".h5";
	  A2AMatrixIo<HADRONS_A2AM_IO_TYPE> io(filename,ioname,nt, N_i, N_j);
	  A2AFewMesonFieldMetadata md;
	  LOG(Message) << "Writing block to " << filename << std::endl;
	  for (auto pmu: mom_[f])
	  {
            md.momentum.push_back(pmu);
	  }
	  md.gamma = gamma;
	  if ( i+j==0 ) io.initFile(md, block);

	  io.saveBlock(mf, f, 0, i, j);
	}
	grid->Barrier();
#endif
      }// i
    }// j
  }// gamma
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_A2AFewMesonField_hpp_
