/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/Modules/MContraction/A2AChromoMagneticOperatorField.hpp

Copyright (C) 2015-2019

Author: Masaaki Tomii <masaaki.tomii@uconn.edu>
*************************************************************************************/
/*  END LEGAL */
#ifndef Hadrons_MContraction_A2AChromoMagneticOperatorField_hpp_
#define Hadrons_MContraction_A2AChromoMagneticOperatorField_hpp_

//#include <Grid/qcd/utils/WilsonLoops.h>
#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/A2AMatrix.hpp>

BEGIN_HADRONS_NAMESPACE
/******************************************************************************
 *                All-to-all extended meson field creation                    *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MContraction)

class A2AChromoMagneticOperatorFieldPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(A2AChromoMagneticOperatorFieldPar,
                                    int, cacheBlock,
				    std::string, parities,
                                    std::string, left,
                                    std::string, right,
				    std::string, gauge,
                                    std::string, output,
				    std::string, ifOrthogs);
};

class A2AChromoMagneticOperatorFieldMetadata: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(A2AChromoMagneticOperatorFieldMetadata,
                                    std::string, meta);
};

template <typename GImpl, typename FImpl>
class TA2AChromoMagneticOperatorField : public Module<A2AChromoMagneticOperatorFieldPar>
{
public:
  //GAUGE_TYPE_ALIASES(GImpl,);
  typedef typename GImpl::GaugeLinkField GaugeMat;
  typedef typename GImpl::SiteGaugeLink SiteGaugeLink;
    FERM_TYPE_ALIASES(FImpl,);
  //typedef typename FImpl::FermionField FermionField;
  typedef typename FImpl::SiteSpinor vobj;
  typedef typename vobj::scalar_object sobj;
  typedef typename vobj::scalar_type scalar_type;
  typedef typename vobj::vector_type vector_type;
public:
    // constructor
    TA2AChromoMagneticOperatorField(const std::string name);
    // destructor
    virtual ~TA2AChromoMagneticOperatorField(void){};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
  std::vector<int> parities_;
  std::vector<int> ifOrthogs_;
};

MODULE_REGISTER(A2AChromoMagneticOperatorField, ARG(TA2AChromoMagneticOperatorField<GIMPL,FIMPL>), MContraction);

/******************************************************************************
*               TA2AChromoMagneticOperatorField implementation                *
******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename GImpl, typename FImpl>
TA2AChromoMagneticOperatorField<GImpl,FImpl>::TA2AChromoMagneticOperatorField(const std::string name)
: Module<A2AChromoMagneticOperatorFieldPar>(name)
{
}

// dependencies/products ///////////////////////////////////////////////////////
// Modification needed??
template <typename GImpl, typename FImpl>
std::vector<std::string> TA2AChromoMagneticOperatorField<GImpl,FImpl>::getInput(void)
{
  std::vector<std::string> in = {par().left, par().right, par().gauge};
  return in;
}

template <typename GImpl, typename FImpl>
std::vector<std::string> TA2AChromoMagneticOperatorField<GImpl,FImpl>::getOutput(void)
{
  std::vector<std::string> out = {};
  return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename GImpl, typename FImpl>
void TA2AChromoMagneticOperatorField<GImpl,FImpl>::setup(void)
{
  parities_ = strToVec<int>(par().parities);
  ifOrthogs_ = strToVec<int>(par().ifOrthogs);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename GImpl, typename FImpl>
void TA2AChromoMagneticOperatorField<GImpl,FImpl>::execute(void)
{
  typedef iSpinColourVector<vector_type> SpinColourVector_v;
  typedef iSpinColourMatrix<vector_type> SpinColourMatrix_v;
  typedef iColourMatrix<vector_type> ColourMatrix_v;
  typedef iSpinMatrix<vector_type> SpinMatrix_v;
  typedef iSinglet<vector_type> Scalar_v;
  typedef iSinglet<scalar_type> Scalar_s;

  auto &left  = envGet(std::vector<FermionField>, par().left);
  auto &right = envGet(std::vector<FermionField>, par().right);
  const auto &U = envGet(GaugeField, par().gauge);

  //std::vector<std::string> gname = strToVec<int>(par().parities);
  GridBase *grid = left[0].Grid();

  int orthogdim = 3;
  int rd=grid->_rdimensions[orthogdim];//2
  int ld=grid->_ldimensions[orthogdim];
  int Nd=grid->_ndimension;
  int Nsimd=grid->Nsimd();
  int e1=    grid->_slice_nblock[orthogdim];//1
  int e2=    grid->_slice_block [orthogdim];//64 must be 4^3
  int stride=grid->_slice_stride[orthogdim];//128
  LOG(Message) << "Computing all-to-all ChromoMagnetic operator fields" << std::endl;
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

  LOG(Message) << "Left: '" << par().left << "' Right: '"
	       << par().right << "'" << std::endl;
  LOG(Message) << "Gauge field: '" << par().gauge << "'" << std::endl;
  LOG(Message) << "Parities:" << std::endl;
  for (auto &p: parities_)
  {
    LOG(Message) << "  " << p << std::endl;
  }
  LOG(Message) << "CMO field size: " << nt << "*" << N_i << "*" << N_j 
	       << " (filesize " << sizeString(nt*N_i*N_j*sizeof(HADRONS_A2AM_IO_TYPE)) << std::endl;

  //GaugeMat F(U.Grid());
  //std::vector<std::vector<ColourMatrix_v> > G;
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

  int ig=0;
  for (auto &ifOrthog: ifOrthogs_){
    std::vector<GaugeMat> G;
    std::vector<Gamma::Algebra> Sigma;
    for(int mu=0;mu<3;mu++){
      if( mu == orthogdim && ifOrthog == 1 ) continue;
      GaugeMat Umu(U.Grid());
      Umu = peekLorentz(U,mu);
      for(int nu=mu+1;nu<4;nu++){
	if( nu == orthogdim && ifOrthog == 1 ) continue;
	if( ifOrthog == 0 ){
	  if( mu!=orthogdim && nu!=orthogdim ) continue;
	}
	GaugeMat Unu(U.Grid());
	Unu = peekLorentz(U,nu);
	GaugeMat Gmunu(U.Grid());
	WilsonLoops<GImpl>::CloverleafMxN(Gmunu,Umu,Unu,mu,nu,1,1);
	G.push_back(Gmunu - adj(Gmunu));

	if( mu == 0 && nu == 1 ){
	  Sigma.push_back(Gamma::Algebra::SigmaXY);
	} else if ( mu == 0 && nu == 2 ) {
	  Sigma.push_back(Gamma::Algebra::SigmaXZ);
	} else if ( mu == 0 && nu == 3 ) {
	  Sigma.push_back(Gamma::Algebra::SigmaXT);
	} else if ( mu == 1 && nu == 2 ) {
	  Sigma.push_back(Gamma::Algebra::SigmaYZ);
	} else if ( mu == 1 && nu == 3 ) {
	  Sigma.push_back(Gamma::Algebra::SigmaYT);
	} else if ( mu == 2 && nu == 3 ) {
	  Sigma.push_back(Gamma::Algebra::SigmaZT);
	}
      }
    }
    assert( Sigma.size() == 3 );

    for (auto &parity: parities_){
      LOG(Message) << "Starting calculation with parity" << parity << std::endl;
      A2AMatrixSet<HADRONS_A2AM_IO_TYPE> cmf(mBuf.data(),1,1,nt,N_i,N_j);
      for ( unsigned int j = 0; j < N_j; j += cacheBlock ){
	int Njj = MIN(N_j-j,cacheBlock);
	std::vector<std::vector<SpinColourVector_v> > rightv(Njj,zeroVec);
	for ( unsigned int jj = 0; jj < Njj; ++jj ){
	  //auto viewme = PeriodicBC::CovShiftForward<vobj, SiteGaugeLink>
	  for ( int munu=0;munu<G.size();++munu){
	    autoView(tmp,right[j+jj],CpuRead);
	    autoView(Gmunu,G[munu],CpuRead);
	    //tmp = Gamma(Sigma[munu]) * tmp;
	    //if ( parity == 1 ) tmp = Gamma(Gamma::Algebra::Gamma5) * tmp;
	    thread_for(r,rd,{
	      int so=r*grid->_ostride[orthogdim];
	      for(int n=0;n<e1;n++)
	      for(int b=0;b<e2;b++){
		int ss = so+n*stride+b;
		SpinColourVector_v  tmp2 = Gamma(Sigma[munu]) * tmp[ss];
		if ( parity == 1 ) tmp2 = Gamma(Gamma::Algebra::Gamma5) * tmp2;
		for(int s=0;s<Ns;++s)
		for(int c=0;c<Nc;++c){
		  rightv[jj][ss]()(s)(c)
		    += Gmunu[ss]()()(c,0) * tmp2()(s)(0)
		    +  Gmunu[ss]()()(c,1) * tmp2()(s)(1)
		    +  Gmunu[ss]()()(c,2) * tmp2()(s)(2);
		}
	      }
	    });
	  }
	}
	//LOG(Message) << "Memory for right vectors allocated" << std::endl;

	for ( unsigned int i = 0; i < N_i; i += cacheBlock ){
	  int Nii = MIN(N_i-i,cacheBlock);

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
		  lvSum[idx]()()() += leftv[i+ii][ss]()(s)(c) * rightv[jj][ss]()(s)(c);
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

	  A2AMatrixSet<Complex> cmfCache(mCacheBuf.data(),1,1,nt,Nii,Njj);
	  int pd = grid->_processors[orthogdim];
	  int pc = grid->_processor_coor[orthogdim];
	  thread_for_collapse(2,lt,ld,{
	    for(int pt=0;pt<pd;pt++){
	      int t = lt + pt*ld;
	      if (pt == pc){
		for(int ii=0;ii<Nii;ii++)
		for(int jj=0;jj<Njj;jj++){
		  int ij_dx = jj+Njj*(ii+Nii*lt);//m+Nmom*i + Nmom*cacheBlock * j + Nmom*lock * Rblock * lt;
		  cmfCache(0,0,t,ii,jj) = lsSum[ij_dx]()()();
		}
	      } else {
		const scalar_type zz(0.0);
		for(int ii=0;ii<Nii;ii++)
		for(int jj=0;jj<Njj;jj++){
		  cmfCache(0,0,t,ii,jj) = zz;
		}
	      }
	    }
	  });
	  grid->GlobalSumVector(&cmfCache(0,0,0,0,0),nt*Nii*Njj);
	  thread_for_collapse( 5,e,1,{
	    for(int s =0;s< 1;s++)
	    for(int t =0;t< nt;t++)
	    for(int ii=0;ii< Nii;ii++)
	    for(int jj=0;jj< Njj;jj++)
	    {
	      cmf(0,0,t,i+ii,j+jj) = cmfCache(0,0,t,ii,jj);
	    }
	  });
	}// i
      }// j
      LOG(Message) << "Done" << std::endl;

      double       blockSize, ioTime;

      //ioTime = -GET_TIMER("IO: write block");
      //START_TIMER("IO: total");
#ifdef HADRONS_A2AM_PARALLEL_IO
      grid->Barrier();
      if ( grid->ThisRank() == 0 ) {
#endif
    // make task list for current node
	std::string ioname = "parity" + std::to_string(parity);
	if ( ifOrthog == 1 ) {
	  ioname = ioname + "_GijSij";
	} else {
	  ioname = ioname + "_GitSit";
	}
	std::string filename = par().output + "." + std::to_string(vm().getTrajectory()) + "/" + ioname + ".h5";
	LOG(Message) << "Writing block to " << filename << std::endl;
	makeFileDir(filename, grid);
	A2AMatrixIo<HADRONS_A2AM_IO_TYPE> io(filename,ioname,nt, N_i, N_j);
	A2AChromoMagneticOperatorFieldMetadata md;
	md.meta = ioname;

	// memory consuming
	io.initFile(md, MAX(N_i,N_j));
	//START_TIMER("IO: write block");
	io.saveBlock(cmf, 0, 0, 0, 0);
	//STOP_TIMER("IO: write block");

#ifdef HADRONS_A2AM_PARALLEL_IO
      }
      grid->Barrier();
#endif
      //STOP_TIMER("IO: total");
      ++ig;
    }// parity
  }// ifOrthog
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_A2AChromoMagneticOperatorField_hpp_
