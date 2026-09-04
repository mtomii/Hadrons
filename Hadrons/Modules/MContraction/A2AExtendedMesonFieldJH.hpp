/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Grid/algorithms/multigrid/GeneralCoarsenedMatrix.hGrid/algorithms/multigrid/GeneralCoarsenedMatrix.hSource file: Hadrons/Modules/MContraction/A2AExtendedMesonFieldJH.hpp

Copyright (C) 2015-2019

Author: Peter Boyle <paboyle@bnl.gov>
Author: Jonas Hildebrand <jonas.hildebrand@uconn.edu>
Author: Masaaki Tomii <masaaki.tomii@uconn.edu>
*************************************************************************************/
/*  END LEGAL */
#ifndef Hadrons_MContraction_A2AExtendedMesonFieldJH_hpp_
#define Hadrons_MContraction_A2AExtendedMesonFieldJH_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/A2AMatrix.hpp>
#include <Grid/qcd/utils/A2Autils.h>

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

class A2AExtendedMesonFieldJHPar: Serializable
{
public:
    // The quark loop arrives one of two mutually exclusive ways:
    //
    //   loop                name of a precomputed PropagatorField - from
    //                       MContraction::A2ALoopNew in this job, or read off
    //                       disk with MIO::LoadProp. This is the production
    //                       path: the module never sees modes, hit counts or
    //                       W representations, and the loop's vector arrays
    //                       are freed before this module runs.
    //
    //   loop_vw1, loop_vw2  names of the loop's A2A vector arrays, built
    //                       in-module. TEMPORARY - expanded W only, plain
    //                       index-for-index pairing, no dense-W handling. It
    //                       exists solely as the reference the loop input is
    //                       verified against (gen_verify_loop.py rung 1) and
    //                       is to be deleted once that passes. It cannot fit
    //                       at production hit counts anyway, since both mode
    //                       arrays must stay resident alongside the external
    //                       legs for this module's whole execution.
    //
    // setup() rejects neither and both.
    GRID_SERIALIZABLE_CLASS_MEMBERS(A2AExtendedMesonFieldJHPar,
                                    int, block,
                                    int, cacheBlock,
				    std::string, types,
                                    std::string, left,
                                    std::string, right,
                                    std::string, loop,
				    std::string, loop_vw1,
				    std::string, loop_vw2,
                                    std::string, output,
                                    std::string, gammas1,
				    std::string, gammas2);
};

class A2AExtendedMesonFieldJHMetadata: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(A2AExtendedMesonFieldJHMetadata,
                                    std::string, gamma1,
				    std::string, gamma2);
};

template <typename FImpl>
class TA2AExtendedMesonFieldJH : public Module<A2AExtendedMesonFieldJHPar>
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
    TA2AExtendedMesonFieldJH(const std::string name);
    // destructor
    virtual ~TA2AExtendedMesonFieldJH(void){};
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

MODULE_REGISTER(A2AExtendedMesonFieldJH, TA2AExtendedMesonFieldJH<FIMPL>, MContraction);

/******************************************************************************
*               TA2AExtendedMesonFieldJH implementation                       *
******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TA2AExtendedMesonFieldJH<FImpl>::TA2AExtendedMesonFieldJH(const std::string name)
: Module<A2AExtendedMesonFieldJHPar>(name)
{
}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TA2AExtendedMesonFieldJH<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().left, par().right};

    // Only the loop source actually in use is declared, so the other one's
    // objects are never kept alive on this module's account. With the loop
    // input that means the loop's mode arrays die with their producer.
    if (!par().loop.empty())
    {
        in.push_back(par().loop);
    }
    else
    {
        in.push_back(par().loop_vw1);
        in.push_back(par().loop_vw2);
    }

    return in;
}

template <typename FImpl>
std::vector<std::string> TA2AExtendedMesonFieldJH<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {};

    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TA2AExtendedMesonFieldJH<FImpl>::setup(void)
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
    bool haveLoop = !par().loop.empty();
    bool haveVecs = !par().loop_vw1.empty() || !par().loop_vw2.empty();

    if (haveLoop && haveVecs)
    {
        HADRONS_ERROR(Argument, "set either 'loop' or 'loop_vw1'/'loop_vw2', "
                                "not both");
    }
    if (!haveLoop && !(!par().loop_vw1.empty() && !par().loop_vw2.empty()))
    {
        HADRONS_ERROR(Argument, "no quark loop: set 'loop', or both of "
                                "'loop_vw1' and 'loop_vw2'");
    }

    // Only the in-module path needs somewhere to build the loop; with the
    // loop input the propagator is owned by whoever produced it.
    if (!haveLoop)
    {
        envCreateLat(PropagatorField, getName() + "_propagatorLoop");
    }
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TA2AExtendedMesonFieldJH<FImpl>::execute(void)
{
  typedef iSpinColourVector<vector_type> SpinColourVector_v;
    auto &left  = envGet(std::vector<FermionField>, par().left);
    auto &right = envGet(std::vector<FermionField>, par().right);

    GridBase *grid = left[0].Grid();
    bool haveLoop = !par().loop.empty();

    LOG(Message) << "Computing all-to-all EXTENDED meson fields" << std::endl;

    int nt         = env().getDim().back();
    int N_i        = left.size();
    int N_j        = right.size();
    int block = par().block;
    int cacheBlock = par().cacheBlock;
    Vector<HADRONS_A2AM_IO_TYPE> mBuf; mBuf.resize(nt*N_i*N_j);

    LOG(Message) << "Left: '" << par().left << "' Right: '"
		 << par().right << "'" << std::endl;
    if (haveLoop)
    {
        LOG(Message) << "Quark loop: '" << par().loop << "' (precomputed)"
                     << std::endl;
    }
    else
    {
        LOG(Message) << "A2AVectors for loop: '" << par().loop_vw1
                     << "' and '" << par().loop_vw2 << "'" << std::endl;
    }
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

    PropagatorField *loopPtr;

    if (haveLoop)
    {
        loopPtr = &envGet(PropagatorField, par().loop);
    }
    else
    {
        // Temporary reference path, see the note on the Par struct: expanded
        // W only, so the mode arrays pair index for index and the whole-array
        // LoopPropagator wrapper is all that is needed.
        auto &loop1 = envGet(std::vector<FermionField>, par().loop_vw1);
        auto &loop2 = envGet(std::vector<FermionField>, par().loop_vw2);

        if (loop1.size() != loop2.size())
        {
            HADRONS_ERROR(Size, "loop_vw1 has " + std::to_string(loop1.size())
                                + " modes, loop_vw2 has "
                                + std::to_string(loop2.size())
                                + "; the in-module path pairs index for index "
                                  "and needs the expanded W representation");
        }
        loopPtr = &envGet(PropagatorField, getName() + "_propagatorLoop");
        startTimer("LoopPropagator");
        Grid::A2AExtendedMesonField<FImpl>::LoopPropagator(*loopPtr, loop1, loop2);
        stopTimer("LoopPropagator");
        LOG(Message) << "Quark loop calculated" << std::endl;
    }

    auto &loop = *loopPtr;
    LOG(Message) << "Quark loop norm2 = " << norm2(loop) << std::endl;

    std::vector<FermionField> loopRight(block, grid);
    PropagatorField tloop(grid);

    std::array<double, 6> sumTimings = {};
    std::array<double, 6> sumBytes   = {};
    //std::array<double, 7> ioTimings  = {};

    for (int &type: types_){

      for (int ig = 0 ; ig < gamma1_.size() ; ++ig ){

	A2AMatrixSet<HADRONS_A2AM_IO_TYPE> emf(mBuf.data(),1,1,nt,N_i,N_j);

	Vector<Gamma::Algebra> gamma1(gamma1_[ig].begin(), gamma1_[ig].end());
	Vector<Gamma::Algebra> gamma2(gamma2_[ig].begin(), gamma2_[ig].end());

	startTimer("Loop contraction");
	tloop = Zero();
	switch (type) {
	case 0: Grid::A2AExtendedMesonField<FImpl>::LoopContractionType0(tloop, loop);                 break;
	case 1: Grid::A2AExtendedMesonField<FImpl>::LoopContractionType1(tloop, loop, gamma1, gamma2); break;
	case 2: Grid::A2AExtendedMesonField<FImpl>::LoopContractionType2(tloop, loop, gamma2);         break;
	case 3: Grid::A2AExtendedMesonField<FImpl>::LoopContractionType3(tloop, loop, gamma1, gamma2); break;
	}
	stopTimer("Loop contraction");
	//LOG(Message) << "tloop contraction done for type " << type << std::endl;
	//LOG(Message) << "Making EMF" << std::endl;
	A2ASpatialSum<SpinColourVector_v> spatial_sum;

	// Result buffers, one per distinct block shape. A block is full or on the
	// tail in each axis independently, so a 2x2 pool indexed by (i on tail, j
	// on tail) covers every case. Each slot is asked for the same dimensions
	// every time it is selected, so the first visit allocates and every later
	// one is a dimension assignment -- Eigen's resize reallocates only when the
	// total element count changes. That replaces one construct/destruct per
	// (i,j) block, which over all the gamma families is thousands of
	// allocations per trajectory of a buffer SumRing overwrites in full anyway.
	//
	// RowMajor is what lets SumRing take its direct device->host path: the
	// gathered panel's [gt][iii][m][jjj] layout and a RowMajor (nt, Nii, 1,
	// Njj) tensor are then the same addresses, so its scatter is skipped and
	// its "scatter" timer stays at zero. ColMajor would put t fastest in
	// memory while the copy-out below walks t outermost, which is both the
	// wrong order for that loop and the reason the direct path could not apply.
	// Element access is layout independent, so the values are unchanged.
	Eigen::Tensor<ComplexD, 4, Eigen::RowMajor> resPool[2][2];

	for ( unsigned int j = 0; j < N_j; j += block ){
	  int Njj = MIN(N_j-j,block);
	  startTimer("LoopRight contraction");
	  for (int jj = 0; jj < Njj; jj++) {
	    switch (type) {
	    case 0: Grid::A2AExtendedMesonField<FImpl>::LoopRightContractionType0(loopRight[jj], tloop, right[j+jj], gamma1, gamma2); break;
	    case 1: Grid::A2AExtendedMesonField<FImpl>::LoopRightContractionType1(loopRight[jj], tloop, right[j+jj]);                 break;
	    case 2: Grid::A2AExtendedMesonField<FImpl>::LoopRightContractionType2(loopRight[jj], tloop, right[j+jj], gamma1);         break;
	    case 3: Grid::A2AExtendedMesonField<FImpl>::LoopRightContractionType3(loopRight[jj], tloop, right[j+jj]);                 break;
	    }
	  }
	  stopTimer("LoopRight contraction");
	  //LOG(Message) << "loopRight packed for j-block " << j/block << " type " << type << std::endl;

	  startTimer("Allocate");
	  spatial_sum.AllocateRight(Njj, grid);
	  stopTimer("Allocate");
	  startTimer("Pack vectors");
	  spatial_sum.PackRight(loopRight, 0, Njj);
	  stopTimer("Pack vectors");

	  for ( unsigned int i = 0; i < N_i; i += block ) {
	    int Nii = MIN(N_i-i,block);

	    startTimer("Allocate");
	    spatial_sum.AllocateLeft(Nii);
	    stopTimer("Allocate");
	    startTimer("Pack vectors");
	    spatial_sum.PackLeftConj(left, i, Nii);
	    stopTimer("Pack vectors");

	    // Rank 4 with a singleton momentum axis: SumRing writes
	    // result[t][i][m][j] for the general nmom case, and EMF carries no
	    // momentum projection.
	    //
	    // No setZero: SumRing writes every element of the tensor on both its
	    // direct and its scatter path, so zeroing first is dead work.
	    auto &emfBlock = resPool[Nii != block][Njj != block];
	    startTimer("Allocate");
	    emfBlock.resize(nt, Nii, 1, Njj);
	    stopTimer("Allocate");

	    startTimer("Sum");
	    spatial_sum.SumRing(emfBlock, cacheBlock, &sumTimings, &sumBytes);
	    stopTimer("Sum");

	    startTimer("Copy out");
	    thread_for_collapse(3, t, nt, {
	      for(int ii=0;ii< Nii;ii++)
	      for(int jj=0;jj< Njj;jj++)
		emf(0,0,(int)t,i+ii,j+jj) = emfBlock((int)t,ii,0,jj);
	    });
	    stopTimer("Copy out");
	}
	//LOG(Message) << "EMF made for j-block " << j/block << " type " << type << std::endl;

	}// i,j
	LOG(Message) << "EMF made for type " << type << "; gamma1: " << nameg1_[ig] << "; gamma2: " << nameg2_[ig] << std::endl;

	std::string ioname = "type" + std::to_string(type) + "_" + nameg1_[ig] + "_" + nameg2_[ig];
	std::string filename = par().output + "." + std::to_string(vm().getTrajectory()) + "/" + ioname + ".h5";
        LOG(Message) << "Writing block to " << filename << std::endl;
        makeFileDir(filename, grid);
        double ioBytes = static_cast<double>(nt) * N_i * N_j * sizeof(HADRONS_A2AM_IO_TYPE);
        startTimer("IO");
        double writeTime = -usecond();
#ifdef HADRONS_A2AM_PARALLEL_IO
        startTimer("Barrier");
        grid->Barrier();
        stopTimer("Barrier");
        LOG(Message) << "HADRONS_A2AM_PARALLEL_IO" << std::endl;
	if ( grid->ThisRank() == 0 ) {
#endif
	  A2AMatrixIo<HADRONS_A2AM_IO_TYPE> io(filename,ioname,nt, N_i, N_j);
	  A2AExtendedMesonFieldJHMetadata md;
	  md.gamma1 = nameg1_[ig];
	  md.gamma2 = nameg2_[ig];
	  startTimer("initFile");
	  io.initFile(md, MAX(N_i,N_j));
	  stopTimer("initFile");
	  //io.saveBlock(emf, 0, 0, 0, 0, &ioTimings);
	  io.saveBlock(emf, 0, 0, 0, 0);
#ifdef HADRONS_A2AM_PARALLEL_IO
	}
	startTimer("Barrier");
	grid->Barrier();
	stopTimer("Barrier");
#endif
        writeTime += usecond();
        stopTimer("IO");
        // writeTime is this rank's own wall time from barrier to barrier
        // (only rank 0's LOG output survives, since Grid_quiesce_nodes
        // suppresses the rest by default).
        if (writeTime > 0.)
            LOG(Message) << "IO type=" << type << " ig=" << ig << ": "
                         << sizeString(ioBytes) << " in " << writeTime
                         << " us local (" << ioBytes / writeTime * 1.e6 / 1024. / 1024.
                         << " MB/s effective)" << std::endl;
      }// ig
    }// type

    // Throughput of the post-GEMM SumRing stages -- bytesMoved[k] and
    // sumTimings[k] accumulate the same way across all (type,ig,i,j) calls
    // and all cacheBlock tiles, so their ratio is the average effective
    // bandwidth of that stage over the whole run, comparable across
    // different cacheBlock choices.
    //
    // The two ring stages report bytes on the wire rather than payload, so
    // their rates are the ones comparable with a link rate; the local stages
    // report the bytes they actually touch. See the SumRing header comment.
    auto gbps = [](double bytes, double us)
    {
        return (us > 0.) ? bytes / us * 1.e6 / 1024. / 1024. / 1024. : 0.;
    };
    LOG(Message) << "Sum detail (us), rank 0:" << std::endl;
    LOG(Message) << "  GEMM            = " << sumTimings[0] << std::endl;
    LOG(Message) << "  device->host    = " << sumTimings[1]
                 << " (" << gbps(sumBytes[1], sumTimings[1]) << " GB/s)" << std::endl;
    LOG(Message) << "  gather to slab  = " << sumTimings[2]
                 << " (" << gbps(sumBytes[2], sumTimings[2]) << " GB/s)" << std::endl;
    LOG(Message) << "  spatial reduce  = " << sumTimings[3]
                 << " (" << gbps(sumBytes[3], sumTimings[3]) << " GB/s wire)" << std::endl;
    LOG(Message) << "  scatter         = " << sumTimings[4]
                 << " (" << gbps(sumBytes[4], sumTimings[4]) << " GB/s)" << std::endl;
    LOG(Message) << "  temporal gather = " << sumTimings[5]
                 << " (" << gbps(sumBytes[5], sumTimings[5]) << " GB/s wire)" << std::endl;
    LOG(Message) << "IO detail (us), rank 0:" << std::endl;
    /*
    LOG(Message) << "  open            = " << ioTimings[0]  << std::endl;
    LOG(Message) << "  push/group      = " << ioTimings[1]  << std::endl;
    LOG(Message) << "  openDataSet     = " << ioTimings[2]  << std::endl;
    LOG(Message) << "  getSpace        = " << ioTimings[3]  << std::endl;
    LOG(Message) << "  selectHyperslab = " << ioTimings[4]  << std::endl;
    LOG(Message) << "  write           = " << ioTimings[5]  << std::endl;
    LOG(Message) << "  close(fsync)    = " << ioTimings[6]  << std::endl;
    */
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_A2AExtendedMesonFieldJH_hpp_
