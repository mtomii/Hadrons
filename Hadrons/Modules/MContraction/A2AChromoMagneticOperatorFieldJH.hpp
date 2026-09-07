/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: Hadrons/Modules/MContraction/A2AChromoMagneticOperatorFieldJH.hpp

Copyright (C) 2015-2019

Author: Peter Boyle <paboyle@bnl.gov>
Author: Jonas Hildebrand <jonas.hildebrand@uconn.edu>
Author: Masaaki Tomii <masaaki.tomii@uconn.edu>
*************************************************************************************/
/*  END LEGAL */
#ifndef Hadrons_MContraction_A2AChromoMagneticOperatorFieldJH_hpp_
#define Hadrons_MContraction_A2AChromoMagneticOperatorFieldJH_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/A2AMatrix.hpp>
#include <Grid/qcd/utils/A2Autils.h>

BEGIN_HADRONS_NAMESPACE
/******************************************************************************
 *          All-to-all chromo-magnetic operator field creation (GPU)          *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MContraction)

class A2AChromoMagneticOperatorFieldJHPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(A2AChromoMagneticOperatorFieldJHPar,
                                    int,         block,
                                    int,         cacheBlock,
                                    std::string, parities,
                                    std::string, left,
                                    std::string, right,
                                    std::string, gauge,
                                    std::string, output,
                                    std::string, ifOrthogs);
};

class A2AChromoMagneticOperatorFieldJHMetadata: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(A2AChromoMagneticOperatorFieldJHMetadata,
                                    std::string, meta);
};

template <typename GImpl, typename FImpl>
class TA2AChromoMagneticOperatorFieldJH : public Module<A2AChromoMagneticOperatorFieldJHPar>
{
public:
  typedef typename GImpl::GaugeLinkField GaugeMat;
  typedef typename GImpl::SiteGaugeLink  SiteGaugeLink;
  FERM_TYPE_ALIASES(FImpl,);
  typedef typename FImpl::SiteSpinor     vobj;
  typedef typename vobj::scalar_object   sobj;
  typedef typename vobj::scalar_type     scalar_type;
  typedef typename vobj::vector_type     vector_type;
public:
    // constructor
    TA2AChromoMagneticOperatorFieldJH(const std::string name);
    // destructor
    virtual ~TA2AChromoMagneticOperatorFieldJH(void){};
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

MODULE_REGISTER(A2AChromoMagneticOperatorFieldJH, ARG(TA2AChromoMagneticOperatorFieldJH<GIMPL,FIMPL>), MContraction);

/******************************************************************************
*               TA2AChromoMagneticOperatorFieldJH implementation                *
******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename GImpl, typename FImpl>
TA2AChromoMagneticOperatorFieldJH<GImpl,FImpl>::TA2AChromoMagneticOperatorFieldJH(const std::string name)
: Module<A2AChromoMagneticOperatorFieldJHPar>(name)
{
}

// dependencies/products ///////////////////////////////////////////////////////
template <typename GImpl, typename FImpl>
std::vector<std::string> TA2AChromoMagneticOperatorFieldJH<GImpl,FImpl>::getInput(void)
{
  std::vector<std::string> in = {par().left, par().right, par().gauge};
  return in;
}

template <typename GImpl, typename FImpl>
std::vector<std::string> TA2AChromoMagneticOperatorFieldJH<GImpl,FImpl>::getOutput(void)
{
  std::vector<std::string> out = {};
  return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename GImpl, typename FImpl>
void TA2AChromoMagneticOperatorFieldJH<GImpl,FImpl>::setup(void)
{
  parities_  = strToVec<int>(par().parities);
  ifOrthogs_ = strToVec<int>(par().ifOrthogs);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename GImpl, typename FImpl>
void TA2AChromoMagneticOperatorFieldJH<GImpl,FImpl>::execute(void)
{
  typedef iSpinColourVector<vector_type> SpinColourVector_v;

  auto &left    = envGet(std::vector<FermionField>, par().left);
  auto &right   = envGet(std::vector<FermionField>, par().right);
  const auto &U = envGet(GaugeField, par().gauge);

  GridBase *grid = left[0].Grid();

  LOG(Message) << "Computing all-to-all ChromoMagnetic operator fields (GPU)" << std::endl;

  int nt         = env().getDim().back();
  int N_i        = left.size();
  int N_j        = right.size();
  int block = par().block;
  int cacheBlock = par().cacheBlock;

  const std::size_t mBufSize =
    static_cast<std::size_t>(nt)
    * static_cast<std::size_t>(N_i)
    * static_cast<std::size_t>(N_j);

  Vector<HADRONS_A2AM_IO_TYPE> mBuf; mBuf.resize(mBufSize);

  LOG(Message) << "Left: '"        << par().left  << "' Right: '"
               << par().right      << "'"          << std::endl;
  LOG(Message) << "Gauge field: '" << par().gauge  << "'"          << std::endl;
  LOG(Message) << "Parities:"      << std::endl;
  for (auto &p: parities_)
    LOG(Message) << "  " << p << std::endl;
  LOG(Message) << "CMO field size: " << nt << "*" << N_i << "*" << N_j
               << " (filesize " << sizeString(mBufSize*sizeof(HADRONS_A2AM_IO_TYPE)) << std::endl;

  std::vector<FermionField> loopRight(block, grid);

  std::array<double, 6> sumTimings = {};
  std::array<double, 6> sumBytes   = {};
  //std::array<double, 7> ioTimings  = {};

  for (auto &ifOrthog: ifOrthogs_) {
    std::vector<GaugeMat>  G;
    Vector<Gamma::Algebra> Sigma;

    startTimer("CMO contraction");
    if (ifOrthog == 0)
      Grid::A2AChromoMagneticOperator<GImpl,FImpl>::CMOContraction0(G, Sigma, U);
    else
      Grid::A2AChromoMagneticOperator<GImpl,FImpl>::CMOContraction1(G, Sigma, U);
    stopTimer("CMO contraction");
    LOG(Message) << "Field strength constructed for ifOrthog=" << ifOrthog << std::endl;

    for (auto &parity: parities_) {
      LOG(Message) << "Starting calculation with ifOrthog=" << ifOrthog
                   << " parity=" << parity << std::endl;
      A2AMatrixSet<HADRONS_A2AM_IO_TYPE> cmf(mBuf.data(), 1, 1, nt, N_i, N_j);
      A2ASpatialSum<SpinColourVector_v> spatial_sum;

      // Result buffers, one per distinct block shape. A block is full or on
      // the tail in each axis independently, so a 2x2 pool indexed by (i on
      // tail, j on tail) covers every case. Each slot is asked for the same
      // dimensions every time it is selected, so the first visit allocates and
      // every later one is a dimension assignment -- Eigen's resize
      // reallocates only when the total element count changes. That replaces
      // one construct/destruct per (i,j) block of a buffer SumRing overwrites
      // in full anyway.
      //
      // RowMajor is what lets SumRing take its direct device->host path: the
      // gathered panel's [gt][iii][m][jjj] layout and a RowMajor (nt, Nii, 1,
      // Njj) tensor are then the same addresses, so its scatter is skipped and
      // its "scatter" timer stays at zero. ColMajor would put t fastest in
      // memory while the copy-out below walks t outermost, which is both the
      // wrong order for that loop and the reason the direct path could not
      // apply. Element access is layout independent, so the values are
      // unchanged.
      Eigen::Tensor<ComplexD, 4, Eigen::RowMajor> resPool[2][2];

      LOG(Message) << "Making CMF" << std::endl;

      for (unsigned int j = 0; j < N_j; j += block) {
        int Njj = MIN(N_j-j, block);

        startTimer("CMOContractRight");
        for (int jj = 0; jj < Njj; jj++)
          Grid::A2AChromoMagneticOperator<GImpl,FImpl>::CMOContractRight(
              loopRight[jj], G, Sigma, right[j+jj], parity);
        stopTimer("CMOContractRight");
        LOG(Message) << "loopRight packed for j-block " << j/block
                     << " ifOrthog=" << ifOrthog << " parity=" << parity << std::endl;

        startTimer("Allocate");
        spatial_sum.AllocateRight(Njj, grid);
        stopTimer("Allocate");
        startTimer("Pack vectors");
        spatial_sum.PackRight(loopRight, 0, Njj);
        stopTimer("Pack vectors");

        for (unsigned int i = 0; i < N_i; i += block) {
          int Nii = MIN(N_i-i, block);

          startTimer("Allocate");
          spatial_sum.AllocateLeft(Nii);
          stopTimer("Allocate");
          startTimer("Pack vectors");
          spatial_sum.PackLeftConj(left, i, Nii);
          stopTimer("Pack vectors");

          // Rank 4 with a singleton momentum axis: SumRing writes
          // result[t][i][m][j] for the general nmom case, and the CMO field
          // carries no momentum projection.
          //
          // No setZero: SumRing writes every element of the tensor on both its
          // direct and its scatter path, so zeroing first is dead work.
          auto &cmfBlock = resPool[Nii != block][Njj != block];
          startTimer("Allocate");
          cmfBlock.resize(nt, Nii, 1, Njj);
          stopTimer("Allocate");

          startTimer("Sum");
          spatial_sum.SumRing(cmfBlock, cacheBlock, &sumTimings, &sumBytes);
          stopTimer("Sum");

          startTimer("Copy out");
          thread_for_collapse(3, t, nt, {
            for (int ii = 0; ii < Nii; ii++)
            for (int jj = 0; jj < Njj; jj++)
              cmf(0,0,(int)t,i+ii,j+jj) = cmfBlock((int)t,ii,0,jj);
          });
          stopTimer("Copy out");

          //LOG(Message) << "CMF made for i-block " << i/block
          //             << " j-block "             << j/block
          //             << " ifOrthog="            << ifOrthog
          //             << " parity="              << parity << std::endl;

        }// i
      }// j

      LOG(Message) << "CMF made for ifOrthog=" << ifOrthog << " parity=" << parity << std::endl;

      std::string ioname = "parity" + std::to_string(parity);
      if (ifOrthog == 1)
        ioname = ioname + "_GijSij";
      else
        ioname = ioname + "_GitSit";
      std::string filename = par().output + "." + std::to_string(vm().getTrajectory())
                           + "/" + ioname + ".h5";
      LOG(Message) << "Writing block to " << filename << std::endl;
      makeFileDir(filename, grid);
      startTimer("IO");
#ifdef HADRONS_A2AM_PARALLEL_IO
      startTimer("Barrier");
      grid->Barrier();
      stopTimer("Barrier");
      if (grid->ThisRank() == 0) {
#endif
      A2AMatrixIo<HADRONS_A2AM_IO_TYPE> io(filename, ioname, nt, N_i, N_j);
      A2AChromoMagneticOperatorFieldJHMetadata md;
      md.meta = ioname;
      startTimer("initFile");
      io.initFile(md, MAX(N_i,N_j));
      stopTimer("initFile");
      //io.saveBlock(cmf, 0, 0, 0, 0, &ioTimings);
      io.saveBlock(cmf, 0, 0, 0, 0);
#ifdef HADRONS_A2AM_PARALLEL_IO
      }
      startTimer("Barrier");
      grid->Barrier();
      stopTimer("Barrier");
#endif
      stopTimer("IO");
    }// parity
  }// ifOrthog

  // Throughput of the post-GEMM SumRing stages -- bytesMoved[k] and
  // sumTimings[k] accumulate the same way across all (ifOrthog,parity,i,j)
  // calls and all cacheBlock tiles, so their ratio is the average effective
  // bandwidth of that stage over the whole run, comparable across different
  // cacheBlock choices.
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

#endif // Hadrons_MContraction_A2AChromoMagneticOperatorFieldJH_hpp_
