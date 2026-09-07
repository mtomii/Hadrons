/*
 * A2ANewMesonField.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2024
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Peter Boyle <paboyle@ph.ed.ac.uk>
 * Author: ferben <ferben@debian.felix.com>
 * Author: paboyle <paboyle@ph.ed.ac.uk>
 * Author: Jonas Hildebrand <jonas.hildebrand@uconn.edu>
 *
 * Hadrons is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * Hadrons is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Hadrons.  If not, see <http://www.gnu.org/licenses/>.
 *
 * See the full license in the file "LICENSE" in the top level distribution
 * directory.
 */

/*  END LEGAL */
#ifndef Hadrons_MContraction_A2ANewMesonField_hpp_
#define Hadrons_MContraction_A2ANewMesonField_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/A2AMatrix.hpp>
#include <Grid/qcd/utils/A2Autils.h>
#include <Grid/algorithms/blas/A2ASpatialSum.h>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *  All-to-all meson field creation. Drives A2ASpatialSum with mode index
 *  blocking exposed on the Hadrons side. Mode and momentum indices are packed
 *  into buffers that are fed into the GEMM call, which batches over timeslice.
 *  Gamma index is the only explicit loop present in the module.
 *
 *  SumRing does the GEMM, and then we execute a spatial ring all reduce to
 *  complete the spatial + spin-colour reduction followed by a purely temporal
 *  gather, constructing the full meson field through a ring rather than
 *  GlobalSumVector.
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MContraction)

class A2ANewMesonFieldPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(A2ANewMesonFieldPar,
                                    int,                     block,
                                    int,                     cacheBlock,
                                    std::string,             left,
                                    std::string,             right,
                                    std::string,             output,
                                    std::string,             gammas,
                                    std::vector<std::string>, mom);
};

class A2ANewMesonFieldMetadata: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(A2ANewMesonFieldMetadata,
                                    std::vector<RealF>, momentum,
                                    Gamma::Algebra,     gamma);
};

template <typename FImpl>
class TA2ANewMesonField : public Module<A2ANewMesonFieldPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
public:
    TA2ANewMesonField(const std::string name);
    virtual ~TA2ANewMesonField(void){};
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    virtual void setup(void);
    virtual void execute(void);
private:
    bool                                      hasPhase_{false};
    std::string                               momphName_;
    std::vector<Gamma::Algebra>               gamma_;
    std::vector<std::vector<Real>>            mom_;
    A2ASpatialSum<typename FImpl::SiteSpinor> spatial_sum_;
};

MODULE_REGISTER(A2ANewMesonField, ARG(TA2ANewMesonField<FIMPL>), MContraction);

/******************************************************************************
 *                  TA2ANewMesonField implementation                       *
 ******************************************************************************/
template <typename FImpl>
TA2ANewMesonField<FImpl>::TA2ANewMesonField(const std::string name)
: Module<A2ANewMesonFieldPar>(name)
, momphName_(name + "_momph")
{}

template <typename FImpl>
std::vector<std::string> TA2ANewMesonField<FImpl>::getInput(void)
{
    return {par().left, par().right};
}

template <typename FImpl>
std::vector<std::string> TA2ANewMesonField<FImpl>::getOutput(void)
{
    return {};
}

template <typename FImpl>
void TA2ANewMesonField<FImpl>::setup(void)
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

template <typename FImpl>
void TA2ANewMesonField<FImpl>::execute(void)
{
    typedef typename FImpl::SiteSpinor      vobj;
    typedef typename vobj::vector_type      vector_type;
    typedef iSpinColourVector<vector_type>  SpinColourVector_v;
    typedef typename SpinColourVector_v::scalar_type scalar_t;

    auto &left  = envGet(std::vector<FermionField>, par().left);
    auto &right = envGet(std::vector<FermionField>, par().right);

    GridBase *grid = envGetGrid(FermionField);

    int nt     = env().getDim().back();
    int N_i    = left.size();
    int N_j    = right.size();
    int ngamma = gamma_.size();
    int nmom   = mom_.size();
    int block  = par().block;
    int cacheBlock = par().cacheBlock;

    const std::size_t mBufSize =
      static_cast<std::size_t>(nt)
      * static_cast<std::size_t>(N_i)
      * static_cast<std::size_t>(N_j);

    LOG(Message) << "Computing all-to-all meson fields" << std::endl;
    LOG(Message) << "Left: '" << par().left << "' Right: '" << par().right << "'" << std::endl;
    LOG(Message) << "Momenta:" << std::endl;
    for (auto &p: mom_)
        LOG(Message) << "  " << p << std::endl;
    LOG(Message) << "Spin bilinears:" << std::endl;
    for (auto &g: gamma_)
        LOG(Message) << "  " << g << std::endl;
    LOG(Message) << "Meson field size: " << nt << "*" << N_i << "*" << N_j
                 << " (filesize " << sizeString(mBufSize*sizeof(HADRONS_A2AM_IO_TYPE))
                 << "/momentum/bilinear)" << std::endl;

    auto &ph = envGet(std::vector<ComplexField>, momphName_);

    if (!hasPhase_)
    {
        startTimer("Momentum phases");
        for (int j = 0; j < nmom; ++j)
        {
            Complex i(0.0, 1.0);
            envGetTmp(ComplexField, coor);
            ph[j] = Zero();
            for (unsigned int mu = 0; mu < mom_[j].size(); mu++)
            {
                LatticeCoordinate(coor, mu);
                ph[j] = ph[j] + (mom_[j][mu]/env().getDim(mu))*coor;
            }
            ph[j] = exp((Real)(2*M_PI)*i*ph[j]);
        }

        hasPhase_ = true;
        stopTimer("Momentum phases");
    }

    auto ionameFn = [this](const int m, const int g)
    {
        std::stringstream ss;
        ss << gamma_[g] << "_";
        for (unsigned int mu = 0; mu < mom_[m].size(); ++mu)
            ss << mom_[m][mu] << ((mu == mom_[m].size() - 1) ? "" : "_");
        return ss.str();
    };

    auto metadataFn = [this](const int m, const int g)
    {
        A2ANewMesonFieldMetadata md;
        for (auto pmu: mom_[m])
            md.momentum.push_back(pmu);
        md.gamma = gamma_[g];
        return md;
    };

    auto filenameFn = [this, &ionameFn](const int m, const int g)
    {
        return par().output + "." + std::to_string(vm().getTrajectory())
               + "/" + ionameFn(m, g) + ".h5";
    };

    // Output buffer: one (nt, Nii, Njj) block at a time.
    Vector<HADRONS_A2AM_IO_TYPE> mBuf;
    mBuf.resize(nt * block * block);

    // Scratch right vectors for GammaRight output (zero-momentum base pack).
    std::vector<FermionField> gammaRight(block, grid);

    // Result buffers, one per distinct block shape. A block is either full or
    // on the tail in each axis independently, so there are at most four shapes
    // and a 2x2 pool indexed by (i on tail, j on tail) covers every case. Each
    // slot is asked for the same dimensions every time it is selected, so the
    // first visit allocates and every later one is a dimension assignment --
    // Eigen's resize reallocates only when the total element count changes.
    //
    // Holding the shapes side by side rather than resizing one buffer is what
    // keeps the kernel from re-zeroing a freshly mapped region on every grow:
    // the tail shape and the full shape otherwise alternate once per jb
    // iteration, and the full shape here is gigabytes. The cost is the sum of
    // the shapes instead of the max, a few GB against the A2A vectors' tens of
    // TB. It also makes `block` a free choice again -- N_i = nLow + Nsc*nHit
    // and N_j = nLow + Nhigh*nHit move with the hit count, so no fixed block
    // divides both, and a block that divides neither now costs two extra
    // buffers rather than anything on the critical path.
    //
    // Dimension order (nt, Nii, nmom, Njj) -- nmom BEFORE N_j -- is the layout
    // SumRing writes, matching its GEMM's [i][m][j] output. RowMajor then makes
    // N_j the fastest dimension, so the IO fill below, which reads at fixed m
    // and walks jj innermost, is contiguous on both sides. Passing the true
    // per-block dimensions rather than one buffer padded to block is also what
    // lets SumRing take its direct device->host path on every block including
    // the tails, which leaves its "scatter" timer at zero.
    Eigen::Tensor<ComplexD, 4, Eigen::RowMajor> resPool[2][2];

    // Every rank creates the output directory itself, rather than relying
    // on makeFileDir's boss-rank-only mkdir. File writes below are spread
    // across ranks (owned by (m*ngamma+g) % nRank), and output can land on
    // physically separate per-node storage (e.g. a node-local NVMe burst
    // buffer) where a directory created on the boss rank's node is simply
    // absent on every other node -- no barrier can fix that, since it isn't
    // a visibility-lag problem, the directory really doesn't exist there.
    // Hadrons::mkdir checks access() first and only costs a few redundant
    // syscalls when the directory already exists, so doing this on every
    // rank is harmless (if a little repetitive) on a shared filesystem too.
    std::string dirBase = par().output + "." + std::to_string(vm().getTrajectory());
    Hadrons::mkdir(dirBase);
    grid->Barrier();

    unsigned int myRank = grid->ThisRank();
    unsigned int nRank  = grid->RankCount();

    // Initialise one HDF5 file per (mom, gamma) pair before the block loops.
    // Each rank initialises only its assigned files; single barrier after.
    for (int m = 0; m < nmom; m++)
    for (int g = 0; g < ngamma; g++)
    {
        if ((unsigned int)(m * ngamma + g) % nRank == myRank)
        {
            std::string ioname   = ionameFn(m, g);
            std::string filename = filenameFn(m, g);
            A2ANewMesonFieldMetadata md = metadataFn(m, g);
            A2AMatrixIo<HADRONS_A2AM_IO_TYPE> io(filename, ioname, nt, N_i, N_j);
            io.initFile(md, block);
        }
    }
    grid->Barrier();

    // Pre-pack flat phase arrays, one absolute phase per momentum -- no
    // difference-encoding needed here since ApplyAllPhaseRight reads a single
    // unphased base pack and writes all nmom copies directly, rather than
    // stepping an in-place buffer through consecutive momenta.
    startTimer("Pack phases");
    std::vector<deviceVector<scalar_t>> ph_flat(nmom);
    for (int m = 0; m < nmom; m++)
        A2ASpatialSum<SpinColourVector_v>::PackPhase(grid, ph[m], ph_flat[m]);
    stopTimer("Pack phases");

    // One-time allocation for the full block size; subsequent pointer rewrites are cheap.
    startTimer("Allocate");
    spatial_sum_.AllocateRight(block, grid, nmom);
    spatial_sum_.AllocateLeft(block);
    stopTimer("Allocate");

    // Loop order (jb, g, ib):
    //   AllocateRight + PackRight + ApplyAllPhaseRight - once per (jb, g)
    //   AllocateLeft  + PackLeftConj + SumRing         - once per (jb, g, ib)

    double                fillTime     = 0.;
    //std::array<double, 7> ioTimings    = {};
    std::array<double, 6> sumTimings   = {};
    std::array<double, 6> sumBytes     = {};

    for (int jb = 0; jb < N_j; jb += block)
    {
        int Njj = std::min(N_j - jb, block);

        startTimer("Allocate");
        spatial_sum_.AllocateRight(Njj, grid, nmom);
        stopTimer("Allocate");

        for (int g = 0; g < ngamma; g++)
        {
            startTimer("GammaRight");
            for (int jj = 0; jj < Njj; jj++)
                A2Autils<FImpl>::GammaRight(gammaRight[jj], gamma_[g], right[jb + jj]);
            stopTimer("GammaRight");

            startTimer("Pack vectors");
            spatial_sum_.PackRight(gammaRight, 0, Njj);
            stopTimer("Pack vectors");

            startTimer("Phase");
            spatial_sum_.ApplyAllPhaseRight(ph_flat);
            stopTimer("Phase");

            for (int ib = 0; ib < N_i; ib += block)
            {
                int Nii = std::min(N_i - ib, block);

                // Pick the pool slot for this block's shape. Allocates on the
                // first visit to each shape, dimension assignment after that,
                // so a nonzero "Allocate" time past the first jb iteration
                // means a shape is being reallocated and the pool is missing.
                auto &all_results = resPool[Nii != block][Njj != block];

                startTimer("Allocate");
                spatial_sum_.AllocateLeft(Nii);
                all_results.resize(nt, Nii, nmom, Njj);
                stopTimer("Allocate");

                startTimer("Pack vectors");
                spatial_sum_.PackLeftConj(left, ib, Nii);
                stopTimer("Pack vectors");

                startTimer("Sum");
                spatial_sum_.SumRing(all_results, cacheBlock, &sumTimings, &sumBytes);
                stopTimer("Sum");

                // Parallel IO: each rank writes its assigned momenta simultaneously.
                // Barrier count drops from 2*nmom to 2 per outer block.
                //
                // Ownership key matches the (m*ngamma+g) % nRank used by the
                // initFile loop above, so the rank that writes a given (m,g)
                // file is always the same rank that created it -- creating on
                // one rank and writing from another would depend on that
                // file being visible from a different rank's node right
                // after the barrier, which a plain MPI_Barrier does not
                // guarantee on a parallel filesystem (client-side metadata
                // caching can lag), and "file not found" races result.
                double ioBytes = static_cast<double>(nmom) * nt * Nii * Njj
                                 * sizeof(HADRONS_A2AM_IO_TYPE);
                startTimer("IO");
                double writeTime = -usecond();
                grid->Barrier();
                for (int m = 0; m < nmom; m++)
                {
                    if ((unsigned int)(m * ngamma + g) % nRank != myRank) continue;

                    A2AMatrixSet<HADRONS_A2AM_IO_TYPE> mf(mBuf.data(), 1, 1, nt, Nii, Njj);
                    double dt = -usecond();
                    thread_for_collapse(3, t, nt, {
                        for (int ii = 0; ii < Nii; ii++)
                        for (int jj = 0; jj < Njj; jj++)
                            mf(0, 0, (int)t, ii, jj) = all_results((int)t, ii, m, jj);
                    });
                    dt += usecond();
                    fillTime += dt;
                    A2AMatrixIo<HADRONS_A2AM_IO_TYPE> io(filenameFn(m, g), ionameFn(m, g),
                                                         nt, N_i, N_j);
                    //io.saveBlock(mf, 0, 0, ib, jb, &ioTimings);
                    io.saveBlock(mf, 0, 0, ib, jb);
                }
                grid->Barrier();
                writeTime += usecond();
                stopTimer("IO");
                // writeTime is this rank's own wall time from barrier to
                // barrier (only rank 0's LOG output survives, since
                // Grid_quiesce_nodes suppresses the rest by default).
                if (writeTime > 0. && ib == jb)
                    LOG(Message) << "IO block i=" << ib << " j=" << jb
                                 << " g=" << gamma_[g] << ": "
                                 << sizeString(ioBytes) << " in "
                                 << writeTime << " us local ("
                                 << ioBytes / writeTime * 1.e6 / 1024. / 1024.
                                 << " MB/s effective)" << std::endl;
            } // ib
        } // g
    } // jb

    // Throughput of the post-GEMM SumRing stages -- bytesMoved[k] and
    // sumTimings[k] accumulate the same way across all (jb,g,ib) calls and
    // all cacheBlock tiles, so their ratio is the average effective
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
    LOG(Message) << "Sum detail (us), rank " << myRank << ":" << std::endl;
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
    LOG(Message) << "IO detail (us), rank " << myRank << ":" << std::endl;
    LOG(Message) << "  fill            = " << fillTime      << std::endl;
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

#endif // Hadrons_MContraction_A2ANewMesonField_hpp_
