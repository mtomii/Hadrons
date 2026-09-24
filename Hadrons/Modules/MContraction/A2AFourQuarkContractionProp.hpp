/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: Hadrons/Modules/MContraction/A2AFourQuarkContractionProp.hpp

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
#ifndef Hadrons_MContraction_A2AFourQuarkContractionProp_hpp_
#define Hadrons_MContraction_A2AFourQuarkContractionProp_hpp_
#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <sstream>
#include <type_traits>

/* Based on A2AFourQuarkContractionMT.hpp. Keep its paired gamma groups,
 * four spin/colour contractions, local SIMD sum and result layout.
 *
 * mat[u*nSpace+p] pairs with prop(x_p, (tmin+u) mod Nt).
 * tmin is an unwrapped physical time, not an index into the short vector.
 * No loop/window copy, extra conjugation, gamma5 or hit normalization.
 * CPU; full standard 4D Grid, MPI_t = SIMD_t = 1. Spatial splitting is allowed.
 * The standard-vector input must already have that Grid's spatial layout.
 */

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *             Four-quark contraction: time-window vector and full field      *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MContraction)
class A2AFourQuarkContractionPropPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(A2AFourQuarkContractionPropPar,
                                    std::string, mat,
                                    std::string, prop,
                                    int,         tmin,
                                    std::string, sctypes,
                                    std::string, gammas1,
                                    std::string, gammas2,
                                    std::string, output);
};

template <typename FImpl>
class TA2AFourQuarkContractionProp: public Module<A2AFourQuarkContractionPropPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    typedef typename FImpl::SiteSpinor vobj;
    typedef typename vobj::vector_type vector_type;
    typedef typename vobj::scalar_type scalar_type;
    typedef iSpinColourMatrix<vector_type> SpinColourMatrix_v;
    typedef iSpinMatrix<vector_type> SpinMatrix_v;
    typedef iColourMatrix<vector_type> ColourMatrix_v;
    typedef iSinglet<vector_type> Scalar_v;
    typedef iSinglet<scalar_type> Scalar_s;
    static_assert(std::is_same<typename PropagatorField::vector_object,
                               SpinColourMatrix_v>::value,
                  "Window matrix and PropagatorField site types must agree");
public:
    // constructor
    TA2AFourQuarkContractionProp(const std::string name);
    // destructor
    virtual ~TA2AFourQuarkContractionProp(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    virtual std::vector<std::string> getOutputFiles(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
    std::vector<std::vector<Gamma::Algebra>> gamma1_, gamma2_;
    std::vector<std::string> nameg1_, nameg2_;
    std::vector<int> types_;
    int nt_{0};
    std::size_t nSpace_{0}, nTime_{0}, nElements_{0}, nCorr_{0}, nWork_{0};
};
MODULE_REGISTER_TMP(A2AFourQuarkContractionProp,
                    TA2AFourQuarkContractionProp<FIMPL>, MContraction);

// Do not reuse CorrelatorResult: the old MT header defines that name already.
class A2AFourQuarkContractionPropResult: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(A2AFourQuarkContractionPropResult,
                                    std::vector<ComplexD>, correlator);
};

/******************************************************************************
 *                 TA2AFourQuarkContractionProp implementation                 *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TA2AFourQuarkContractionProp<FImpl>::
TA2AFourQuarkContractionProp(const std::string name)
: Module<A2AFourQuarkContractionPropPar>(name)
{}
// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TA2AFourQuarkContractionProp<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().mat, par().prop};

    return in;
}

template <typename FImpl>
std::vector<std::string> TA2AFourQuarkContractionProp<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {};

    return out;
}

template <typename FImpl>
std::vector<std::string> TA2AFourQuarkContractionProp<FImpl>::getOutputFiles(void)
{
    std::vector<std::string> out = {par().output + "/"
        + ModuleBase::resultFilename(getName(), vm().getTrajectory())};

    return out;
}
// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TA2AFourQuarkContractionProp<FImpl>::setup(void)
{
    if (par().mat.empty() || par().prop.empty() || par().output.empty())
    {
        HADRONS_ERROR(Argument, "A2AFourQuarkContractionProp: empty mat, prop or output (use . for the current directory)");
    }
    gamma1_.clear();
    gamma2_.clear();
    types_.clear();
    nameg1_ = strToVec<std::string>(par().gammas1);
    nameg2_ = strToVec<std::string>(par().gammas2);
    if (nameg1_.empty() || nameg1_.size() != nameg2_.size())
    {
        HADRONS_ERROR(Argument, "A2AFourQuarkContractionProp: gamma lists must be nonempty and have equal lengths");
    }
    for (std::size_t ig = 0; ig < nameg1_.size(); ++ig)
    {
        std::vector<Gamma::Algebra> vec;
        if (nameg1_[ig] == "GammaMU")
        {
            vec = {Gamma::Algebra::GammaX, Gamma::Algebra::GammaY,
                   Gamma::Algebra::GammaZ, Gamma::Algebra::GammaT};
        }
        else if (nameg1_[ig] == "GammaMUGamma5")
        {
            vec = {Gamma::Algebra::GammaXGamma5, Gamma::Algebra::GammaYGamma5,
                   Gamma::Algebra::GammaZGamma5, Gamma::Algebra::GammaTGamma5};
        }
        else
        {
            vec = strToVec<Gamma::Algebra>(nameg1_[ig]);
        }
        gamma1_.push_back(vec);
        vec.clear();
        if (nameg2_[ig] == "GammaMU")
        {
            vec = {Gamma::Algebra::GammaX, Gamma::Algebra::GammaY,
                   Gamma::Algebra::GammaZ, Gamma::Algebra::GammaT};
        }
        else if (nameg2_[ig] == "GammaMUGamma5")
        {
            vec = {Gamma::Algebra::GammaXGamma5, Gamma::Algebra::GammaYGamma5,
                   Gamma::Algebra::GammaZGamma5, Gamma::Algebra::GammaTGamma5};
        }
        else
        {
            vec = strToVec<Gamma::Algebra>(nameg2_[ig]);
        }
        gamma2_.push_back(vec);
        if (gamma1_[ig].empty() || gamma1_[ig].size() != gamma2_[ig].size())
        {
            HADRONS_ERROR(Argument, "A2AFourQuarkContractionProp: expanded gamma groups have different lengths");
        }
        for (std::size_t igg = 0; igg < gamma1_[ig].size(); ++igg)
        {
            if (gamma1_[ig][igg] == Gamma::Algebra::undef ||
                gamma2_[ig][igg] == Gamma::Algebra::undef)
            {
                HADRONS_ERROR(Argument, "A2AFourQuarkContractionProp: unknown gamma in pair '"
                              + nameg1_[ig] + " " + nameg2_[ig] + "'");
            }
        }
        for (std::size_t j = 0; j < ig; ++j)
        {
            if (nameg1_[j] == nameg1_[ig] && nameg2_[j] == nameg2_[ig])
            {
                HADRONS_ERROR(Argument, "A2AFourQuarkContractionProp: repeated gamma pair would duplicate a result name");
            }
        }
    }
    // Parse strictly: do not silently accept a prefix such as '1' in '1abc'.
    std::istringstream typeStream(par().sctypes);
    while (typeStream >> std::ws && !typeStream.eof())
    {
        int value;
        if (!(typeStream >> value) || value < 0 || value > 3)
        {
            HADRONS_ERROR(Argument, "A2AFourQuarkContractionProp: sctypes must contain integers from 0 to 3");
        }
        if (std::find(types_.begin(), types_.end(), value) != types_.end())
        {
            HADRONS_ERROR(Argument, "A2AFourQuarkContractionProp: repeated sctype would duplicate a result name");
        }
        types_.push_back(value);
    }
    if (types_.empty())
    {
        HADRONS_ERROR(Argument, "A2AFourQuarkContractionProp: empty sctypes");
    }

    const auto &mat = envGet(std::vector<SpinColourMatrix_v>, par().mat);
    const auto &prop = envGet(PropagatorField, par().prop);
    GridBase *grid = prop.Grid();
    if (grid->_ndimension != 4 || grid->_isCheckerBoarded ||
        env().getObjectLs(par().prop) != 1 || env().getObjectLs(par().mat) != 1)
    {
        HADRONS_ERROR(Argument, "A2AFourQuarkContractionProp requires full, non-checkerboard 4D inputs");
    }
    if (grid->_processors[Tp] != 1 || grid->_simd_layout[Tp] != 1)
    {
        HADRONS_ERROR(Argument, "A2AFourQuarkContractionProp requires MPI_t = SIMD_t = 1");
    }
    if (grid != envGetGrid(PropagatorField))
    {
        HADRONS_ERROR(Size, "A2AFourQuarkContractionProp: prop must use the standard Environment Grid");
    }
    nt_ = grid->_gdimensions[Tp];
    if (nt_ <= 0 || grid->_slice_nblock[Tp] <= 0 || grid->_slice_block[Tp] <= 0 ||
        grid->_slice_stride[Tp] <= 0 || grid->_ostride[Tp] <= 0)
    {
        HADRONS_ERROR(Size, "A2AFourQuarkContractionProp: invalid spatial slice layout");
    }
    const std::size_t e1 = grid->_slice_nblock[Tp];
    const std::size_t e2 = grid->_slice_block[Tp];
    if (e1 > std::numeric_limits<std::size_t>::max()/e2)
    {
        HADRONS_ERROR(Size, "A2AFourQuarkContractionProp: spatial size overflow");
    }
    nSpace_ = e1*e2;
    if (nSpace_ > std::numeric_limits<std::size_t>::max()/std::size_t(nt_) ||
        nSpace_*std::size_t(nt_) != std::size_t(grid->oSites()) ||
        mat.empty() || mat.size() % nSpace_ != 0)
    {
        HADRONS_ERROR(Size, "A2AFourQuarkContractionProp: mat size is not a positive number of standard spatial slices");
    }
    nElements_ = mat.size();
    nTime_ = nElements_/nSpace_;
    if (nTime_ > std::size_t(nt_))
    {
        HADRONS_ERROR(Size, "A2AFourQuarkContractionProp: window is longer than Nt");
    }
    const std::size_t maxSize = std::numeric_limits<std::size_t>::max();
    if (gamma1_.size() > maxSize/types_.size())
    {
        HADRONS_ERROR(Size, "A2AFourQuarkContractionProp: correlator count overflow");
    }
    nCorr_ = gamma1_.size()*types_.size();
    if (nCorr_ > maxSize/nTime_)
    {
        HADRONS_ERROR(Size, "A2AFourQuarkContractionProp: work count overflow");
    }
    nWork_ = nCorr_*nTime_;
    if (nWork_ > maxSize/(sizeof(Scalar_v) + sizeof(ComplexD)))
    {
        HADRONS_ERROR(Size, "A2AFourQuarkContractionProp: correlator byte size overflow");
    }
    LOG(Message) << "A2AFourQuarkContractionProp: " << nTime_
                 << " time slices starting at unwrapped tmin=" << par().tmin
                 << ", " << nCorr_ << " selected correlators; no loop copy"
                 << std::endl;
}
// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TA2AFourQuarkContractionProp<FImpl>::execute(void)
{
    const auto &mat = envGet(std::vector<SpinColourMatrix_v>, par().mat);
    const auto &prop = envGet(PropagatorField, par().prop);
    GridBase *grid = prop.Grid();
    assert(mat.size() == nElements_);
    if (mat.size() != nElements_ || grid != envGetGrid(PropagatorField) ||
        grid->_gdimensions[Tp] != nt_)
    {
        HADRONS_ERROR(Size, "A2AFourQuarkContractionProp: input shape/Grid changed after setup");
    }
    const int Nsimd = grid->Nsimd();
    const std::size_t e2 = grid->_slice_block[Tp];
    const std::size_t stride = grid->_slice_stride[Tp];
    std::vector<std::size_t> timeOffset(nTime_);
    for (std::size_t it = 0; it < nTime_; ++it)
    {
        const std::int64_t t = std::int64_t(par().tmin) + std::int64_t(it);
        const int tx = int((t % nt_ + nt_) % nt_);
        timeOffset[it] = std::size_t(tx)*grid->_ostride[Tp];
    }
    Scalar_v Cv0 = Zero();
    std::vector<Scalar_v> corr_v0(nTime_, Cv0);
    std::vector<std::vector<Scalar_v>> corr_v(nCorr_, corr_v0);

    startTimer("Four-quark contraction");
    {
        autoView(prop_v, prop, CpuRead);
        thread_for(ittg, nWork_, {
            const std::size_t ig = ittg % gamma1_.size();
            const std::size_t itt = ittg / gamma1_.size();
            const std::size_t isct = itt % types_.size();
            const std::size_t it = itt / types_.size();
            const std::size_t itg = ig + gamma1_.size()*isct;
            const int sctype = types_[isct];
            const auto &gvec1 = gamma1_[ig];
            const auto &gvec2 = gamma2_[ig];
            for (std::size_t ix3d = 0; ix3d < nSpace_; ++ix3d)
            {
                const std::size_t ix1 = ix3d + nSpace_*it;
                const std::size_t ix2 = timeOffset[it] + (ix3d/e2)*stride + ix3d%e2;
                for (std::size_t igg = 0; igg < gvec1.size(); ++igg)
                {
                    SpinColourMatrix_v WM1 = mat[ix1] * Gamma(gvec1[igg]);
                    SpinColourMatrix_v WM2 = prop_v[ix2] * Gamma(gvec2[igg]);
                    Scalar_v val = Zero();
                    // Preserve the four component contractions of the MT module.
                    if (sctype == 0)
                    {
                        val = trace(WM1) * trace(WM2);
                    }
                    else if (sctype == 1)
                    {
                        for (int s1=0; s1<Ns; ++s1)
                        for (int s2=0; s2<Ns; ++s2)
                        for (int c1=0; c1<Nc; ++c1)
                        for (int c2=0; c2<Nc; ++c2)
                        {
                            val()()() += WM1()(s1,s2)(c1,c2) * WM2()(s2,s1)(c2,c1);
                        }
                    }
                    else if (sctype == 2)
                    {
                        ColourMatrix_v CM1 = Zero();
                        ColourMatrix_v CM2 = Zero();
                        for (int s1=0; s1<Ns; ++s1)
                        for (int c1=0; c1<Nc; ++c1)
                        for (int c2=0; c2<Nc; ++c2)
                        {
                            CM1()()(c1,c2) += WM1()(s1,s1)(c1,c2);
                            CM2()()(c1,c2) += WM2()(s1,s1)(c1,c2);
                        }
                        for (int c1=0; c1<Nc; ++c1)
                        for (int c2=0; c2<Nc; ++c2)
                        {
                            val()()() += CM1()()(c1,c2) * CM2()()(c2,c1);
                        }
                    }
                    else if (sctype == 3)
                    {
                        SpinMatrix_v SM1 = Zero();
                        SpinMatrix_v SM2 = Zero();
                        for (int s1=0; s1<Ns; ++s1)
                        for (int s2=0; s2<Ns; ++s2)
                        for (int c1=0; c1<Nc; ++c1)
                        {
                            SM1()(s1,s2)() += WM1()(s1,s2)(c1,c1);
                            SM2()(s1,s2)() += WM2()(s1,s2)(c1,c1);
                        }
                        for (int s1=0; s1<Ns; ++s1)
                        for (int s2=0; s2<Ns; ++s2)
                        {
                            val()()() += SM1()(s1,s2)() * SM2()(s2,s1)();
                        }
                    }
                    corr_v[itg][it]()()() += val()()();
                }
            }
        });
    }
    stopTimer("Four-quark contraction");

    startTimer("Spatial reduction");
    const ComplexD C0(0., 0.);
    std::vector<ComplexD> corr0(nTime_, C0);
    std::vector<std::vector<ComplexD>> corr(nCorr_, corr0);
    thread_for(ittg, nWork_, {
        const std::size_t ig = ittg % gamma1_.size();
        const std::size_t itt = ittg / gamma1_.size();
        const std::size_t isct = itt % types_.size();
        const std::size_t it = itt / types_.size();
        const std::size_t itg = ig + gamma1_.size()*isct;
        ExtractBuffer<Scalar_s> extracted(Nsimd);
        extract(corr_v[itg][it], extracted);
        for (int isimd=0; isimd<Nsimd; ++isimd)
        {
            corr[itg][it] = corr[itg][it] + extracted[isimd];
        }
    });
    // Collective, outside OpenMP and outside the boss-only output block.
    for (std::size_t i = 0; i < nWork_; ++i)
    {
        const std::size_t it = i % nTime_;
        const std::size_t itg = i / nTime_;
        grid->GlobalSum(corr[itg][it]);
    }
    stopTimer("Spatial reduction");

    const std::string filename = par().output + "/"
        + ModuleBase::resultFilename(getName(), vm().getTrajectory());
    LOG(Message) << "Saving correlator to '" << filename << "'" << std::endl;
    startTimer("Write correlators");
    if (grid->IsBoss())
    {
        makeFileDir(filename);
        ResultWriter writer(filename);
        for (std::size_t itg = 0; itg < nCorr_; ++itg)
        {
            A2AFourQuarkContractionPropResult out;
            const std::size_t ig = itg % gamma1_.size();
            const std::size_t isct = itg / gamma1_.size();
            out.correlator = corr[itg];
            const std::string dataSet = getName() + "_" + nameg1_[ig] + "_"
                + nameg2_[ig] + "_sort" + std::to_string(types_[isct]);
            write(writer, dataSet, out);
        }
    }
    stopTimer("Write correlators");
}

END_MODULE_NAMESPACE
END_HADRONS_NAMESPACE
#endif // Hadrons_MContraction_A2AFourQuarkContractionProp_hpp_
