/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: Hadrons/Modules/MContraction/A2ALoopNew.hpp

Copyright (C) 2015-2026

Author: Jonas Hildebrand <jonas.hildebrand@uconn.edu>

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
#ifndef Hadrons_MContraction_A2ALoopNew_hpp_
#define Hadrons_MContraction_A2ALoopNew_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Grid/qcd/utils/A2Autils.h>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *  Closed quark loop propagator from A2A vectors:                            *
 *                                                                            *
 *      loop = sum_k outerProduct(V[k], W[k])                                 *
 *                                                                            *
 *  Superset of MContraction::A2ALoop: it handles both the fully expanded W   *
 *  representation and the dense (time-undiluted) one, blocks the mode sum so  *
 *  device residency is bounded independently of the mode count, and is meant  *
 *  to be run in its own job that writes the loop to disk for later consumers  *
 *  (e.g. MContraction::A2AExtendedMesonField) rather than building it inline. *
 *                                                                            *
 *  Array layout, matching MIO::LoadCombinedA2AVecsV / LoadCombinedA2AVecsW:   *
 *                                                                            *
 *      V : nLow low modes, then per hit nt*Nsc high modes, index t*Nsc + sc   *
 *      W : nLow low modes, then per hit either                                *
 *            Nsc dense slots (index sc)                    [dense]            *
 *            nt*Nsc expanded modes (index t*Nsc + sc)      [expanded]         *
 *                                                                            *
 *  Which W representation is in use is deduced from the array sizes and       *
 *  logged; the two differ by a factor nt in the high-mode block, so there is  *
 *  no ambiguity. nLow must be given because it cannot be inferred alone.      *
 *                                                                            *
 *  Expanded W is handled by plain index-for-index pairing over the whole      *
 *  array. Dense W pairs low modes the same way, then for each hit gathers the *
 *  surviving timeslice of each expanded V mode into Nsc dense fields (see     *
 *  A2AExtendedMesonField::GatherDenseTimeslices) and contracts those against  *
 *  that hit's Nsc dense W slots. The gather buffer is reused across hits, so  *
 *  nothing here scales with the hit count.                                    *
 *                                                                            *
 *  No normalization is applied: by convention the 1/nHit hit-average factor   *
 *  lives on the V side and is already baked in by the loader.                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MContraction)

class A2ALoopNewPar: Serializable
{
public:
    // nLow:  size of the low-mode block common to both arrays; 0 for flavors
    //        without low modes (strange, charm).
    // block: how many modes are contracted per kernel call, i.e. how many
    //        field views are open on the device at once. Bounds device memory
    //        independently of the mode count; the high-mode phase in the dense
    //        case is only Nsc modes and ignores it.
    GRID_SERIALIZABLE_CLASS_MEMBERS(A2ALoopNewPar,
                                    std::string,  left,
                                    std::string,  right,
                                    unsigned int, nLow,
                                    unsigned int, block);
};

template <typename FImpl>
class TA2ALoopNew: public Module<A2ALoopNewPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
public:
    // constructor
    TA2ALoopNew(const std::string name);
    // destructor
    virtual ~TA2ALoopNew(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
    // geometry deduced in setup, reused in execute
    bool denseW_;
    int  nt_, nsc_, nHit_;
};

MODULE_REGISTER_TMP(A2ALoopNew, TA2ALoopNew<FIMPL>, MContraction);

/******************************************************************************
 *                       TA2ALoopNew implementation                           *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TA2ALoopNew<FImpl>::TA2ALoopNew(const std::string name)
: Module<A2ALoopNewPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TA2ALoopNew<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().left, par().right};

    return in;
}

template <typename FImpl>
std::vector<std::string> TA2ALoopNew<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};

    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TA2ALoopNew<FImpl>::setup(void)
{
    auto &left  = envGet(std::vector<FermionField>, par().left);
    auto &right = envGet(std::vector<FermionField>, par().right);

    int nLow = (int)par().nLow;

    nt_  = env().getDim().back();
    nsc_ = Ns*FImpl::Dimension;

    if (nLow > (int)left.size() || nLow > (int)right.size())
    {
        HADRONS_ERROR(Argument, "nLow larger than the vector arrays");
    }

    // The two W representations differ by a factor nt in the high-mode block,
    // so the sizes alone say which one this is.
    denseW_ = (left.size() != right.size());
    if (denseW_)
    {
        int highW = (int)right.size() - nLow;

        if (highW <= 0 || highW % nsc_ != 0)
        {
            HADRONS_ERROR(Size, "dense W high-mode block is not a multiple of "
                                + std::to_string(nsc_));
        }
        nHit_ = highW/nsc_;
        if ((int)left.size() != nLow + nHit_*nt_*nsc_)
        {
            HADRONS_ERROR(Size, "V has " + std::to_string(left.size())
                                + " modes, expected "
                                + std::to_string(nLow + nHit_*nt_*nsc_)
                                + " for nLow=" + std::to_string(nLow)
                                + " nHit=" + std::to_string(nHit_)
                                + " nt=" + std::to_string(nt_)
                                + " Nsc=" + std::to_string(nsc_));
        }

        // Gather buffer, reused across hits: Nsc fields regardless of nHit.
        // Built empty and filled in place so the vector (count, value)
        // constructor does not copy-broadcast a temporary into every slot.
        auto grid = envGetGrid(FermionField);

        envTmp(std::vector<FermionField>, "vtilde", 1, 0, grid);
        envGetTmp(std::vector<FermionField>, vtilde);

        vtilde.reserve(nsc_);
        for (int sc = 0; sc < nsc_; ++sc)
        {
            vtilde.emplace_back(grid);
        }
    }
    else
    {
        int high = (int)left.size() - nLow;

        if (high < 0 || (nt_*nsc_ > 0 && high % (nt_*nsc_) != 0))
        {
            HADRONS_ERROR(Size, "expanded high-mode block is not a multiple of "
                                + std::to_string(nt_*nsc_));
        }
        nHit_ = high/(nt_*nsc_);
    }

    envCreateLat(PropagatorField, getName());
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TA2ALoopNew<FImpl>::execute(void)
{
    // Qualified: MODULE_REGISTER declares a class A2AExtendedMesonField in
    // this same namespace, which would otherwise shadow the Grid template.
    typedef Grid::A2AExtendedMesonField<FImpl> EMF;

    auto &left  = envGet(std::vector<FermionField>, par().left);
    auto &right = envGet(std::vector<FermionField>, par().right);
    auto &loop  = envGet(PropagatorField, getName());

    int nLow  = (int)par().nLow;
    int block = (int)par().block;

    if (block <= 0)
    {
        HADRONS_ERROR(Argument, "block must be positive");
    }

    LOG(Message) << "Computing A2A loop propagator from '" << par().left
                 << "' and '" << par().right << "'" << std::endl;
    LOG(Message) << "W representation: " << (denseW_ ? "dense" : "expanded")
                 << " (V " << left.size() << " modes, W " << right.size()
                 << " modes)" << std::endl;
    LOG(Message) << "nLow " << nLow << ", nHit " << nHit_ << ", nt " << nt_
                 << ", Nsc " << nsc_ << ", block " << block << std::endl;

    // Every phase below accumulates, so the sum starts here and nowhere else.
    loop = Zero();

    if (!denseW_)
    {
        // Expanded W pairs index for index over the whole array, low and high
        // alike, so there is nothing to compress and one blocked sweep does it.
        int N = (int)left.size();

        startTimer("All modes");
        for (int i = 0; i < N; i += block)
        {
            EMF::LoopPropagator(loop, left, i, right, i, std::min(block, N - i));
        }
        stopTimer("All modes");
    }
    else
    {
        envGetTmp(std::vector<FermionField>, vtilde);

        // Low modes carry no dilution, so they pair directly in both
        // representations.
        startTimer("Low modes");
        for (int i = 0; i < nLow; i += block)
        {
            EMF::LoopPropagator(loop, left, i, right, i, std::min(block, nLow - i));
        }
        stopTimer("Low modes");

        // One hit at a time: compress that hit's expanded V block down to Nsc
        // fields, contract them against the hit's Nsc dense W slots, then reuse
        // the buffer for the next hit.
        for (int hit = 0; hit < nHit_; ++hit)
        {
            startTimer("Gather");
            EMF::GatherDenseTimeslices(vtilde, left, nLow + hit*nt_*nsc_,
                                       0, nt_*nsc_, nsc_);
            stopTimer("Gather");

            startTimer("High modes");
            EMF::LoopPropagator(loop, vtilde, 0, right, nLow + hit*nsc_, nsc_);
            stopTimer("High modes");
        }
    }

    LOG(Message) << "Loop propagator norm2 = " << norm2(loop) << std::endl;
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_A2ALoopNew_hpp_
