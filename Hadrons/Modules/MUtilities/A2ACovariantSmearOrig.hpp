#ifndef Hadrons_MUtilities_A2ACovariantSmearOrig_hpp_
#define Hadrons_MUtilities_A2ACovariantSmearOrig_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Grid/qcd/utils/CovariantSmearing.h>
#include <Hadrons/A2AVectors.hpp>

BEGIN_HADRONS_NAMESPACE

/* Comments by Masaaki
 This module Looks overwriting the input A2AVectors instead of Outputting as new object

 'auto &a2aTmp = envGet(FermionField, getName());'
 vs
 'auto &a2aVec = envGet(std::vector<FermionField>, par().a2aVectors);'
 indicates inconsistent format of the output

 The module makes sense when we mean to output smeared vectors and to run another job to use them later on as it can save memory
*/
/******************************************************************************
 *                         A2ACovariantSmearOrig                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MUtilities)

class A2ACovariantSmearOrigPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(A2ACovariantSmearOrigPar,
                                    std::string, a2aVectors,
                                    std::string, gauge,
                                    double, alpha,
                                    unsigned int, N,
                                    unsigned int, orthog_axis,
                                    std::string, output,
                                    bool, multiFile);
};

template <typename GImpl, typename FImpl>
class TA2ACovariantSmearOrig: public Module<A2ACovariantSmearOrigPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
public:
    // constructor
    TA2ACovariantSmearOrig(const std::string name);
    // destructor
    virtual ~TA2ACovariantSmearOrig(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(A2ACovariantSmearOrig, ARG(TA2ACovariantSmearOrig<GIMPL, FIMPL>), MUtilities);

/******************************************************************************
 *                 TA2ACovariantSmearOrig implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename GImpl, typename FImpl>
TA2ACovariantSmearOrig<GImpl, FImpl>::TA2ACovariantSmearOrig(const std::string name)
: Module<A2ACovariantSmearOrigPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename GImpl, typename FImpl>
std::vector<std::string> TA2ACovariantSmearOrig<GImpl, FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().a2aVectors, par().gauge};
    
    return in;
}

template <typename GImpl, typename FImpl>
std::vector<std::string> TA2ACovariantSmearOrig<GImpl, FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename GImpl, typename FImpl>
void TA2ACovariantSmearOrig<GImpl, FImpl>::setup(void)
{
	envCreateLat(FermionField, getName());    
}

// execution ///////////////////////////////////////////////////////////////////
template <typename GImpl, typename FImpl>
void TA2ACovariantSmearOrig<GImpl, FImpl>::execute(void)
{
	auto &a2aVec = envGet(std::vector<FermionField>, par().a2aVectors);
	unsigned int Ni = a2aVec.size();

	const auto &U = envGet(GaugeField, par().gauge);

	std::vector<LatticeColourMatrix> Umu(Nd,env().getGrid());
	for(int i=0;i<Nd;i++)
	{
        Umu[i] = PeekIndex<LorentzIndex>(U,i);
	}
	
    auto &a2aTmp = envGet(FermionField, getName());
	
	RealD width = par().alpha;
	unsigned int iterations = par().N;
	unsigned int orthog_axis = par().orthog_axis;

    LOG(Message) << "Starting A2A covariant smearing." << std::endl;
    LOG(Message) << par().a2aVectors << " has size " << Ni << std::endl;
    
	CovariantSmearing<GImpl> CovSmear;
	
	for (int i = 0; i < Ni; i++)
	{
		a2aTmp = a2aVec[i];

        CovSmear.template GaussianSmear<FermionField>(Umu, a2aTmp, width, iterations, orthog_axis);
		a2aVec[i] = a2aTmp;
	}

    // I/O if necessary
    if (!par().output.empty())
    {
        startTimer("I/O");
        A2AVectorsIo::write(par().output, a2aVec, par().multiFile, vm().getTrajectory());
        stopTimer("I/O");
    }
    
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MUtilities_A2ACovariantSmearOrig_hpp_
