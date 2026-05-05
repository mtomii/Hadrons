#ifndef Hadrons_MUtilities_A2ACovariantSmear_hpp_
#define Hadrons_MUtilities_A2ACovariantSmear_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Grid/qcd/utils/CovariantSmearing.h>
#include <Hadrons/A2AVectors.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         A2ACovariantSmear                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MUtilities)

class A2ACovariantSmearPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(A2ACovariantSmearPar,
                                    std::string, a2aVectors,
                                    std::string, gauge,
                                    double, alpha,
                                    unsigned int, N,
                                    unsigned int, orthog_axis,
                                    std::string, output,
                                    bool, multiFile);
};

template <typename GImpl, typename FImpl>
class TA2ACovariantSmear: public Module<A2ACovariantSmearPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
public:
    // constructor
    TA2ACovariantSmear(const std::string name);
    // destructor
    virtual ~TA2ACovariantSmear(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(A2ACovariantSmear, ARG(TA2ACovariantSmear<GIMPL, FIMPL>), MUtilities);

/******************************************************************************
 *                 TA2ACovariantSmear implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename GImpl, typename FImpl>
TA2ACovariantSmear<GImpl, FImpl>::TA2ACovariantSmear(const std::string name)
: Module<A2ACovariantSmearPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename GImpl, typename FImpl>
std::vector<std::string> TA2ACovariantSmear<GImpl, FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().a2aVectors, par().gauge};
    
    return in;
}

template <typename GImpl, typename FImpl>
std::vector<std::string> TA2ACovariantSmear<GImpl, FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename GImpl, typename FImpl>
void TA2ACovariantSmear<GImpl, FImpl>::setup(void)
{
  auto &a2aVec = envGet(std::vector<FermionField>, par().a2aVectors);
  envCreate(std::vector<FermionField>, getName(), 1, a2aVec.size(),
	    envGetGrid(FermionField));
}

// execution ///////////////////////////////////////////////////////////////////
template <typename GImpl, typename FImpl>
void TA2ACovariantSmear<GImpl, FImpl>::execute(void)
{
  const auto &a2aVec = envGet(std::vector<FermionField>, par().a2aVectors);
  unsigned int Ni = a2aVec.size();

  const auto &U = envGet(GaugeField, par().gauge);

  std::vector<LatticeColourMatrix> Umu(Nd,env().getGrid());
  for(int i=0;i<Nd;i++)
  {
    Umu[i] = PeekIndex<LorentzIndex>(U,i);
  }
	
  auto &a2aSmr = envGet(std::vector<FermionField>, getName());
	
  RealD width = par().alpha;
  unsigned int iterations = par().N;
  unsigned int orthog_axis = par().orthog_axis;

  LOG(Message) << "Starting A2A covariant smearing." << std::endl;
  LOG(Message) << par().a2aVectors << " has size " << Ni << std::endl;
    
  CovariantSmearing<GImpl> CovSmear;
	
  for (int i = 0; i < Ni; i++)
  {
    a2aSmr[i] = a2aVec[i];

    CovSmear.template GaussianSmear<FermionField>(Umu, a2aSmr[i], width, iterations, orthog_axis);
  }

  // I/O if necessary,
  // !!! Probably not needed as smearing is cheap
  // !!!  but consider binned IO if implement
  if (!par().output.empty())
  {
    startTimer("I/O");
    A2AVectorsIo::write(par().output, a2aSmr, par().multiFile, vm().getTrajectory());
    stopTimer("I/O");
  }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MUtilities_A2ACovariantSmear_hpp_
