#ifndef Hadrons_MUtilities_EigenPackLCLazyDecompress_hpp_
#define Hadrons_MUtilities_EigenPackLCLazyDecompress_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/EigenPack.hpp>
#include <Hadrons/EigenPackLazy.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         EigenPackLCLazyDecompress                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MUtilities)

class EigenPackLCLazyDecompressPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(EigenPackLCLazyDecompressPar,
                                    std::string,  epack,
                                    std::string,  blockSize,
				    std::string, action,
                                    unsigned int, size,
                                    unsigned int, Ls);
};


template <typename FImpl, int nBasis, typename CoarseFImpl=FImpl, typename FImplIo=FImpl, typename CoarseFImplIo=CoarseFImpl>
class TEigenPackLCLazyDecompress: public Module<EigenPackLCLazyDecompressPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    typedef LazilyDecompressedFermionEigenPack<FImpl, nBasis, CoarseFImpl, CoarseFImplIo> LazyPack;
    typedef CoarseFermionEigenPack<CoarseFImpl, nBasis, CoarseFImplIo> CoarsePack;
    typedef typename CoarsePack::Field Field;
    // constructor
    TEigenPackLCLazyDecompress(const std::string name);
    // destructor
    virtual ~TEigenPackLCLazyDecompress(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    virtual DependencyMap getObjectDependencies(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(EigenPackLCLazyDecompress, ARG(TEigenPackLCLazyDecompress<FIMPL, HADRONS_DEFAULT_LANCZOS_NBASIS>), MUtilities);
MODULE_REGISTER_TMP(EigenPackLCLazyDecompressF, ARG(TEigenPackLCLazyDecompress<FIMPLF, HADRONS_DEFAULT_LANCZOS_NBASIS>), MUtilities);
MODULE_REGISTER_TMP(EigenPackLCLazyDecompressIo32, ARG(TEigenPackLCLazyDecompress<FIMPL, HADRONS_DEFAULT_LANCZOS_NBASIS, FIMPLF>), MUtilities);
MODULE_REGISTER_TMP(EigenPackLCLazyDecompress250, ARG(TEigenPackLCLazyDecompress<FIMPL, 250>), MUtilities);
MODULE_REGISTER_TMP(EigenPackLCLazyDecompress250Io32, ARG(TEigenPackLCLazyDecompress<FIMPL, 250, FIMPLF>), MUtilities);
MODULE_REGISTER_TMP(EigenPackLCLazyDecompress1000, ARG(TEigenPackLCLazyDecompress<FIMPL, 1000>), MUtilities);
MODULE_REGISTER_TMP(EigenPackLCLazyDecompress1000Io32, ARG(TEigenPackLCLazyDecompress<FIMPL, 1000, FIMPLF>), MUtilities);
MODULE_REGISTER_TMP(ZEigenPackLCLazyDecompress, ARG(TEigenPackLCLazyDecompress<ZFIMPL, HADRONS_DEFAULT_LANCZOS_NBASIS>), MUtilities);
MODULE_REGISTER_TMP(ZEigenPackLCLazyDecompress250, ARG(TEigenPackLCLazyDecompress<ZFIMPL, 250>), MUtilities);
MODULE_REGISTER_TMP(ZEigenPackLCLazyDecompress250Io32, ARG(TEigenPackLCLazyDecompress<ZFIMPL, 250, ZFIMPLF>), MUtilities);
MODULE_REGISTER_TMP(ZEigenPackLCLazyDecompress1000, ARG(TEigenPackLCLazyDecompress<ZFIMPL, 1000>), MUtilities);
MODULE_REGISTER_TMP(ZEigenPackLCLazyDecompress1000Io32, ARG(TEigenPackLCLazyDecompress<ZFIMPL, 1000, ZFIMPLF>), MUtilities);

/******************************************************************************
 *                 TEigenPackLCLazyDecompress implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl, int nBasis, typename CoarseFImpl, typename FImplIo, typename CoarseFImplIo>
TEigenPackLCLazyDecompress<FImpl, nBasis, CoarseFImpl, FImplIo, CoarseFImplIo>::TEigenPackLCLazyDecompress(const std::string name)
: Module<EigenPackLCLazyDecompressPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl, int nBasis, typename CoarseFImpl, typename FImplIo, typename CoarseFImplIo>
std::vector<std::string> TEigenPackLCLazyDecompress<FImpl, nBasis, CoarseFImpl, FImplIo, CoarseFImplIo>::getInput(void)
{
    std::vector<std::string> in;
    
    return in;
}

template <typename FImpl, int nBasis, typename CoarseFImpl, typename FImplIo, typename CoarseFImplIo>
std::vector<std::string> TEigenPackLCLazyDecompress<FImpl, nBasis, CoarseFImpl, FImplIo, CoarseFImplIo>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}


template <typename FImpl, int nBasis, typename CoarseFImpl, typename FImplIo, typename CoarseFImplIo>
DependencyMap TEigenPackLCLazyDecompress<FImpl, nBasis, CoarseFImpl, FImplIo, CoarseFImplIo>::getObjectDependencies(void)
{
    DependencyMap dep;
    
    dep.insert({par().epack, getName()});

    return dep;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl, int nBasis, typename CoarseFImpl, typename FImplIo, typename CoarseFImplIo>
void TEigenPackLCLazyDecompress<FImpl, nBasis, CoarseFImpl, FImplIo, CoarseFImplIo>::setup(void)
{
  std::cout << "START SETUP" << std::endl;
    auto blockSize   = strToVec<int>(par().blockSize);
  std::cout << "AAA" << std::endl;
    auto& coarsePack = envGet(CoarsePack, par().epack);
  std::cout << "BBB" << std::endl;
  //GridBase* grid   = envGetGrid(CoarsePack::Field, par().Ls);
  //GridBase* grid = coarsePack.evec[0].Grid();
  //GridBase* grid        = getGrid<Field>(true, par().Ls);
  auto        &action    = envGet(FMat, par().action);
  GridBase* frbGrid = action.FermionRedBlackGrid();
  std::cout << "CCC" << std::endl;
    envCreate(LazyPack, getName(), par().Ls, coarsePack, par().size, frbGrid);
    //envCreate(LazyPack, getName(), 1, coarsePack, par().size, grid); commented by MT
  std::cout << "END SETUP" << std::endl;
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl, int nBasis, typename CoarseFImpl, typename FImplIo, typename CoarseFImplIo>
void TEigenPackLCLazyDecompress<FImpl, nBasis, CoarseFImpl, FImplIo, CoarseFImplIo>::execute(void)
{
  std::cout << "START EXE" << std::endl;
    const LazyPack& lpack = envGet(LazyPack, getName());
  std::cout << "AAA" << std::endl;
    for (size_t i=0; i < par().size; ++i)
    {
      if ( i==0 ) std::cout << "BBB" << std::endl;
        LOG(Message) << "Getting vector " << i << std::endl;
	if ( i==0 ) std::cout << "CCC" << std::endl;
        //auto res = lpack.evec[i];
	if ( i==0 ) std::cout << "DDD" << std::endl;
        LOG(Message) << "Norm " << i << ": " << norm2(lpack.evec[i]) << std::endl;
        //LOG(Message) << "Norm " << i << ": " << norm2(res) << std::endl;
    }
  std::cout << "END EXE" << std::endl;
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MUtilities_EigenPackLCLazyDecompress_hpp_
