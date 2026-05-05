#ifndef Hadrons_EigenPackLazy_hpp_
#define Hadrons_EigenPackLazy_hpp_

#include <Hadrons/Global.hpp>

BEGIN_HADRONS_NAMESPACE

template <typename FineField, typename CoarsePack>
class LazilyDecompressedEigenPack
{
private:
    class LazyDecompressionVector
    {
    private:
        CoarsePack& coarsePack;
        size_t    _size;
        GridBase* grid;
    
    public:
        LazyDecompressionVector(CoarsePack& coarsePack, unsigned int size, GridBase* grid)
            : coarsePack{coarsePack}
            , _size     {size}
            , grid      {grid}
        {}

        FineField operator[](size_t i) const
        {
            typename CoarsePack::CoarseField& coarseEvec = coarsePack.evecCoarse[i];
	    GridBase* gridEvec = coarsePack.evec[0].Grid();
            typename CoarsePack::Field decompressedField(gridEvec);
            LOG(Message) << "Block-Promote " << i << std::endl;
	    {
	      GridBase* grid1 = decompressedField.Grid();
	      GridBase* grid2 = coarsePack.evec[0].Grid();
	      GridBase* grid3 = coarseEvec.Grid();
	      std::cout << "Grid0: " << grid  << std::endl;
	      std::cout << "Grid1: " << grid1 << std::endl;
	      std::cout << "Grid2: " << grid2 << std::endl;
	      std::cout << "Grid3: " << grid3 << std::endl;
	    }
            blockPromote(coarsePack.evecCoarse[i], decompressedField, coarsePack.evec);
            if constexpr (!std::is_same_v<FineField, typename CoarsePack::Field>)
            {
                LOG(Message) << "Precision change... " << i << std::endl;
                FineField out(grid);
                precisionChange(decompressedField, out);
                LOG(Message) << "Precision change successful!" << std::endl;
                return out;
            }
            else
            {
                return decompressedField;
            }
        }

        size_t size()
        {
            return this->_size;
        }
    };

private:
    CoarsePack& coarsePack;
public:
    LazyDecompressionVector evec;
    std::vector<RealD>& eval;
    
    LazilyDecompressedEigenPack(void)          = default;
    virtual ~LazilyDecompressedEigenPack(void) = default;

    LazilyDecompressedEigenPack(CoarsePack& coarsePack, unsigned int size, GridBase* grid)
        : coarsePack{coarsePack}
        , evec      {coarsePack, size, grid}
        , eval      {coarsePack.eval}
    {}
};

template <typename FImplFine, int nBasis, typename FImplCoarse=FImplFine, typename FImplCoarseIO=FImplCoarse>
using LazilyDecompressedFermionEigenPack = LazilyDecompressedEigenPack<
    typename FImplFine::FermionField,
    CoarseFermionEigenPack<FImplCoarse, nBasis, FImplCoarseIO>>;

END_HADRONS_NAMESPACE

#endif // Hadrons_EigenPackLazy_hpp_
