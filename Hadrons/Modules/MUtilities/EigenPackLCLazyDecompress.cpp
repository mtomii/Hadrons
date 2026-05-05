#include <Hadrons/Modules/MUtilities/EigenPackLCLazyDecompress.hpp>

using namespace Grid;
using namespace Hadrons;
using namespace MUtilities;

template class HADRONS_NAMESPACE::MUtilities::TEigenPackLCLazyDecompress<FIMPL, HADRONS_DEFAULT_LANCZOS_NBASIS>;
template class HADRONS_NAMESPACE::MUtilities::TEigenPackLCLazyDecompress<FIMPLF, HADRONS_DEFAULT_LANCZOS_NBASIS>;
template class HADRONS_NAMESPACE::MUtilities::TEigenPackLCLazyDecompress<FIMPL, HADRONS_DEFAULT_LANCZOS_NBASIS, FIMPLF>;
template class HADRONS_NAMESPACE::MUtilities::TEigenPackLCLazyDecompress<FIMPL,  120>;
template class HADRONS_NAMESPACE::MUtilities::TEigenPackLCLazyDecompress<FIMPLF, 120>;
template class HADRONS_NAMESPACE::MUtilities::TEigenPackLCLazyDecompress<FIMPL,  120, FIMPLF>;
template class HADRONS_NAMESPACE::MUtilities::TEigenPackLCLazyDecompress<FIMPL, 250>;
template class HADRONS_NAMESPACE::MUtilities::TEigenPackLCLazyDecompress<FIMPLF, 250>;
template class HADRONS_NAMESPACE::MUtilities::TEigenPackLCLazyDecompress<FIMPL, 250, FIMPLF>;
template class HADRONS_NAMESPACE::MUtilities::TEigenPackLCLazyDecompress<FIMPL, 1000>;
template class HADRONS_NAMESPACE::MUtilities::TEigenPackLCLazyDecompress<FIMPLF, 1000>;
template class HADRONS_NAMESPACE::MUtilities::TEigenPackLCLazyDecompress<FIMPL, 1000, FIMPLF>;
template class HADRONS_NAMESPACE::MUtilities::TEigenPackLCLazyDecompress<ZFIMPL, HADRONS_DEFAULT_LANCZOS_NBASIS>;
template class HADRONS_NAMESPACE::MUtilities::TEigenPackLCLazyDecompress<ZFIMPL, 250>;
template class HADRONS_NAMESPACE::MUtilities::TEigenPackLCLazyDecompress<ZFIMPLF, 250>;
template class HADRONS_NAMESPACE::MUtilities::TEigenPackLCLazyDecompress<ZFIMPL, 250, ZFIMPLF>;
template class HADRONS_NAMESPACE::MUtilities::TEigenPackLCLazyDecompress<ZFIMPL, 1000>;
template class HADRONS_NAMESPACE::MUtilities::TEigenPackLCLazyDecompress<ZFIMPLF, 1000>;
template class HADRONS_NAMESPACE::MUtilities::TEigenPackLCLazyDecompress<ZFIMPL, 1000, ZFIMPLF>;
