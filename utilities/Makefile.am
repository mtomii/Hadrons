AM_LDFLAGS += -L$(top_builddir)/Hadrons
AM_CXXFLAGS += -I$(top_srcdir)

bin_PROGRAMS = HadronsXmlRun HadronsXmlValidate 

EXTRA_PROGRAMS = \
  BatchDeflationBenchmark    \
  HadronsContractor          \
  HadronsContractorBenchmark \
  HadronsXmlRun              \
  HadronsXmlValidate         \
  HadronsFermionEP64To32     \
  K2pipiContractor           \
  BubbleContractor           \
  FlexibleContractor         \
  HadronsFermionEP64To32 

CLEANFILES = $(EXTRA_PROGRAMS)

BatchDeflationBenchmark_SOURCES = BatchDeflationBenchmark.cpp
BatchDeflationBenchmark_LDADD = -lHadrons -lGrid

HadronsXmlRun_SOURCES = HadronsXmlRun.cpp
HadronsXmlRun_LDADD   = -lHadrons -lGrid

HadronsXmlValidate_SOURCES = HadronsXmlValidate.cpp
HadronsXmlValidate_LDADD   = -lHadrons -lGrid

HadronsFermionEP64To32_SOURCES  = EigenPackCast.cpp
HadronsFermionEP64To32_CXXFLAGS = $(AM_CXXFLAGS) -DFIN=WilsonImplD::FermionField -DFOUT=WilsonImplF::FermionField
HadronsFermionEP64To32_LDADD    = -lHadrons -lGrid

HadronsContractor_SOURCES = Contractor.cpp
HadronsContractor_LDADD   = -lHadrons -lGrid

FlexibleContractor_SOURCES = FlexibleContractor.cpp
FlexibleContractor_LDADD   = -lHadrons -lGrid

K2pipiContractor_SOURCES = K2pipiContractor.cpp
K2pipiContractor_LDADD   = -lHadrons -lGrid

BubbleContractor_SOURCES = BubbleContractor.cpp
BubbleContractor_LDADD   = -lHadrons -lGrid

HadronsContractorBenchmark_SOURCES = ContractorBenchmark.cpp
HadronsContractorBenchmark_LDADD   = -lHadrons -lGrid
