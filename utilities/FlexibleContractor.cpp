/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/Utilities/FlexibleContractor.cc

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
#include <Hadrons/Global.hpp>
#include <Hadrons/A2AMatrix.hpp>
#include <Hadrons/DiskVector.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/TimerArray.hpp>

using namespace Grid;
using namespace Hadrons;

#define TIME_MOD(t) (((t) + par.global.nt) % par.global.nt)

namespace FlexibleContractor
{
  class TrajRange: Serializable
  {
  public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(TrajRange,
				    unsigned int, start,
				    unsigned int, end,
				    unsigned int, step);
  };
    
  class GlobalPar: Serializable
  {
  public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(GlobalPar,
				    TrajRange, trajCounter,
				    unsigned int, nt,
				    std::string, diskVectorDir,
				    std::string, output);
  };

  class A2AMatrixPar: Serializable
  {
  public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(A2AMatrixPar,
				    std::string, file,
				    std::string, dataset,
				    unsigned int, cacheSize,
				    std::string, name);
  };

  class TermsPar: Serializable
  {
  public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(TermsPar,
				    std::string, term,
				    bool, tdps,
				    std::string, time);
  };

  class ProductPar: Serializable
  {
  public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(ProductPar,
				    std::vector<FlexibleContractor::TermsPar>, terms,
				    std::string, translations,
				    unsigned int, min_t,
				    unsigned int, max_t,
				    bool, translationAverage);
  };

  class CorrelatorResult: Serializable
  {
  public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(CorrelatorResult,
				    std::vector<FlexibleContractor::A2AMatrixPar>,  a2aMatrix,
				    ProductPar, contraction,
				    std::vector<unsigned int>, times,
				    std::vector<ComplexD>, correlator);
  };
}

struct FlexibleContractorPar
{
    FlexibleContractor::GlobalPar                  global;
    std::vector<FlexibleContractor::A2AMatrixPar>  a2aMatrix;
    std::vector<FlexibleContractor::ProductPar>    product;
};

void makeTimeSeq(std::vector<std::vector<unsigned int>> &timeSeq, 
                 const std::vector<std::set<unsigned int>> &times,
                 std::vector<unsigned int> &current,
                 const unsigned int depth)
{
    if (depth > 0)
    {
        for (auto t: times[times.size() - depth])
        {
            current[times.size() - depth] = t;
            makeTimeSeq(timeSeq, times, current, depth - 1);
        }
    }
    else
    {
        timeSeq.push_back(current);
    }
}

void makeTimeSeq(std::vector<std::vector<unsigned int>> &timeSeq, 
                 const std::vector<std::set<unsigned int>> &times)
{
    std::vector<unsigned int> current(times.size());

    makeTimeSeq(timeSeq, times, current, times.size());
}

std::set<unsigned int> parseTimeRange(const std::string str, const unsigned int nt)
{
  std::regex               rex("([0-9]+)|(([0-9]+)\\.\\.([0-9]+))");
  std::smatch              sm;
  std::vector<std::string> rstr = strToVec<std::string>(str);
  std::set<unsigned int>   tSet;

  for (auto &s: rstr)
  {
    std::regex_match(s, sm, rex);
    if (sm[1].matched)
    {
      unsigned int t;
            
      t = std::stoi(sm[1].str());
      if (t >= nt)
      {
	HADRONS_ERROR(Range, "time out of range (from expression '" + str + "')");
      }
      tSet.insert(t);
    }
    else if (sm[2].matched)
    {
      unsigned int ta, tb;

      ta = std::stoi(sm[3].str());
      tb = std::stoi(sm[4].str());
      if ((ta >= nt) or (tb >= nt))
      {
	HADRONS_ERROR(Range, "time out of range (from expression '" + str + "')");
      }
      for (unsigned int ti = ta; ti <= tb; ++ti)
      {
	tSet.insert(ti);
      }
    }
  }

  return tSet;
}

struct Sec
{
  Sec(const double usec)
  {
    seconds = usec/1.0e6;
  }
    
  double seconds;
};

inline std::ostream & operator<< (std::ostream& s, const Sec &&sec)
{
  s << std::setw(10) << sec.seconds << " sec";

  return s;
}

struct Flops
{
  Flops(const double flops, const double fusec)
  {
    gFlopsPerSec = flops/fusec/1.0e3;
  }
    
  double gFlopsPerSec;
};

inline std::ostream & operator<< (std::ostream& s, const Flops &&f)
{
  s << std::setw(10) << f.gFlopsPerSec << " GFlop/s";

  return s;
}

struct Bytes
{
  Bytes(const double bytes, const double busec)
  {
    gBytesPerSec = bytes/busec*1.0e6/1024/1024/1024;
  }
    
  double gBytesPerSec;
};

inline std::ostream & operator<< (std::ostream& s, const Bytes &&b)
{
  s << std::setw(10) << b.gBytesPerSec << " GB/s";

  return s;
}

int main(int argc, char* argv[])
{
  // parse command line
  std::string   parFilename;

  if (argc != 2)
  {
    std::cerr << "usage: " << argv[0] << " <parameter file>";
    std::cerr << std::endl;
        
    return EXIT_FAILURE;
  }
  parFilename = argv[1];

  // parse parameter file
  FlexibleContractorPar par;
  //    unsigned int  nMat,nCont;
  XmlReader     reader(parFilename);

  read(reader, "global",    par.global);
  read(reader, "a2aMatrix", par.a2aMatrix);
  read(reader, "product",   par.product);
  //    nMat  = par.a2aMatrix.size();
  //    nCont = par.product.size();

  // create diskvectors
  std::map<std::string, EigenDiskVector<ComplexD>> a2aMat;
  //    unsigned int                                     cacheSize;

  for (auto &p: par.a2aMatrix)
  {
    std::string dirName = par.global.diskVectorDir + "/" + p.name;

    a2aMat.emplace(p.name, EigenDiskVector<ComplexD>(dirName, par.global.nt, p.cacheSize));
  }

  // trajectory loop
  for (unsigned int traj = par.global.trajCounter.start; 
       traj < par.global.trajCounter.end; traj += par.global.trajCounter.step)
  {
    std::cout << ":::::::: Trajectory " << traj << std::endl;

    // load data
    for (auto &p: par.a2aMatrix)
    {
      std::string filename = p.file;
      double      t;
      //	    double  size;

      tokenReplace(filename, "traj", traj);
      std::cout << "======== Loading '" << filename << "'" << std::endl;

      A2AMatrixIo<HADRONS_A2AM_IO_TYPE> a2aIo(filename, p.dataset, par.global.nt);

      a2aIo.load(a2aMat.at(p.name), &t);
      std::cout << "Read " << a2aIo.getSize() << " bytes in " << t/1.0e6 
		<< " sec, " << a2aIo.getSize()/t*1.0e6/1024/1024 << " MB/s" << std::endl;
    }

    // contract
    //EigenDiskVector<ComplexD>::Matrix buf;

    for (auto &p: par.product)
    {
      std::vector<FlexibleContractor::TermsPar>                  terms = p.terms;
      std::vector<std::set<unsigned int>>    times;
      std::vector<std::vector<unsigned int>> timeSeq;
      std::set<unsigned int>                 translations;
      A2AMatrix<ComplexD>                    prod0, prod, buf, tmp;
      TimerArray                             tAr;
      double                                 fusec, busec, flops, bytes;
      double                                 fusecp, busecp, flopsp, bytesp;
      std::vector<A2AMatrixTr<ComplexD>>     lastTerm(par.global.nt);
      //	    double  tusec;
      FlexibleContractor::CorrelatorResult   result;

      if ( terms[0].tdps ) {
	std::cout << "You have to specify tdps=False for the 1st meson field of each product" << std::endl;
	assert(!terms[0].tdps);
      }

      tAr.startTimer("Total");
      std::cout << "======== Contraction tr(";
      for (unsigned int g = 0; g < terms.size(); ++g)
      {
	std::cout << terms[g].term << ((g == terms.size() - 1) ? ')' : '*');
	times.push_back(parseTimeRange(terms[g].time, par.global.nt));
      }
      std::cout << std::endl;
      for (auto &m: par.a2aMatrix)
      {
	if (std::find(result.a2aMatrix.begin(), result.a2aMatrix.end(), m) == result.a2aMatrix.end())
	{
	  result.a2aMatrix.push_back(m);
	  tokenReplace(result.a2aMatrix.back().file, "traj", traj);
	}
      }
      result.contraction = p;
      result.correlator.resize(par.global.nt, 0.);

      translations = parseTimeRange(p.translations, par.global.nt);
      makeTimeSeq(timeSeq, times);

      std::cout << "* Caching transposed last term " << std::endl;
      for (unsigned int t = 0; t < par.global.nt; ++t)
      {
	tAr.startTimer("Disk vector overhead");
	const A2AMatrix<ComplexD> &ref = a2aMat.at(terms.back().term)[t];
	tAr.stopTimer("Disk vector overhead");

	tAr.startTimer("Transpose caching");
	lastTerm[t].resize(ref.rows(), ref.cols());
	thread_for( j,ref.cols(),{
	  for (unsigned int i = 0; i < ref.rows(); ++i)
	  {
	    lastTerm[t](i, j) = ref(i, j);
	  }
	});
	tAr.stopTimer("Transpose caching");
      }
      bytes = par.global.nt*lastTerm[0].rows()*lastTerm[0].cols()*sizeof(ComplexD);
      std::cout << Sec(tAr.getDTimer("Transpose caching")) << " " 
		<< Bytes(bytes, tAr.getDTimer("Transpose caching")) << std::endl;

      std::string              fileStem = "", filename, dataSet;
      std::vector<FlexibleContractor::TermsPar> termsv = result.contraction.terms;
      for (unsigned int i = 0; i < termsv.size() ; i++)
      {
	fileStem += termsv[i].term + "_";
	if ( termsv[i].tdps ) {
	  fileStem += "x";
	} else {
	  fileStem += "y";
	}
	fileStem += termsv[i].time;
	if ( i < termsv.size() - 1 ) fileStem += "_";
      }
      const std::string dir = par.global.output;
      makeFileDir(dir);
      filename = dir + "/" + ModuleBase::resultFilename(fileStem, traj);
      //filename = dir + "/" + RESULT_FILE_NAME(fileStem, traj);
      std::cout << "Saving correlator to '" << filename << "'" << std::endl;
      ResultWriter writer(filename);

      for (unsigned int i = 0; i < timeSeq.size(); ++i)
      {
	unsigned int dti = 0;
	auto         &t = timeSeq[i];

	result.times = t;
	for (unsigned int tLast = 0; tLast < par.global.nt; ++tLast)
	{
	  result.correlator[tLast] = 0.;
	}
	for (auto &dt: translations)
	{
	  std::cout << "* Step " << i*translations.size() + dti + 1
		    << "/" << timeSeq.size()*translations.size()
		    << " -- positions= " << t << ", dt= " << dt << std::endl;
	  flops  = 0.;
	  bytes  = 0.;
	  fusec  = tAr.getDTimer("A*B algebra prec");
	  busec  = tAr.getDTimer("A*B total prec");
	  tAr.startTimer("Linear algebra");
	  unsigned int nlast = terms.size()-1;
	  tAr.startTimer("Disk vector overhead");
	  prod = a2aMat.at(terms[0].term)[TIME_MOD(t[0] + dt)];
	  tAr.stopTimer("Disk vector overhead");
	  for (unsigned int j = 1; j < terms.size(); ++j)
	  {
	    if ( terms[j].tdps ) break;
	    --nlast;
	    tAr.startTimer("Disk vector overhead");
	    const A2AMatrix<ComplexD> &ref = a2aMat.at(terms[j].term)[TIME_MOD(t[j] + dt)];
	    tAr.stopTimer("Disk vector overhead");

	    tAr.startTimer("A*B total prec");
	    tAr.startTimer("A*B algebra prec");
	    A2AContraction::mul(tmp, prod, ref);
	    tAr.stopTimer("A*B algebra prec");
	    flops += A2AContraction::mulFlops(prod, ref);
	    prod   = tmp;
	    tAr.stopTimer("A*B total prec");
	    bytes += 3.*tmp.rows()*tmp.cols()*sizeof(ComplexD);
	  }
	  if (terms.size() - nlast > 1 )
	  {
	    std::cout << std::setw(8) << "t-independent products"
		      << Sec(tAr.getDTimer("A*B total prec") - busec) << " "
		      << Flops(flops, tAr.getDTimer("A*B algebra prec") - fusec) << " "
		      << Bytes(bytes, tAr.getDTimer("A*B total prec") - busec) << std::endl;
	  }
	  std::cout << std::setw(8) << "traces";
	  flops  = 0.;
	  bytes  = 0.;
	  fusec  = tAr.getDTimer("tr(A*B)");
	  busec  = tAr.getDTimer("tr(A*B)");
	  flopsp  = 0.;
	  bytesp  = 0.;
	  fusecp  = tAr.getDTimer("A*B algebra rest");
	  busecp  = tAr.getDTimer("A*B total rest");
	  if ( nlast > 1 ) {
	    tAr.startTimer("Copying A2AMatrix");
	    prod0 = prod;
	    tAr.stopTimer("Copying A2AMatrix");
	    for (unsigned int tLast = p.min_t+dt; tLast <= p.max_t+dt; ++tLast)
	    {
	      tAr.startTimer("Copying A2AMatrix");
	      prod = prod0;
	      tAr.stopTimer("Copying A2AMatrix");
	      for (unsigned int j = terms.size() - nlast ; j < terms.size() - 1 ; ++j) {
		tAr.startTimer("Disk vector overhead");
		unsigned int tidx = TIME_MOD(t[j]+dt);
		if ( terms[j].tdps ) tidx = TIME_MOD(t[j] + tLast);
		const A2AMatrix<ComplexD> &ref = a2aMat.at(terms[j].term)[tidx];
		tAr.stopTimer("Disk vector overhead");

		tAr.startTimer("A*B total rest");
		tAr.startTimer("A*B algebra rest");
		A2AContraction::mul(tmp, prod, ref);
		tAr.stopTimer("A*B algebra rest");
		flopsp += A2AContraction::mulFlops(prod, ref);
		prod   = tmp;
		tAr.stopTimer("A*B total rest");
		bytesp += 3.*tmp.rows()*tmp.cols()*sizeof(ComplexD);
	      }
	      tAr.startTimer("tr(A*B)");
	      unsigned int tidx = TIME_MOD(t[terms.size()-1]+dt);
	      if ( terms[terms.size()-1].tdps )
		tidx = TIME_MOD(t[terms.size()-1] + tLast);
	      A2AContraction::accTrMul(result.correlator[TIME_MOD(tLast - dt)], prod, lastTerm[tidx]);
	      tAr.stopTimer("tr(A*B)");
	      flops += A2AContraction::accTrMulFlops(prod, lastTerm[tidx]);
	      bytes += 2.*prod.rows()*prod.cols()*sizeof(ComplexD);
	    }
	    std::cout << Sec(tAr.getDTimer("A*B total rest") - busecp) << " "
		      << Flops(flopsp, tAr.getDTimer("A*B algebra rest") - fusecp) << " "
		      << Bytes(bytesp, tAr.getDTimer("A*B total rest") - busecp) << std::endl;
	  } else {
	    for (unsigned int tLast = p.min_t+dt; tLast <= p.max_t+dt; ++tLast)
	    {
	      tAr.startTimer("tr(A*B)");
	      unsigned int tidx = TIME_MOD(t[terms.size()-1]+tLast);
	      assert( terms[terms.size()-1].tdps );
	      A2AContraction::accTrMul(result.correlator[TIME_MOD(tLast - dt)], prod, lastTerm[tidx]);
	      tAr.stopTimer("tr(A*B)");
	      flops += A2AContraction::accTrMulFlops(prod, lastTerm[tidx]);
	      bytes += 2.*prod.rows()*prod.cols()*sizeof(ComplexD);
	    }
	  }
	  tAr.stopTimer("Linear algebra");
	  std::cout << Sec(tAr.getDTimer("tr(A*B)") - busec) << " "
		    << Flops(flops, tAr.getDTimer("tr(A*B)") - fusec) << " " 
		    << Bytes(bytes, tAr.getDTimer("tr(A*B)") - busec) << std::endl;
	  if (!p.translationAverage)
	  {
	    dataSet = fileStem + "_dt_" + std::to_string(dt);
	    tAr.startTimer("Write");
	    write(writer, dataSet, result);
	    tAr.stopTimer("Write");
	    //saveCorrelator(result, par.global.output, dt, traj);
	    for (unsigned int tLast = 0; tLast < par.global.nt; ++tLast)
	    {
	      result.correlator[tLast] = 0.;
	    }
	  }
	  dti++;
	} // dt
	if (p.translationAverage)
	{
	  for (unsigned int tLast = 0; tLast < par.global.nt; ++tLast)
	  {
	    result.correlator[tLast] /= translations.size();
	  }
	  tAr.startTimer("Write");
	  write(writer, fileStem, result);
	  tAr.stopTimer("Write");
	  //saveCorrelator(result, par.global.output, 0, traj);
	}
      }
      tAr.stopTimer("Total");
      printTimeProfile(tAr.getTimings(), tAr.getTimer("Total"));
    }
  }
    
  return EXIT_SUCCESS;
}
