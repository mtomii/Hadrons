/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/Modules/MIO/LoadBinary.hpp

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

/* This file is modified from Hadrons/Modules/MIO/LoadBinary.hpp in Grid */
#ifndef Hadrons_MIO_WriteProp_hpp_
#define Hadrons_MIO_WriteProp_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                       Load a binary configurations                         *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MIO)

class WritePropPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(WritePropPar,
                                    std::string, prop,
                                    std::string, file,
                                    std::string, format);
};

template <typename Impl>
class TWriteProp: public Module<WritePropPar>
{
public:
    typedef typename Impl::PropagatorField        Field;
    typedef typename Impl::Simd                   Simd;
    typedef typename Field::vector_object         vobj;
    typedef typename vobj::scalar_object          sobj;
    typedef typename sobj::DoublePrecision        sobj_double;
    typedef BinarySimpleMunger<sobj, sobj_double> Munger;
public:
    // constructor
    TWriteProp(const std::string name);
    // destructor
    virtual ~TWriteProp(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(WriteProp, TWriteProp<FIMPL>, MIO);

/******************************************************************************
 *                         TWriteProp implementation                         *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename Impl>
TWriteProp<Impl>::TWriteProp(const std::string name)
: Module<WritePropPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename Impl>
std::vector<std::string> TWriteProp<Impl>::getInput(void)
{
    std::vector<std::string> in = {par().prop};
    
    return in;
}

template <typename Impl>
std::vector<std::string> TWriteProp<Impl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename Impl>
void TWriteProp<Impl>::setup(void)
{ }

// execution ///////////////////////////////////////////////////////////////////
template <typename Impl>
void TWriteProp<Impl>::execute(void)
{
    Munger      munge;
    uint32_t    nersc_csum, scidac_csuma, scidac_csumb;
    auto        &U = envGet(Field, par().prop);
    std::string filename = par().file + "."
                           + std::to_string(vm().getTrajectory());

    LOG(Message) << "Writing " << par().format 
                 << " binary configuration from file '" << filename
                 << "'" << std::endl;
    BinaryIO::writeLatticeObject<vobj, sobj_double>(U, filename, munge, 0, 
                                                   par().format, nersc_csum,
                                                   scidac_csuma, scidac_csumb);
    LOG(Message) << "Checksums:" << std::endl;
    LOG(Message) << "  NERSC    " << nersc_csum << std::endl;
    LOG(Message) << "  SciDAC A " << scidac_csuma << std::endl;
    LOG(Message) << "  SciDAC B " << scidac_csumb << std::endl;
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MIO_WriteProp_hpp_
