/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: Hadrons/TimeWindowIo.hpp

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
#ifndef Hadrons_TimeWindowIo_hpp_
#define Hadrons_TimeWindowIo_hpp_
#include <Hadrons/Global.hpp>

BEGIN_HADRONS_NAMESPACE

// Shared by the time-window saver and its reader. Geometry and window LENGTH
// are automatic; the original times/source times are deliberately NOT stored.
// The binary field has dimensions (Nx, Ny, Nz, nTime), with time coordinate
// u indexing the original vector's slices in their unchanged input order.
class PropagatorTimeWindowMetadata: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(PropagatorTimeWindowMetadata,
                                    std::string,     format,
                                    unsigned int,    version,
                                    std::vector<int>, fullDimensions,
                                    unsigned int,    nTime,
                                    unsigned int,    trajectory,
                                    std::string,     sourceObject);
};

END_HADRONS_NAMESPACE
#endif // Hadrons_TimeWindowIo_hpp_
