/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/Modules/MIO/SaveBinnedA2AVecs.cc

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
#include <Hadrons/Modules/MIO/SaveBinnedA2AVecs.hpp>

using namespace Grid;
using namespace Hadrons;
using namespace MIO;

template class Grid::Hadrons::MIO::TSaveBinnedA2AVecs<FIMPL, 48>;
template class Grid::Hadrons::MIO::TSaveBinnedA2AVecs<FIMPL, 64>;
template class Grid::Hadrons::MIO::TSaveBinnedA2AVecs<FIMPL, 96>;
template class Grid::Hadrons::MIO::TSaveBinnedA2AVecs<FIMPL, 144>;
template class Grid::Hadrons::MIO::TSaveBinnedA2AVecs<FIMPL, 149>;
template class Grid::Hadrons::MIO::TSaveBinnedA2AVecs<FIMPL, 173>;
template class Grid::Hadrons::MIO::TSaveBinnedA2AVecs<FIMPL, 192>;
template class Grid::Hadrons::MIO::TSaveBinnedA2AVecs<FIMPL, 196>;
template class Grid::Hadrons::MIO::TSaveBinnedA2AVecs<FIMPL, 197>;
template class Grid::Hadrons::MIO::TSaveBinnedA2AVecs<FIMPL, 216>;
template class Grid::Hadrons::MIO::TSaveBinnedA2AVecs<FIMPL, 298>;
