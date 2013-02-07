/*
 *    vmg - a versatile multigrid solver
 *    Copyright (C) 2012 Institute for Numerical Simulation, University of Bonn
 *
 *  vmg is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  vmg is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/*
 * interface_sinus.cpp
 *
 *  Created on: 31.03.2011
 *      Author: Julian Iseringhausen
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <cmath>
#include <iostream>

#include "grid/grid.hpp"
#include "grid/multigrid.hpp"

#include "interface_sinus.hpp"

using namespace VMG;
using VMGInterfaces::InterfaceSinus;

void InterfaceSinus::ImportRightHandSide(Multigrid& multigrid)
{
  Index i;
  Vector pos;

  Grid& grid = multigrid(multigrid.MaxLevel());
  grid.ClearBoundary();

  const Index begin_local = grid.Global().LocalBegin() - grid.Local().HaloSize1();

  for (i.X()=grid.Local().Begin().X(); i.X()<grid.Local().End().X(); ++i.X())
    for (i.Y()=grid.Local().Begin().Y(); i.Y()<grid.Local().End().Y(); ++i.Y())
      for (i.Z()=grid.Local().Begin().Z(); i.Z()<grid.Local().End().Z(); ++i.Z()) {
        pos = grid.Extent().MeshWidth() * static_cast<Vector>(begin_local + i);
	grid(i) = 3.0 * sine_factor * sine_factor * std::sin(sine_factor * pos.X()) * std::sin(sine_factor * pos.Y()) * std::sin(sine_factor * pos.Z());
      }
}

void InterfaceSinus::ExportSolution(Grid& grid)
{
}
