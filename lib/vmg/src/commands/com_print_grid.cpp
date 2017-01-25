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

/**
 * @file   com_print_grid.cpp
 * @author Julian Iseringhausen <isering@ins.uni-bonn.de>
 * @date   Mon Apr 18 12:43:54 2011
 *
 * @brief  Writes a given grid including boundary and
 *         halo to a VTKStructuredPoints file in
 *         output/date_time
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <sstream>

#include "base/command.hpp"
#include "comm/comm.hpp"
#include "grid/grid.hpp"
#include "grid/multigrid.hpp"

using namespace VMG;

class VMGCommandPrintGrid : public Command
{
public:
  Request Run(Command::argument_vector arguments)
  {
    MPE_EVENT_BEGIN()

    std::ostringstream buffer;

    Multigrid& grid = *MG::GetFactory().Get(arguments[0])->Cast<Multigrid>();

    buffer << "Level " << grid.Level() << " ";

    switch (grid().Global().BoundaryType())
      {
      case LocallyRefined:
	buffer << "Locally Refined ";
	break;
      case GlobalMax:
	buffer << "Global Max ";
	break;
      case GlobalCoarsened:
	buffer << "Global Coarsened ";
	break;
      case EmptyGrid:
	buffer << "Empty Grid ";
	break;
      case BTUndefined:
        buffer << "Boundary Type Undefined";
        break;
      }

    buffer << arguments[0];

    MG::GetComm()->PrintGrid(grid(), buffer.str().c_str());

    MPE_EVENT_END()

    return Continue;
  }

  static const char* Name() {return "PrintGrid";}
  static int Arguments() {return 1;}
};

CREATE_INITIALIZER(VMGCommandPrintGrid)
