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
 * @file   com_copy_boundary.cpp
 * @author Julian Iseringhausen <isering@ins.uni-bonn.de>
 * @date   Mon Apr 18 12:33:14 2011
 *
 * @brief  Copies the boundary values from one grid to another
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "base/command.hpp"
#include "grid/grid.hpp"
#include "grid/multigrid.hpp"

using namespace VMG;

class VMGCommandCopyBoundary : public Command
{
public:
  Request Run(Command::argument_vector arguments)
  {
    MPE_EVENT_BEGIN()

    Grid& from = (*MG::GetFactory().Get(arguments[0])->Cast<Multigrid>())();
    Grid& to = (*MG::GetFactory().Get(arguments[1])->Cast<Multigrid>())();

    to.SetBoundary(from);

    MPE_EVENT_END()

    return Continue;
  }

  static const char* Name() {return "CopyBoundary";}
  static int Arguments() {return 2;}
};

CREATE_INITIALIZER(VMGCommandCopyBoundary)

