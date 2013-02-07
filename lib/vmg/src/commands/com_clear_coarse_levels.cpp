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
 * @file   com_clear_coarse_levels.cpp
 * @author Julian Iseringhausen <isering@ins.uni-bonn.de>
 * @date   Mon Apr 18 12:32:12 2011
 *
 * @brief  Overwrites all grid values of all grids except
 *         the finest one with zeros.
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "base/command.hpp"
#include "base/factory.hpp"
#include "grid/multigrid.hpp"

using namespace VMG;

class VMGCommandClearCoarseLevels : public Command
{
public:
  Request Run(Command::argument_vector arguments)
  {
    MPE_EVENT_BEGIN()

    MG::GetFactory().Get(arguments[0])->Cast<Multigrid>()->ClearAllCoarseLevels();

    MPE_EVENT_END()

    return Continue;
  }

  static const char* Name() {return "ClearCoarseLevels";}
  static int Arguments() {return 1;}
};

CREATE_INITIALIZER(VMGCommandClearCoarseLevels)

