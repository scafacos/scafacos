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
 * @file   com_check_iteration_counter.cpp
 * @author Julian Iseringhausen <isering@ins.uni-bonn.de>
 * @date   Mon Apr 18 12:29:39 2011
 *
 * @brief  Increases the iteration counter and checks if
 *         the maximum number of iterations is reached.
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "base/command.hpp"
#include "base/object.hpp"
#include "mg.hpp"

using namespace VMG;

class VMGCommandCheckIterationCounter : public Command
{
public:
  VMGCommandCheckIterationCounter()
  {
    new ObjectStorage<int>("ITERATION", 0);
  }

  Request Run(Command::argument_vector arguments)
  {
    MPE_EVENT_BEGIN()

    const int& max_iteration = MG::GetFactory().GetObjectStorageVal<int>("MAX_ITERATION");
    int& iteration = MG::GetFactory().GetObjectStorageVal<int>("ITERATION");

    MPE_EVENT_END()

    if (++iteration >= max_iteration)
      return StopCycleLater;
    else
      return Continue;
  }

  static const char* Name() {return "CheckIterationCounter";}
  static int Arguments() {return 0;}
};

CREATE_INITIALIZER(VMGCommandCheckIterationCounter)
