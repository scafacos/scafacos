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
 * @file   com_print_running_time.cpp
 * @author Julian Iseringhausen <isering@ins.uni-bonn.de>
 * @date   Mon Nov 21 11:01:58 2011
 *
 * @brief  Print the running time of various subsystems of vmg.
 *         Needs vmg to be configured with --enable-debug-output
 *         and --enable-debug-measure-time
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <cstdio>

#include "base/command.hpp"
#include "base/timer.hpp"

using namespace VMG;

class VMGCommandPrintRunningTime : public Command
{
public:
  Request Run(Command::argument_vector arguments)
  {
    MPE_EVENT_BEGIN()

      Timer::Print();

    MPE_EVENT_END()

    return Continue;
  }

  static const char* Name() {return "PrintRunningTime";}
  static int Arguments() {return 0;}
};

CREATE_INITIALIZER(VMGCommandPrintRunningTime)
