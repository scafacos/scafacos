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
 * @file   com_execute_full_cycle_loop.cpp
 * @author Julian Iseringhausen <isering@ins.uni-bonn.de>
 * @date   Mon Apr 18 12:35:37 2011
 *
 * @brief  Executes a set of cycles. foo_INIT and foo_FINALIZE
 *         get executed once and foo_LOOP in a loop.
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "base/command.hpp"
#include "base/command_list.hpp"

using namespace VMG;

class VMGCommandExecuteFullCycleLoop : public Command
{
public:
  Request Run(Command::argument_vector arguments)
  {
    MPE_EVENT_BEGIN()

    std::string str_init = arguments[0] + "_INIT";
    std::string str_loop = arguments[0] + "_LOOP";
    std::string str_finalize = arguments[0] + "_FINALIZE";

    MG::GetFactory().Get(str_init)->Cast<CommandList>()->ExecuteList();

    while (MG::GetFactory().Get(str_loop)->Cast<CommandList>()->ExecuteList() == Continue);

    MG::GetFactory().Get(str_finalize)->Cast<CommandList>()->ExecuteList();

    MPE_EVENT_END()

    return Continue;
  }

  static const char* Name() {return "ExecuteFullCycleLoop";}
  static int Arguments() {return 1;}
};

CREATE_INITIALIZER(VMGCommandExecuteFullCycleLoop)
