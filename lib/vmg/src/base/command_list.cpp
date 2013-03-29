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
 * @file   command_list.cpp
 * @author Julian Iseringhausen <isering@ins.uni-bonn.de>
 * @date   Tue Apr  5 20:15:06 2011
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef DEBUG_BARRIER
#ifdef HAVE_MPI
#include <mpi.h>
#ifdef HAVE_MARMOT
#include <enhancempicalls.h>
#include <sourceinfompicalls.h>
#endif
#endif
#endif

#include <cstdio>

#include "base/command.hpp"
#include "base/command_factory.hpp"
#include "base/command_list.hpp"
#include "base/timer.hpp"
#include "base/defs.hpp"
#include "mg.hpp"

using namespace VMG;

Request CommandList::ExecuteList()
{
  Request request;
  Request final_request = (commands.size() == 0 ? StopCycleNow : Continue);

  for (CommandList::iterator iter=commands.begin(); iter!=commands.end(); ++iter) {

#ifdef DEBUG
    const int num_args = (iter->second.size() > 1 ? iter->second.size() : (iter->second[0] == "" ? 0 : 1));
    MG::GetCommands().CheckNumberOfArguments(iter->first, num_args);
#endif

#ifdef DEBUG_BARRIER
    Comm& comm = *MG::GetComm();
    comm.Barrier();
    comm.PrintOnce(Debug, "Command \"%s\" start", iter->first.c_str());
#endif

    Timer::Start(iter->first);
    request = MG::GetCommands().Get(iter->first)->Run(iter->second);
    Timer::Stop(iter->first);

#ifdef DEBUG_BARRIER
    comm.Barrier();
    comm.PrintOnce(Debug, "Command \"%s\" done", iter->first.c_str());
#endif

    if (request == StopCycleLater)
      final_request = StopCycleNow;
    else if (request == StopCycleNow) {
      final_request = StopCycleNow;
      break;
    }
  }

  return final_request;
}

void CommandList::AddCommand(std::string command, std::string arguments)
{
  std::vector<std::string> argument_list;
  size_t pos;

  do {
    pos = arguments.find(':');
    argument_list.push_back(arguments.substr(0, pos));
    arguments.erase(0, pos+1);
  }while (pos != std::string::npos);

  commands.push_back(std::pair<std::string, std::vector<std::string> >(command, argument_list));
}

void CommandList::DeleteCommand(const CommandList::iterator& iter)
{
  commands.erase(iter);
}

void CommandList::Print()
{
  for (CommandList::iterator iter=commands.begin(); iter!=commands.end(); ++iter)
    printf("%s\n", (*iter).first.c_str());
}

void CommandList::Clear()
{
  commands.clear();
}
