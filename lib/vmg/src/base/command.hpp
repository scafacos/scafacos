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
 * @file   command.hpp
 * @author Julian Iseringhausen <isering@ins.uni-bonn.de>
 * @date   Tue Apr  5 19:21:37 2011
 *
 * @brief  Base class for commands that can be added to several
 *         command lists. Commands in the lists COMMANDLIST_INIT,
 *         COMMAND_LIST_LOOP and COMMANDLIST_FINALIZE will be
 *         executed automatically. However the user may create
 *         other command lists and execute them manually.
 *
 */

#ifndef COMMAND_HPP_
#define COMMAND_HPP_

#include <string>
#include <vector>

#include "base/defs.hpp"
#include "base/proxy.hpp"
#include "comm/comm.hpp"
#include "mg.hpp"

#ifdef HAVE_MPE
#include <mpe_log.h>
#define MPE_EVENT_BEGIN() static int eventID_begin;\
  static int eventID_end;\
  if (eventID_begin == 0 && eventID_end == 0) {\
    MPE_Log_get_state_eventIDs(&eventID_begin, &eventID_end);\
    if (MG::GetComm()->GlobalRank() == 0)\
      MPE_Describe_state(eventID_begin, eventID_end, Name(), "pink");	\
  }\
  MPE_Log_event(eventID_begin, 0, NULL);
#define MPE_EVENT_END() MPE_Log_event(eventID_end, 0, NULL);
#else
#define MPE_EVENT_BEGIN()
#define MPE_EVENT_END()
#endif /* HAVE_MPE */

#define CREATE_INITIALIZER(a) void Initialize##a() {		\
    VMG::MG::GetCommands().Register(new VMG::CommandProxy<a>);	\
  }

namespace VMG
{

class Command
{
public:
  typedef std::vector<std::string> argument_vector; ///< Holds arguments for execution

  virtual ~Command() {}

  /**
   * This function will be executed when command flow in a certain command
   * list reaches this command.
   *
   * @param arguments Holds arguments for execution as strings. When using multiple arguments,
   *                  ':' serves as the separator.
   *
   * @return Certain commands may want to change the following command flow of the command list
   *         (e.g. stopping criteria). This can be achieved by returning either "StopNow", which
   *         stops execution right after this command or by returning "StopLater", which will stop
   *         execution of a loop after the remaining commands have also been executed.
   */
  virtual Request Run(Command::argument_vector arguments) = 0;
};

}

#endif /* COMMAND_HPP_ */
