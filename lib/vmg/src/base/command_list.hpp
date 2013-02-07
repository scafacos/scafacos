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
 * @file   command_list.hpp
 * @author Julian Iseringhausen <isering@ins.uni-bonn.de>
 * @date   Tue Apr  5 20:28:02 2011
 *
 * @brief  This class holds a list of commands. These commands
 *         can be executed by calling ExecuteList.
 *
 */

#ifndef COMMAND_LIST_HPP_
#define COMMAND_LIST_HPP_

#include <list>
#include <string>
#include <vector>

#include "base/defs.hpp"
#include "base/object.hpp"

namespace VMG
{

class Command;

class CommandList : public Object
{
public:
  typedef std::list< std::pair<std::string, std::vector<std::string> > >::iterator iterator;

  Request ExecuteList(); ///< Execute all commands in this list.

  void AddCommand(std::string command, std::string arguments = ""); ///< Add a command to the back of the list.
  void DeleteCommand(const CommandList::iterator& iter);          ///< Remove a command from the list.

  void Print(); ///< Print all commands in list.

  void Clear(); ///< Remove all commands from list.

private:
  std::list< std::pair<std::string, std::vector< std::string> > > commands;
};

}

#endif /* COMMAND_LIST_HPP_ */
