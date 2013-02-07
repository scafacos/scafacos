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
 * @file   command_factory.hpp
 * @author Julian Iseringhausen <isering@ins.uni-bonn.de>
 * @date   Thu May 19 16:48:40 2011
 *
 * @brief  Class to store possible commands.
 *
 */

#ifndef COMMAND_FACTORY_HPP_
#define COMMAND_FACTORY_HPP_

#include <map>
#include <string>

namespace VMG
{

class Command;
class CommandProxyBase;

class CommandFactory
{
public:
  CommandFactory();
  virtual ~CommandFactory();

  void Register(CommandProxyBase* command); ///< Registers a command
  Command* Get(const std::string& id);             ///< Returns a command
  void Delete(const std::string& id);              ///< Deletes a command

  bool CheckNumberOfArguments(const std::string& id, const int& num_arguments);

  void PrintAvailableCommands(); ///< Prints the names of all registered commands

private:
  std::map<std::string, CommandProxyBase*> proxy_map;
  std::map<std::string, Command*> instance_map;

};

}

#endif /* COMMAND_FACTORY_HPP_ */

