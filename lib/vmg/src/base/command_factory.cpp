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
 * @file   command_factory.cpp
 * @author Julian Iseringhausen <isering@ins.uni-bonn.de>
 * @date   Thu May 19 16:48:40 2011
 *
 * @brief  Class to store possible commands.
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <cstdio>

#include "base/command.hpp"
#include "base/command_factory.hpp"
#include "base/defs.hpp"
#include "base/proxy.hpp"

using namespace VMG;

CommandFactory::CommandFactory()
{
}

CommandFactory::~CommandFactory()
{
  for (std::map<std::string, CommandProxyBase*>::iterator iter=proxy_map.begin(); iter!=proxy_map.end(); ++iter)
    delete iter->second;

  for (std::map<std::string, Command*>::iterator iter=instance_map.begin(); iter!=instance_map.end(); ++iter)
    delete iter->second;
}

void CommandFactory::Register(CommandProxyBase* object)
{
  proxy_map.insert(std::make_pair(object->Name(), object));
}

Command* CommandFactory::Get(const std::string& id)
{
  std::map<std::string, Command*>::iterator iter_instance = instance_map.find(id);
  if (iter_instance != instance_map.end())
    return iter_instance->second;

  std::map<std::string, CommandProxyBase*>::iterator iter_proxy = proxy_map.find(id);
  if (iter_proxy != proxy_map.end()) {
    Command *command = iter_proxy->second->CreateCommand();
    instance_map.insert(std::make_pair(id, command));
    return command;
  }

  assert(0 == "This command is not registered.");

  return NULL;
}

void CommandFactory::Delete(const std::string& id)
{
  std::map<std::string, CommandProxyBase*>::iterator iter_proxy = proxy_map.find(id);

  if (iter_proxy != proxy_map.end()) {
    delete iter_proxy->second;
    proxy_map.erase(iter_proxy);
  }

  std::map<std::string, Command*>::iterator iter_instance = instance_map.find(id);

  if (iter_instance != instance_map.end()) {
    delete iter_instance->second;
    instance_map.erase(iter_instance);
  }
}


bool CommandFactory::CheckNumberOfArguments(const std::string& id, const int& num_arguments)
{
  std::map<std::string, CommandProxyBase*>::iterator iter = proxy_map.find(id);

  assert(iter != proxy_map.end());
  assert(iter->second->Arguments() == num_arguments);

  return iter->second->Arguments() == num_arguments;
}


void CommandFactory::PrintAvailableCommands()
{
  for (std::map<std::string, CommandProxyBase*>::iterator iter=proxy_map.begin(); iter!=proxy_map.end(); ++iter)
    printf("%s\n", iter->second->Name());
}
