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
 * @file   proxy.hpp
 * @author Julian Iseringhausen <isering@ins.uni-bonn.de>
 * @date   Mon Apr 18 12:23:43 2011
 *
 * @brief  Header file for the classes VMG::CommandProxyBase and
 *         VMG::CommandProxy, used register commands at the factory.
 *
 */

#ifndef COMMAND_PROXY_HPP_
#define COMMAND_PROXY_HPP_

namespace VMG
{

class Command;

class CommandProxyBase
{
public:
  CommandProxyBase();
  virtual Command* CreateCommand() const = 0;

  virtual const char* Name() const = 0;
  virtual int Arguments() const = 0;
};

template <class T>
class CommandProxy : public CommandProxyBase
{
public:
  Command* CreateCommand() const
  {
    return new T;
  }

  const char* Name() const
  {
    return T::Name();
  }

  int Arguments() const
  {
    return T::Arguments();
  }
};

}

#endif /* COMMAND_PROXY_HPP_ */
