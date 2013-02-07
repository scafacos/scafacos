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
 * @file   cycle.hpp
 * @author Julian Iseringhausen <isering@ins.uni-bonn.de>
 * @date   Wed Jun 27 14:45:29 2012
 *
 * @brief  Base class for multigrid cycles.
 *
 */

#ifndef CYCLE_HPP_
#define CYCLE_HPP_

#include "base/object.hpp"

namespace VMG
{

class CommandList;

class Cycle : public Object
{
public:
  Cycle(const int& gamma, bool register_ = true) :
    Object("CYCLE", register_),
    gamma(gamma)
  {}

  virtual ~Cycle() {}

  virtual void Generate() = 0;

protected:
  void AddCycle(CommandList& list, int numLevels, int gamma);
  void AddCycleDebug(CommandList& list, int numLevels, int gamma);

  int gamma;
};

}

#endif /* CYCLE_HPP_ */
