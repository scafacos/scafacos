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
 * @file   cycle_fas_dirichlet_debug.hpp
 * @author Julian Iseringhausen <isering@ins.uni-bonn.de>
 * @date   Wed Jun 27 14:48:42 2012
 *
 * @brief  Generates a FAS cycle for
 *         Dirichlet boundary conditions
 *         with debug output.
 *
 */

#ifndef CYCLE_FAS_DIRICHLET_DEBUG_HPP_
#define CYCLE_FAS_DIRICHLET_DEBUG_HPP_

#include "cycles/cycle.hpp"

namespace VMG
{

class CycleFASDirichletDebug : public Cycle
{
public:
  CycleFASDirichletDebug(const int& gamma) :
    Cycle(gamma)
  {}

  virtual ~CycleFASDirichletDebug() {}

  void Generate();
};

}

#endif /* CYCLE_FAS_DIRICHLET_DEBUG_HPP_ */
