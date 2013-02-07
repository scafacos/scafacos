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
 * @file   cycle_cs_dirichlet.hpp
 * @author Julian Iseringhausen <isering@ins.uni-bonn.de>
 * @date   Wed Jun 27 14:48:42 2012
 *
 * @brief  Generates a Correction Scheme cycle for
 *         Dirichlet boundary conditions.
 *
 */

#ifndef CYCLE_CS_DIRICHLET_HPP_
#define CYCLE_CS_DIRICHLET_HPP_

#include "cycles/cycle.hpp"

namespace VMG
{

class CycleCSDirichlet : public Cycle
{
public:
  CycleCSDirichlet(const int& gamma, bool register_ = true) :
    Cycle(gamma, register_)
  {}

  virtual ~CycleCSDirichlet() {}

  void Generate();
};

}

#endif /* CYCLE_CS_DIRICHLET_HPP_ */
