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
 * @file   smoother.cpp
 * @author Julian Iseringhausen <isering@ins.uni-bonn.de>
 * @date   Mon Apr 18 13:08:57 2011
 *
 * @brief  VMG::Smoother
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "base/discretization.hpp"
#include "comm/comm.hpp"
#include "smoother/smoother.hpp"
#include "mg.hpp"

using namespace VMG;

void Smoother::Run(Multigrid& sol, Multigrid& rhs, int steps)
{
  for (int i=0; i<steps; i++) {

    if (sol().Global().BoundaryType() == LocallyRefined)
      MG::GetDiscretization()->SetInnerBoundary(sol(), rhs(), sol(sol.Level()-1));

    this->Compute(sol(), rhs());

  }
}
