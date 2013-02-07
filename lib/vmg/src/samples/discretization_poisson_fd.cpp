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
 * @file   discretization_poisson_fd_collatz.cpp
 * @author Julian Iseringhausen <isering@ins.uni-bonn.de>
 * @date   Mon Apr 18 13:03:47 2011
 *
 * @brief  Discretization of the poisson equation
 *         using the Collatz Mehrstellen Ansatz.
 *         Discretization error: O(h^4)
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "grid/grid.hpp"
#include "samples/discretization_poisson_fd.hpp"
#include "mg.hpp"

using namespace VMG;

DiscretizationPoissonFD::DiscretizationPoissonFD(const int& order) :
  Discretization(order)
{
  switch (order)
    {
    case 2:
      stencil.SetDiag(6.0);
      stencil.push_back(-1,  0,  0, -1.0);
      stencil.push_back( 1,  0,  0, -1.0);
      stencil.push_back( 0, -1,  0, -1.0);
      stencil.push_back( 0,  1,  0, -1.0);
      stencil.push_back( 0,  0, -1, -1.0);
      stencil.push_back( 0,  0,  1, -1.0);
      break;
    case 4:
      stencil.SetDiag(24.0/6.0);
      stencil.push_back(-1,  0,  0, -2.0/6.0);
      stencil.push_back( 1,  0,  0, -2.0/6.0);
      stencil.push_back( 0, -1,  0, -2.0/6.0);
      stencil.push_back( 0,  1,  0, -2.0/6.0);
      stencil.push_back( 0,  0, -1, -2.0/6.0);
      stencil.push_back( 0,  0,  1, -2.0/6.0);
      stencil.push_back(-1, -1,  0, -1.0/6.0);
      stencil.push_back(-1,  1,  0, -1.0/6.0);
      stencil.push_back( 1, -1,  0, -1.0/6.0);
      stencil.push_back( 1,  1,  0, -1.0/6.0);
      stencil.push_back(-1,  0, -1, -1.0/6.0);
      stencil.push_back(-1,  0,  1, -1.0/6.0);
      stencil.push_back( 1,  0, -1, -1.0/6.0);
      stencil.push_back( 1,  0,  1, -1.0/6.0);
      stencil.push_back( 0, -1, -1, -1.0/6.0);
      stencil.push_back( 0, -1,  1, -1.0/6.0);
      stencil.push_back( 0,  1, -1, -1.0/6.0);
      stencil.push_back( 0,  1,  1, -1.0/6.0);
      break;
    default:
      assert(0 != "vmg choose discretization order 2 or 4");
      break;
    }
}

void DiscretizationPoissonFD::ModifyRightHandSide()
{
  if (order == 4) {

    Grid& rhs = MG::GetRhsMaxLevel();

    Stencil stencil(6.0/12.0);
    stencil.push_back(-1,  0,  0, 1.0/12.0);
    stencil.push_back( 1,  0,  0, 1.0/12.0);
    stencil.push_back( 0, -1,  0, 1.0/12.0);
    stencil.push_back( 0,  1,  0, 1.0/12.0);
    stencil.push_back( 0,  0, -1, 1.0/12.0);
    stencil.push_back( 0,  0,  1, 1.0/12.0);

    stencil.Apply(rhs);

  }
}
