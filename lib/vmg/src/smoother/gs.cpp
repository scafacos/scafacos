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
 * @file   gs.cpp
 * @author Julian Iseringhausen <isering@ins.uni-bonn.de>
 * @date   Mon Apr 18 13:07:44 2011
 *
 * @brief  Gauss-Seidel method
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <cassert>

#include "base/discretization.hpp"
#include "base/stencil.hpp"
#include "comm/comm.hpp"
#include "grid/multigrid.hpp"
#include "smoother/gs.hpp"
#include "mg.hpp"

using namespace VMG;

void GaussSeidel::Compute(Grid& sol, Grid& rhs)
{
#ifdef DEBUG_MATRIX_CHECKS
  sol.IsConsistent();
  rhs.IsConsistent();
  sol.IsCompatible(rhs);
#endif

  Grid::iterator grid_iter;
  Stencil::iterator stencil_iter;
  vmg_float temp;

  const Stencil& A = MG::GetDiscretization()->GetStencil();
  const vmg_float prefactor_inv = 1.0 / MG::GetDiscretization()->OperatorPrefactor(sol);
  const vmg_float diag_inv = 1.0 / A.GetDiag();

  MG::GetComm()->CommToGhosts(sol);

  for (grid_iter = rhs.Iterators().Local().Begin(); grid_iter != rhs.Iterators().Local().End(); ++grid_iter) {

    temp = prefactor_inv * rhs.GetVal(*grid_iter);

    for (stencil_iter=A.begin(); stencil_iter!=A.end(); ++stencil_iter)
      temp -= stencil_iter->Val() * sol.GetVal(*grid_iter + stencil_iter->Disp());

    sol(*grid_iter) = temp * diag_inv;

  }
}
