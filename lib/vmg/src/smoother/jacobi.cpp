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
 * @file   jacobi.hpp
 * @author Julian Iseringhausen <isering@ins.uni-bonn.de>
 * @date   Fri Apr 27 19:50:36 2012
 *
 * @brief  Jacobi smoother
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#include "base/discretization.hpp"
#include "base/stencil.hpp"
#include "comm/comm.hpp"
#include "grid/grid.hpp"
#include "grid/tempgrid.hpp"
#include "smoother/jacobi.hpp"
#include "mg.hpp"

using namespace VMG;

static inline void ComputePartial(Grid& sol_new, const Grid& sol_old, const Grid& rhs,
				  const Stencil& mat, const GridIteratorSet& bounds,
				  const vmg_float& prefactor, const vmg_float& diag_inv)
{
  Grid::iterator grid_iter;
  Stencil::iterator stencil_iter;

  for (grid_iter=bounds.Begin(); grid_iter!=bounds.End(); ++grid_iter) {
    sol_new(*grid_iter) = prefactor * rhs.GetVal(*grid_iter);
    for (stencil_iter=mat.begin(); stencil_iter!=mat.end(); ++stencil_iter)
      sol_new(*grid_iter) -= stencil_iter->Val() * sol_old.GetVal(*grid_iter + stencil_iter->Disp());
    sol_new(*grid_iter) *= diag_inv;
  }
}

void Jacobi::Compute(Grid& sol, Grid& rhs)
{
  Comm& comm = *MG::GetComm();
  TempGrid& sol_new = *MG::GetTempGrid();

  const Stencil& mat = MG::GetDiscretization()->GetStencil();
  const vmg_float prefactor_inv = 1.0 / MG::GetDiscretization()->OperatorPrefactor(sol);
  const vmg_float diag_inv = 1.0 / mat.GetDiag();

  sol_new.SetProperties(sol);

  comm.CommToGhostsAsyncStart(sol);
  ComputePartial(sol_new, sol, rhs, mat, sol.Iterators().InnerLocalGrid(), prefactor_inv, diag_inv);
  comm.CommToGhostsAsyncFinish(sol);

  for (int i=0; i<3; ++i) {
    ComputePartial(sol_new, sol, rhs, mat, sol.Iterators().NearBoundary1()[i], prefactor_inv, diag_inv);
    ComputePartial(sol_new, sol, rhs, mat, sol.Iterators().NearBoundary2()[i], prefactor_inv, diag_inv);
  }

  sol.SetGrid(sol_new);
}
