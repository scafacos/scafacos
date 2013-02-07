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
 * @file   level_operator_fas.cpp
 * @author Julian Iseringhausen <isering@ins.uni-bonn.de>
 * @date   Mon Apr 18 13:00:23 2011
 *
 * @brief  VMG::LevelOperatorFAS
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "base/discretization.hpp"
#include "base/index.hpp"
#include "comm/comm.hpp"
#include "grid/multigrid.hpp"
#include "grid/tempgrid.hpp"
#include "level/level_operator_fas.hpp"
#include "mg.hpp"

using namespace VMG;

void LevelOperatorFAS::Restrict(Multigrid& sol, Multigrid& rhs)
{
  Grid::iterator iter_f, iter_c;
  Stencil::iterator stencil_iter;

  Comm& comm = *MG::GetComm();

  Grid& sol_f = sol(sol.Level());
  Grid& rhs_f = rhs(rhs.Level());

  Grid& sol_c_dist = sol(sol.Level()-1);
  Grid& rhs_c_dist = rhs(rhs.Level()-1);

  Grid& sol_c_undist = comm.GetCoarserGrid(sol);
  Grid& rhs_c_undist = comm.GetCoarserGrid(rhs);

  Index begin_f = rhs_f.Local().Begin();
  Index end_f = rhs_f.Local().End();

  if (rhs_c_undist.Global().BoundaryType() == GlobalCoarsened) {
    begin_f += rhs_f.Local().BoundarySize1();
    end_f -= rhs_f.Local().BoundarySize2();
  }

  const GridIteratorSet bounds_f(begin_f, end_f);
  const GridIteratorSet bounds_c(rhs_c_undist.Local().FinerBegin(), rhs_c_undist.Local().FinerEnd());

  const Stencil& op_res = OperatorRestrict();
  const Stencil& op_sol = OperatorSol();
  const Stencil& op_pde = MG::GetDiscretization()->GetStencil();
  const vmg_float prefactor = MG::GetDiscretization()->OperatorPrefactor(sol_c_undist);

  // Communication step
  comm.CommToGhosts(sol_f);
  comm.CommSubgrid(sol_c_dist, sol_c_undist, 0);

  MG::GetDiscretization()->SetInnerBoundary(sol_f, rhs_f, sol_c_undist);

  // Compute residual
  VMG::TempGrid *res_grid = MG::GetTempGrid();
  res_grid->SetProperties(rhs_f);
  res_grid->ImportFromResidual(sol_f, rhs_f);

  for (iter_f=bounds_f.Begin(), iter_c=bounds_c.Begin(); iter_c!=bounds_c.End(); iter_f+=2, ++iter_c)
    sol_c_undist(*iter_c) = op_sol.Apply(sol_f, *iter_f);

  for (iter_f=bounds_f.Begin(), iter_c=bounds_c.Begin(); iter_c!=bounds_c.End(); iter_f+=2, ++iter_c)
    rhs_c_undist(*iter_c) = op_res.Apply(*res_grid, *iter_f);

  comm.CommSubgrid(sol_c_undist, sol_c_dist, 1);
  comm.CommToGhosts(sol_c_dist);

  comm.CommSubgrid(rhs_c_undist, rhs_c_dist, 1);

  VMG::TempGrid* sol_old = this->GetTempGrid(sol_c_dist.Level());
  sol_old->SetProperties(sol_c_dist);
  sol_old->SetGrid(sol_c_dist);

  for (iter_c=rhs_c_dist.Iterators().Local().Begin(); iter_c!=rhs_c_dist.Iterators().Local().End(); ++iter_c)
    rhs_c_dist(*iter_c) += prefactor * op_pde.Apply(sol_c_dist, *iter_c);

  comm.CommToGhosts(rhs_c_dist);

  sol.ToCoarserLevel();
  rhs.ToCoarserLevel();
}

void LevelOperatorFAS::Prolongate(Multigrid& sol, Multigrid& rhs)
{
  Grid::iterator iter_f, iter_c;
  Stencil::iterator stencil_iter;
  vmg_float val;

  Comm& comm = *MG::GetComm();

  Grid& sol_c = sol(sol.Level());
  Grid& sol_f_dist = sol(sol.Level()+1);
  Grid& rhs_f_dist = rhs(rhs.Level()+1);
  Grid& sol_f_undist = comm.GetFinerGrid(sol);
  Grid& rhs_f_undist = comm.GetFinerGrid(rhs);

  const Stencil& op = OperatorProlongate();

  TempGrid *sol_old = this->GetTempGrid(sol_c.Level());

  Index begin_f = sol_f_undist.Local().Begin();
  Index end_f = sol_f_undist.Local().End();

  if (sol_c.Global().BoundaryType() == GlobalCoarsened) {
    begin_f += rhs_f_undist.Local().BoundarySize1();
    end_f -= rhs_f_undist.Local().BoundarySize2();
  }

  const GridIteratorSet bounds_f(begin_f, end_f);
  const GridIteratorSet bounds_c(sol_c.Local().FinerBegin(), sol_c.Local().FinerEnd());

  sol_f_undist.ClearHalo();

  for (iter_f=bounds_f.Begin(), iter_c=bounds_c.Begin(); iter_f!=bounds_f.End(); iter_f+=2, ++iter_c) {
    val = sol_c.GetVal(*iter_c) - sol_old->GetVal(*iter_c);
    sol_f_undist(*iter_f) += op.GetDiag() * val;
    for (stencil_iter = op.begin(); stencil_iter != op.end(); ++stencil_iter)
	sol_f_undist(*iter_f + stencil_iter->Disp()) += stencil_iter->Val() * val;
  }

  comm.CommFromGhosts(sol_f_undist);
  comm.CommSubgrid(sol_f_undist, sol_f_dist, 1);

  if (sol_f_dist.Global().BoundaryType() == LocallyRefined)
    MG::GetDiscretization()->SetInnerBoundary(sol_f_dist, rhs_f_dist, sol_c);

  sol.ToFinerLevel();
  rhs.ToFinerLevel();
}
