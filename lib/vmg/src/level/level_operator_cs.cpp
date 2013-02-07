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
 * @file   level_operator_cs.cpp
 * @author Julian Iseringhausen <isering@ins.uni-bonn.de>
 * @date   Mon Apr 18 12:59:46 2011
 *
 * @brief  VMG::LevelOperatorCS
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "base/index.hpp"
#include "comm/comm.hpp"
#include "grid/grid_iterator.hpp"
#include "grid/multigrid.hpp"
#include "grid/tempgrid.hpp"
#include "level/level_operator_cs.hpp"
#include "mg.hpp"

using namespace VMG;

void LevelOperatorCS::Restrict(Multigrid& sol, Multigrid& rhs)
{
  Grid::iterator iter_f, iter_c;

  Comm& comm = *MG::GetComm();

  Grid& sol_f = sol(sol.Level());
  Grid& rhs_f = rhs(rhs.Level());
  Grid& rhs_c_dist = rhs(rhs.Level()-1);
  Grid& rhs_c_undist = comm.GetCoarserGrid(rhs);

  const Stencil& op = OperatorRestrict();

  const Index begin_f_global = rhs_f.Global().LocalSize().Product() > 0
    ? rhs_f.Global().LocalBegin() + rhs_f.Global().LocalBegin() % 2
    : rhs_f.Global().LocalBegin();
  const Index end_f_global = rhs_f.Global().LocalSize().Product() > 0
    ? rhs_f.Global().LocalEnd() - (rhs_f.Global().LocalEnd() - 1) % 2
    : rhs_f.Global().LocalBegin();

  Index begin_f = begin_f_global - rhs_f.Global().LocalBegin() + rhs_f.Local().Begin();
  Index end_f = end_f_global - rhs_f.Global().LocalBegin() + rhs_f.Local().Begin();

  /* Modify fine begin/end to align the points on both levels correctly */
  if (rhs_c_undist.Global().BoundaryType() == GlobalCoarsened) {
    begin_f += rhs_f.Local().BoundarySize1();
    end_f -= rhs_f.Local().BoundarySize2();
  }

  const Index begin_c = rhs_c_undist.Local().Begin() + begin_f_global / 2 - rhs_c_undist.Global().LocalBegin();
  const Index end_c = rhs_c_undist.Local().Begin() + end_f_global / 2 - rhs_c_undist.Global().LocalBegin() + 1;

  const GridIteratorSet bounds_f(begin_f, end_f);
  const GridIteratorSet bounds_c(begin_c, end_c);

  // Compute residual
  TempGrid *temp = MG::GetTempGrid();
  temp->SetProperties(sol_f);
  temp->ImportFromResidual(sol_f, rhs_f);
  comm.CommToGhosts(*temp);

  for (iter_f=bounds_f.Begin(), iter_c=bounds_c.Begin(); iter_f!=bounds_f.End(); iter_f+=2, ++iter_c)
    rhs_c_undist(*iter_c) = op.Apply(*temp, *iter_f);

  comm.CommSubgrid(rhs_c_undist, rhs_c_dist, 1);

  if (rhs_c_dist.Global().BoundaryType() == GlobalCoarsened)
    rhs_c_dist.ClearBoundary();

  sol.ToCoarserLevel();
  rhs.ToCoarserLevel();
}

void LevelOperatorCS::Prolongate(Multigrid& sol, Multigrid& rhs)
{
  Grid::iterator iter_f, iter_c;
  Stencil::iterator stencil_iter;
  vmg_float val;

  Comm& comm = *MG::GetComm();

  Grid& sol_c = sol(sol.Level());
  Grid& sol_f_dist = sol(sol.Level()+1);
  Grid& sol_f_undist = comm.GetFinerGrid(sol);
  Grid& rhs_f_undist = comm.GetFinerGrid(rhs);

  const Stencil& op = OperatorProlongate();

  Index begin_f = sol_f_undist.Local().Begin() + 2*sol_c.Global().LocalBegin() - sol_f_undist.Global().LocalBegin();
  Index end_f = sol_f_undist.Local().End();

  if (sol_c.Global().BoundaryType() == GlobalCoarsened) {
    begin_f += rhs_f_undist.Local().BoundarySize1();
    end_f -= rhs_f_undist.Local().BoundarySize2();
  }

  const GridIteratorSet bounds_f(begin_f, end_f);
  const GridIteratorSet bounds_c(sol_c.Local().FinerBegin(), sol_c.Local().FinerEnd());

  sol_f_undist.ClearHalo();

  comm.CommSubgrid(sol_f_dist, sol_f_undist, 0);

  for (iter_f=bounds_f.Begin(), iter_c=bounds_c.Begin(); iter_c!=bounds_c.End(); iter_f+=2, ++iter_c) {
    val = sol_c.GetVal(*iter_c);
    sol_f_undist(*iter_f) += op.GetDiag() * val;
    for (stencil_iter = op.begin(); stencil_iter != op.end(); ++stencil_iter)
      sol_f_undist(iter_f->X() + stencil_iter->Disp().X(),
		   iter_f->Y() + stencil_iter->Disp().Y(),
		   iter_f->Z() + stencil_iter->Disp().Z()) += stencil_iter->Val() * val;
  }

  comm.CommFromGhosts(sol_f_undist);
  comm.CommSubgrid(sol_f_undist, sol_f_dist, 1);

  sol.ToFinerLevel();
  rhs.ToFinerLevel();
}
