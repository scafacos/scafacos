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
 * @file   stencil.cpp
 * @author Julian Iseringhausen <isering@ins.uni-bonn.de>
 * @date   Mon Apr 18 12:25:05 2011
 *
 * @brief  VMG::Stencil
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#include "base/stencil.hpp"
#include "comm/comm.hpp"
#include "grid/grid.hpp"
#include "grid/tempgrid.hpp"
#include "mg.hpp"

using namespace VMG;

void Stencil::Apply(Grid& grid) const
{
  Grid::iterator grid_iter;
  Stencil::iterator stencil_iter;

  TempGrid& temp = *MG::GetTempGrid();
  temp.SetProperties(grid);

  MG::GetComm()->CommToGhosts(grid);

  for (grid_iter=grid.Iterators().Local().Begin(); grid_iter!=grid.Iterators().Local().End(); ++grid_iter) {
    temp(*grid_iter) = diag * grid.GetVal(*grid_iter);
    for (stencil_iter=disp.begin(); stencil_iter!=disp.end(); ++stencil_iter)
      temp(*grid_iter) += stencil_iter->Val() * grid.GetVal(*grid_iter + stencil_iter->Disp());
  }

  grid.SetGrid(temp);
}
