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
 * @file   solver_regular.cpp
 * @author Julian Iseringhausen <isering@ins.uni-bonn.de>
 * @date   Mon Apr 18 13:11:32 2011
 *
 * @brief  VMG::SolverRegular
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>

#include "base/discretization.hpp"
#include "base/stencil.hpp"
#include "grid/multigrid.hpp"
#include "solver/solver_regular.hpp"
#include "mg.hpp"

using namespace VMG;

// TODO: Implement global communication here
// TODO: Implement this more efficiently
// TODO: This has to be reviewed for parallelization

void SolverRegular::AssembleMatrix(const Grid& rhs)
{
  Grid::iterator grid_iter;
  Stencil::iterator stencil_iter;
  int mat_index, mat_index2;
  vmg_float prefactor_inv = 1.0 / MG::GetDiscretization()->OperatorPrefactor(rhs);
  const Stencil& A = MG::GetDiscretization()->GetStencil();

#ifdef DEBUG_MATRIX_CHECKS
  rhs.IsConsistent();
#endif

  this->Realloc(rhs.Global().GlobalSize().Product());

  for (grid_iter = rhs.Iterators().Local().Begin(); grid_iter != rhs.Iterators().Local().End(); ++grid_iter) {

    mat_index = rhs.GlobalLinearIndex(*grid_iter + rhs.Global().LocalBegin());

    assert(mat_index >= 0 && mat_index<this->Size());

    this->Sol(mat_index) = 0.0;
    this->Rhs(mat_index) = prefactor_inv * rhs.GetVal(*grid_iter);

    for (int l=0; l<this->Size(); l++)
      this->Mat(mat_index, l) = 0.0;

    this->Mat(mat_index, mat_index) = A.GetDiag();

    for (stencil_iter = A.begin(); stencil_iter != A.end(); ++stencil_iter) {

      mat_index2 = rhs.GlobalLinearIndex(*grid_iter + rhs.Global().LocalBegin() + stencil_iter->Disp());

      assert(mat_index2 >= 0 && mat_index2<this->Size());

      this->Mat(mat_index, mat_index2) += stencil_iter->Val();

    }
  }

  for (int i=0; i<3; ++i) {

    for (grid_iter = rhs.Iterators().Boundary1()[i].Begin(); grid_iter != rhs.Iterators().Boundary1()[i].End(); ++grid_iter) {

      mat_index = rhs.GlobalLinearIndex(*grid_iter + rhs.Global().LocalBegin());

      assert(mat_index >= 0 && mat_index<this->Size());

      this->Sol(mat_index) = this->Rhs(mat_index) = rhs.GetVal(*grid_iter);

      for (int l=0; l<this->Size(); l++)
	this->Mat(mat_index, l) = 0.0;

      this->Mat(mat_index, mat_index) = 1.0;

    }

    for (grid_iter = rhs.Iterators().Boundary2()[i].Begin(); grid_iter != rhs.Iterators().Boundary2()[i].End(); ++grid_iter) {

      mat_index = rhs.GlobalLinearIndex(*grid_iter + rhs.Global().LocalBegin());

      assert(mat_index >= 0 && mat_index<this->Size());

      this->Sol(mat_index) = this->Rhs(mat_index) = rhs.GetVal(*grid_iter);

      for (int l=0; l<this->Size(); l++)
	this->Mat(mat_index, l) = 0.0;

      this->Mat(mat_index, mat_index) = 1.0;

    }

  }

}

void SolverRegular::ExportSol(Grid& sol, Grid& rhs)
{
  int index;
  Index offset;

  for (int i=0; i<3; ++i)
    offset[i] = (sol.Local().HaloEnd1()[i] > 0 ? 1 : 0);

  for (Grid::iterator iter = sol.Iterators().CompleteGrid().Begin(); iter != sol.Iterators().CompleteGrid().End(); ++iter) {
    index = sol.GlobalLinearIndex(sol.Global().LocalBegin() + *iter - offset);
    sol(*iter) = this->Sol(index);
  }

#ifdef DEBUG_MATRIX_CHECKS
  sol.IsConsistent();
#endif
}
