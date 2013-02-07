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
 * @file   solver.cpp
 * @author Julian Iseringhausen <isering@ins.uni-bonn.de>
 * @date   Mon Apr 18 13:10:55 2011
 *
 * @brief  VMG::Solver
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "solver/solver.hpp"

using namespace VMG;

void Solver::Run(Grid& sol, Grid& rhs)
{
#ifdef DEBUG_MATRIX_CHECKS
  sol.IsConsistent();
  rhs.IsConsistent();
#endif

#ifdef DEBUG
  sol.IsCompatible(rhs);
#endif

  if (rhs.Global().LocalSize().Product() > 0) {
    this->Realloc(rhs);
    this->AssembleMatrix(rhs);
    this->Compute();
    this->ExportSol(sol, rhs);
  }
}

void Solver::Realloc(int n)
{
  //Reallocate memory if necessary
  this->size = n;

  if (static_cast<int>(A.size()) < n*n)
    A.resize(n*n);

  if (static_cast<int>(b.size()) < n)
    b.resize(n);

  if (static_cast<int>(x.size()) < n)
    x.resize(n);
}

void Solver::Realloc(Grid& sol)
{
  this->Realloc(sol.Global().GlobalSize().Product());
}
