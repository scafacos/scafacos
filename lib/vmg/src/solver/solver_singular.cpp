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
 * @file   solver_singular.cpp
 * @author Julian Iseringhausen <isering@ins.uni-bonn.de>
 * @date   Mon Apr 18 13:12:02 2011
 *
 * @brief  VMG::SolverSingular
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <cmath>
#include <cassert>
#include <iostream>
#include <limits>

#include "base/discretization.hpp"
#include "base/stencil.hpp"
#include "comm/comm.hpp"
#include "grid/multigrid.hpp"
#include "solver/solver_singular.hpp"
#include "mg.hpp"

using namespace VMG;

//TODO: Implement global MPI communication here

void SolverSingular::AssembleMatrix(const Grid& rhs)
{
  Grid::iterator grid_iter;
  Stencil::iterator stencil_iter;
  Index i;
  int index, index2;
  vmg_float row_sum;

  Comm* comm = MG::GetComm();

  const vmg_float prefactor = MG::GetDiscretization()->OperatorPrefactor(rhs);
  const Stencil& A = MG::GetDiscretization()->GetStencil();

  // Make sure that arrays are big enough to hold expanded system of equations
  this->Realloc(rhs.Global().GlobalSize().Product() + 1);

  for (grid_iter = rhs.Iterators().Local().Begin(); grid_iter != rhs.Iterators().Local().End(); ++grid_iter) {

    // Compute 1-dimensional index from 3-dimensional grid
    index = rhs.GlobalLinearIndex(*grid_iter - rhs.Local().Begin() + rhs.Global().LocalBegin());

    // Check if we computed the index correctly
    assert(index >= 0 && index < this->Size()-1);

    // Set solution and right hand side vectors
    this->Sol(index) = 0.0;
    this->Rhs(index) = rhs.GetVal(*grid_iter);

    // Initialize matrix with zeros and then set entries according to the stencil
    for (int l=0; l<this->Size(); l++)
      this->Mat(index,l) = 0.0;

    this->Mat(index,index) = prefactor * A.GetDiag();

    for (stencil_iter = A.begin(); stencil_iter != A.end(); ++stencil_iter) {

      i = *grid_iter - rhs.Local().Begin() + rhs.Global().LocalBegin() + stencil_iter->Disp();

      for (int j=0; j<3; ++j)
	if (comm->BoundaryConditions()[j] == Periodic) {
	  if (i[j] < 0)
	    i[j] += rhs.Global().GlobalSize()[j];
	  else if (i[j] >= rhs.Global().GlobalSize()[j])
	    i[j] -= rhs.Global().GlobalSize()[j];
	}

      // Compute global 1-dimensional index
      index2 = rhs.GlobalLinearIndex(i);

      // Set matrix entry
      this->Mat(index,index2) += prefactor * stencil_iter->Val();
    }
  }

  // Check if matrix has zero row sum (i.e. (1,1,...,1) is an Eigenvector to the Eigenvalue 0)
  row_sum = A.GetDiag();
  for (Stencil::iterator iter=A.begin(); iter!=A.end(); iter++)
    row_sum += iter->Val();

  if (std::abs(row_sum) <= (A.size()+1) * std::numeric_limits<vmg_float>::epsilon()) {

    // Expand equation system in order to make the system regular.
    // The last entry of the solution vector will hold a correction to the right hand side,
    // ensuring that the discrete compatibility condition holds. (Compare Trottenberg)
    for (int i=0; i<this->Size()-1; i++)
      this->Mat(this->Size()-1, i) = this->Mat(i, this->Size()-1) = 1.0;

    this->Mat(this->Size()-1, this->Size()-1) = 0.0;
    this->Sol(this->Size()-1) = 0.0;
    this->Rhs(this->Size()-1) = 0.0;

  }else {
    //TODO: Implement this
    assert(0 == "At the first glance your stencil does not seem to be singular. Try SolverRegular instead.");
  }
}

void SolverSingular::ExportSol(Grid& sol, Grid& rhs)
{
  int index;
  const vmg_float correction = this->Sol(this->Size()-1);

  for (int i=0; i<sol.Local().Size().X(); i++)
    for (int j=0; j<sol.Local().Size().Y(); j++)
      for (int k=0; k<sol.Local().Size().Z(); k++) {

      	// Compute global 1-dimensional index
      	index = sol.GlobalLinearIndex(sol.Global().LocalBegin().X()+i,
				      sol.Global().LocalBegin().Y()+j,
				      sol.Global().LocalBegin().Z()+k);

      	// Set solution
      	sol(sol.Local().Begin().X()+i, sol.Local().Begin().Y()+j, sol.Local().Begin().Z()+k) = this->Sol(index) - correction;

      }
}
