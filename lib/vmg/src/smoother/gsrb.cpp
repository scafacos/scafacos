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
 * @file   gsrb.cpp
 * @author Julian Iseringhausen <isering@ins.uni-bonn.de>
 * @date   Mon Apr 18 13:08:20 2011
 *
 * @brief  Gauss-Seidel Red Black method
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
#include "smoother/gsrb.hpp"

using namespace VMG;

static inline void ComputePartial(Grid& sol, Grid& rhs, const Stencil& mat,
				  const Index& begin, const Index& end,
				  const vmg_float& prefactor, const vmg_float& diag_inv,
				  const int& off)
{
  int i,j,k;
  vmg_float temp;
  Stencil::iterator iter;

  for (i=begin.X(); i<end.X(); ++i)
    for (j=begin.Y(); j<end.Y(); ++j) {

      int z_begin = begin.Z() + (i + j + begin.Z() + off) % 2;

#ifdef DEBUG
      int off_sum = MG::GetComm()->LevelSum(rhs, z_begin - begin.Z());
      assert(z_begin - begin.Z() == 0 || z_begin - begin.Z() == 1);
      assert(off_sum == 0 || off_sum == MG::GetComm()->Size(rhs));
#endif /* DEBUG */

      for (k=z_begin; k<end.Z(); k+=2) {

	temp = prefactor * rhs.GetVal(i,j,k);

	for (iter=mat.begin(); iter!=mat.end(); ++iter)
	  temp -= iter->Val() * sol.GetVal(i+iter->Disp().X(),
					   j+iter->Disp().Y(),
					   k+iter->Disp().Z());

	sol(i,j,k) = temp * diag_inv;

      }
    }
}

void GaussSeidelRB::Compute(Grid& sol, Grid& rhs)
{
  const Stencil& mat = MG::GetDiscretization()->GetStencil();
  const vmg_float prefactor_inv = 1.0 / MG::GetDiscretization()->OperatorPrefactor(sol);
  const vmg_float diag_inv = 1.0 / mat.GetDiag();
  const int off = rhs.Global().LocalBegin().Sum() - rhs.Local().HaloSize1().Sum();
  const LocalIndices& local = rhs.Local();
  Comm& comm = *MG::GetComm();

  /*
   * Compute first halfstep
   */

  // Start asynchronous communication
  comm.CommToGhostsAsyncStart(sol);

  // Smooth part not depending on ghost cells
  ComputePartial(sol, rhs, mat,
		 local.Begin()+1, local.End()-1,
		 prefactor_inv, diag_inv, off+1);

  // Finish asynchronous communication
  comm.CommToGhostsAsyncFinish(sol);

  /*
   * Smooth near boundary cells
   */

  ComputePartial(sol, rhs, mat,
		 local.Begin(),
		 Index(local.Begin().X()+1, local.End().Y(), local.End().Z()),
		 prefactor_inv, diag_inv, off+1);

  ComputePartial(sol, rhs, mat,
		 Index(local.End().X()-1, local.Begin().Y(), local.Begin().Z()),
		 local.End(),
		 prefactor_inv, diag_inv, off+1);

  ComputePartial(sol, rhs, mat,
		 Index(local.Begin().X()+1, local.Begin().Y(), local.Begin().Z()),
		 Index(local.End().X()-1, local.Begin().Y()+1, local.End().Z()),
		 prefactor_inv, diag_inv, off+1);

  ComputePartial(sol, rhs, mat,
		 Index(local.Begin().X()+1, local.End().Y()-1, local.Begin().Z()),
		 Index(local.End().X()-1, local.End().Y(), local.End().Z()),
		 prefactor_inv, diag_inv, off+1);

  ComputePartial(sol, rhs, mat,
		 Index(local.Begin().X()+1, local.Begin().Y()+1, local.Begin().Z()),
		 Index(local.End().X()-1, local.End().Y()-1, local.Begin().Z()+1),
		 prefactor_inv, diag_inv, off+1);

  ComputePartial(sol, rhs, mat,
		 Index(local.Begin().X()+1, local.Begin().Y()+1, local.End().Z()-1),
		 Index(local.End().X()-1, local.End().Y()-1, local.End().Z()),
		 prefactor_inv, diag_inv, off+1);

  /*
   * Compute second halfstep
   */

  // Start asynchronous communication
  comm.CommToGhostsAsyncStart(sol);

  // Smooth part not depending on ghost cells
  ComputePartial(sol, rhs, mat,
		 local.Begin()+1, local.End()-1,
		 prefactor_inv, diag_inv, off);

  // Finish asynchronous communication
  comm.CommToGhostsAsyncFinish(sol);

  /*
   * Smooth near boundary cells
   */

  ComputePartial(sol, rhs, mat,
		 local.Begin(),
		 Index(local.Begin().X()+1, local.End().Y(), local.End().Z()),
		 prefactor_inv, diag_inv, off);

  ComputePartial(sol, rhs, mat,
		 Index(local.End().X()-1, local.Begin().Y(), local.Begin().Z()),
		 local.End(),
		 prefactor_inv, diag_inv, off);

  ComputePartial(sol, rhs, mat,
		 Index(local.Begin().X()+1, local.Begin().Y(), local.Begin().Z()),
		 Index(local.End().X()-1, local.Begin().Y()+1, local.End().Z()),
		 prefactor_inv, diag_inv, off);

  ComputePartial(sol, rhs, mat,
		 Index(local.Begin().X()+1, local.End().Y()-1, local.Begin().Z()),
		 Index(local.End().X()-1, local.End().Y(), local.End().Z()),
		 prefactor_inv, diag_inv, off);

  ComputePartial(sol, rhs, mat,
		 Index(local.Begin().X()+1, local.Begin().Y()+1, local.Begin().Z()),
		 Index(local.End().X()-1, local.End().Y()-1, local.Begin().Z()+1),
		 prefactor_inv, diag_inv, off);

  ComputePartial(sol, rhs, mat,
		 Index(local.Begin().X()+1, local.Begin().Y()+1, local.End().Z()-1),
		 Index(local.End().X()-1, local.End().Y()-1, local.End().Z()),
		 prefactor_inv, diag_inv, off);
}
