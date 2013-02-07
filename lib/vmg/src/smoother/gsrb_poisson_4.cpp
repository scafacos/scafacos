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
 * @file   gsrb_poisson_4.cpp
 * @author Julian Iseringhausen <isering@ins.uni-bonn.de>
 * @date   Fri May 11 18:30:20 2012
 *
 * @brief  Gauss-Seidel Red Black method, specialized to
 *         the Poisson equation. Performance improved by
 *         explicit loop unrolling.
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#include "base/helper.hpp"
#include "comm/comm.hpp"
#include "grid/grid.hpp"
#include "smoother/gsrb_poisson_4.hpp"
#include "mg.hpp"

using namespace VMG;

static inline void ComputePartial(Grid& sol, Grid& rhs,
				  const Index& begin, const Index& end,
				  const vmg_float& prefactor, const int& off)
{
  const vmg_float fac_1 = 1.0 / 12.0;
  const vmg_float fac_2 = 1.0 / 24.0;

  for (int i=begin.X(); i<end.X(); ++i)
    for (int j=begin.Y(); j<end.Y(); ++j)
      for (int k=begin.Z() + (i + j + begin.Z() + off) % 2; k<end.Z(); k+=2)
	sol(i,j,k) = prefactor * rhs.GetVal(i,j,k) + fac_1 * (sol.GetVal(i-1,j  ,k  ) +
							      sol.GetVal(i+1,j  ,k  ) +
							      sol.GetVal(i  ,j-1,k  ) +
							      sol.GetVal(i  ,j+1,k  ) +
							      sol.GetVal(i  ,j  ,k-1) +
							      sol.GetVal(i  ,j  ,k+1))
	                                           + fac_2 * (sol.GetVal(i-1,j-1,k  ) +
							      sol.GetVal(i-1,j+1,k  ) +
							      sol.GetVal(i+1,j-1,k  ) +
							      sol.GetVal(i+1,j+1,k  ) +
							      sol.GetVal(i-1,j  ,k-1) +
							      sol.GetVal(i-1,j  ,k+1) +
							      sol.GetVal(i+1,j  ,k-1) +
							      sol.GetVal(i+1,j  ,k+1) +
							      sol.GetVal(i  ,j-1,k-1) +
							      sol.GetVal(i  ,j-1,k+1) +
							      sol.GetVal(i  ,j+1,k-1) +
							      sol.GetVal(i  ,j+1,k+1));
}

void GaussSeidelRBPoisson4::Compute(Grid& sol, Grid& rhs)
{
  const vmg_float prefactor_inv = Helper::pow_2(sol.Extent().MeshWidth().Max()) / 4.0;
  const int off = rhs.Global().LocalBegin().Sum() - rhs.Local().HaloSize1().Sum();
  const LocalIndices& local = rhs.Local();
  Comm& comm = *MG::GetComm();

  /*
   * Compute first halfstep
   */

  // Start asynchronous communication
  comm.CommToGhostsAsyncStart(sol);

  // Smooth part not depending on ghost cells
  ComputePartial(sol, rhs,
		 local.Begin()+1, local.End()-1,
		 prefactor_inv, off+1);

  // Finish asynchronous communication
  comm.CommToGhostsAsyncFinish(sol);

  /*
   * Smooth near boundary cells
   */

  ComputePartial(sol, rhs,
		 local.Begin(),
		 Index(local.Begin().X()+1, local.End().Y(), local.End().Z()),
		 prefactor_inv, off+1);

  ComputePartial(sol, rhs,
		 Index(local.End().X()-1, local.Begin().Y(), local.Begin().Z()),
		 local.End(),
		 prefactor_inv, off+1);

  ComputePartial(sol, rhs,
		 Index(local.Begin().X()+1, local.Begin().Y(), local.Begin().Z()),
		 Index(local.End().X()-1, local.Begin().Y()+1, local.End().Z()),
		 prefactor_inv, off+1);

  ComputePartial(sol, rhs,
		 Index(local.Begin().X()+1, local.End().Y()-1, local.Begin().Z()),
		 Index(local.End().X()-1, local.End().Y(), local.End().Z()),
		 prefactor_inv, off+1);

  ComputePartial(sol, rhs,
		 Index(local.Begin().X()+1, local.Begin().Y()+1, local.Begin().Z()),
		 Index(local.End().X()-1, local.End().Y()-1, local.Begin().Z()+1),
		 prefactor_inv, off+1);

  ComputePartial(sol, rhs,
		 Index(local.Begin().X()+1, local.Begin().Y()+1, local.End().Z()-1),
		 Index(local.End().X()-1, local.End().Y()-1, local.End().Z()),
		 prefactor_inv, off+1);

  /*
   * Compute second halfstep
   */

  // Start asynchronous communication
  comm.CommToGhostsAsyncStart(sol);

  // Smooth part not depending on ghost cells
  ComputePartial(sol, rhs,
		 local.Begin()+1, local.End()-1,
		 prefactor_inv, off);

  // Finish asynchronous communication
  comm.CommToGhostsAsyncFinish(sol);

  /*
   * Smooth near boundary cells
   */

  ComputePartial(sol, rhs,
		 local.Begin(),
		 Index(local.Begin().X()+1, local.End().Y(), local.End().Z()),
		 prefactor_inv, off);

  ComputePartial(sol, rhs,
		 Index(local.End().X()-1, local.Begin().Y(), local.Begin().Z()),
		 local.End(),
		 prefactor_inv, off);

  ComputePartial(sol, rhs,
		 Index(local.Begin().X()+1, local.Begin().Y(), local.Begin().Z()),
		 Index(local.End().X()-1, local.Begin().Y()+1, local.End().Z()),
		 prefactor_inv, off);

  ComputePartial(sol, rhs,
		 Index(local.Begin().X()+1, local.End().Y()-1, local.Begin().Z()),
		 Index(local.End().X()-1, local.End().Y(), local.End().Z()),
		 prefactor_inv, off);

  ComputePartial(sol, rhs,
		 Index(local.Begin().X()+1, local.Begin().Y()+1, local.Begin().Z()),
		 Index(local.End().X()-1, local.End().Y()-1, local.Begin().Z()+1),
		 prefactor_inv, off);

  ComputePartial(sol, rhs,
		 Index(local.Begin().X()+1, local.Begin().Y()+1, local.End().Z()-1),
		 Index(local.End().X()-1, local.End().Y()-1, local.End().Z()),
		 prefactor_inv, off);
}
