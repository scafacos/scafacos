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
 * @file   domain_decomposition_mpi.cpp
 * @author Julian Iseringhausen <isering@ins.uni-bonn.de>
 * @date   Mon Jun 27 12:53:50 2011
 *
 * @brief  Computes a domain decomposition which separates
 *         the finest grid equally for all processes.
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "base/interface.hpp"
#include "comm/comm.hpp"
#include "comm/domain_decomposition_mpi.hpp"
#include "grid/grid.hpp"
#include "grid/multigrid.hpp"

using namespace VMG;

void DomainDecompositionMPI::Compute(Comm* comm, const Interface* interface, std::vector<GlobalIndices>& global)
{
  GlobalIndices global_l;
  Index remainder, procs;
  Index last_procs = comm->GlobalProcs();

  global.clear();

  for (unsigned int i=0; i<interface->Global().size(); ++i) {

    /*
     * Inherit global properties from interface
     */
    global_l.GlobalFinerBegin() = interface->Global()[i].GlobalFinerBegin();
    global_l.GlobalFinerEnd() = interface->Global()[i].GlobalFinerEnd();
    global_l.GlobalFinerSize() = interface->Global()[i].GlobalFinerSize();
    global_l.FinestAbsBegin() = interface->Global()[i].FinestAbsBegin();
    global_l.FinestAbsEnd() = interface->Global()[i].FinestAbsEnd();
    global_l.FinestAbsSize() = interface->Global()[i].FinestAbsSize();
    global_l.GlobalSize() = interface->Global()[i].GlobalSize();
    global_l.BoundaryType() = interface->Global()[i].BoundaryType();

    if (IsActive(comm, global_l.GlobalSize(), procs)) {

      if (i == 0) {

	global_l.LocalSize() = global_l.GlobalSize() / procs;

	remainder = global_l.GlobalSize() % procs;
	for (int j=0; j<3; ++j)
	  if (comm->GlobalPos()[j] < remainder[j])
	    ++(global_l.LocalSize()[j]);

	global_l.LocalBegin() = comm->GlobalPos() * global_l.LocalSize();

	for (int j=0; j<3; ++j)
	  if (comm->GlobalPos()[j] >= remainder[j])
	    global_l.LocalBegin()[j] += remainder[j];

	global_l.LocalEnd() = global_l.LocalBegin() + global_l.LocalSize();

	global_l.LocalFinerBegin() = 0;
	global_l.LocalFinerEnd() = 0;
	global_l.LocalFinerSize() = 0;

      }else {

	for (int j=0; j<3; ++j) {

	  if (procs[j] == last_procs[j]) {

	    if (global.back().LocalBegin()[j] == 0)
	      global_l.LocalBegin()[j] = 0;
	    else
	      global_l.LocalBegin()[j] = global.back().LocalBegin()[j] / 2 + global_l.GlobalFinerBegin()[j];

	    if (global.back().LocalEnd()[j] == global.back().GlobalSize()[j])
	      global_l.LocalEnd()[j] = global_l.GlobalSize()[j];
	    else
	      global_l.LocalEnd()[j] = global.back().LocalEnd()[j] / 2 + global_l.GlobalFinerBegin()[j];

	    global_l.LocalSize()[j] = global_l.LocalEnd()[j] - global_l.LocalBegin()[j];

	  }else {

	    global_l.LocalSize()[j] = global_l.GlobalSize()[j] / procs[j];

	    remainder[j] = global_l.GlobalSize()[j] % procs[j];
	    if (comm->GlobalPos()[j] < remainder[j])
	      ++(global_l.LocalSize()[j]);

	    global_l.LocalBegin()[j] = comm->GlobalPos()[j] * global_l.LocalSize()[j];

	    if (comm->GlobalPos()[j] >= remainder[j])
	      global_l.LocalBegin()[j] += remainder[j];

	    global_l.LocalEnd()[j] = global_l.LocalBegin()[j] + global_l.LocalSize()[j];

	  }
	}

	global_l.LocalFinerBegin() = global_l.LocalBegin().Clamp(global_l.GlobalFinerBegin(), global_l.GlobalFinerEnd());
	global_l.LocalFinerEnd() = global_l.LocalEnd().Clamp(global_l.GlobalFinerBegin(), global_l.GlobalFinerEnd());
	global_l.LocalFinerSize() = global_l.LocalFinerEnd() - global_l.LocalFinerBegin();

      }

    }else {

      global_l.LocalBegin() = 0;
      global_l.LocalEnd() = 0;
      global_l.LocalSize() = 0;
      global_l.LocalFinerBegin() = 0;
      global_l.LocalFinerEnd() = 0;
      global_l.LocalFinerSize() = 0;

    }

    last_procs = procs;

    global.push_back(global_l);

  }
}

bool DomainDecompositionMPI::IsActive(Comm* comm, const Index& size_global, Index& procs)
{
  bool is_active = true;
  const int points_min = 5;

  procs = size_global / points_min + 1;

  for (int i=0; i<3; ++i) {
    procs[i] = std::min(procs[i], comm->GlobalProcs()[i]);
    is_active &= comm->GlobalPos()[i] < procs[i];
  }

  return is_active;
}

void DomainDecompositionMPI::FineToCoarse(Comm* comm, int& begin, int& end, int levels)
{
  int last_point = end - 1;

  for (int i=0; i<levels; ++i) {

    if (begin % 2 == 0)
      begin /= 2;
    else
      begin = (begin+1) / 2;

    if (last_point % 2 == 0)
      last_point /= 2;
    else
      last_point = (last_point-1) / 2;

  }

end = last_point + 1;
}
