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
 * @file   global_grid_partitioning.hpp
 * @author Julian Iseringhausen <isering@ins.uni-bonn.de>
 * @date   Thu May 19 14:02:46 2011
 *
 * @brief  Class to store global grid partitioning.
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <cassert>

#include "grid/global_grid_partitioning.hpp"

using namespace VMG;

GlobalGridPartitioning::GlobalGridPartitioning(Index pos, Index procs, Index size, int points_min)
{
  Index i;
  Index index_beg, index_end;
  Index local_size;

  const Index active_procs = size / points_min;
  const Index size_per_proc = size / active_procs;
  const Index remainder = size % active_procs;

  for (i.X() = 0; i.X() < active_procs.X(); ++i.X())
    for (i.Y() = 0; i.Y() < active_procs.Y(); ++i.Y())
      for (i.Z() = 0; i.Z() < active_procs.Z(); ++i.Z()) {

	local_size = size_per_proc;

	for (int j=0; j<3; ++j)
	  if (pos[j] < remainder[j])
	    ++(local_size[j]);

	index_beg = i * local_size;

	for (int j=0; j<3; ++j)
	  if (i[j] >= remainder[j])
	    index_beg[j] += remainder[j];

	index_end = index_beg + local_size;

	begin.insert(std::make_pair(i, index_beg));
	end.insert(std::make_pair(i, index_end));

      }

}

const Index& GlobalGridPartitioning::GlobalBegin(const Index& pos) const
{
  std::map<Index, Index>::const_iterator iter = begin.find(pos);

  assert(iter != begin.end());

  return iter->second;
}

const Index& GlobalGridPartitioning::GlobalEnd(const Index& pos) const
{
  std::map<Index, Index>::const_iterator iter = end.find(pos);

  assert(iter != end.end());

  return iter->second;
}
