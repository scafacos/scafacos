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
 * @file   grid_iterator_suite.cpp
 * @author Julian Iseringhausen <isering@ins.uni-bonn.de>
 * @date   Thu Apr 21 18:52:43 2011
 *
 * @brief  Some useful grid iterators.
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "grid/grid_iterator_suite.hpp"
#include "grid/grid_properties.hpp"

using namespace VMG;

void GridIteratorSuite::SetSubgrids(const LocalIndices& local_)
{
  local.SetBounds(local_.Begin(), local_.End());
  complete_grid.SetBounds(0, local_.SizeTotal());
  inner_local_grid.SetBounds(local_.Begin()+local_.HaloSize1(), local_.End()-local_.HaloSize2());

  halo_1.X().SetBounds(Index(local_.HaloBegin1().X(), 0, 0),
		       Index(local_.HaloEnd1().X(), local_.SizeTotal().Y(), local_.SizeTotal().Z()));

  halo_1.Y().SetBounds(Index(local_.Begin().X(), local_.HaloBegin1().Y(), 0),
		       Index(local_.End().X(), local_.HaloEnd1().Y(), local_.SizeTotal().Z()));

  halo_1.Z().SetBounds(Index(local_.Begin().X(), local_.Begin().Y(), local_.HaloBegin1().Z()),
		       Index(local_.End().X(), local_.End().Y(), local_.HaloEnd1().Z()));

  halo_2.X().SetBounds(Index(local_.HaloBegin2().X(), 0, 0),
		       Index(local_.HaloEnd2().X(), local_.SizeTotal().Y(), local_.SizeTotal().Z()));

  halo_2.Y().SetBounds(Index(local_.Begin().X(), local_.HaloBegin2().Y(), 0),
		       Index(local_.End().X(), local_.HaloEnd2().Y(), local_.SizeTotal().Z()));

  halo_2.Z().SetBounds(Index(local_.Begin().X(), local_.Begin().Y(), local_.HaloBegin2().Z()),
		       Index(local_.End().X(), local_.End().Y(), local_.HaloEnd2().Z()));

  boundary_1.X().SetBounds(Index(local_.BoundaryBegin1().X(), 0, 0),
			   Index(local_.BoundaryEnd1().X(), local_.SizeTotal().Y(), local_.SizeTotal().Z()));

  boundary_1.Y().SetBounds(Index(local_.Begin().X(), local_.BoundaryBegin1().Y(), 0),
			   Index(local_.End().X(), local_.BoundaryEnd1().Y(), local_.SizeTotal().Z()));

  boundary_1.Z().SetBounds(Index(local_.Begin().X(), local_.Begin().Y(), local_.BoundaryBegin1().Z()),
			   Index(local_.End().X(), local_.End().Y(), local_.BoundaryEnd1().Z()));

  boundary_2.X().SetBounds(Index(local_.BoundaryBegin2().X(), 0, 0),
			   Index(local_.BoundaryEnd2().X(), local_.SizeTotal().Y(), local_.SizeTotal().Z()));


  boundary_2.Y().SetBounds(Index(local_.Begin().X(), local_.BoundaryBegin2().Y(), 0),
			   Index(local_.End().X(), local_.BoundaryEnd2().Y(), local_.SizeTotal().Z()));

  boundary_2.Z().SetBounds(Index(local_.Begin().X(), local_.Begin().Y(), local_.BoundaryBegin2().Z()),
			   Index(local_.End().X(), local_.End().Y(), local_.BoundaryEnd2().Z()));

  near_boundary_1.X().SetBounds(Index(local_.Begin().X(), 0, 0),
				Index(local_.Begin().X()+local_.HaloSize1().X(),
				      local_.SizeTotal().Y(),
				      local_.SizeTotal().Z()));

  near_boundary_1.Y().SetBounds(Index(local_.Begin().X(), local_.Begin().Y(), 0),
				Index(local_.End().X(),
				      local_.Begin().Y()+local_.HaloSize1().Y(),
				      local_.SizeTotal().Z()));

  near_boundary_1.Z().SetBounds(Index(local_.Begin().X(), local_.Begin().Y(), local_.Begin().Z()),
				Index(local_.End().X(),
				      local_.End().Y(),
				      local_.Begin().Z()+local_.HaloSize1().Z()));

  near_boundary_2.X().SetBounds(Index(local_.End().X()-local_.HaloSize2().X(),
				      0,
				      0),
				Index(local_.End().X(), local_.SizeTotal().Y(), local_.SizeTotal().Z()));

  near_boundary_2.Y().SetBounds(Index(local_.Begin().X(),
				      local_.End().Y()-local_.HaloSize2().Y(),
				      0),
				Index(local_.End().X(), local_.End().Y(), local_.SizeTotal().Z()));

  near_boundary_2.Z().SetBounds(Index(local_.Begin().X(),
				      local_.Begin().Y(),
				      local_.End().Z()-local_.HaloSize2().Z()),
				local_.End());
}
