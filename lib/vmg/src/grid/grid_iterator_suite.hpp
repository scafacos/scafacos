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
 * @file   grid_iterator_suite.hpp
 * @author Julian Iseringhausen <isering@ins.uni-bonn.de>
 * @date   Thu Apr 21 18:48:09 2011
 *
 * @brief  Some useful grid iterators.
 *
 */

#ifndef GRID_ITERATOR_SUITE_HPP_
#define GRID_ITERATOR_SUITE_HPP_

#include "grid/grid_iterator_set.hpp"
#include "grid/grid_iterator.hpp"
#include "grid/grid_properties.hpp"

namespace VMG
{

class GridIteratorSuite
{
public:
  GridIteratorSuite()
  {
    LocalIndices local;
    SetSubgrids(local);
  }

  GridIteratorSuite(const LocalIndices& local)
  {
    SetSubgrids(local);
  }

  GridIteratorSuite(const GridIteratorSuite& other) :
    local(other.local),
    complete_grid(other.complete_grid),
    inner_local_grid(other.inner_local_grid),
    halo_1(other.halo_1),
    halo_2(other.halo_2),
    boundary_1(other.boundary_1),
    boundary_2(other.boundary_2),
    near_boundary_1(other.near_boundary_1),
    near_boundary_2(other.near_boundary_2)
  {
  }

  virtual ~GridIteratorSuite() {}

  void SetSubgrids(const LocalIndices& local);

  const GridIteratorSet& Local() const {return local;}
  const GridIteratorSet& CompleteGrid() const {return complete_grid;}
  const GridIteratorSet& InnerLocalGrid() const {return inner_local_grid;}

  const GridIteratorSet3& Halo1() const {return halo_1;}
  const GridIteratorSet3& Halo2() const {return halo_2;}

  const GridIteratorSet3& Boundary1() const {return boundary_1;}
  const GridIteratorSet3& Boundary2() const {return boundary_2;}

  const GridIteratorSet3& NearBoundary1() const {return near_boundary_1;}
  const GridIteratorSet3& NearBoundary2() const {return near_boundary_2;}

private:
  GridIteratorSet local;
  GridIteratorSet complete_grid;
  GridIteratorSet inner_local_grid;

  GridIteratorSet3 halo_1, halo_2;
  GridIteratorSet3 boundary_1, boundary_2;
  GridIteratorSet3 near_boundary_1, near_boundary_2;
};

}

#endif /* GRID_ITERATOR_SUITE_HPP_ */
