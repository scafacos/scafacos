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
 * @file   tempgrid.hpp
 * @author Julian Iseringhausen <isering@ins.uni-bonn.de>
 * @date   Mon Apr 18 12:55:18 2011
 *
 * @brief  VMG::TempGrid
 *
 */

#ifndef TEMPGRID_HPP_
#define TEMPGRID_HPP_

#include "grid/grid.hpp"

namespace VMG
{

class Index;
class Vector;

class TempGrid : public Grid
{
public:
  TempGrid();
  TempGrid(const Grid& rhs);
  TempGrid(const GlobalIndices& global, const LocalIndices& local, const SpatialExtent& extent);
  TempGrid(const Index& size, const Index& halo_size,
           const Vector& spatial_begin, const Vector& spatial_end);
  virtual ~TempGrid();

  void SetProperties(const Grid& rhs);
  void SetProperties(const GlobalIndices& global, const LocalIndices& local, const SpatialExtent& extent);
  void SetProperties(const Index& size, const Index& halo_size,
                     const Vector& spatial_begin, const Vector& spatial_end);

  void SetPropertiesToFiner(const Grid& grid, const Boundary& boundary);
  void SetPropertiesToCoarser(const Grid& grid, const Boundary& boundary);

  void ImportFromResidual(Grid& sol, Grid& rhs);

private:
  void Allocate();

  int size_max;
};

}

#endif /* TEMPGRID_HPP_ */
