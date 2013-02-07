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
 * @file   grid.hpp
 * @author Julian Iseringhausen <isering@ins.uni-bonn.de>
 * @date   Mon Apr 18 12:53:45 2011
 *
 * @brief  VMG::Grid
 *
 */

#ifndef GRID_HPP_
#define GRID_HPP_

#include "base/object.hpp"
#include "grid/grid_index_translations.hpp"
#include "grid/grid_iterator.hpp"
#include "grid/grid_iterator_suite.hpp"
#include "grid/grid_properties.hpp"

namespace VMG
{

class Comm;
class Multigrid;
class Stencil;

class Grid : public Object
{
public:
  typedef GridIterator iterator;

  Grid(int level_ = 0, Multigrid* father_ = NULL) :
    index_translations(this),
    level(level_),
    father(father_)
  {
    grid = NULL;
  }

  Grid(const GlobalIndices& global_,
       const LocalIndices& local_,
       const SpatialExtent& extent_,
       int level_ = 0,
       Multigrid* father_ = NULL) :
    index_translations(this),
    level(level_),
    global(global_),
    local(local_),
    extent(extent_),
    iterators(local_),
    father(father_)
  {
    InitGrid();
  }

  Grid(const Grid& rhs) :
    index_translations(rhs.Indexing()),
    level(rhs.Level()),
    global(rhs.Global()),
    local(rhs.Local()),
    extent(rhs.Extent()),
    iterators(rhs.Iterators()),
    father(rhs.Father())
  {
    InitGrid();
    SetGrid(rhs);
  }

  virtual ~Grid();

  Grid& operator=(const Grid& rhs);

  GlobalIndices& Global() {return global;}
  LocalIndices& Local() {return local;}
  SpatialExtent& Extent() {return extent;}

  const GlobalIndices& Global() const {return global;}
  const LocalIndices& Local() const {return local;}
  const SpatialExtent& Extent() const {return extent;}

  GridIteratorSuite& Iterators() {return iterators;}
  const GridIteratorSuite& Iterators() const {return iterators;}

  void Clear();         ///< Overwrites all grid points on current level with zeros
  void ClearInner();
  void ClearHalo();     ///< Overwrites all halo points on current level with zeros
  void ClearBoundary(); ///< Overwrites all boundary points on current level with zeros

  vmg_float& operator()(int x, int y, int z);  ///< Returns a reference to the requested gridpoint.
  vmg_float& operator()(const Index& index);

  const vmg_float& GetVal(int x, int y, int z) const; ///< Returns the value of a requested gridpoint.
  const vmg_float& GetVal(const Index& index) const;

  void ForceDiscreteCompatibilityCondition();
  void SetAverageToZero();

  void SetGrid(const Grid& rhs);     ///< Overwrite current grid with values from another grid
  void SetBoundary(const Grid& rhs); ///< Overwrite boundary with values from rhs

  void AddGrid(const Grid& rhs);             ///< Add values of another grid
  void SubtractGrid(const Grid& rhs);        ///< Subtract values of another grid
  void MultiplyScalar(const vmg_float& scalar);   ///< Multiply grid values with scalar
  void ApplyStencil(const Stencil& stencil); ///< Apply stencil to grid

  int GlobalLinearIndex(int x, int y, int z) const;  ///< Compute a unique 1-dimensional global index
  int GlobalLinearIndex(const Index& index) const;

  bool IsCompatible(const Grid& rhs) const; ///< Check if two grids share compatible settings
  bool IsConsistent() const;                  ///< Check grid for nan and inf

  Multigrid* Father() const {return father;}

  virtual vmg_float DebugKnownSolution(Vector x) const {return 0.0;}
  vmg_float DebugKnownSolution(Index i) const {return DebugKnownSolution(extent.Begin() + i * extent.MeshWidth());}

  const GridIndexTranslations& Indexing() const {return index_translations;}

  const int& Level() const {return level;}

  bool IsActive() const {return Local().Size().Product() > 0;}

private:
  void InitGrid();

  GridIndexTranslations index_translations;

protected:
  int level;

  GlobalIndices global;
  LocalIndices local;
  SpatialExtent extent;

  GridIteratorSuite iterators;

  vmg_float *grid;

  Multigrid* father;
};

inline vmg_float& Grid::operator()(int x, int y, int z)
{
  return grid[z + local.SizeTotal().Z() * (y + local.SizeTotal().Y() * x)];
}

inline vmg_float& Grid::operator()(const Index& index)
{
  return this->operator()(index.X(), index.Y(), index.Z());
}

inline const vmg_float& Grid::GetVal(int x, int y, int z) const
{
  return grid[z + local.SizeTotal().Z() * (y + local.SizeTotal().Y() * x)];
}

inline const vmg_float& Grid::GetVal(const Index& index) const
{
  return this->GetVal(index.X(), index.Y(), index.Z());
}

inline int Grid::GlobalLinearIndex(int x, int y, int z) const
{
  return z + global.GlobalSize().Z() * (y + global.GlobalSize().Y() * x);
}

inline int Grid::GlobalLinearIndex(const Index& index) const
{
  return GlobalLinearIndex(index.X(), index.Y(), index.Z());
}

}

#endif /* GRID_HPP_ */
