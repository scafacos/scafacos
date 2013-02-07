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
 * @file   is_grid.hpp
 * @author Julian Iseringhausen <isering@ins.uni-bonn.de>
 * @date   Mon Apr 18 12:53:45 2011
 *
 * @brief  Grid-like class that holds arbitrary data.
 *
 */

#ifndef IS_GRID_HPP_
#define IS_GRID_HPP_

#include "base/object.hpp"
#include "grid/grid.hpp"
#include "grid/grid_index_translations.hpp"
#include "grid/grid_iterator.hpp"
#include "grid/grid_iterator_suite.hpp"
#include "grid/grid_properties.hpp"

namespace VMG
{

class Comm;
class Grid;
class Multigrid;
class Stencil;

template <class T>
class IsGrid
{
public:
  typedef GridIterator iterator;

  IsGrid(int level = 0) :
    level(level)
  {
    grid = NULL;
  }

  IsGrid(const GlobalIndices& global,
	 const LocalIndices& local,
	 const SpatialExtent& extent,
	 int level = 0) :
    level(level),
    global(global),
    local(local),
    extent(extent),
    iterators(local)
  {
    InitIsGrid();
  }

  IsGrid(const IsGrid& rhs) :
    level(rhs.Level()),
    global(rhs.Global()),
    local(rhs.Local()),
    extent(rhs.Extent()),
    iterators(rhs.Iterators())
  {
    InitIsGrid();
    SetGrid(rhs);
  }

  IsGrid(const Grid& rhs) :
    level(rhs.Level()),
    global(rhs.Global()),
    local(rhs.Local()),
    extent(rhs.Extent()),
    iterators(rhs.Iterators())
  {
    InitIsGrid();
  }

  virtual ~IsGrid();

  void SetGridSize(const Grid& grid);
  void SetGridSize(const GlobalIndices& global,
		   const LocalIndices& local,
		   const SpatialExtent& extent);

  IsGrid& operator=(const IsGrid& rhs);

  const GlobalIndices& Global() const {return global;}
  const LocalIndices& Local() const {return local;}
  const SpatialExtent& Extent() const {return extent;}

  GlobalIndices& Global() {return global;}
  LocalIndices& Local() {return local;}
  SpatialExtent& Extent() {return extent;}

  GridIteratorSuite& Iterators() {return iterators;}
  const GridIteratorSuite& Iterators() const {return iterators;}

  const vmg_float& MeshWidth() const {return Extent().MeshWidth().X();} ///< Mesh width of current level

  T& operator()(int x, int y, int z);  ///< Returns a reference to the requested gridpoint.
  T& operator()(const Index& index);

  const T& GetVal(int x, int y, int z) const; ///< Returns the value of a requested gridpoint.
  const T& GetVal(const Index& index) const;

  void SetGrid(const IsGrid& rhs);     ///< Overwrite current grid with values from another grid
  void SetBoundary(const IsGrid& rhs); ///< Overwrite boundary with values from rhs

  int GlobalLinearIndex(int x, int y, int z) const;  ///< Compute a unique 1-dimensional global index
  int GlobalLinearIndex(const Index& index) const;

  const int& Level() const {return level;}

  bool IsActive() const {return Local().Size().Product() > 0;}

private:
  void InitIsGrid();

protected:
  int level;

  GlobalIndices global;
  LocalIndices local;
  SpatialExtent extent;

  GridIteratorSuite iterators;

  T* grid;

  static vmg_float correction;
};

template <class T>
inline T& IsGrid<T>::operator()(int x, int y, int z)
{
  return grid[z + local.SizeTotal().Z() * (y + local.SizeTotal().Y() * x)];
}

template <class T>
inline T& IsGrid<T>::operator()(const Index& index)
{
  return this->operator()(index.X(), index.Y(), index.Z());
}

template <class T>
inline const T& IsGrid<T>::GetVal(int x, int y, int z) const
{
  return grid[z + local.SizeTotal().Z() * (y + local.SizeTotal().Y() * x)];
}

template <class T>
inline const T& IsGrid<T>::GetVal(const Index& index) const
{
  return this->GetVal(index.X(), index.Y(), index.Z());
}

template <class T>
inline int IsGrid<T>::GlobalLinearIndex(int x, int y, int z) const
{
  return z + global.GlobalSize().Z() * (y + global.GlobalSize().Y() * x);
}

template <class T>
inline int IsGrid<T>::GlobalLinearIndex(const Index& index) const
{
  return index.Z() + global.GlobalSize().Z() * (index.Y() + global.GlobalSize().Y() * index.X());
}

template <class T>
void IsGrid<T>::InitIsGrid()
{
  grid = new T[local.SizeTotal().Product()];
}

template <class T>
IsGrid<T>::~IsGrid()
{
  delete [] grid;
}

template <class T>
void IsGrid<T>::SetGridSize(const Grid& rhs)
{
  global = rhs.Global();
  local = rhs.Local();
  extent = rhs.Extent();
  iterators = rhs.Iterators();
  level = rhs.Level();

  delete [] grid;
  grid = new T[local.SizeTotal().Product()];
}

template <class T>
void IsGrid<T>::SetGridSize(const GlobalIndices& global_,
			    const LocalIndices& local_,
			    const SpatialExtent& extent_)
{
  global = global_;
  local = local_;
  extent = extent_;
  iterators.SetSubgrids(local);
  level = 0;
  delete [] grid;
  grid = new T[local.SizeTotal().Product()];
}

template <class T>
IsGrid<T>& IsGrid<T>::operator=(const IsGrid<T>& rhs)
{
  if (this != &rhs) {

    global = rhs.Global();
    local = rhs.Local();
    extent = rhs.Extent();
    iterators = rhs.Iterators();
    level = rhs.Level();

    delete [] grid;
    grid = new T[local.SizeTotal().Product()];

    SetGrid(rhs);

  }

  return *this;
}

template <class T>
void IsGrid<T>::SetGrid(const IsGrid<T>& rhs)
{
  for (typename IsGrid<T>::iterator iter = Iterators().CompleteGrid().Begin(); iter != Iterators().CompleteGrid().End(); ++iter)
    (*this)(*iter) = rhs.GetVal(*iter);
}

template <class T>
void IsGrid<T>::SetBoundary(const IsGrid<T>& rhs)
{
  typename IsGrid<T>::iterator iter;

  for (int i=0; i<3; ++i) {

    for (iter = Iterators().Boundary1()[i].Begin(); iter != Iterators().Boundary1()[i].End(); ++iter)
      (*this)(*iter) = rhs.GetVal(*iter);

    for (iter = Iterators().Boundary2()[i].Begin(); iter != Iterators().Boundary2()[i].End(); ++iter)
      (*this)(*iter) = rhs.GetVal(*iter);

  }
}


}

#endif /* GRID_HPP_ */
