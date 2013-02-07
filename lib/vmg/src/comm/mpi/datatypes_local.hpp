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

#ifndef DATATYPES_LOCAL_HPP_
#define DATATYPES_LOCAL_HPP_

#include <vector>

#include "base/index.hpp"
#include "comm/mpi/datatype.hpp"

namespace VMG
{

namespace MPI
{

class DatatypesLocal
{
public:
  template <class T>
  DatatypesLocal(const T& grid, const MPI_Comm& comm, const bool& alloc)
  {
    InitDatatypesLocal(grid, comm, alloc);
  }

  ~DatatypesLocal() {}

  std::vector<Datatype>& Halo() {return _halo;}
  std::vector<Datatype>& NB() {return _nb;}

  const std::vector<Datatype>& Halo() const {return _halo;}
  const std::vector<Datatype>& NB() const {return _nb;}

  const std::vector<Index>& Offset() const {return _offset;}

private:
  template <class T>
  void InitDatatypesLocal(const T& grid, const MPI_Comm& comm, const bool& alloc_buffer);

  std::vector<Datatype> _halo, _nb;
  std::vector<Index> _offset;
};

namespace
{

inline int to_1d(const Index& i)
{
  return i.Z()+3*(i.Y()+3*i.X());
}

inline bool _is_valid(const Index& coord, const Index& dims, const Index& periods)
{
  return (periods[0] || (coord[0] >= 0 && coord[0] < dims[0])) &&
         (periods[1] || (coord[1] >= 0 && coord[1] < dims[1])) &&
         (periods[2] || (coord[2] >= 0 && coord[2] < dims[2]));
}

}

template <class T>
void DatatypesLocal::InitDatatypesLocal(const T& grid, const MPI_Comm& comm, const bool& alloc_buffer)
{
  if (comm != MPI_COMM_NULL) {

    int index, ranks[27];
    Index dims, periods, coords;
    Index sizes, subsizes, starts;
    Index offset, i;
    const LocalIndices& l = grid.Local();

    MPI_Cart_get(comm, 3, dims.vec(), periods.vec(), coords.vec());

    for (i.X()=-1; i.X()<=1; ++i.X())
      for (i.Y()=-1; i.Y()<=1; ++i.Y())
	for (i.Z()=-1; i.Z()<=1; ++i.Z())
	  if (_is_valid(coords + i, dims, periods))
	    MPI_Cart_rank(comm, (coords + i).vec(), &ranks[to_1d(i+1)]);

    sizes = l.SizeTotal();

    /* -1 0 0 */
    offset = Index(-1,0,0);
    if (_is_valid(coords + offset, dims, periods)) {
      index = to_1d(offset+1);
      subsizes = Index(l.HaloSize1().X(), l.Size().Y(), l.Size().Z());
      starts = Index(0, l.Begin().Y(), l.Begin().Z());
      _halo.push_back(VMG::MPI::Datatype(sizes, subsizes, starts, ranks[index], 0, 1, alloc_buffer));
      starts = l.Begin();
      _nb.push_back(VMG::MPI::Datatype(sizes, subsizes, starts, ranks[index], 0, 1, alloc_buffer));
      _offset.push_back(offset);
    }

    /* 1 0 0 */
    offset = Index(1,0,0);
    if (_is_valid(coords + offset, dims, periods)) {
      index = to_1d(offset+1);
      subsizes = Index(l.HaloSize2().X(), l.Size().Y(), l.Size().Z());
      starts = Index(l.End().X(), l.Begin().Y(), l.Begin().Z());
      _halo.push_back(VMG::MPI::Datatype(sizes, subsizes, starts, ranks[index], 1, 0, alloc_buffer));
      starts = Index(l.End().X()-l.HaloSize2().X(), l.Begin().Y(), l.Begin().Z());
      _nb.push_back(VMG::MPI::Datatype(sizes, subsizes, starts, ranks[index], 1, 0, alloc_buffer));
      _offset.push_back(offset);
    }

    /* 0 -1 0 */
    offset = Index(0,-1,0);
    if (_is_valid(coords + offset, dims, periods)) {
      index = to_1d(offset+1);
      subsizes = Index(l.Size().X(), l.HaloSize1().Y(), l.Size().Z());
      starts = Index(l.Begin().X(), 0, l.Begin().Z());
      _halo.push_back(VMG::MPI::Datatype(sizes, subsizes, starts, ranks[index], 2, 3, alloc_buffer));
      starts = l.Begin();
      _nb.push_back(VMG::MPI::Datatype(sizes, subsizes, starts, ranks[index], 2, 3, alloc_buffer));
      _offset.push_back(offset);
    }

    /* 0 1 0 */
    offset = Index(0,1,0);
    if (_is_valid(coords + offset, dims, periods)) {
      index = to_1d(offset+1);
      subsizes = Index(l.Size().X(), l.HaloSize2().Y(), l.Size().Z());
      starts = Index(l.Begin().X(), l.End().Y(), l.Begin().Z());
      _halo.push_back(VMG::MPI::Datatype(sizes, subsizes, starts, ranks[index], 3, 2, alloc_buffer));
      starts = Index(l.Begin().X(), l.End().Y()-l.HaloSize2().Y(), l.Begin().Z());
      _nb.push_back(VMG::MPI::Datatype(sizes, subsizes, starts, ranks[index], 3, 2, alloc_buffer));
      _offset.push_back(offset);
    }

    /* 0 0 -1 */
    offset = Index(0,0,-1);
    if (_is_valid(coords + offset, dims, periods)) {
      index = to_1d(offset+1);
      subsizes = Index(l.Size().X(), l.Size().Y(), l.HaloSize1().Z());
      starts = Index(l.Begin().X(), l.Begin().Y(), 0);
      _halo.push_back(VMG::MPI::Datatype(sizes, subsizes, starts, ranks[index], 4, 5, alloc_buffer));
      starts = l.Begin();
      _nb.push_back(VMG::MPI::Datatype(sizes, subsizes, starts, ranks[index], 4, 5, alloc_buffer));
      _offset.push_back(offset);
    }

    /* 0 0 1 */
    offset = Index(0,0,1);
    if (_is_valid(coords + offset, dims, periods)) {
      index = to_1d(offset+1);
      subsizes = Index(l.Size().X(), l.Size().Y(), l.HaloSize2().Z());
      starts = Index(l.Begin().X(), l.Begin().Y(), l.End().Z());
      _halo.push_back(VMG::MPI::Datatype(sizes, subsizes, starts, ranks[index], 5, 4, alloc_buffer));
      starts = Index(l.Begin().X(), l.Begin().Y(), l.End().Z()-l.HaloSize2().Z());
      _nb.push_back(VMG::MPI::Datatype(sizes, subsizes, starts, ranks[index], 5, 4, alloc_buffer));
      _offset.push_back(offset);
    }

    /* -1 -1 0 */
    offset = Index(-1,-1,0);
    if (_is_valid(coords + offset, dims, periods)) {
      index = to_1d(offset+1);
      subsizes = Index(l.HaloSize1().X(), l.HaloSize1().Y(), l.Size().Z());
      starts = Index(0, 0, l.Begin().Z());
      _halo.push_back(VMG::MPI::Datatype(sizes, subsizes, starts, ranks[index], 6, 7, alloc_buffer));
      starts = l.Begin();
      _nb.push_back(VMG::MPI::Datatype(sizes, subsizes, starts, ranks[index], 6, 7, alloc_buffer));
      _offset.push_back(offset);
    }

    /* -1 1 0 */
    offset = Index(-1,1,0);
    if (_is_valid(coords + offset, dims, periods)) {
      index = to_1d(offset+1);
      subsizes = Index(l.HaloSize1().X(), l.HaloSize2().Y(), l.Size().Z());
      starts = Index(0, l.End().Y(), l.Begin().Z());
      _halo.push_back(VMG::MPI::Datatype(sizes, subsizes, starts, ranks[index], 8, 9, alloc_buffer));
      starts = Index(l.Begin().X(), l.End().Y()-l.HaloSize2().Y(), l.Begin().Z());
      _nb.push_back(VMG::MPI::Datatype(sizes, subsizes, starts, ranks[index], 8, 9, alloc_buffer));
      _offset.push_back(offset);
    }

    /* 1 -1 0 */
    offset = Index(1,-1,0);
    if (_is_valid(coords + offset, dims, periods)) {
      index = to_1d(offset+1);
      subsizes = Index(l.HaloSize2().X(), l.HaloSize1().Y(), l.Size().Z());
      starts = Index(l.End().X(), 0, l.Begin().Z());
      _halo.push_back(VMG::MPI::Datatype(sizes, subsizes, starts, ranks[index], 9, 8, alloc_buffer));
      starts = Index(l.End().X()-l.HaloSize2().X(), l.Begin().Y(), l.Begin().Z());
      _nb.push_back(VMG::MPI::Datatype(sizes, subsizes, starts, ranks[index], 9, 8, alloc_buffer));
      _offset.push_back(offset);
    }

    /* 1 1 0 */
    offset = Index(1,1,0);
    if (_is_valid(coords + offset, dims, periods)) {
      index = to_1d(offset+1);
      subsizes = Index(l.HaloSize2().X(), l.HaloSize2().Y(), l.Size().Z());
      starts = Index(l.End().X(), l.End().Y(), l.Begin().Z());
      _halo.push_back(VMG::MPI::Datatype(sizes, subsizes, starts, ranks[index], 7, 6, alloc_buffer));
      starts = Index(l.End().X()-l.HaloSize2().X(), l.End().Y()-l.HaloSize2().Y(), l.Begin().Z());
      _nb.push_back(VMG::MPI::Datatype(sizes, subsizes, starts, ranks[index], 7, 6, alloc_buffer));
      _offset.push_back(offset);
    }

    /* -1 0 -1 */
    offset = Index(-1,0,-1);
    if (_is_valid(coords + offset, dims, periods)) {
      index = to_1d(offset+1);
      subsizes = Index(l.HaloSize1().X(), l.Size().Y(), l.HaloSize1().Z());
      starts = Index(0, l.Begin().Y(), 0);
      _halo.push_back(VMG::MPI::Datatype(sizes, subsizes, starts, ranks[index], 10, 11, alloc_buffer));
      starts = l.Begin();
      _nb.push_back(VMG::MPI::Datatype(sizes, subsizes, starts, ranks[index], 10, 11, alloc_buffer));
      _offset.push_back(offset);
    }

    /* -1 0 1 */
    offset = Index(-1,0,1);
    if (_is_valid(coords + offset, dims, periods)) {
      index = to_1d(offset+1);
      subsizes = Index(l.HaloSize1().X(), l.Size().Y(), l.HaloSize2().Z());
      starts = Index(0, l.Begin().Y(), l.End().Z());
      _halo.push_back(VMG::MPI::Datatype(sizes, subsizes, starts, ranks[index], 12, 13, alloc_buffer));
      starts = Index(l.Begin().X(), l.Begin().Y(), l.End().Z()-l.HaloSize2().Z());
      _nb.push_back(VMG::MPI::Datatype(sizes, subsizes, starts, ranks[index], 12, 13, alloc_buffer));
      _offset.push_back(offset);
    }

    /* 1 0 -1 */
    offset = Index(1,0,-1);
    if (_is_valid(coords + offset, dims, periods)) {
      index = to_1d(offset+1);
      subsizes = Index(l.HaloSize2().X(), l.Size().Y(), l.HaloSize1().Z());
      starts = Index(l.End().X(), l.Begin().Y(), 0);
      _halo.push_back(VMG::MPI::Datatype(sizes, subsizes, starts, ranks[index], 13, 12, alloc_buffer));
      starts = Index(l.End().X()-l.HaloSize2().X(), l.Begin().Y(), l.Begin().Z());
      _nb.push_back(VMG::MPI::Datatype(sizes, subsizes, starts, ranks[index], 13, 12, alloc_buffer));
      _offset.push_back(offset);
    }

    /* 1 0 1 */
    offset = Index(1,0,1);
    if (_is_valid(coords + offset, dims, periods)) {
      index = to_1d(offset+1);
      subsizes = Index(l.HaloSize2().X(), l.Size().Y(), l.HaloSize2().Z());
      starts = Index(l.End().X(), l.Begin().Y(), l.End().Z());
      _halo.push_back(VMG::MPI::Datatype(sizes, subsizes, starts, ranks[index], 11, 10, alloc_buffer));
      starts = Index(l.End().X()-l.HaloSize2().X(), l.Begin().Y(), l.End().Z()-l.HaloSize2().Z());
      _nb.push_back(VMG::MPI::Datatype(sizes, subsizes, starts, ranks[index], 11, 10, alloc_buffer));
      _offset.push_back(offset);
    }

    /* 0 -1 -1 */
    offset = Index(0,-1,-1);
     if (_is_valid(coords + offset, dims, periods)) {
       index = to_1d(offset+1);
       subsizes = Index(l.Size().X(), l.HaloSize1().Y(), l.HaloSize1().Z());
       starts = Index(l.Begin().X(), 0, 0);
       _halo.push_back(VMG::MPI::Datatype(sizes, subsizes, starts, ranks[index], 14, 15, alloc_buffer));
       starts = l.Begin();
       _nb.push_back(VMG::MPI::Datatype(sizes, subsizes, starts, ranks[index], 14, 15, alloc_buffer));
       _offset.push_back(offset);
     }

    /* 0 -1 1 */
     offset = Index(0,-1,1);
     if (_is_valid(coords + offset, dims, periods)) {
       index = to_1d(offset+1);
       subsizes = Index(l.Size().X(), l.HaloSize1().Y(), l.HaloSize2().Z());
       starts = Index(l.Begin().X(), 0, l.End().Z());
       _halo.push_back(VMG::MPI::Datatype(sizes, subsizes, starts, ranks[index], 16, 17, alloc_buffer));
       starts = Index(l.Begin().X(), l.Begin().Y(), l.End().Z()-l.HaloSize2().Z());
       _nb.push_back(VMG::MPI::Datatype(sizes, subsizes, starts, ranks[index], 16, 17, alloc_buffer));
       _offset.push_back(offset);
     }

    /* 0 1 -1 */
     offset = Index(0,1,-1);
     if (_is_valid(coords + offset, dims, periods)) {
       index = to_1d(offset+1);
       subsizes = Index(l.Size().X(), l.HaloSize2().Y(), l.HaloSize1().Z());
       starts = Index(l.Begin().X(), l.End().Y(), 0);
       _halo.push_back(VMG::MPI::Datatype(sizes, subsizes, starts, ranks[index], 17, 16, alloc_buffer));
       starts = Index(l.Begin().X(), l.End().Y()-l.HaloSize2().Y(), l.Begin().Z());
       _nb.push_back(VMG::MPI::Datatype(sizes, subsizes, starts, ranks[index], 17, 16, alloc_buffer));
       _offset.push_back(offset);
     }

    /* 0 1 1 */
     offset = Index(0,1,1);
     if (_is_valid(coords + offset, dims, periods)) {
       index = to_1d(offset+1);
       subsizes = Index(l.Size().X(), l.HaloSize2().Y(), l.HaloSize2().Z());
       starts = Index(l.Begin().X(), l.End().Y(), l.End().Z());
       _halo.push_back(VMG::MPI::Datatype(sizes, subsizes, starts, ranks[index], 15, 14, alloc_buffer));
       starts = Index(l.Begin().X(), l.End().Y()-l.HaloSize2().Y(), l.End().Z()-l.HaloSize2().Z());
       _nb.push_back(VMG::MPI::Datatype(sizes, subsizes, starts, ranks[index], 15, 14, alloc_buffer));
       _offset.push_back(offset);
     }

    /* -1 -1 -1 */
     offset = Index(-1,-1,-1);
     if (_is_valid(coords + offset, dims, periods)) {
       index = to_1d(offset+1);
       subsizes = l.HaloSize1();
       starts = 0;
       _halo.push_back(VMG::MPI::Datatype(sizes, subsizes, starts, ranks[index], 18, 19, alloc_buffer));
       starts = l.Begin();
       _nb.push_back(VMG::MPI::Datatype(sizes, subsizes, starts, ranks[index], 18, 19, alloc_buffer));
       _offset.push_back(offset);
     }

    /* -1 -1 1 */
     offset = Index(-1,-1,1);
     if (_is_valid(coords + offset, dims, periods)) {
       index = to_1d(offset+1);
       subsizes = Index(l.HaloSize1().X(), l.HaloSize1().Y(), l.HaloSize2().Z());
       starts = Index(0, 0, l.End().Z());;
       _halo.push_back(VMG::MPI::Datatype(sizes, subsizes, starts, ranks[index], 20, 21, alloc_buffer));
       starts = Index(l.Begin().X(), l.Begin().Y(), l.End().Z()-l.HaloSize2().Z());
       _nb.push_back(VMG::MPI::Datatype(sizes, subsizes, starts, ranks[index], 20, 21, alloc_buffer));
       _offset.push_back(offset);
     }

    /* -1 1 -1 */
     offset = Index(-1,1,-1);
     if (_is_valid(coords + offset, dims, periods)) {
       index = to_1d(offset+1);
       subsizes = Index(l.HaloSize1().X(), l.HaloSize2().Y(), l.HaloSize1().Z());
       starts = Index(0, l.End().Y(), 0);
       _halo.push_back(VMG::MPI::Datatype(sizes, subsizes, starts, ranks[index], 22, 23, alloc_buffer));
       starts = Index(l.Begin().X(), l.End().Y()-l.HaloSize2().Y(), l.Begin().Z());
       _nb.push_back(VMG::MPI::Datatype(sizes, subsizes, starts, ranks[index], 22, 23, alloc_buffer));
       _offset.push_back(offset);
     }

    /* 1 -1 -1 */
     offset = Index(1,-1,-1);
     if (_is_valid(coords + offset, dims, periods)) {
       index = to_1d(offset+1);
       subsizes = Index(l.HaloSize2().X(), l.HaloSize1().Y(), l.HaloSize1().Z());
       starts = Index(l.End().X(), 0, 0);
       _halo.push_back(VMG::MPI::Datatype(sizes, subsizes, starts, ranks[index], 24, 25, alloc_buffer));
       starts = Index(l.End().X()-l.HaloSize2().X(), l.Begin().Y(), l.Begin().Z());
       _nb.push_back(VMG::MPI::Datatype(sizes, subsizes, starts, ranks[index], 24, 25, alloc_buffer));
       _offset.push_back(offset);
     }

    /* -1 1 1 */
     offset = Index(-1,1,1);
     if (_is_valid(coords + offset, dims, periods)) {
       index = to_1d(offset+1);
       subsizes = Index(l.HaloSize1().X(), l.HaloSize2().Y(), l.HaloSize2().Z());
       starts = Index(0, l.End().Y(), l.End().Z());
       _halo.push_back(VMG::MPI::Datatype(sizes, subsizes, starts, ranks[index], 25, 24, alloc_buffer));
       starts = Index(l.Begin().X(), l.End().Y()-l.HaloSize2().Y(), l.End().Z()-l.HaloSize2().Z());
       _nb.push_back(VMG::MPI::Datatype(sizes, subsizes, starts, ranks[index], 25, 24, alloc_buffer));
       _offset.push_back(offset);
     }

    /* 1 -1 1 */
     offset = Index(1,-1,1);
     if (_is_valid(coords + offset, dims, periods)) {
       index = to_1d(offset+1);
       subsizes = Index(l.HaloSize2().X(), l.HaloSize1().Y(), l.HaloSize2().Z());
       starts = Index(l.End().X(), 0, l.End().Z());
       _halo.push_back(VMG::MPI::Datatype(sizes, subsizes, starts, ranks[index], 23, 22, alloc_buffer));
       starts = Index(l.End().X()-l.HaloSize2().X(), l.Begin().Y(), l.End().Z()-l.HaloSize2().Z());
       _nb.push_back(VMG::MPI::Datatype(sizes, subsizes, starts, ranks[index], 23, 22, alloc_buffer));
       _offset.push_back(offset);
     }

    /* 1 1 -1 */
     offset = Index(1,1,-1);
     if (_is_valid(coords + offset, dims, periods)) {
       index = to_1d(offset+1);
       subsizes = Index(l.HaloSize2().X(), l.HaloSize2().Y(), l.HaloSize1().Z());
       starts = Index(l.End().X(), l.End().Y(), 0);
       _halo.push_back(VMG::MPI::Datatype(sizes, subsizes, starts, ranks[index], 21, 20, alloc_buffer));
       starts = Index(l.End().X()-l.HaloSize2().X(), l.End().Y()-l.HaloSize2().Y(), l.Begin().Z());
       _nb.push_back(VMG::MPI::Datatype(sizes, subsizes, starts, ranks[index], 21, 20, alloc_buffer));
       _offset.push_back(offset);
     }

    /* 1 1 1 */
     offset = Index(1,1,1);
     if (_is_valid(coords + offset, dims, periods)) {
       index = to_1d(offset+1);
       subsizes = l.HaloSize2();
       starts = l.End();
       _halo.push_back(VMG::MPI::Datatype(sizes, subsizes, starts, ranks[index], 19, 18, alloc_buffer));
       starts = l.End()-l.HaloSize2();
       _nb.push_back(VMG::MPI::Datatype(sizes, subsizes, starts, ranks[index], 19, 18, alloc_buffer));
       _offset.push_back(offset);
     }

  }
}

}

}

#endif /* DATATYPES_LOCAL_HPP_ */
