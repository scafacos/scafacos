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

#ifndef DATATYPE_HPP_
#define DATATYPE_HPP_

#include <vector>

#include "base/index.hpp"
#include "grid/grid.hpp"

namespace VMG
{

namespace MPI
{

class Datatype
{
public:
  Datatype() :
    _sizes(0),
    _subsizes(0),
    _starts(0),
    _rank(-1),
    _tag_send(0),
    _tag_recv(0),
    _type(MPI_DATATYPE_NULL),
    _alloc_buffer(false)
  {}

  Datatype(Index sizes, Index subsizes, Index starts, const int& rank,
	   const int& tag_send, const int& tag_receive,
           const bool& alloc_buffer) :
    _sizes(sizes),
    _subsizes(subsizes),
    _starts(starts),
    _rank(rank),
    _tag_send(tag_send),
    _tag_recv(tag_receive),
    _alloc_buffer(alloc_buffer)
  {
    InitDatatype();
  }

  Datatype(const GridIteratorSet& bounds, const Grid& grid, const int& rank,
	   const int& tag_send, const int& tag_receive,
           const bool& alloc_buffer) :
    _sizes(grid.Local().SizeTotal()),
    _subsizes(bounds.Begin().GetEnd() - bounds.Begin().GetBegin()),
    _starts(bounds.Begin().GetBegin()),
    _rank(rank),
    _tag_send(tag_send),
    _tag_recv(tag_receive),
    _alloc_buffer(alloc_buffer)
  {
    InitDatatype();
  }

  Datatype(const Datatype& other) :
    _sizes(other._sizes),
    _subsizes(other._subsizes),
    _starts(other._starts),
    _rank(other._rank),
    _tag_send(other._tag_send),
    _tag_recv(other._tag_recv),
    _alloc_buffer(other._alloc_buffer)
  {
    InitDatatype();
  }

  ~Datatype()
  {
    if (_type != MPI_DATATYPE_NULL)
      MPI_Type_free(&_type);
  }

  void Set(const GridIteratorSet& bounds, const Grid& grid, const int& rank,
	   const int& tag_send, const int& tag_receive);

  void Set(const Index& sizes, const Index& subsizes, const Index& starts, const int& rank,
	   const int& tag_send, const int& tag_receive);

  const Index& Sizes() const {return _sizes;}
  const Index& Subsizes() const {return _subsizes;}
  const Index& Starts() const {return _starts;}

  const int& Rank() const {return _rank;}
  const int& TagSend() const {return _tag_send;}
  const int& TagReceive() const {return _tag_recv;}

  const MPI_Datatype& Type() const {return _type;}

  std::vector<vmg_float>& Buffer() {return _buffer;}
  const std::vector<vmg_float>& Buffer() const {return _buffer;}

  void Send(Grid& grid, const int& tag, const MPI_Comm& comm) const;
  void Isend(Grid& grid, const int& tag, const MPI_Comm& comm, MPI_Request& request) const;
  void Recv(Grid& grid, const int& tag, const MPI_Comm& comm) const;
  void Irecv(Grid& grid, const int& tag, const MPI_Comm& comm, MPI_Request& request) const;

  void SendBuffered(const Grid& grid, const int& tag, const MPI_Comm& comm);
  void IsendBuffered(const Grid& grid, const int& tag, const MPI_Comm& comm, MPI_Request& request);
  void RecvBuffered(const int& tag, const MPI_Comm& comm);
  void IrecvBuffered(const int& tag, const MPI_Comm& comm, MPI_Request& request);

  void GridReplace(Grid& grid) const;
  void GridSum(Grid& grid) const;

  bool Feasible() const
  {
    return _sizes.Product() > 0 && _subsizes.Product() > 0;
  }

private:
  void InitDatatype();

  Index _sizes, _subsizes, _starts;
  int _rank, _tag_send, _tag_recv;
  MPI_Datatype _type;
  std::vector<vmg_float> _buffer;
  bool _alloc_buffer;
};

}

}

#endif /* DATATYPE_HPP_ */
