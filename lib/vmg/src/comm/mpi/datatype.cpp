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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifndef HAVE_MPI
#error MPI is needed to compile VMG::MPI::Datatype
#endif

#include <mpi.h>
#ifdef HAVE_MARMOT
#include <enhancempicalls.h>
#include <sourceinfompicalls.h>
#endif

#include <cstring>

#include "comm/mpi/datatype.hpp"

using namespace VMG;

void VMG::MPI::Datatype::Send(Grid& grid, const int& tag, const MPI_Comm& comm) const
{
  if (Feasible())
    MPI_Send(&grid(0), 1, _type, _rank, _tag_send+tag, comm);
}

void VMG::MPI::Datatype::Isend(Grid& grid, const int& tag, const MPI_Comm& comm, MPI_Request& request) const
{
  if (Feasible())
    MPI_Isend(&grid(0), 1, _type, _rank, _tag_send+tag, comm, &request);
}

void VMG::MPI::Datatype::Recv(Grid& grid, const int& tag, const MPI_Comm& comm) const
{
  if (Feasible())
    MPI_Recv(&grid(0), 1 ,_type, _rank, _tag_recv+tag, comm, MPI_STATUS_IGNORE);
}

void VMG::MPI::Datatype::Irecv(Grid& grid, const int& tag, const MPI_Comm& comm, MPI_Request& request) const
{
  if (Feasible())
    MPI_Irecv(&grid(0), 1, _type, _rank, _tag_recv+tag, comm, &request);
}

void VMG::MPI::Datatype::SendBuffered(const Grid& grid, const int& tag, const MPI_Comm& comm)
{
  if (Feasible()) {

    Index i;
    int c = 0;
    const Index end = _starts + _subsizes;
    const size_t memcpy_size = _subsizes.Z() * sizeof(vmg_float);

    for (i.X()=_starts.X(); i.X()<end.X();  ++i.X())
      for (i.Y()=_starts.Y(); i.Y()<end.Y(); ++i.Y()) {
	std::memcpy(&_buffer[c], &grid.GetVal(i.X(), i.Y(), _starts.Z()), memcpy_size);
	c += _subsizes.Z();
      }

    MPI_Send(&_buffer.front(), _buffer.size(), MPI_DOUBLE, _rank, _tag_send+tag, comm);
  }
}

void VMG::MPI::Datatype::IsendBuffered(const Grid& grid, const int& tag, const MPI_Comm& comm, MPI_Request& request)
{
  if (Feasible()) {

    Index i;
    unsigned int c = 0;
    const Index end = _starts + _subsizes;
    const size_t memcpy_size = _subsizes.Z() * sizeof(vmg_float);

    for (i.X()=_starts.X(); i.X()<end.X();  ++i.X())
      for (i.Y()=_starts.Y(); i.Y()<end.Y(); ++i.Y()) {
	std::memcpy(&_buffer[c], &grid.GetVal(i.X(), i.Y(), _starts.Z()), memcpy_size);
	c += _subsizes.Z();
      }

    assert(c == _buffer.size());

    MPI_Isend(&_buffer.front(), _buffer.size(), MPI_DOUBLE, _rank, _tag_send+tag, comm, &request);
  }
}

void VMG::MPI::Datatype::RecvBuffered(const int& tag, const MPI_Comm& comm)
{
  if (Feasible())
    MPI_Recv(&_buffer.front(), _buffer.size(), MPI_DOUBLE, _rank, _tag_recv+tag, comm, MPI_STATUS_IGNORE);
}

void VMG::MPI::Datatype::IrecvBuffered(const int& tag, const MPI_Comm& comm, MPI_Request& request)
{
  if (Feasible())
    MPI_Irecv(&_buffer.front(), _buffer.size(), MPI_DOUBLE, _rank, _tag_recv+tag, comm, &request);
}

void VMG::MPI::Datatype::GridReplace(Grid& grid) const
{
  if (Feasible()) {

    Index i;
    unsigned int c = 0;
    const Index end = _starts + _subsizes;
    const size_t memcpy_size = _subsizes.Z() * sizeof(vmg_float);

    for (i.X()=_starts.X(); i.X()<end.X(); ++i.X())
      for (i.Y()=_starts.Y(); i.Y()<end.Y(); ++i.Y()) {
	std::memcpy(&grid(i.X(), i.Y(), _starts.Z()), &_buffer[c], memcpy_size);
	c += _subsizes.Z();
      }

    assert(c == _buffer.size());
  }
}

void VMG::MPI::Datatype::GridSum(Grid& grid) const
{
  if (Feasible()) {

    Index i;
    const Index end = _starts + _subsizes;
    std::vector<vmg_float>::const_iterator iter = _buffer.begin();

    for (i.X()=_starts.X(); i.X()<end.X(); ++i.X())
      for (i.Y()=_starts.Y(); i.Y()<end.Y(); ++i.Y())
	for (i.Z()=_starts.Z(); i.Z()<end.Z(); ++i.Z())
	  grid(i) += *iter++;

    assert(iter == _buffer.end());
  }
}

void VMG::MPI::Datatype::Set(const GridIteratorSet& bounds, const Grid& grid, const int& rank,
			     const int& tag_send, const int& tag_receive)
{
  _sizes = grid.Local().SizeTotal();
  _subsizes = bounds.Begin().GetEnd() - bounds.Begin().GetBegin();
  _starts = bounds.Begin().GetBegin();
  _rank = rank;
  _tag_send = tag_send;
  _tag_recv = tag_receive;

  if (_type != MPI_DATATYPE_NULL)
    MPI_Type_free(&_type);

  InitDatatype();
}

void VMG::MPI::Datatype::Set(const Index& sizes, const Index& subsizes, const Index& starts, const int& rank,
			     const int& tag_send, const int& tag_receive)
{
  _sizes = sizes;
  _subsizes = subsizes;
  _starts = starts;
  _rank = rank;
  _tag_send = tag_send;
  _tag_recv = tag_receive;

  if (_type != MPI_DATATYPE_NULL)
    MPI_Type_free(&_type);

  InitDatatype();
}

void VMG::MPI::Datatype::InitDatatype()
{
  if (Feasible()) {
    MPI_Type_create_subarray(3, _sizes.vec(), _subsizes.vec(), _starts.vec(), MPI_ORDER_C, MPI_DOUBLE, &_type);
    MPI_Type_commit(&_type);
    if (_alloc_buffer)
      _buffer.resize(_subsizes.Product());
  }else  {
    _type = MPI_DATATYPE_NULL;
  }
}
