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
 * @file   datatypes_global.hpp
 * @author Julian Iseringhausen <isering@ins.uni-bonn.de>
 * @date   Mon Jan 2 18:45:22 2012
 *
 * @brief  Stores some CommSubgrid related information in order to minimize
 *         the communication in this routine.
 *
 */

#ifndef DATATYPES_GLOBAL_HPP_
#define DATATYPES_GLOBAL_HPP_

#include <vector>

#include "comm/mpi/datatype.hpp"

namespace VMG
{

namespace MPI
{

class DatatypesGlobal
{
public:
  DatatypesGlobal() {}
  ~DatatypesGlobal() {}

  typedef std::vector<Datatype>::iterator iterator;
  typedef std::vector<Datatype>::const_iterator const_iterator;
  typedef std::vector<Datatype>::reverse_iterator reverse_iterator;
  typedef std::vector<Datatype>::const_reverse_iterator const_reverse_iterator;

  std::vector<Datatype>& Send() {return send;}
  const std::vector<Datatype>& Send() const {return send;}

  std::vector<Datatype>& Receive() {return recv;}
  const std::vector<Datatype>& Receive() const {return recv;}

private:
  std::vector<Datatype> send, recv;
};

}

}

#endif /* DATATYPES_GLOBAL_HPP_ */
