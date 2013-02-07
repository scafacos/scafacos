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
 * @file   domain_decomposition_mpi.hpp
 * @author Julian Iseringhausen <isering@ins.uni-bonn.de>
 * @date   Mon Jun 27 12:53:50 2011
 *
 * @brief  Computes a domain decomposition which separates
 *         the finest grid equally for all processes.
 *
 */

#ifndef DOMAIN_DECOMPOSITION_MPI_HPP_
#define DOMAIN_DECOMPOSITION_MPI_HPP_

#include "comm/domain_decomposition.hpp"

namespace VMG
{

class Index;

class DomainDecompositionMPI : public DomainDecomposition
{
public:
  void Compute(Comm* comm, const Interface* interface, std::vector<GlobalIndices>& global);

private:
  bool IsActive(Comm* comm, const Index& size_global, Index& procs);
  void FineToCoarse(Comm* comm, int& begin, int& end, int levels);
};

}

#endif /* DOMAIN_DECOMPOSITION_MPI_HPP_ */
