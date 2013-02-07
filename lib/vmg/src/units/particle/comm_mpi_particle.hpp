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

#ifndef COMM_MPI_PARTICLE_HPP_
#define COMM_MPI_PARTICLE_HPP_

#include <list>

#include "comm/comm_mpi.hpp"
#include "units/particle/particle.hpp"

namespace VMG
{

namespace Particle
{

class LinkedCellList;

class CommMPI : public VMG::CommMPI
{
public:
  CommMPI(const Boundary& boundary, DomainDecomposition* domain_dec, const MPI_Comm& mpi_comm) :
    VMG::CommMPI(boundary, domain_dec, mpi_comm)
  {}

  CommMPI(const Boundary& boundary, DomainDecomposition* domain_dec) :
    VMG::CommMPI(boundary, domain_dec)
  {}

  virtual ~CommMPI() {}

  void CommParticles(const Grid& grid, std::list<Particle>& particles);
  void CommParticlesBack(std::list<Particle>& particles);
  void CommLCListToGhosts(LinkedCellList& lc);
};

}

}

#endif /* COMM_MPI_PARTICLE_HPP_ */
