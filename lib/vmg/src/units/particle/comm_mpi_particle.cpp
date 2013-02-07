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
 * @file   comm_mpi_particle.cpp
 * @author Julian Iseringhausen <isering@ins.uni-bonn.de>
 * @date   Tue May 8 15:27:06 2012
 *
 * @brief  Class for MPI-based particle-related communication.
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef HAVE_MPI

#include <mpi.h>
#ifdef HAVE_MARMOT
#include <enhancempicalls.h>
#include <sourceinfompicalls.h>
#endif

#include "units/particle/comm_mpi_particle.hpp"
#include "units/particle/linked_cell_list.hpp"
#include "units/particle/particle.hpp"

using namespace VMG;

void Particle::CommMPI::CommParticles(const Grid& grid, std::list<Particle>& particles)
{
  Factory& factory = MG::GetFactory();

  const vmg_int& num_particles_local = factory.GetObjectStorageVal<vmg_int>("PARTICLE_NUM_LOCAL");
  vmg_float* x = factory.GetObjectStorageArray<vmg_float>("PARTICLE_POS_ARRAY");
  vmg_float* q = factory.GetObjectStorageArray<vmg_float>("PARTICLE_CHARGE_ARRAY");

  int rank, size;
  MPI_Comm_rank(comm_global, &rank);
  MPI_Comm_size(comm_global, &size);

#ifndef VMG_ONE_SIDED
  vmg_int* receiver;
  if (!factory.TestObject("PARTICLE_RECEIVER_ARRAY"))
    new ObjectStorageArray<vmg_int>("PARTICLE_RECEIVER_ARRAY", size);
  receiver = factory.GetObjectStorageArray<vmg_int>("PARTICLE_RECEIVER_ARRAY");
#endif

  Index index;
  std::vector<int> global_extent(6*size);
  std::vector<int> send_sizes(size);
  std::vector<int> recv_sizes(size);
  std::vector<Index> begin_remote(size);
  std::vector<Index> end_remote(size);
  std::vector< std::vector<vmg_float> > send_buffer_x(size);
  std::vector< std::vector<vmg_float> > send_buffer_q(size);
  std::vector< std::vector<vmg_int> > send_buffer_ind(size);
  std::vector< std::vector<vmg_float> > recv_buffer_x(size);
  std::vector< std::vector<vmg_float> > recv_buffer_q(size);
  std::vector< std::vector<vmg_int> > recv_buffer_ind(size);

  std::memcpy(&global_extent[6*rank], grid.Global().LocalBegin().vec(), 3*sizeof(int));
  std::memcpy(&global_extent[6*rank+3], grid.Global().LocalEnd().vec(), 3*sizeof(int));

  MPI_Allgather(MPI_IN_PLACE, 6, MPI_INT, &global_extent.front(), 6, MPI_INT, comm_global);

  for (int i=0; i<size; ++i) {
    begin_remote[i] = static_cast<Index>(&global_extent[6*i]);
    end_remote[i] = static_cast<Index>(&global_extent[6*i+3]);
  }

  for (int i=0; i<num_particles_local; ++i) {
    index = static_cast<Index>((Vector(&x[3*i]) - grid.Extent().Begin()) / grid.Extent().MeshWidth());
    for (int j=0; j<size; ++j)
      if (index.IsInBounds(begin_remote[j], end_remote[j])) {
	send_buffer_x[j].push_back(x[3*i+0]);
	send_buffer_x[j].push_back(x[3*i+1]);
	send_buffer_x[j].push_back(x[3*i+2]);
	send_buffer_q[j].push_back(q[i]);
	send_buffer_ind[j].push_back(i);
	break;
      }
  }

  /*
   * Communicate which process gets how many particles
   */
  for (int i=0; i<size; ++i)
    send_sizes[i] = send_buffer_q[i].size();

  MPI_Alltoall(&send_sizes.front(), 1, MPI_INT, &recv_sizes.front(), 1, MPI_INT, comm_global);

  assert(RequestsPending() == 0);

  /*
   * Send particles
   */
  for (int i=0; i<size; ++i) {

    if (!send_buffer_q[i].empty()) {
      MPI_Isend(&send_buffer_x[i].front(), send_buffer_x[i].size(), MPI_DOUBLE, i, 0, comm_global, &Request());
      MPI_Isend(&send_buffer_q[i].front(), send_buffer_q[i].size(), MPI_DOUBLE, i, 1, comm_global, &Request());
      MPI_Isend(&send_buffer_ind[i].front(), send_buffer_ind[i].size(), MPI_INT, i, 2, comm_global, &Request());
    }

#ifndef VMG_ONE_SIDED
    receiver[i] = send_buffer_q[i].size();
#endif
  }

  /*
   * Receive particles
   */
  for (int i=0; i<size; ++i) {

    if (recv_sizes[i] > 0) {

      recv_buffer_x[i].resize(3*recv_sizes[i]);
      recv_buffer_q[i].resize(recv_sizes[i]);
      recv_buffer_ind[i].resize(recv_sizes[i]);

      MPI_Irecv(&recv_buffer_x[i].front(), 3*recv_sizes[i], MPI_DOUBLE, i, 0, comm_global, &Request());
      MPI_Irecv(&recv_buffer_q[i].front(), recv_sizes[i], MPI_DOUBLE, i, 1, comm_global, &Request());
      MPI_Irecv(&recv_buffer_ind[i].front(), recv_sizes[i], MPI_INT, i, 2, comm_global, &Request());

    }

  }

  WaitAll();

  particles.clear();

  for (int i=0; i<size; ++i)
    for (int j=0; j<recv_sizes[i]; ++j)
      particles.push_back(Particle(&recv_buffer_x[i][3*j], recv_buffer_q[i][j], 0.0, 0.0, i, recv_buffer_ind[i][j]));
}

void Particle::CommMPI::CommParticlesBack(std::list<Particle>& particles)
{
  std::list<Particle>::iterator iter;

#ifdef VMG_ONE_SIDED
  if (!win_created) {
    vmg_float* p = MG::GetFactory().GetObjectStorageArray<vmg_float>("PARTICLE_POTENTIAL_ARRAY");
    const vmg_int& num_particles_local = MG::GetFactory().GetObjectStorageVal<vmg_int>("PARTICLE_NUM_LOCAL");
    MPI_Win_create(p, num_particles_local*sizeof(vmg_float), sizeof(vmg_float), info, comm_global, &win);
    win_created = true;
  }

  MPI_Win_fence(MPI_MODE_NOPRECEDE, win);

  for (iter=particles.begin(); iter!=particles.end(); ++iter)
    MPI_Put(&iter->Pot(), 1, MPI_DOUBLE, iter->Rank(), iter->Index(), 1, MPI_DOUBLE, win);

  MPI_Win_fence(MPI_MODE_NOSTORE | MPI_MODE_NOSUCCEED, win);
#else
  int rank, size;
  MPI_Comm_rank(comm_global, &rank);
  MPI_Comm_size(comm_global, &size);

  std::vector< std::vector<vmg_float> > send_buffer_float(size);
  std::vector< std::vector<vmg_float> > recv_buffer_float(size);
  std::vector< std::vector<vmg_int> > send_buffer_index(size);
  std::vector< std::vector<vmg_int> > recv_buffer_index(size);

  vmg_int* size_receive = MG::GetFactory().GetObjectStorageArray<vmg_int>("PARTICLE_RECEIVER_ARRAY");
  vmg_float* p = MG::GetFactory().GetObjectStorageArray<vmg_float>("PARTICLE_POTENTIAL_ARRAY");
  vmg_float* f = MG::GetFactory().GetObjectStorageArray<vmg_float>("PARTICLE_FIELD_ARRAY");

  // Build send buffer
  for (iter=particles.begin(); iter!=particles.end(); ++iter) {
    send_buffer_float[iter->Rank()].push_back(iter->Pot());
    send_buffer_float[iter->Rank()].push_back(iter->Field()[0]);
    send_buffer_float[iter->Rank()].push_back(iter->Field()[1]);
    send_buffer_float[iter->Rank()].push_back(iter->Field()[2]);
    send_buffer_index[iter->Rank()].push_back(iter->Index());
  }

  // Send potentials
  for (int i=0; i<size; ++i) {
    if (!send_buffer_float[i].empty()) {
      MPI_Isend(&send_buffer_float[i].front(), send_buffer_float[i].size(), MPI_DOUBLE, i, 699+rank, comm_global, &Request());
      MPI_Isend(&send_buffer_index[i].front(), send_buffer_index[i].size(), MPI_INT, i, 32111+rank, comm_global, &Request());
    }
  }

  //Receive potentials
  for (int i=0; i<size; ++i) {
    if (size_receive[i] > 0) {
      recv_buffer_float[i].resize(4*size_receive[i]);
      recv_buffer_index[i].resize(size_receive[i]);
      MPI_Irecv(&recv_buffer_float[i].front(), 4*size_receive[i], MPI_DOUBLE, i, 699+i, comm_global, &Request());
      MPI_Irecv(&recv_buffer_index[i].front(), size_receive[i], MPI_INT, i, 32111+i, comm_global, &Request());
    }
  }

  WaitAll();

  // Add potential values
  for (int i=0; i<size; ++i)
    for (vmg_int j=0; j<size_receive[i]; ++j) {
      p[recv_buffer_index[i][j]] = recv_buffer_float[i][4*j];
      std::memcpy(&f[3*recv_buffer_index[i][j]], &recv_buffer_float[i][4*j+1], 3*sizeof(vmg_float));
    }
#endif

}

void Particle::CommMPI::CommLCListToGhosts(LinkedCellList& lc)
{
  VMG::MPI::DatatypesLocal types(lc, comm_global, false);
  std::vector<int> send_size(types.NB().size());
  vmg_int recv_size;
  std::list<Particle*>::iterator iter;
  Index ind;
  Vector offset;

  const Vector halo_length = lc.Local().HaloSize1() * lc.Extent().MeshWidth();

  lc.ClearHalo();

  for (unsigned int i=0; i<types.NB().size(); ++i)
    if (types.NB()[i].Feasible()) {

      for (int j=0; j<3; ++j)
        if ((types.Offset()[i][j] < 0 && lc.Global().LocalBegin()[j] == 0) ||
            (types.Offset()[i][j] > 0 && lc.Global().LocalEnd()[j] == lc.Global().GlobalSize()[j]))
          offset[j] = -1.0 * types.Offset()[i][j] * lc.Extent().Size()[j];
        else
          offset[j] = 0.0;

      for (ind.X() = types.NB()[i].Starts().X(); ind.X() < types.NB()[i].Starts().X()+types.NB()[i].Subsizes().X(); ++ind.X())
	for (ind.Y() = types.NB()[i].Starts().Y(); ind.Y() < types.NB()[i].Starts().Y()+types.NB()[i].Subsizes().Y(); ++ind.Y())
	  for (ind.Z() = types.NB()[i].Starts().Z(); ind.Z() < types.NB()[i].Starts().Z()+types.NB()[i].Subsizes().Z(); ++ind.Z())
	    for (iter=lc(ind).begin(); iter!=lc(ind).end(); ++iter) {

              for (int j=0; j<3; ++j)
                types.NB()[i].Buffer().push_back((*iter)->Pos()[j] + offset[j]);
              types.NB()[i].Buffer().push_back((*iter)->Charge());

              assert(lc.Extent().Begin().IsComponentwiseLessOrEqual((*iter)->Pos()));
              assert(lc.Extent().End().IsComponentwiseGreaterOrEqual((*iter)->Pos()));
              assert(lc.Extent().Begin().IsComponentwiseLessOrEqual((*iter)->Pos() + offset + halo_length));
              assert(lc.Extent().End().IsComponentwiseGreaterOrEqual((*iter)->Pos() + offset - halo_length));
	    }

      send_size[i] = types.NB()[i].Buffer().size();
      MPI_Isend(&send_size[i], 1, MPI_INT, types.NB()[i].Rank(), 2048+types.NB()[i].TagSend(), comm_global, &Request());

      if (send_size[i] > 0)
	MPI_Isend(&types.NB()[i].Buffer().front(), send_size[i], MPI_DOUBLE,
		  types.NB()[i].Rank(), 4096+types.NB()[i].TagSend(),
		  comm_global, &Request());
    }

  for (unsigned int i=0; i<types.Halo().size(); ++i)
    if (types.Halo()[i].Feasible()) {
      MPI_Recv(&recv_size, 1, MPI_INT, types.Halo()[i].Rank(), 2048+types.Halo()[i].TagReceive(), comm_global, MPI_STATUS_IGNORE);
      if (recv_size > 0) {
	types.Halo()[i].Buffer().resize(recv_size);
	MPI_Irecv(&types.Halo()[i].Buffer().front(), recv_size, MPI_DOUBLE,
		  types.Halo()[i].Rank(), 4096+types.Halo()[i].TagReceive(),
		  comm_global, &Request());
      }
    }

  WaitAll();

  for (unsigned int i=0; i<types.Halo().size(); ++i)
    for (unsigned int j=0; j<types.Halo()[i].Buffer().size(); j+=4)
      lc.AddParticleToHalo(&types.Halo()[i].Buffer()[j], types.Halo()[i].Buffer()[j+3]);
}

#endif /* HAVE_MPI */
