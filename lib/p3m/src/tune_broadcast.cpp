/*
  Copyright (C) 2013 Olaf Lenz
  
  This file is part of ScaFaCoS.
  
  ScaFaCoS is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  ScaFaCoS is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>. 
*/
#include "tune_broadcast.hpp"

#include "utils.hpp"
#include <cstdio>
#include <stdexcept>
#include "error_estimate.hpp"
#include "timing.hpp"

namespace P3M {
  /***************************************************/
  /* TYPES AND CONSTANTS */
  /***************************************************/

  /***************************************************/
  /* FORWARD DECLARATIONS OF INTERNAL FUNCTIONS */
  /***************************************************/

  void
  tune_broadcast_params(data_struct *d);

  void
  tune_receive_params(data_struct *d);


  /***************************************************/
  /* IMPLEMENTATION */
  /***************************************************/
  void
  tune_broadcast_command(data_struct *d, fcs_int command) {
    /* First send the command */
    P3M_DEBUG_LOCAL(printf("       %2d: Broadcasting command %d.\n", \
                           d->comm.rank, command));
    MPI_Bcast(&command, 1, FCS_MPI_INT, 0, d->comm.mpicomm);

    /* Now send the parameters, depending on the command */
    switch (command) {
    case CMD_FINISHED:
    case CMD_TIMING:
    case CMD_COMPUTE_ERROR_ESTIMATE:
      tune_broadcast_params(d);
      return;
    case CMD_FAILED: 
      return;
    }
  }

  void
  tune_broadcast_slave(data_struct *d, fcs_int num_particles,
                       fcs_float *positions, fcs_float *charges) {
    P3M_DEBUG(printf( "tune_broadcast_slave() started...\n"));
    if (d->comm.rank == 0)
      throw std::logic_error("Internal error: tune_broadcast_slave " 
                             "should not be called on master!");

    for (;;) {
      /* Receive the command */
      fcs_int command;
      P3M_DEBUG_LOCAL(printf("      %2d: Waiting to receive command.\n", \
                             d->comm.rank));
      MPI_Bcast(&command, 1, FCS_MPI_INT, 0, d->comm.mpicomm);
      P3M_DEBUG_LOCAL(printf("      %2d: Received command %d.\n", \
                             d->comm.rank, command));

      switch (command) {
      case CMD_COMPUTE_ERROR_ESTIMATE:
        tune_receive_params(d);
        k_space_error(d);
        break;
      case CMD_TIMING:
        tune_receive_params(d);
        timing(d, num_particles, positions, charges);
        break;
      case CMD_FINISHED:
        tune_receive_params(d);
        return;
      case CMD_FAILED: {
        char msg[255];
        sprintf(msg, 
                "Cannot achieve required accuracy (p3m_tolerance_field=" FFLOATE \
                ") for given parameters.", 
                d->tolerance_field);
        throw std::logic_error(msg);
      }
      }
    }
    P3M_DEBUG(printf( "tune_broadcast_slave() finished.\n"));
  }

  void
  tune_broadcast_params(data_struct *d) {
    // broadcast parameters and mark as final
    fcs_int int_buffer[4];
    fcs_float float_buffer[5];
  
    // pack int data
    int_buffer[0] = d->grid[0];
    int_buffer[1] = d->grid[1];
    int_buffer[2] = d->grid[2];
    int_buffer[3] = d->cao;
    MPI_Bcast(int_buffer, 4, FCS_MPI_INT, 0, d->comm.mpicomm);
  
    // pack float data
    float_buffer[0] = d->alpha;
    float_buffer[1] = d->r_cut;
    float_buffer[2] = d->error;
    float_buffer[3] = d->rs_error;
    float_buffer[4] = d->ks_error;
    MPI_Bcast(float_buffer, 5, FCS_MPI_FLOAT, 0, d->comm.mpicomm);
  }

  void
  tune_receive_params(data_struct *d) {
    fcs_int int_buffer[4];
    fcs_float float_buffer[5];

    MPI_Bcast(int_buffer, 4, FCS_MPI_INT, 0, d->comm.mpicomm);
    d->grid[0] = int_buffer[0];
    d->grid[1] = int_buffer[1];
    d->grid[2] = int_buffer[2];
    d->cao = int_buffer[3];
  
    // unpack float data
    MPI_Bcast(float_buffer, 5, FCS_MPI_FLOAT, 0,  d->comm.mpicomm);
    d->alpha = float_buffer[0];
    d->r_cut = float_buffer[1];
    d->error = float_buffer[2];
    d->rs_error = float_buffer[3];
    d->ks_error = float_buffer[4];
  }
}
