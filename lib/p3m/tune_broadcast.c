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
#include "tune_broadcast.h"
#include "error_estimate.h"
#include "utils.h"
#include <stdio.h>

/***************************************************/
/* TYPES AND CONSTANTS */
/***************************************************/

/***************************************************/
/* FORWARD DECLARATIONS OF INTERNAL FUNCTIONS */
/***************************************************/

static void
ifcs_p3m_tune_broadcast_params(ifcs_p3m_data_struct *d);

static void
ifcs_p3m_tune_receive_params(ifcs_p3m_data_struct *d);


/***************************************************/
/* IMPLEMENTATION */
/***************************************************/
void
ifcs_p3m_tune_broadcast_command
(ifcs_p3m_data_struct *d, fcs_int command) {
  /* First send the command */
  P3M_INFO_LOCAL(printf("       %2d: Broadcasting command %d.\n", d->comm.rank, command));
  MPI_Bcast(&command, 1, FCS_MPI_INT, 0, d->comm.mpicomm);

  /* Now send the parameters, depending on the command */
  switch (command) {
  case FINISHED:
  case TEST_RUN:
  case COMPUTE_ERROR_ESTIMATE:
    ifcs_p3m_tune_broadcast_params(d);
    return;
  case FAILED: 
    return;
  }
}

FCSResult
ifcs_p3m_tune_broadcast_slave
(ifcs_p3m_data_struct *d, fcs_int num_particles, fcs_int max_particles,
 fcs_float *positions, fcs_float *charges) {
  const char* fnc_name = "ifcs_p3m_tune_broadcast_slave";

  P3M_INFO(printf( "ifcs_p3m_tune_broadcast_slave() started...\n"));
  if (d->comm.rank == 0) {
    return fcsResult_create
      (FCS_LOGICAL_ERROR, fnc_name, 
       "Internal error: Function should not be called on master node.");
  }

  for (;;) {
    /* Receive the command */
    fcs_int command;
    P3M_INFO_LOCAL(printf("      %2d: Waiting to receive command.\n", d->comm.rank));
    MPI_Bcast(&command, 1, FCS_MPI_INT, 0, d->comm.mpicomm);
    P3M_INFO_LOCAL(printf("      %2d: Received command %d.\n", d->comm.rank, command));

    switch (command) {
    case COMPUTE_ERROR_ESTIMATE:
      ifcs_p3m_tune_receive_params(d);
      ifcs_p3m_k_space_error(d);
      break;
    case TEST_RUN:
      ifcs_p3m_tune_receive_params(d);
      //ifcs_p3m_test_run();
      break;
    case FINISHED:
      ifcs_p3m_tune_receive_params(d);
      return NULL;
    case FAILED: {
      char msg[255];
      sprintf(msg, 
              "Cannot achieve required accuracy (p3m_tolerance_field=%" 
              FCS_LMOD_FLOAT 
              "e) for given parameters.", 
              d->tolerance_field);
      return fcsResult_create(FCS_LOGICAL_ERROR, fnc_name, msg);
    }
    }
  }
  P3M_INFO(printf( "ifcs_p3m_tune_broadcast_slave() finished.\n"));
}

static void
ifcs_p3m_tune_broadcast_params(ifcs_p3m_data_struct *d) {
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

static void
ifcs_p3m_tune_receive_params(ifcs_p3m_data_struct *d) {
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
