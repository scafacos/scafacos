/*
  Copyright (C) 2010,2011 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 
    Max-Planck-Institute for Polymer Research, Theory Group
  
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

#ifndef _P3M_ERROR_ESTIMATE_H
#define _P3M_ERROR_ESTIMATE_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <mpi.h>
#include "types.h"

/** Calculates the rms error estimate in the force (as described in
    the book of Hockney and Eastwood (Eqn. 8.23) for a system of N
    randomly distributed particles.
    \param N        	number of charged particles in the system
    \param sum_q2   	sum of square of charges in the system
    \param box_l    	system size
    \param r_cut 	cutoff
    \param grid     	number of grid points in the different directions
    \param alpha	ewald splitting parameter
    \param cao		charge assignment order
    \param err		(out) total error estimate
    \param rs_error	(out) real-space error estimate
    \param ks_error	(out) k-space error estimate
*/
void
ifcs_p3m_compute_error_estimate(fcs_int N, 
				      fcs_float sum_q2, 
				      fcs_float box_l[3], 
				      fcs_float r_cut, 
				      fcs_int grid[3], 
				      fcs_float alpha, 
				      fcs_int cao,
				      fcs_float *err,
				      fcs_float *rs_error,
				      fcs_float *ks_error,
				      MPI_Comm comm);

/** Determines a value for alpha that achieves the wanted_error, if at
    all possible. Also returns the achieved errors with these
    parameters. Check whether wanted_error > achieved_error to see
    whether the required error can actually be met.
    \param N        	number of charged particles in the system
    \param sum_q2   	sum of square of charges in the system
    \param box_l    	system size
    \param r_cut 	cutoff
    \param grid     	number of grid points in the different directions
    \param cao		charge assignment order
    \param wanted_err	the wanted error
    \param alpha	(out) ewald splitting parameter
    \param err		(out) total error estimate
    \param rs_error	(out) real-space error estimate
    \param ks_error	(out) k-space error estimate
    \param comm		MPI communicator
*/
void
ifcs_p3m_determine_good_alpha(fcs_int N, 
                              fcs_float sum_q2, 
                              fcs_float box_l[3], 
                              fcs_float r_cut, 
                              fcs_int grid[3], 
                              fcs_int cao,
                              fcs_float wanted_error,
                              fcs_float *alpha, 
                              fcs_float *err,
                              fcs_float *rs_error,
                              fcs_float *ks_error,
                              MPI_Comm comm);

/** Broadcast the final parameters at the end of tuning. Also ends the
    slave loops.
    \param r_cut 	(out) cutoff
    \param grid     	(out) number of grid points in the different directions
    \param alpha	(out) ewald splitting parameter
    \param cao		(out) charge assignment order
    \param error	(out) total error estimate
    \param rs_error	(out) real-space error estimate
    \param ks_error	(out) k-space error estimate
    \param comm		MPI communicator
*/
void
ifcs_p3m_param_broadcast(fcs_float r_cut, fcs_int grid[3], 
                         fcs_float alpha, fcs_int cao,
                         fcs_float error, 
                         fcs_float rs_error, fcs_float ks_error,
                         MPI_Comm comm);


/** Run the loop on the slaves that waits for requests from one of the
    other functions. The loop will after when
    ifcs_p3m_param_broadcast has been called on the
    master. Afterwards, r_cut, grid, alpha, cao, error, rs_error and
    ks_error will be set to their respective values on all nodes.
    \param N        	number of charged particles in the system
    \param sum_q2   	sum of square of charges in the system
    \param box_l    	system size
    \param r_cut 	(out) cutoff
    \param grid     	(out) number of grid points in the different directions
    \param alpha	(out) ewald splitting parameter
    \param cao		(out) charge assignment order
    \param error	(out) total error estimate
    \param rs_error	(out) real-space error estimate
    \param ks_error	(out) k-space error estimate
    \param comm		MPI communicator
*/
void
ifcs_p3m_param_broadcast_slave(fcs_int N, fcs_float sum_q2,
                               fcs_float box_l[3],
                               fcs_float *r_cut, fcs_int grid[3],
                               fcs_float *alpha, fcs_int *cao,
                               fcs_float *error, 
                               fcs_float *rs_error, fcs_float *ks_error,
                               MPI_Comm comm);

#endif
