/*
 * Copyright (C) 2011-2013 Michael Pippig
 * Copyright (C) 2011 Sebastian Banert
 *
 * This file is part of ScaFaCoS.
 * 
 * ScaFaCoS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ScaFaCoS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *	
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _P2NFFT_RUN_H
#define _P2NFFT_RUN_H
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "FCSResult.h"

/** @brief
 *  @param rd
 *  @param num_particles
 *  @param positions
 *  @param charges
 *  @returns
 */
FCSResult ifcs_p2nfft_run(
    void *rd, fcs_int num_particles, fcs_int max_local_num_particles,
    fcs_float *positions, fcs_float *charges,
    fcs_float *potentials, fcs_float *field);

void ifcs_p2nfft_set_max_particle_move(void *rd, fcs_float max_particle_move);
void ifcs_p2nfft_set_resort(void *rd, fcs_int resort);
void ifcs_p2nfft_get_resort(void *rd, fcs_int *resort);
void ifcs_p2nfft_get_resort_availability(void *rd, fcs_int *availability);
void ifcs_p2nfft_get_resort_particles(void *rd, fcs_int *resort_particles);
void ifcs_p2nfft_resort_ints(void *rd, fcs_int *src, fcs_int *dst, fcs_int n, MPI_Comm comm);
void ifcs_p2nfft_resort_floats(void *rd, fcs_float *src, fcs_float *dst, fcs_int n, MPI_Comm comm);
void ifcs_p2nfft_resort_bytes(void *rd, void *src, void *dst, fcs_int n, MPI_Comm comm);

#endif
