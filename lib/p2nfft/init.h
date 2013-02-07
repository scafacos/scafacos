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

#ifndef _P2NFFT_INIT_H
#define _P2NFFT_INIT_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <mpi.h>
#include <FCSResult.h>

/** @brief Initialize all structures, parameters and arrays needed for
 *         the P2NFFT algorithm and set their default values.
 *  @param[inout] rd If @p *rd is equal to @c NULL, some memory will be
 *         allocated for the P2NFFT data structure, otherwise it is
 *         assumed that the memory is already allocated. Afterwards,
 *         in either of these cases, @p *rd is a pointer to a valid 
 *         P2NFFT data structure.
 *  @param[in] comm The MPI communicator used for the P2NFFT
 *         operation.
 */
FCSResult ifcs_p2nfft_init(void **rd, MPI_Comm comm);

/** @brief Clean up all memory allocations done by ifcs_p2nfft_init().
 *  @param[in] rd A pointer to a valid P2NFFT data structure.
 */
void ifcs_p2nfft_destroy(void *rd);

#endif
