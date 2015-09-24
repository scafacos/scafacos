/*
 * Copyright (C) 2015 Michael Pippig
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

#ifndef _P2NFFT_NEARFIELD_P_H
#define _P2NFFT_NEARFIELD_P_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

fcs_float ifcs_p2nfft_compute_self_potential(
    const void* param);
fcs_float ifcs_p2nfft_compute_near_potential(
    const void* param, fcs_float dist);
fcs_float ifcs_p2nfft_compute_near_field(
    const void* param, fcs_float dist);
void ifcs_p2nfft_compute_near_field_and_potential(
    const void* param, fcs_float dist, 
    fcs_float *potential, fcs_float *field);

#endif
