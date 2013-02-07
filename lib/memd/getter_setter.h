/*
 Copyright (C) 2010/2011/2012 Florian Fahrenberger
 
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

#ifndef _MEMD_GETTER_SETTER_H
#define _MEMD_GETTER_SETTER_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "FCSResult.h"

FCSResult memd_set_init_flag(void* rawdata, fcs_int flagvalue);

FCSResult memd_set_speed_of_light(void* rawdata, fcs_float lightspeed);

FCSResult memd_set_time_step(void* rawdata, fcs_float timestep);

FCSResult memd_set_permittivity(void* rawdata, fcs_float epsilon);

FCSResult memd_set_temperature(void* rawdata, fcs_float temperature);

FCSResult memd_set_bjerrum_length(void* rawdata, fcs_float bjerrum);

FCSResult memd_set_total_number_of_particles(void* rawdata, fcs_int number_of_particles);

FCSResult memd_set_local_number_of_particles(void* rawdata, fcs_int number_of_particles);

FCSResult memd_set_box_size(void* rawdata, fcs_float length_x, fcs_float length_y, fcs_float length_z);

FCSResult memd_set_mesh_size_1D(void* rawdata, fcs_int mesh_size);

fcs_float memd_get_time_step(void* rawdata);

fcs_int memd_needs_retuning(void* rawdata, fcs_int local_particles, fcs_float* positions, fcs_float* charges);

FCSResult memd_tune_method(void* rawdata, fcs_int local_particles, fcs_float* positions, fcs_float* charges);

#endif