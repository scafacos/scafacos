/*
  Copyright (C) 2011-2012 Rene Halver

  This file is part of ScaFaCoS.

  ScaFaCoS is free software: you can redistribute it and/or modify
  it under the terms of the GNU Lesser Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  ScaFaCoS is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU Lesser Public License for more details.

  You should have received a copy of the GNU Lesser Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>. 
*/



#ifndef FCS_MEMD_P_INCLUDED
#define FCS_MEMD_P_INCLUDED

#include "FCSDefinitions.h"
#include "FCSResult_p.h"
#include "FCSInterface_p.h"

/**
 * @file fcs_memd_p.h
 * @brief file containing all memd specific functions
 * @author Florian Fahrenberger <florian.fahrenberger@icp.uni-stuttgart.de>
 */

typedef struct fcs_memd_parameters_t *fcs_memd_parameters;

/**
 * @brief Function to set all memd parameters with a single call
 *


 * @param handle FCS-object that contains the parameters
 * @param max_level The maximum level of the algorithm, i.e. n_gridpoints = 2^max_level.
 * @param periodic Periodicity.
 * @param spline_degree The degree of the interpolating B-Splines.
 * @param max_iterations The maximum number of multigrid iterations.
 * @param smoothing_steps Number of pre/postsmoothing steps on each level.
 * @param gamma Cycle-number.
 * @param precision Desired precision.
 *
 * @return FCSResult-object containing the return value.
 */
FCSResult fcs_memd_setup(FCS handle,
                         fcs_float box_size,
                         fcs_float timestep,
                         fcs_int local_number_of_particles,
                         fcs_int mesh_size,
                         fcs_float lightspeed,
                         fcs_float temperature,
                         fcs_float permittivity
                         );

/**
 * @brief Set the periodicity for memd.
 *
 * @param handle FCS-object that contains the parameter
 * @param periodic Periodicity
 *
 * @return FCSResult-object containing the return state
 */
FCSResult fcs_fcs_memd_set_box_size(FCS handle, fcs_float length_x, fcs_float length_y, fcs_float length_z);
FCSResult fcs_fcs_memd_set_time_step(FCS handle, fcs_float timestep);
FCSResult fcs_fcs_memd_set_total_number_of_particles(FCS handle, fcs_int number_of_particles);
FCSResult fcs_fcs_memd_set_local_number_of_particles(FCS handle, fcs_int number_of_particles);

FCSResult fcs_fcs_memd_set_init_flag(FCS handle, fcs_int flagvalue);

FCSResult fcs_fcs_memd_set_mesh_size_1D(FCS handle, fcs_int mesh_size);
FCSResult fcs_fcs_memd_set_speed_of_light(FCS handle, fcs_float lightspeed);

FCSResult fcs_fcs_memd_set_permittivity(FCS handle, fcs_float epsilon);
FCSResult fcs_fcs_memd_set_temperature(FCS handle, fcs_float temperature);
FCSResult fcs_fcs_memd_set_bjerrum_length(FCS handle, fcs_float bjerrum);


#endif
