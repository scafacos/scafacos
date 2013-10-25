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

#include <stdio.h>

#include "getter_setter.h"
#include "communication.h"
#include "data_types.h"
#include "helper_functions.h"
#include "FCSCommon.h"
#include <math.h>

FCSResult fcs_memd_set_init_flag(void* rawdata, fcs_int flagvalue)
{
    const char* fnc_name = "fcs_memd_set_init_flag";

    memd_struct* memd = (memd_struct*) rawdata;
    memd->init_flag = flagvalue;
    return NULL;
}

FCSResult fcs_memd_set_speed_of_light(void* rawdata, fcs_float lightspeed)
{
    const char* fnc_name = "fcs_memd_set_speed_of_light";
    memd_struct* memd = (memd_struct*) rawdata;

    fcs_float f_mass = 1.0/lightspeed;
    fcs_float invsqrt_f_mass = lightspeed * lightspeed;

    memd->parameters.f_mass = f_mass;
    memd->parameters.invsqrt_f_mass = invsqrt_f_mass;
    
    return NULL;
}

FCSResult fcs_memd_set_time_step(void* rawdata, fcs_float timestep)
{
    const char* fnc_name = "fcs_memd_set_time_step";
    memd_struct* memd = (memd_struct*) rawdata;

    memd->parameters.time_step = timestep;
    return NULL;
}

FCSResult fcs_memd_set_permittivity(void* rawdata, fcs_float epsilon)
{
    const char* fnc_name = "fcs_memd_set_permittivity";
    memd_struct* memd = (memd_struct*) rawdata;
    
    if ( ( memd->parameters.temperature > 0.0 ) && ( memd->parameters.bjerrum > 0.0 ) ){
        return fcsResult_create(FCS_LOGICAL_ERROR, fnc_name, "You can only set 2 of the parameters bjerrum, temperature and epsilon.");
    }
    else if ( memd->parameters.temperature > 0.0 ) {
        memd->parameters.permittivity = epsilon;
        memd->parameters.bjerrum = 1.0 / memd->parameters.permittivity / memd->parameters.temperature;
        memd->parameters.prefactor = sqrt( 4.0 * M_PI / memd->parameters.permittivity );
        memd->parameters.pref2 = 1.0 / memd->parameters.permittivity ;
    }
    else if ( memd->parameters.bjerrum > 0.0 ) {
        memd->parameters.permittivity = epsilon;
        memd->parameters.temperature = 1.0 / memd->parameters.permittivity / memd->parameters.bjerrum;
        memd->parameters.prefactor = sqrt( 4.0 * M_PI / memd->parameters.permittivity );
        memd->parameters.pref2 = 1.0 / memd->parameters.permittivity ;
    }
    else memd->parameters.permittivity = epsilon;
    return NULL;
}

FCSResult fcs_memd_set_temperature(void* rawdata, fcs_float temperature)
{
    const char* fnc_name = "fcs_memd_set_temperature";
    memd_struct* memd = (memd_struct*) rawdata;
    
    if ( ( memd->parameters.permittivity > 0.0 ) && ( memd->parameters.bjerrum > 0.0 ) ){
        return fcsResult_create(FCS_LOGICAL_ERROR, fnc_name, "You can only set 2 of the parameters bjerrum, temperature and epsilon.");
    }
    else if ( memd->parameters.permittivity > 0.0 ) {
        memd->parameters.temperature = temperature;
        memd->parameters.bjerrum = 1.0 / memd->parameters.permittivity / memd->parameters.temperature;
        memd->parameters.prefactor = sqrt( 4.0 * M_PI / memd->parameters.permittivity );
        memd->parameters.pref2 = 1.0 / memd->parameters.permittivity ;
    }
    else if ( memd->parameters.bjerrum > 0.0 ) {
        memd->parameters.temperature = temperature;
        memd->parameters.permittivity = 1.0 / memd->parameters.temperature / memd->parameters.bjerrum;
        memd->parameters.prefactor = sqrt( 4.0 * M_PI / memd->parameters.permittivity );
        memd->parameters.pref2 = 1.0 / memd->parameters.permittivity ;
    }
    else memd->parameters.temperature = temperature;
    return NULL;
}

FCSResult fcs_memd_set_bjerrum_length(void* rawdata, fcs_float bjerrum)
{
    const char* fnc_name = "fcs_memd_set_bjerrum_length";
    memd_struct* memd = (memd_struct*) rawdata;
    
    if ( ( memd->parameters.permittivity > 0.0 ) && ( memd->parameters.temperature > 0.0 ) ){
        return fcsResult_create(FCS_LOGICAL_ERROR, fnc_name,
                                "You can only set 2 of the parameters bjerrum, temperature and epsilon.");
    }
    else if ( memd->parameters.permittivity > 0.0 ) {
        memd->parameters.bjerrum = bjerrum;
        memd->parameters.temperature = 1.0 / memd->parameters.permittivity / memd->parameters.bjerrum;
        memd->parameters.prefactor = sqrt( 4.0 * M_PI / memd->parameters.permittivity );
        memd->parameters.pref2 = 1.0 / memd->parameters.permittivity ;
        return NULL;
    }
    else if ( memd->parameters.temperature > 0.0 ) {
        memd->parameters.bjerrum = bjerrum;
        memd->parameters.permittivity = 1.0 / memd->parameters.temperature / memd->parameters.bjerrum;
        memd->parameters.prefactor = sqrt( 4.0 * M_PI / memd->parameters.permittivity );
        memd->parameters.pref2 = 1.0 / memd->parameters.permittivity ;
        return NULL;
    }
    else memd->parameters.bjerrum = bjerrum;
    return NULL;
}

FCSResult fcs_memd_set_total_number_of_particles(void* rawdata, fcs_int number_of_particles)
{
    const char* fnc_name = "fcs_memd_set_total_number_of_particles";
    memd_struct* memd = (memd_struct*) rawdata;
    memd->parameters.n_part_total = number_of_particles;
    return NULL;
}

FCSResult fcs_memd_set_local_number_of_particles(void* rawdata, fcs_int number_of_particles)
{
    const char* fnc_name = "fcs_memd_set_local_number_of_particles";
    memd_struct* memd = (memd_struct*) rawdata;
    memd->parameters.n_part = number_of_particles; 
    return NULL;
}

FCSResult fcs_memd_set_box_size(void* rawdata, fcs_float length_x, fcs_float length_y, fcs_float length_z)
{
    const char* fnc_name = "fcs_memd_set_box_size";
    memd_struct* memd = (memd_struct*) rawdata;

    /* Set boxlength parameter in memd struct */
    if ( ( ! fcs_float_is_equal(length_x, length_y) ) || ( ! fcs_float_is_equal(length_y, length_z) ) ) {
        return fcsResult_create(FCS_LOGICAL_ERROR, fnc_name,
                                "Only cubic boxes are allowed at this stage.");
    }
    else {
        memd->parameters.box_length[0] = length_x;
        memd->parameters.box_length[1] = length_y;
        memd->parameters.box_length[2] = length_z;
        return NULL;
    }

    /* set all dependent parameters */
    if ( (memd->parameters.mesh>0) && (memd->parameters.box_length[0]>0.0) ) {
        memd->parameters.inva  = (fcs_float) memd->parameters.mesh/memd->parameters.box_length[0];
        memd->parameters.a     = 1.0/memd->parameters.inva;
    } else {
        int k;
        //        fprintf(stdout, "box_l: %f\n", memd->parameters.box_length[0]); fflush(stdout);
        memd->parameters.mesh=32;
        if (memd->parameters.box_length[0]<ROUND_ERROR_PREC) {
            FOR3D(k) memd->parameters.box_length[k] = 10.0;
        }
        memd->parameters.inva  = (fcs_float) memd->parameters.mesh/memd->parameters.box_length[0];
        memd->parameters.a     = 1.0/memd->parameters.inva;
    }
    
}

FCSResult fcs_memd_set_mesh_size_1D(void* rawdata, fcs_int mesh_size)
{
    const char* fnc_name = "fcs_memd_set_mesh_size_1D";
    memd_struct* memd = (memd_struct*) rawdata;
    memd->parameters.mesh = mesh_size;
    return NULL;
}

fcs_float fcs_memd_get_time_step(void* rawdata)
{
    const char* fnc_name = "fcs_memd_set_time_step";
    memd_struct* memd = (memd_struct*) rawdata;
    return memd->parameters.time_step;
}

fcs_int fcs_memd_needs_retuning(void* rawdata, fcs_int local_particles, fcs_float* positions, fcs_float* charges)
{
    const char* fnc_name = "fcs_memd_needs_retuning";
    memd_struct* memd = (memd_struct*) rawdata;
    return 0;
}

FCSResult fcs_memd_tune_method(void* rawdata, fcs_int local_particles, fcs_float* positions, fcs_float* charges)
{
    const char* fnc_name = "fcs_memd_tune_method";
    memd_struct* memd = (memd_struct*) rawdata;
    fcs_int n_part_total = memd->parameters.n_part_total;

    if (n_part_total < 1) {
        n_part_total = ifcs_memd_count_total_charges(memd, local_particles);
        memd->parameters.n_part_total = n_part_total;
    }
    
    fcs_float boxlength = memd->parameters.box_length[0];
    
    fcs_int meshsize1D = 16;
    fcs_float lightspeed = 1.0;
    fcs_float timestep = 0.01;
    
    fcs_float avg_dist = pow((boxlength*boxlength*boxlength) / (fcs_float)n_part_total, 0.33333);
    fcs_float dists_per_box = boxlength/avg_dist;

    fcs_float epsilon = 1.0;
    fcs_float temperature = 1.0;
    
    
    meshsize1D = ifcs_memd_get_next_higher_power_of_two(dists_per_box);
    
    /* The factor 0.01 is "<<" for the stability criterion */
    lightspeed = 0.01 * (boxlength/meshsize1D) / timestep;
    
    fcs_memd_set_speed_of_light(rawdata, lightspeed);
    fcs_memd_set_mesh_size_1D(rawdata, meshsize1D);
    fcs_memd_set_time_step(rawdata, timestep);
//    printf("mesh_size_1D: %d\n", meshsize1D);
    
    fcs_memd_set_permittivity(rawdata, epsilon);
    fcs_memd_set_temperature(rawdata, temperature);
    
    fcs_memd_setup_local_lattice(memd);
    
    return NULL;    
}
