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
#include "data_types.h"
#include "FCSCommon.h"
#include <math.h>

FCSResult memd_set_init_flag(void* rawdata, fcs_int flagvalue)
{
    const char* fnc_name = "memd_set_init_flag";

    memd_struct* memd = (memd_struct*) rawdata;
    memd->init_flag = flagvalue;
    return NULL;
}

FCSResult memd_set_speed_of_light(void* rawdata, fcs_float lightspeed)
{
    const char* fnc_name = "memd_set_speed_of_light";
    memd_struct* memd = (memd_struct*) rawdata;

    fcs_float f_mass = 1.0/lightspeed;
    fcs_float invsqrt_f_mass = lightspeed * lightspeed;

    memd->parameters.f_mass = f_mass;
    memd->parameters.invsqrt_f_mass = invsqrt_f_mass;
    
    return NULL;
}

FCSResult memd_set_time_step(void* rawdata, fcs_float timestep)
{
    const char* fnc_name = "memd_set_time_step";
    memd_struct* memd = (memd_struct*) rawdata;

    memd->parameters.time_step = timestep;
    return NULL;
}

FCSResult memd_set_permittivity(void* rawdata, fcs_float epsilon)
{
    const char* fnc_name = "memd_set_permittivity";
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

FCSResult memd_set_temperature(void* rawdata, fcs_float temperature)
{
    const char* fnc_name = "memd_set_temperature";
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

FCSResult memd_set_bjerrum_length(void* rawdata, fcs_float bjerrum)
{
    const char* fnc_name = "memd_set_bjerrum_length";
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

FCSResult memd_set_total_number_of_particles(void* rawdata, fcs_int number_of_particles)
{
    const char* fnc_name = "memd_set_total_number_of_particles";
    memd_struct* memd = (memd_struct*) rawdata;
    memd->parameters.n_part_total = number_of_particles;
    return NULL;
}

FCSResult memd_set_local_number_of_particles(void* rawdata, fcs_int number_of_particles)
{
    const char* fnc_name = "memd_set_local_number_of_particles";
    memd_struct* memd = (memd_struct*) rawdata;
    memd->parameters.n_part = number_of_particles; 
    return NULL;
}

FCSResult memd_set_box_size(void* rawdata, fcs_float length_x, fcs_float length_y, fcs_float length_z)
{
    const char* fnc_name = "memd_set_box_size";
    memd_struct* memd = (memd_struct*) rawdata;

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
}

FCSResult memd_set_mesh_size_1D(void* rawdata, fcs_int mesh_size)
{
    const char* fnc_name = "memd_set_mesh_size_1D";
    memd_struct* memd = (memd_struct*) rawdata;
    memd->parameters.mesh = mesh_size;
    return NULL;
}

fcs_float memd_get_time_step(void* rawdata)
{
    const char* fnc_name = "memd_set_time_step";
    memd_struct* memd = (memd_struct*) rawdata;
    return memd->parameters.time_step;
}

fcs_int memd_needs_retuning(void* rawdata, fcs_int local_particles, fcs_float* positions, fcs_float* charges)
{
    const char* fnc_name = "memd_needs_retuning";
    memd_struct* memd = (memd_struct*) rawdata;
    return 0;
}

FCSResult memd_tune_method(void* rawdata, fcs_int local_particles, fcs_float* positions, fcs_float* charges)
{
    const char* fnc_name = "memd_tune_method";
    memd_struct* memd = (memd_struct*) rawdata;
    return NULL;    
}
