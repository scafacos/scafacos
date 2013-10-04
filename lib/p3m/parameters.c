/*
  Copyright (C) 2011,2012 Olaf Lenz
  
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
#include "parameters.h"
#include "types.h"
#include "FCSCommon.h"
#include <stdio.h>

void ifcs_p3m_set_near_field_flag(void *rd, fcs_int flag) {
    ifcs_p3m_data_struct *d = (ifcs_p3m_data_struct*) rd;
    d->near_field_flag = flag;
}

void ifcs_p3m_set_box_a(void* rd, fcs_float a) {
    ifcs_p3m_data_struct *d = (ifcs_p3m_data_struct*) rd;
    if (!fcs_float_is_equal(a, d->box_l[0]))
        d->needs_retune = 1;
    d->box_l[0] = a;
}

void ifcs_p3m_set_box_b(void *rd, fcs_float b) {
    ifcs_p3m_data_struct *d = (ifcs_p3m_data_struct*) rd;
    if (!fcs_float_is_equal(b, d->box_l[1]))
        d->needs_retune = 1;
    d->box_l[1] = b;
}

void ifcs_p3m_set_box_c(void *rd, fcs_float c) {
    ifcs_p3m_data_struct *d = (ifcs_p3m_data_struct*) rd;
    if (!fcs_float_is_equal(c, d->box_l[2]))
        d->needs_retune = 1;
    d->box_l[2] = c;
}


//TODO I assume these functions should be placed somewhere else

fcs_float unit_volume(fcs_float alpha, fcs_float beta, fcs_float gamma) {
    return (fcs_float)fcs_sqrt((1 - cos(alpha) * cos(alpha) - cos(beta) * cos(beta) - cos(gamma) * cos(gamma) + 2 * cos(alpha) * cos(gamma) * cos(beta)));
}

void tricFROMcart(void* rd, fcs_float *vec) {
    ifcs_p3m_data_struct *d = (ifcs_p3m_data_struct*) rd;
    vec[0] = vec[0] / d->box_a - cos(d->box_gamma) / (d->box_a * sin(d->box_gamma)) * vec[1]
            + (cos(d->box_alpha) * cos(d->box_gamma) - cos(d->box_beta))
            / (d->box_a * unit_volume(d->box_alpha, d->box_beta, d->box_gamma) * sin(d->box_gamma))
            * vec[2];
    vec[1] = 1 / (d->box_b * sin(d->box_gamma)) * vec[1]
            + (cos(d->box_beta) * cos(d->box_gamma) - cos(d->box_alpha)) / sin(d->box_gamma) / d->box_b
            / unit_volume(d->box_alpha, d->box_beta, d->box_gamma) * vec[2];
    vec[2] = 1 / (d->box_c * unit_volume(d->box_alpha, d->box_beta, d->box_gamma)) * sin(d->box_gamma)
            * vec[2];
}
fcs_float angle_between_vectors(fcs_float *vec_a, fcs_float *vec_b) {
    
    return acos(
            (vec_a[0] * vec_b[0] + vec_a[1] * vec_b[1] + vec_a[2] * vec_b[2])
            / (fcs_norm(vec_a) * fcs_norm(vec_b)));
}


void cartFROMtric(void* rd, fcs_float *vec) {//@todo move those calculations to paramters.c to avoid errors (data struct unkown and request for member in something not a struct)
    ifcs_p3m_data_struct *d = (ifcs_p3m_data_struct*) rd;
    vec[0] = d->box_a * vec[0] + d->box_b * cos(d->box_gamma) * vec[1]
            + d->box_c * fcs_cos(d->box_beta) * vec[2]; //a*m+b*cos(gamma)*n+c*cos(beta)*l
    vec[1] = d->box_b * fcs_sin(d->box_gamma) * vec[1]
            + d->box_c * (fcs_cos(d->box_alpha) - cos(d->box_beta) * cos(d->box_gamma)) / sin(d->box_gamma)
            * vec[2]; //b*sin(gamma)*n+c*(cos(alpha)-cos(beta)*cos(gamma))/(sin(gamma))*l
    vec[2] = d->box_c * unit_volume(d->box_alpha, d->box_beta, d->box_gamma) / sin(d->box_gamma) * vec[2]; //c*v/sin(gamma)*l
}


void ifcs_p3m_set_box_geometry(void *rd, fcs_float *a, fcs_float *b, fcs_float *c) {
    ifcs_p3m_data_struct *d = (ifcs_p3m_data_struct*) rd;
    // if (!fcs_float_is_equal(a, d->box_l[0])||!fcs_float_is_equal(a, d->box_l[0])||!fcs_float_is_equal(c, d->box_l[2])) //TODO check if d_boxlength still fits (first i need to know what the box lengths will be)
    //d->needs_retune = 1;

        //from here on: angles and lengths are set
        d->box_alpha = angle_between_vectors(b, c); //angle between box vector b and c
        d->box_beta = angle_between_vectors(a, c);
        d->box_gamma = angle_between_vectors(a, b); //angle between box vector a and b;

        //TODO what size should the  
        d->box_a = fcs_norm(a);
        d->box_b = fcs_norm(b);
        d->box_c = fcs_norm(c);

}

void ifcs_p3m_set_r_cut(void *rd, fcs_float r_cut) {
    ifcs_p3m_data_struct *d = (ifcs_p3m_data_struct*) rd;
    if (!fcs_float_is_equal(r_cut, d->r_cut))
        d->needs_retune = 1;
    d->r_cut = r_cut;
    d->tune_r_cut = 0;
}

void ifcs_p3m_set_r_cut_tune(void *rd) {
    ifcs_p3m_data_struct *d = (ifcs_p3m_data_struct*) rd;
    d->needs_retune = 1;
    d->tune_r_cut = 1;
}

void ifcs_p3m_get_r_cut(void *rd, fcs_float *r_cut) {
    ifcs_p3m_data_struct *d = (ifcs_p3m_data_struct*) rd;
    *r_cut = d->r_cut;
}

void ifcs_p3m_set_alpha(void *rd, fcs_float alpha) {
    ifcs_p3m_data_struct *d = (ifcs_p3m_data_struct*) rd;
    if (!fcs_float_is_equal(alpha, d->alpha))
        d->needs_retune = 1;
    d->alpha = alpha;
    d->tune_alpha = 0;
}

void ifcs_p3m_set_alpha_tune(void *rd) {
    ifcs_p3m_data_struct *d = (ifcs_p3m_data_struct*) rd;
    d->needs_retune = 1;
    d->tune_alpha = 1;
}

void ifcs_p3m_get_alpha(void *rd, fcs_float *alpha) {
    ifcs_p3m_data_struct *d = (ifcs_p3m_data_struct*) rd;
    *alpha = d->alpha;
}

void ifcs_p3m_set_grid(void *rd, fcs_int grid) {
    ifcs_p3m_data_struct *d = (ifcs_p3m_data_struct*) rd;
    if (!(grid == d->grid[0] && grid == d->grid[1] && grid == d->grid[2]))
        d->needs_retune = 1;
    d->grid[0] = grid;
    d->grid[1] = grid;
    d->grid[2] = grid;
    d->tune_grid = 0;
}

void ifcs_p3m_set_grid_tune(void *rd) {
    ifcs_p3m_data_struct *d = (ifcs_p3m_data_struct*) rd;
    d->needs_retune = 1;
    d->tune_grid = 1;
}

void ifcs_p3m_get_grid(void *rd, fcs_int *grid) {
    ifcs_p3m_data_struct *d = (ifcs_p3m_data_struct*) rd;
    grid[0] = d->grid[0];
    grid[1] = d->grid[1];
    grid[2] = d->grid[2];
}

void ifcs_p3m_set_cao(void *rd, fcs_int cao) {
    ifcs_p3m_data_struct *d = (ifcs_p3m_data_struct*) rd;
    if (cao != d->cao)
        d->needs_retune = 1;
    d->cao = cao;
    d->tune_cao = 0;
}

void ifcs_p3m_set_cao_tune(void *rd) {
    ifcs_p3m_data_struct *d = (ifcs_p3m_data_struct*) rd;
    d->needs_retune = 1;
    d->tune_cao = 1;
}

void ifcs_p3m_get_cao(void *rd, fcs_int *cao) {
    ifcs_p3m_data_struct *d = (ifcs_p3m_data_struct*) rd;
    *cao = d->cao;
}

void ifcs_p3m_set_tolerance_field(void *rd, fcs_float tolerance_field) {
    ifcs_p3m_data_struct *d = (ifcs_p3m_data_struct*) rd;
    if (!fcs_float_is_equal(tolerance_field, d->tolerance_field))
        d->needs_retune = 1;
    d->tolerance_field = tolerance_field;
}

void ifcs_p3m_set_tolerance_field_tune(void *rd) {
    ifcs_p3m_data_struct *d = (ifcs_p3m_data_struct*) rd;
    d->needs_retune = 1;
    d->tolerance_field = P3M_DEFAULT_TOLERANCE_FIELD;
}

void ifcs_p3m_get_tolerance_field(void *rd, fcs_float* tolerance_field) {
    ifcs_p3m_data_struct *d = (ifcs_p3m_data_struct*) rd;
    *tolerance_field = d->tolerance_field;
}

void ifcs_p3m_require_total_energy(void *rd, fcs_int flag) {
    ifcs_p3m_data_struct *d = (ifcs_p3m_data_struct*) rd;
    d->require_total_energy = flag;
}

FCSResult ifcs_p3m_get_total_energy(void *rd, fcs_float *total_energy) {
    const char* fnc_name = "ifcs_p3m_get_total_energy";
    ifcs_p3m_data_struct *d = (ifcs_p3m_data_struct*) rd;
    if (d->require_total_energy) {
        *total_energy = d->total_energy;
        return NULL;
    } else
        return
        fcsResult_create
            (FCS_LOGICAL_ERROR, fnc_name,
            "Trying to get total energy, but computation was not requested.");
}

void ifcs_p3m_require_timings(void *rd, fcs_int flag) {
    ifcs_p3m_data_struct *d = (ifcs_p3m_data_struct*) rd;
    d->require_timings = flag;
}

FCSResult
ifcs_p3m_get_timings(void *rd,
        double *timing,
        double *timing_near_field,
        double *timing_far_field) {
    const char* fnc_name = "ifcs_p3m_get_timings";
    ifcs_p3m_data_struct *d = (ifcs_p3m_data_struct*) rd;

    if (!d->require_timings)
        return
        fcsResult_create
            (FCS_LOGICAL_ERROR, fnc_name,
            "Trying to get timings, but timings were not requested.");

    *timing = d->timings[TIMING];
    *timing_near_field = d->timings[TIMING_NEAR];
    *timing_far_field = d->timings[TIMING_FAR];

    return NULL;
}
