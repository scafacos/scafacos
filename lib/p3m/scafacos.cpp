/*
 Copyright (C) 2014 Olaf Lenz

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
/* Implements the functions that translate between the scafacos
 C-interface and the C++-P3M interface. */
#include "scafacos.h"
#include "src/Solver.hpp"
#include "src/utils.hpp"
#include <stdexcept>

extern "C" {

using namespace P3M;

FCSResult ifcs_p3m_init(void **rd, MPI_Comm communicator) {
	Solver *d;

	try {
		if (*rd == NULL) {
			d = new Solver(communicator);

			/* store the new pointer in rd */
			*rd = d;
		} else
			d = static_cast<Solver*>(*rd);
	} catch (std::exception &e) {
		return fcs_result_create(FCS_ERROR_LOGICAL_ERROR, "ifcs_p3m_init",
				e.what());
	}
	return FCS_RESULT_SUCCESS;
}

void ifcs_p3m_destroy(void *rd) {
	if (rd != NULL) {
		Solver *d = static_cast<Solver*>(rd);
		delete d;
	}
}

fcs_int ifcs_p3m_check_triclinic_box(fcs_float ay, fcs_float az, fcs_float bz) {
    if (ay == 0 && az == 0 && bz == 0) return 1;
    return 0;
}

void ifcs_p3m_set_near_field_flag(void *rd, fcs_int flag) {
	Solver *d = static_cast<Solver *>(rd);
	d->near_field_flag = flag;
    }

FCSResult ifcs_p3m_set_triclinic_flag(void *rd) {
        P3M_DEBUG_LOCAL(printf("P3M triclinic: triclinic box detected.\n"));
#if defined(P3M_INTERLACE) && defined(P3M_AD)
        Solver *d = static_cast<Solver *> (rd);
        d->isTriclinic = true;
        return FCS_RESULT_SUCCESS;
#else
        const char* fnc_name = "ifcs_p3m_set_triclinic_flag";
        return fcs_result_create(FCS_ERROR_NOT_IMPLEMENTED, fnc_name,
                "p3m triclinic is only implemented with interlaced AD so far.");
#endif
}

void ifcs_p3m_set_box_a(void* rd, fcs_float a) {
	Solver *d = static_cast<Solver *>(rd);
	if (!float_is_equal(a, d->box_l[0]))
		d->needs_retune = 1;
	d->box_l[0] = a;
}

void ifcs_p3m_set_box_b(void *rd, fcs_float b) {
	Solver *d = static_cast<Solver *>(rd);
	if (!float_is_equal(b, d->box_l[1]))
		d->needs_retune = 1;
	d->box_l[1] = b;
}

void ifcs_p3m_set_box_c(void *rd, fcs_float c) {
	Solver *d = static_cast<Solver *>(rd);
	if (!float_is_equal(c, d->box_l[2]))
		d->needs_retune = 1;
	d->box_l[2] = c;
}

void ifcs_p3m_set_box_geometry(void *rd, const fcs_float *a, const fcs_float *b, const fcs_float *c) {
    Solver *d = static_cast<Solver *>(rd);

    for (p3m_int i = 0; i < 3; i++) {
        d->box_vectors[0][i] = a[i];
        d->box_vectors[1][i] = b[i];
        d->box_vectors[2][i] = c[i];
    }
    
    d->volume = d->box_vectors[0][0]*d->box_vectors[1][1]*d->box_vectors[2][2];
}

void ifcs_p3m_set_r_cut(void *rd, fcs_float r_cut) {
	Solver *d = static_cast<Solver *>(rd);
	if (!float_is_equal(r_cut, d->r_cut))
		d->needs_retune = 1;
	d->r_cut = r_cut;
	d->tune_r_cut = 0;
}

void ifcs_p3m_set_r_cut_tune(void *rd) {
	Solver *d = static_cast<Solver *>(rd);
	d->needs_retune = 1;
	d->tune_r_cut = 1;
}

void ifcs_p3m_get_r_cut(void *rd, fcs_float *r_cut) {
	Solver *d = static_cast<Solver *>(rd);
	*r_cut = d->r_cut;
}

void ifcs_p3m_set_alpha(void *rd, fcs_float alpha) {
	Solver *d = static_cast<Solver *>(rd);
	if (!float_is_equal(alpha, d->alpha))
		d->needs_retune = 1;
	d->alpha = alpha;
	d->tune_alpha = 0;
}

void ifcs_p3m_set_alpha_tune(void *rd) {
	Solver *d = static_cast<Solver *>(rd);
	d->needs_retune = 1;
	d->tune_alpha = 1;
}

void ifcs_p3m_get_alpha(void *rd, fcs_float *alpha) {
	Solver *d = static_cast<Solver *>(rd);
	*alpha = d->alpha;
}

void ifcs_p3m_get_near_params(void* rd, fcs_float *alpha, fcs_float *offset) {
  Solver *d = static_cast<Solver *>(rd);
  *alpha = d->alpha;
  *offset = d->shiftGaussians?(1-erf(d->r_cut * d->alpha)) / d->r_cut:0.0; //todo: is it better to add the offset as a property of a Solver object?
  printf("get offset %e", *offset);
}

void ifcs_p3m_set_grid(void *rd, fcs_int grid) {
	Solver *d = static_cast<Solver *>(rd);
	if (!(grid == d->grid[0] && grid == d->grid[1] && grid == d->grid[2]))
		d->needs_retune = 1;
	d->grid[0] = grid;
	d->grid[1] = grid;
	d->grid[2] = grid;
	d->tune_grid = 0;
}

void ifcs_p3m_set_grid_tune(void *rd) {
	Solver *d = static_cast<Solver *>(rd);
	d->needs_retune = 1;
	d->tune_grid = 1;
}

void ifcs_p3m_get_grid(void *rd, fcs_int *grid) {
	Solver *d = static_cast<Solver *>(rd);
	grid[0] = d->grid[0];
	grid[1] = d->grid[1];
	grid[2] = d->grid[2];
}

void ifcs_p3m_set_cao(void *rd, fcs_int cao) {
	Solver *d = static_cast<Solver *>(rd);
	if (cao != d->cao)
		d->needs_retune = 1;
	d->cao = cao;
	d->tune_cao = 0;
}

void ifcs_p3m_set_cao_tune(void *rd) {
	Solver *d = static_cast<Solver *>(rd);
	d->needs_retune = 1;
	d->tune_cao = 1;
}

void ifcs_p3m_get_cao(void *rd, fcs_int *cao) {
	Solver *d = static_cast<Solver *>(rd);
	*cao = d->cao;
}

void ifcs_p3m_set_tolerance_field(void *rd, fcs_float tolerance_field) {
	Solver *d = static_cast<Solver *>(rd);
	if (!float_is_equal(tolerance_field, d->tolerance_field))
		d->needs_retune = 1;
	d->tolerance_field = tolerance_field;
}

void ifcs_p3m_set_tolerance_field_tune(void *rd) {
	Solver *d = static_cast<Solver *>(rd);
	d->needs_retune = 1;
	d->tolerance_field = P3M_DEFAULT_TOLERANCE_FIELD;
}

void ifcs_p3m_get_tolerance_field(void *rd, fcs_float* tolerance_field) {
	Solver *d = static_cast<Solver *>(rd);
	*tolerance_field = d->tolerance_field;
}

void ifcs_p3m_require_total_energy(void *rd, fcs_int flag) {
	Solver *d = static_cast<Solver *>(rd);
	d->require_total_energy = flag;
}

FCSResult ifcs_p3m_get_total_energy(void *rd, fcs_float *total_energy) {
	const char* fnc_name = "ifcs_p3m_get_total_energy";
	Solver *d = static_cast<Solver *>(rd);
	if (d->require_total_energy) {
		*total_energy = d->total_energy;
		return FCS_RESULT_SUCCESS;
	} else
		return fcs_result_create(FCS_ERROR_LOGICAL_ERROR, fnc_name,
				"Trying to get total energy, but computation was not requested.");
}

void ifcs_p3m_require_timings(void *rd, fcs_int flag) {
	Solver *d = static_cast<Solver *>(rd);
	switch (flag) {
	case 0:
		d->require_timings = NONE;
		break;
	case 1:
		d->require_timings = ESTIMATE_ALL;
		break;
	case 2:
		d->require_timings = ESTIMATE_FFT;
		break;
	case 3:
		d->require_timings = ESTIMATE_ASSIGNMENT;
		break;
	case 4:
		d->require_timings = FULL;
		break;
	}
}

FCSResult ifcs_p3m_get_timings(void *rd, double *timing,
		double *timing_near_field, double *timing_far_field) {
	const char* fnc_name = "ifcs_p3m_get_timings";
	Solver *d = static_cast<Solver *>(rd);

	if (d->require_timings == NONE)
		return fcs_result_create(FCS_ERROR_LOGICAL_ERROR, fnc_name,
				"Trying to get timings, but timings were not requested.");

	*timing = d->timings[P3M::Solver::TIMING];
	*timing_near_field = d->timings[P3M::Solver::TIMING_NEAR];
	*timing_far_field = d->timings[P3M::Solver::TIMING_FAR];

	return NULL;
}

FCSResult ifcs_p3m_tune(void* rd,
fcs_int num_particles,
fcs_int max_particles,
fcs_float *positions,
fcs_float *charges) {
	Solver *d = static_cast<Solver *>(rd);

	try {
		d->tune(num_particles, positions, charges);
	} catch (std::exception &e) {
		return fcs_result_create(FCS_ERROR_LOGICAL_ERROR, "ifcs_p3m_init",
				e.what());
	}

	return FCS_RESULT_SUCCESS;
}

void ifcs_p3m_set_potential_shift(void* rd, fcs_int flag){
    Solver *d = static_cast<Solver *>(rd);
    d->shiftGaussians=(flag!=0);
}

void ifcs_p3m_get_potential_shift(void* rd, fcs_int *flag){
    Solver *d = static_cast<Solver *>(rd);
    *flag = d->shiftGaussians?1:0;
}
}
