/*
 *    vmg - a versatile multigrid solver
 *    Copyright (C) 2012 Institute for Numerical Simulation, University of Bonn
 *
 *  vmg is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  vmg is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * @file   interface_fcs.h
 * @author Julian Iseringhausen <isering@ins.uni-bonn.de>
 * @date   Tue Apr 12 17:39:27 2011
 *
 * @brief  Scafacos C interface.
 *
 *
 */

#ifdef __cplusplus
extern "C" {
#endif

void VMG_fcs_setup(fcs_int level, fcs_int* periodic, fcs_int max_iter,
		   fcs_int smoothing_steps, fcs_int cycle_type, fcs_float precision,
		   fcs_float* box_offset, fcs_float box_size,
		   fcs_int near_field_cells, fcs_int interpolation_degree,
                   fcs_int discretization_order, MPI_Comm mpi_comm);

int VMG_fcs_check();

void VMG_fcs_run(fcs_float* pos, fcs_float* charge, fcs_float* potential, fcs_float* f, fcs_int num_particles_local);

void VMG_fcs_print_timer(void);

void VMG_fcs_destroy(void);

#ifdef __cplusplus
} /* extern "C" */
#endif
