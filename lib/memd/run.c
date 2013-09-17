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
#include <stdlib.h>

#include "run.h"
#include "data_types.h"
#include <math.h>
#include "FCSResult.h"
#include "forces.h"
#include "getter_setter.h"
#include "initial_solution.h"
#include "helper_functions.h"
#include "communication.h"
#include "init.h"


void ifcs_memd_assign_charges(memd_struct* memd, fcs_int local_num_real_particles, fcs_float* local_positions, fcs_float* local_charges, fcs_float* local_fields){
    int k, cell_shift[3], cell_id, part_id, cell_part_id;
    
    for (part_id=0; part_id<local_num_real_particles; part_id++) {
        fcs_float part_position[3] = {local_positions[part_id*3], local_positions[part_id*3+1], local_positions[part_id*3+2]};
        fcs_float part_charge = local_charges[part_id];
        FOR3D(k) cell_shift[k] = (int) floor(
                    (part_position[k] - memd->lparams.left_down_position[k])
                                             / memd->parameters.a );
        printf("lparams.dim: %d %d %d\n", memd->lparams.dim[0], memd->lparams.dim[1], memd->lparams.dim[2]);
        cell_id = ifcs_memd_get_linear_index(cell_shift[0], cell_shift[1], cell_shift[2], memd->lparams.dim);
        printf("cell_id: %d\n", cell_id);
/*
        fprintf(stdout, "cell_shift: %d %d %d\n", cell_shift[0], cell_shift[1], cell_shift[2]); fflush(stdout);
        fprintf(stdout, "position: %f %f %f\n", part_position[0], part_position[1], part_position[2]); fflush(stdout);
        fprintf(stdout, "box_length: %f\n", memd->parameters.box_length[0]); fflush(stdout);
        fprintf(stdout, "lparams.dim: %d\n", memd->lparams.dim[2]); fflush(stdout);
        fprintf(stdout, "parameters.a: %f\n", memd->parameters.a); fflush(stdout);
        fprintf(stdout, "CellID: %d\n", cell_id); fflush(stdout);
*/
        memd->local_cells.cell[cell_id]->n = memd->local_cells.cell[cell_id]->n + 1;
        cell_part_id = memd->local_cells.cell[cell_id]->n - 1;
        memd->local_cells.cell[cell_id]->part[cell_part_id].r = part_position;
        memd->local_cells.cell[cell_id]->part[cell_part_id].q = part_charge;
        memd->local_cells.cell[cell_id]->part[cell_part_id].f = &local_fields[part_id*3];
    }
    
}

void ifcs_memd_run(void* rawdata, fcs_int num_particles, fcs_int max_num_particles, fcs_float *positions, fcs_float *charges, fcs_float *fields, fcs_float *potentials){

    memd_struct* memd = (memd_struct*) rawdata;;
    /*
    if (rawdata == NULL) {
        memd = calloc(1,sizeof(memd_struct));
        *rawdata = memd;
        printf("New handle created!\n"); fflush(stdout);
    } else {
        memd = (memd_struct*) rawdata;
        //fcs_memd_setup_communicator(memd, memd->mpiparams.communicator);
    }
     */

    /* decompose system */
    fcs_int local_num_particles;
    fcs_int local_num_real_particles;
    fcs_int local_num_ghost_particles;
    fcs_float *local_positions, *local_ghost_positions;
    fcs_float *local_charges, *local_ghost_charges;
    fcs_gridsort_index_t *local_indices, *local_ghost_indices;
    fcs_gridsort_t gridsort;
    fcs_float box_base[3] = {0.0, 0.0, 0.0 };
    fcs_float box_a[3] = {memd->parameters.box_length[0], 0.0, 0.0 };
    fcs_float box_b[3] = {0.0, memd->parameters.box_length[1], 0.0 };
    fcs_float box_c[3] = {0.0, 0.0, memd->parameters.box_length[2] };
    
    fcs_gridsort_create(&gridsort);
    fcs_gridsort_set_system(&gridsort, box_base, box_a, box_b, box_c, NULL);
    fcs_gridsort_set_particles(&gridsort, num_particles, max_num_particles, positions, charges);
    fcs_gridsort_sort_forward(&gridsort, 0.0, memd->mpiparams.communicator);
    fcs_gridsort_separate_ghosts(&gridsort, &local_num_real_particles, &local_num_ghost_particles);
    fcs_gridsort_get_sorted_particles(&gridsort, &local_num_particles, NULL, NULL, NULL, NULL);
    fcs_gridsort_get_real_particles(&gridsort, &local_num_real_particles, &local_positions, &local_charges, &local_indices);
    fcs_gridsort_get_ghost_particles(&gridsort, &local_num_ghost_particles, &local_ghost_positions, &local_ghost_charges, &local_ghost_indices);

    
    /* allocate local fields and potentials */
    fcs_float *local_fields = NULL;
    fcs_float *local_potentials = NULL;
    if (fields != NULL)
        local_fields = malloc(sizeof(fcs_float)*3*local_num_real_particles);
    if (potentials != NULL || memd->total_energy_flag)
        local_potentials = malloc(sizeof(fcs_float)*local_num_real_particles);

    if (memd->init_flag) {
        ifcs_memd_init(rawdata, memd->mpiparams.communicator);
        fcs_memd_setup_local_lattice(memd);
        /* charge assignment */
        ifcs_memd_assign_charges(memd, local_num_real_particles, local_positions, local_charges, local_fields);
        /* enforce electric field onto the Born-Oppenheimer surface */
        ifcs_memd_calc_init_e_field(memd);
    } else {
        fcs_memd_setup_local_lattice(memd);
        ifcs_memd_assign_charges(memd, local_num_real_particles, local_positions, local_charges, local_fields);
        printf("charges assigned.\n");
    }
    
    printf("Charge transfer done\n"); fflush(stdout);

    fcs_float timestep = fcs_memd_get_time_step(rawdata);
    printf("Timestep: %f\n", timestep); fflush(stdout);
    fcs_memd_propagate_B_field(memd, (timestep/2.0) );
    fcs_memd_calc_forces(memd);
    fcs_memd_propagate_B_field(memd, (timestep/2.0) );
}
