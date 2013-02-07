/*
  Copyright (C) 2011,2012 Olaf Lenz
  
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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#define FCS_NEAR_FIELD

#include <stdio.h>

#include <stdlib.h>
#include "fcs.h"
void assert_fcs(FCSResult r)
{
  if (r) {
    fcsResult_printResult(r);
    MPI_Finalize();
    exit(-1);
  }
}

int main(int argc, char **argv) {
  const char* method = "p3m";
  const char* datafile = "../inp_data/p3m/p3m_wall.dat";
  fcs_int periodicity[3] = { 1, 1, 1 };
  fcs_float box_l[3] = { 10.0, 10.0, 10.0 };
  fcs_float offset[3] = {0.0, 0.0, 0.0};
  fcs_float p3m_tolerance_field = 1.e-3;
  fcs_int total_particles = 300;

#ifndef FCS_NEAR_FIELD
  fcs_int i,j;
  fcs_float r_cut;
#endif


  /***************************************************/
  /* SET UP MPI */
  MPI_Init(&argc, &argv);

  /* Create cartesian communicator */
  MPI_Comm comm = MPI_COMM_WORLD;
  /* fcs_int node_grid[3] = {0, 0, 0}; */
  /* fcs_int node[3]; */
  int comm_rank, comm_size;

  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
  if (comm_rank == 0) {
    fprintf(stderr, "----------------\n");
    fprintf(stderr, "Running p3m test\n");
    fprintf(stderr, "----------------\n");
    fprintf(stderr, "Setting up MPI...\n");
    fprintf(stderr, "  Using %d tasks.\n", comm_size);
  }

  /***************************************************/
  /* READ POSITIONS */
  fcs_int pid;

  fcs_float charges[300];
  fcs_float positions[900];
  fcs_float fields[900];
  fcs_float forces[900];
  fcs_float reference_forces[900];
  fcs_float potentials[300];
  fcs_float energies[300];
  fcs_float total_energy;
  fcs_int n_particles;

  if (comm_rank == 0) {
    n_particles = total_particles;
    fprintf(stderr, "Reading %s...\n", datafile);
    FILE *data = fopen(datafile, "r");
    if (!data) {
      fprintf(stderr, "ERROR: Can't read %s!", datafile);
      perror("ERROR");
      exit(1);
    }
    fcs_float charge_sum = 0.0;
    for (pid = 0; pid < n_particles; pid++) {
      fscanf(data, "%" FCS_CONV_FLOAT "f %" FCS_CONV_FLOAT "f %" FCS_CONV_FLOAT "f", &positions[3*pid], &positions[3*pid+1], &positions[3*pid+2]);
      fscanf(data, "%" FCS_CONV_FLOAT "f", &charges[pid]);
      fscanf(data, "%" FCS_CONV_FLOAT "f %" FCS_CONV_FLOAT "f %" FCS_CONV_FLOAT "f",
	     &reference_forces[3*pid],
	     &reference_forces[3*pid+1],
	     &reference_forces[3*pid+2]);
      /* fprintf(stderr, "%d: pos=(%lf, %lf, %lf), q=%lf, field=(%lf, %lf, %lf)\n", */
      /* 	   pid, positions[3*pid], positions[3*pid+1], positions[3*pid+2], */
      /* 	   charges[pid], */
      /* 	   reference_forces[3*pid], */
      /* 	   reference_forces[3*pid+1], */
      /* 	   reference_forces[3*pid+2]); */
      charge_sum += charges[pid];
    }
    fclose(data);
    /* fprintf(stderr, "  charge_sum=%lf\n", charge_sum); */
  } else n_particles = 0;

  /***************************************************/
  /* INIT SCAFACOS */

  FCS handle = NULL;
  FCSResult result = NULL;

  if (comm_rank == 0)
    fprintf(stderr, "Initializing p3m...\n");
  result = fcs_init(&handle, method, comm);
  assert_fcs(result);

  if (comm_rank == 0)
    fprintf(stderr, "  setting parameters...\n");

  fcs_float box_a[3] = { 0.0, 0.0, 0.0 };
  fcs_float box_b[3] = { 0.0, 0.0, 0.0 };
  fcs_float box_c[3] = { 0.0, 0.0, 0.0 };
  box_a[0] = box_l[0];
  box_b[1] = box_l[1];
  box_c[2] = box_l[2];

#ifdef FCS_NEAR_FIELD
  result = fcs_common_set(handle, 1, box_a, box_b, box_c, 
			  offset, periodicity, total_particles);
#else
  result = fcs_common_set(handle, 0, box_a, box_b, box_c, 
			  offset, periodicity, total_particles);
#endif

  if (result != NULL) {
    fcsResult_printResult(result);
    MPI_Abort(comm, 1);
  }

  assert_fcs(result);

  /***************************************************/
  /* TUNING */
  fcs_p3m_set_tolerance_field(handle, p3m_tolerance_field);

  /* Known reference values for this system */
  fcs_float p3m_r_cut = 1.001;
  fcs_int p3m_grid = 64;
  fcs_int p3m_cao = 7;
  fcs_float p3m_alpha = 2.70746;

  fcs_p3m_set_r_cut(handle, p3m_r_cut);
  fcs_p3m_set_grid(handle, p3m_grid);
  fcs_p3m_set_cao(handle, p3m_cao);
  fcs_p3m_set_alpha(handle, p3m_alpha);

  if (comm_rank == 0)
    fprintf(stderr, "Tuning p3m to accuracy %" FCS_LMOD_FLOAT "e...\n", p3m_tolerance_field);
  MPI_Barrier(comm);
  result = fcs_tune(handle, n_particles, n_particles, positions, charges);
  if (result != NULL) {
    fcsResult_printResult(result);
    MPI_Abort(comm, 1);
  }
  if (comm_rank == 0)
    fprintf(stderr, "Finished tuning.\n");
  MPI_Barrier(comm);
  
  /***************************************************/
  /* RUN */
  if (comm_rank == 0)
    fprintf(stderr, "Running p3m...\n");

  /* fcs_int grid_size; */
  /* fcs_p3m_get_grid_size(handle, &grid_size); */
  /* fcs_float *grid = malloc(grid_size*sizeof(fcs_float)); */
  /* fcs_p3m_set_grid_ptr(handle, grid); */

  /* fcs_float virial[9]; */
  /* result = fcs_require_virial(handle, 1); */

  result = fcs_p3m_require_total_energy(handle, 1);
  result = fcs_run(handle, n_particles, n_particles, 
		   positions, charges, fields, potentials);
  assert_fcs(result);
  result = fcs_p3m_get_total_energy(handle, &total_energy);
  assert_fcs(result);

#ifndef FCS_NEAR_FIELD
  if (comm_rank == 0) {
    fcs_p3m_get_r_cut(handle, &r_cut);

    fprintf(stderr, "Computing near fields (n*n-loop) for r_cut=%g...\n", r_cut);

    fcs_p3m_near_parameters_t near_params;
    fcs_p3m_get_near_parameters(handle, &near_params);
    
    /* compute near fields */
    for (i = 0; i < n_particles; i++) {
      for (j = 0; j < i; j++) {
	fcs_float d[3];
	
	d[0] = positions[3*j] - positions[3*i];
	while (d[0] < -0.5*box_a[0]) d[0] += box_a[0];
	while (d[0] > 0.5*box_a[0]) d[0] -= box_a[0];
      
	d[1] = positions[3*j+1] - positions[3*i+1];
	while (d[1] < -0.5*box_b[1]) d[1] += box_b[1];
	while (d[1] > 0.5*box_b[1]) d[1] -= box_b[1];
      
	d[2] = positions[3*j+2] - positions[3*i+2];
	while (d[2] < -0.5*box_c[2]) d[2] += box_c[2];
	while (d[2] > 0.5*box_c[2]) d[2] -= box_c[2];
      
	fcs_float length = sqrt(d[0]*d[0]+d[1]*d[1]+d[2]*d[2]);
	if (length < r_cut) {
	  fcs_float field;
	  fcs_float potential;
	  fcs_p3m_compute_near(near_params, length, &potential, &field);
	  fcs_float field_L = field / length;
	  fields[3*i] += field_L*charges[j]*d[0];
	  fields[3*i+1] += field_L*charges[j]*d[1];
	  fields[3*i+2] += field_L*charges[j]*d[2];
	  fields[3*j] += -field_L*charges[i]*d[0];
	  fields[3*j+1] += -field_L*charges[i]*d[1];
	  fields[3*j+2] += -field_L*charges[i]*d[2];
	  potentials[i] += potential*charges[j];
	  potentials[j] += potential*charges[i];
	}
      }
    }
  }
#endif

  if (comm_rank == 0) {

    /***************************************************/
    /* COMPUTE FORCES AND ENERGIES FROM FIELDS AND POTENTIALS */
    for (pid = 0; pid < n_particles; pid++) {
      forces[3*pid] = charges[pid] * fields[3*pid];
      forces[3*pid+1] = charges[pid] * fields[3*pid+1];
      forces[3*pid+2] = charges[pid] * fields[3*pid+2];
      energies[pid] = 0.5 * charges[pid] * potentials[pid];
    }

    fcs_float sqr_sum = 0.0;
    fcs_float sum_energy = 0.0;
    /***************************************************/
    /* COMPARE RESULTS TO REFERENCE */
    for (pid = 0; pid < n_particles; pid++) {
      sum_energy += energies[pid];

      fcs_float d0 = forces[3*pid] - reference_forces[3*pid];
      fcs_float d1 = forces[3*pid+1] - reference_forces[3*pid+1];
      fcs_float d2 = forces[3*pid+2] - reference_forces[3*pid+2];
      fcs_float sqr_err = d0*d0+d1*d1+d2*d2;
      sqr_sum += sqr_err;

      if (sqrt(sqr_err) > 2.0*p3m_tolerance_field) {
        fprintf(stderr, "pid=%d error=%f q=%" FCS_LMOD_FLOAT "f pos=(%+9.5" FCS_LMOD_FLOAT "f, %+9.5" FCS_LMOD_FLOAT "f, %+9.5" FCS_LMOD_FLOAT "f)\n",
		pid, sqrt(sqr_err), charges[pid],  
		positions[pid*3], positions[pid*3+1], positions[pid*3+2]);
        fprintf(stderr, "  reference_force=(%+9.5" FCS_LMOD_FLOAT "f, %+9.5" FCS_LMOD_FLOAT "f, %+9.5" FCS_LMOD_FLOAT "f)\n", 
		reference_forces[pid*3], reference_forces[pid*3+1], reference_forces[pid*3+2]);
        fprintf(stderr, "            force=(%+9.5" FCS_LMOD_FLOAT "f, %+9.5" FCS_LMOD_FLOAT "f, %+9.5" FCS_LMOD_FLOAT "f)\n", 
		forces[pid*3], forces[pid*3+1], forces[pid*3+2]);
        fprintf(stderr, "            field=(%" FCS_LMOD_FLOAT "f, %" FCS_LMOD_FLOAT "f, %" FCS_LMOD_FLOAT "f)\n", 
		fields[pid*3], fields[pid*3+1], fields[pid*3+2]);
      }
    }

    fprintf(stderr, "total_energy=%" FCS_LMOD_FLOAT "f\n", total_energy);
    fprintf(stderr, "sum_energy=%" FCS_LMOD_FLOAT "f\n", sum_energy);
    fprintf(stderr, "rms_error=%e\n", sqrt(sqr_sum / (fcs_float)n_particles));

    fprintf(stderr, "Finalizing...\n");

  }
  fcs_destroy(handle);

  MPI_Finalize();
  if (comm_rank == 0)
    fprintf(stderr, "Done.\n");

  return 0;
}
