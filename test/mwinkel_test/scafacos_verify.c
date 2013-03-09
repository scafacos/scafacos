/*************************************************************
 *
 * ScaFaCoS-Test programme
 *
 * Copyright (C) 2011 Mathias Winkel, Lukas Arnold
 *
 * This file is part of ScaFaCoS.
 *
 * ScaFaCoS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ScaFaCoS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser Public License for more details.
 *
 * You should have received a copy of the GNU Lesser Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 *************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <fcs_config.h>
#include <fcs.h>

#include <mpi.h>

#include "setup.h"

#define checkres(res) {if(res){fcsResult_printResult(res);fcsResult_destroy(res);exit(-1);}}
#define true 1
#define false 0
#define min(a,b) ( a>b ? b : a )


const char *fcs_method[NUMMETHODS] = {"DIRECT", "FMM", "MEMD", "P2NFFT", "P3M", "PEPC", "PP3MG", "VMG"};
const char *test_setup[NUMSETUPS]  = {"TRIVIAL", "NACL", "SIO2"};


// default values taken from interface test
#if FCS_ENABLE_PP3MG_PMG
    int PP3MG_dims[3];
    fcs_int PP3MG_cells_x = 64;/*64;*/
    fcs_int PP3MG_cells_y = 64;/*64;*/
    fcs_int PP3MG_cells_z = 64;/*64;*/
    fcs_int PP3MG_ghost_cells = 8;/*8*/
    fcs_int PP3MG_pol_degree = 4;/*4*/
    fcs_float PP3MG_err_bound = 1e-7;
    fcs_int PP3MG_max_iterations = 10;/*4*/
    fcs_int PP3MG_max_particles = 1000;/*10000;*/
    fcs_int PP3MG_periodic = 0;
#endif
#if FCS_ENABLE_P3M
    fcs_float p3m_required_accuracy = 1.e-3;
#endif
#if FCS_ENABLE_PEPC
    fcs_float PEPC_epsilon  =  0.0;
    fcs_float PEPC_theta    =  0.3;
    fcs_int PEPC_debuglevel = -1;
#endif
#if FCS_ENABLE_FMM
    fcs_int FMM_absrel = FCS_FMM_CUSTOM_RELATIVE;
    fcs_int FMM_dipole_correction = FCS_FMM_STANDARD_DIPOLE_CORRECTION;
    fcs_float FMM_deltaE = 1e-2;
#endif




FCSResult setup_methodspecific(FCS handle, int my_rank, int num_ranks)
{
  switch (fcs_get_method(handle))
    {
#if FCS_ENABLE_PP3MG_PMG
    case FCS_PP3MG:
        PP3MG_dims[0] = PP3MG_dims[1] = PP3MG_dims[2] = 0;
        if (num_ranks > 1)
            MPI_Dims_create(num_ranks,3,PP3MG_dims);
        else
            PP3MG_dims[0] = PP3MG_dims[1] = PP3MG_dims[2] = 1;

        return fcs_PP3MG_setup(handle, PP3MG_dims, PP3MG_cells_x, PP3MG_cells_y, PP3MG_cells_z, PP3MG_ghost_cells, PP3MG_max_iterations, PP3MG_max_particles, PP3MG_periodic, PP3MG_pol_degree, PP3MG_err_bound);
#endif
#if FCS_ENABLE_P3M
    case FCS_P3M:
        fcs_P3M_set_required_accuracy(handle, p3m_required_accuracy);
        return NULL;
#endif
#if FCS_ENABLE_PEPC
    case FCS_PEPC:
        return fcs_PEPC_setup(handle, PEPC_epsilon, PEPC_theta, PEPC_debuglevel);
#endif
#if FCS_ENABLE_VMG
    case FCS_VMG:
        return fcs_VMG_setup(handle, 6,0,6,25,3,1,1.0e-6);
#endif
#if FCS_ENABLE_FMM
    case FCS_FMM:
        return fcs_FMM_setup(handle, FMM_absrel, FMM_deltaE, FMM_dipole_correction);
#endif
    default:
        return NULL;
    }
}






void readparams(int argc, char** argv, int *method, int* nparticles, int* setup_id, int per[3])
{
    int i;

    if (argc > 1) *method     = atoi(argv[1]);
    if (argc > 2) *nparticles = atoi(argv[2]);
    if (argc > 3) *setup_id   = atoi(argv[3]);
    if (argc > 4)
        for (i=0;i<3;i++)
          per[i] = (argv[4][i] == '1') ? 1 : 0;

    if (*method   > NUMMETHODS-1 || *method < 0)   *method   = 0;
    if (*setup_id > NUMSETUPS -1 || *setup_id < 0) *setup_id = 0;
}


void printparams(char* commandline, int num_ranks, int my_rank, int method, int ntotal, int setup_id, int per[3])
{
  int i;

  if(my_rank!=0) return;

  printf("\n ***** ScaFaCoS Test ***** \n\n");
  printf("Call with '%s [METHOD [NUMPARTICLES [SETUP_ID [PERIODICITY]]]]'\n", commandline);
  printf("     e.g. '%s 5 10 0 011'\n", commandline);
  printf("Parameters:\n");
  printf("     NUMPARTICLES - total particle number (default %d)\n", NUMPARTICLES_DEFAULT);
  printf("     SETUP_ID     - one of the following integer values (default %d, current highlighted)\n", SETUP_DEFAULT);
  for (i=0;i<NUMSETUPS;i++)
  {
    if (i==setup_id) printf("     ==>"); else printf("        ");

    printf("  %2d  - %8s\n", i, test_setup[i]);
  }

  printf("     PERIODICITY  - bitmask for periodicity in three spatial dimensions (default '000')\n");
  printf("     METHOD       - one of the following integer values (default %d, current highlighted)\n", METHOD_DEFAULT);

  for (i=0;i<NUMMETHODS;i++)
  {
    if (i==method) printf("     -->"); else printf("        ");

    printf("  %2d  - %8s\n", i, fcs_method[i]);
  }

  printf("\nRunning on %d MPI ranks with in total %d particles\n", num_ranks, ntotal);
}




int main(int argc, char** argv)
{
  int num_ranks, my_rank;
  int method = METHOD_DEFAULT;
  int ntotal = NUMPARTICLES_DEFAULT;
  int per[3] = {false, false, false};
  fcs_float boxa[3]   = {1., 0., 0.};
  fcs_float boxb[3]   = {0., 1., 0.};
  fcs_float boxc[3]   = {0., 0., 1.};
  fcs_float virial[9];
  double energy_local  = 0.;
  double energy_global = 0.;
  double trvirial      = 0.;
  double trvirial2     = 0.;
  double trvirial2l    = 0.;
  int setup_id         = SETUP_DEFAULT;
  MPI_Comm comm;


  int nparts, max_nparts;
  double *pos, *qs, *es, *pot;

  FCSResult fcs_res;
  FCS fcs_handle;

  MPI_Init(&argc, &argv);

  MPI_Comm_size(MPI_COMM_WORLD, &num_ranks);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  readparams(argc, argv, &method, &ntotal, &setup_id, per);
  // Particle setup
  max_nparts = MAXPARTS(ntotal/num_ranks + 1);

  pos = (fcs_float*) malloc(max_nparts*3*sizeof(fcs_float));
  qs  = (fcs_float*) malloc(max_nparts*1*sizeof(fcs_float));
  es  = (fcs_float*) malloc(max_nparts*3*sizeof(fcs_float));
  pot = (fcs_float*) malloc(max_nparts*1*sizeof(fcs_float));

  setup_particles(&my_rank, &num_ranks, &setup_id, &ntotal, &nparts, &max_nparts, pos, qs, es, pot, boxa, boxb, boxc, &comm);

  printparams(argv[0], num_ranks, my_rank, method, ntotal, setup_id, per);

  ////////////////////////////////////////////////////////////////////////////////////////////////
  // ScaFaCoS-Initialization
  fcs_res = fcs_init( &fcs_handle, fcs_method[method], comm);
  checkres(fcs_res);
  // ScaFaCoS: Generic Parameter Setup
  fcs_res = fcs_set_common( fcs_handle, true, boxa, boxb, boxc, per, ntotal );
  checkres(fcs_res);
  // ScaFaCoS: Method-specific Parameter Setup
  setup_methodspecific( fcs_handle, my_rank, num_ranks);
  checkres(fcs_res);
  // ScaFaCos: tell the library that we are interested in the virial tensor
  fcs_res = fcs_require_virial( fcs_handle, 1);
  checkres(fcs_res);
  // ScaFaCoS: Method-specific Internal Tuning
  fcs_res = fcs_tune( fcs_handle, nparts, max_nparts, pos, qs );
  checkres(fcs_res);
  // ScaFaCoS: Output of internal Parameters
  if (my_rank == 0) fcs_printContent(fcs_handle);
  // ScaFaCoS: Invocation of Coulomb Solver
  fcs_res = fcs_run( fcs_handle, nparts, max_nparts, pos, qs, es, pot);
  checkres(fcs_res);
  // ScaFaCos: read virial from results
  fcs_res = fcs_get_virial( fcs_handle, virial);
  checkres(fcs_res);
  // ScaFaCoS: Free any allocated objects
  fcs_res = fcs_destroy( fcs_handle );
  checkres(fcs_res);
  ////////////////////////////////////////////////////////////////////////////////////////////////

  // Evaluation of results
  {
      int i;
      for (i=0; i<nparts; i++)
      {
          energy_local +=   pot[i]*qs[i]/2.;
          trvirial2l   += ( pos[3*i+0]*es[3*i+0]
                          + pos[3*i+1]*es[3*i+1]
                          + pos[3*i+2]*es[3*i+2] ) * qs[i];
      }

      energy_local /= 1.*ntotal;

      MPI_Allreduce(&energy_local, &energy_global, 1, MPI_DOUBLE, MPI_SUM, comm);

      for (i=0; i<3; i++)
          trvirial += virial[4*i];

      trvirial  /= (double)ntotal;
      trvirial2l/= (double)ntotal;

      MPI_Allreduce(&trvirial2l, &trvirial2, 1, MPI_DOUBLE, MPI_SUM, comm);
  }
  // Output of results
  {
    int i;
    if (my_rank ==0) {
      printf("\n\n******* debug output: energy and trace of virial ******\n\n");
      printf("energy/particle:          %14.8g\n", energy_global);
      printf("trace of virial/particle: %14.8g\n", trvirial);
      printf("dito (from forces):       %14.8g\n", trvirial2);

      printf("\n\n******* debug output: particles and fields ************\n\n");
      printf("[Rank] Part. | (      x     ,     y      ,     z      )       q       |  (     Ex     ,     Ey     ,     Ez     )     pot     \n");
      printf("-------------+--------------------------------------------------------+-------------------------------------------------------\n");
      fflush(stdout);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    for(i=0; i < min(nparts,MAX_PART_PRINT/num_ranks); i++){
      printf("[%4.4d] %5.5d | (%12.3g,%12.3g,%12.3g) %12.3g  |  (%12.3g,%12.3g,%12.3g) %12.3g\n",
                 my_rank, i, pos[3*i+0], pos[3*i+1], pos[3*i+2], qs[i],
                 es[3*i+0], es[3*i+1], es[3*i+2], pot[i]);
    }

  }
  
  {
    MPI_Finalize();
    return 0;
  }

}
