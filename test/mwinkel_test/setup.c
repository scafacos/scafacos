#include <fcs.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "setup.h"
#include <mpi.h>

void setup_trivial(int* my_rank, int* num_ranks, fcs_int* ntotal, fcs_int* nparts, fcs_int* max_nparts,
                      fcs_float* pos, fcs_float* qs, fcs_float* es, fcs_float* pot, fcs_float boxa[3], fcs_float boxb[3], fcs_float boxc[3],
                      MPI_Comm* comm);

void setup_crystal(int* my_rank, int* num_ranks, fcs_int* ntotal, fcs_int* nparts, fcs_int* max_nparts,
                      fcs_float* pos, fcs_float* qs, fcs_float* es, fcs_float* pot, fcs_float boxa[3], fcs_float boxb[3], fcs_float boxc[3],
                      MPI_Comm* comm,
                      fcs_float cell_sz[3], fcs_float cell_chg[], fcs_float cell_pos[][3], int natoms);

void setup_nacl(int* my_rank, int* num_ranks, fcs_int* ntotal, fcs_int* nparts, fcs_int* max_nparts,
                      fcs_float* pos, fcs_float* qs, fcs_float* es, fcs_float* pot, fcs_float boxa[3], fcs_float boxb[3], fcs_float boxc[3],
                      MPI_Comm* comm);

void setup_sio2(int* my_rank, int* num_ranks, fcs_int* ntotal, fcs_int* nparts, fcs_int* max_nparts,
                      fcs_float* pos, fcs_float* qs, fcs_float* es, fcs_float* pot, fcs_float boxa[3], fcs_float boxb[3], fcs_float boxc[3],
                      MPI_Comm* comm);


void setup_particles(int* my_rank, int* num_ranks, int* id, fcs_int* ntotal, fcs_int* nparts, fcs_int* max_nparts,
             fcs_float* pos, fcs_float* qs, fcs_float* es, fcs_float* pot, fcs_float boxa[3], fcs_float boxb[3], fcs_float boxc[3],
             MPI_Comm* comm)
{
  switch(*id)
  {
  case (SETUP_TRIVIAL):
          setup_trivial(my_rank, num_ranks, ntotal, nparts, max_nparts, pos, qs, es, pot, boxa, boxb, boxc, comm);
          break;
  case (SETUP_NACL):
          setup_nacl(my_rank, num_ranks, ntotal, nparts, max_nparts, pos, qs, es, pot, boxa, boxb, boxc, comm);
          break;
  case (SETUP_SiO2):
          setup_sio2(my_rank, num_ranks, ntotal, nparts, max_nparts, pos, qs, es, pot, boxa, boxb, boxc, comm);
          break;
  default:
      printf("Unknown setup ID: %d.", id);
      exit(1);
  }

}

void setup_particles_(int* my_rank, int* num_ranks, int* id, fcs_int* ntotal, fcs_int* nparts, fcs_int* max_nparts,
             fcs_float* pos, fcs_float* qs, fcs_float* es, fcs_float* pot, fcs_float boxa[3], fcs_float boxb[3], fcs_float boxc[3],
             MPI_Comm* comm)
{
    setup_particles(my_rank, num_ranks, id, ntotal, nparts, max_nparts,
                 pos, qs, es, pot, boxa, boxb, boxc,
                 comm);
}

void setup_trivial(int* my_rank, int* num_ranks, fcs_int* ntotal, fcs_int* nparts, fcs_int* max_nparts,
                      fcs_float* pos, fcs_float* qs, fcs_float* es, fcs_float* pot, fcs_float boxa[3], fcs_float boxb[3], fcs_float boxc[3],
                      MPI_Comm* comm)
{
  int i;

  *nparts = *ntotal / *num_ranks;
  if (*my_rank < *ntotal%*num_ranks) *nparts++;

  for(i=0; i < *nparts; i++){
      pos[3*i + 0] = ((fcs_float)i       ) / ((fcs_float)*nparts);
      pos[3*i + 1] = ((fcs_float)*my_rank) / ((fcs_float)*num_ranks);
      pos[3*i + 2] = sin(((fcs_float)i) / ((fcs_float)*nparts) + ((fcs_float)*my_rank) / ((fcs_float)*num_ranks) ) / 2. + 0.5 ;

    if (*my_rank == 0 && i == 0) {
      qs[i] = -(*ntotal-1);
    } else {
      qs[i] = +1.0;
    }

    es[3*i + 0] = es[3*i + 1] = es[3*i + 2] = pot[i] = -12.3;
  }

  for(i=0;i<3;i++)
      boxa[i] = boxb[i] = boxc[i] = 0.;

  boxa[0] = 1.;
  boxb[1] = 1.;
  boxc[2] = 1.;

  *comm = MPI_COMM_WORLD;
}


void setup_crystal(int* my_rank, int* num_ranks, fcs_int* ntotal, fcs_int* nparts, fcs_int* max_nparts,
                      fcs_float* pos, fcs_float* qs, fcs_float* es, fcs_float* pot, fcs_float boxa[3], fcs_float boxb[3], fcs_float boxc[3],
                      MPI_Comm* comm,
                      fcs_float cell_sz[3], fcs_float cell_chg[], fcs_float cell_pos[][3], int natoms)
{
    int nx, ny, nz;
    int i, j, k;
    int coord[3], dims[3] = {0,1,1}, period[3] = {1,1,1};

    MPI_Dims_create(*num_ranks, 3, dims);
    MPI_Cart_create(MPI_COMM_WORLD, 3, dims, period, 1, comm);
    MPI_Comm_rank(*comm, my_rank);
    MPI_Cart_coords(*comm, *my_rank, 3, coord);

    nx = ny = nz = (int)(pow(1.* *ntotal/natoms,1./3.));

    for(i=0;i<3;i++)
        boxa[i] = boxb[i] = boxc[i] = 0.;

    boxa[0] = nx*cell_sz[0];
    boxb[1] = ny*cell_sz[1];
    boxc[2] = nz*cell_sz[2];

    *ntotal = nx * ny * nz * natoms;
    *nparts = *ntotal / *num_ranks;

    /* create crystal */
    { int i, j, k, l, chunk, cnt;
      chunk = nx / *num_ranks;
      cnt   = 0;
      for (i=coord[0]*chunk; i<(coord[0]+1)*chunk; ++i)
        for (j=0; j<ny; ++j)
          for (k=0; k<nz; ++k)
            for (l=0; l<natoms; ++l) {
              pos[3*cnt  ] = cell_pos[l][0] + i*cell_sz[0];
              pos[3*cnt+1] = cell_pos[l][1] + j*cell_sz[1];
              pos[3*cnt+2] = cell_pos[l][2] + k*cell_sz[2];
              qs [cnt++  ] = cell_chg[l];
      }
      printf("%d %d %d\n", *my_rank, *nparts, cnt);
    }
}



#define NUM_NACL 8

void setup_nacl(int* my_rank, int* num_ranks, fcs_int* ntotal, fcs_int* nparts, fcs_int* max_nparts,
                      fcs_float* pos, fcs_float* qs, fcs_float* es, fcs_float* pot, fcs_float boxa[3], fcs_float boxb[3], fcs_float boxc[3],
                      MPI_Comm* comm)
{
  fcs_float cell_sz [3]   = { 2.0, 2.0, 2.0 };
  fcs_float cell_chg[NUM_NACL] = { 1.0, -1.0, -1.0, 1.0, -1.0, 1.0, 1.0, -1.0 };
  fcs_float cell_pos[NUM_NACL][3] = {
    { 0.5, 0.5, 0.5 }, { 0.5, 0.5, 1.5 }, { 0.5, 1.5, 0.5 }, { 0.5, 1.5, 1.5 },
    { 1.5, 0.5, 0.5 }, { 1.5, 0.5, 1.5 }, { 1.5, 1.5, 0.5 }, { 1.5, 1.5, 1.5 }   };

  setup_crystal(my_rank, num_ranks, ntotal, nparts, max_nparts, pos, qs, es, pot, boxa, boxb, boxc, comm, cell_sz, cell_chg, cell_pos, NUM_NACL);
}


#define NUM_SIO2 18

void setup_sio2(int* my_rank, int* num_ranks, fcs_int* ntotal, fcs_int* nparts, fcs_int* max_nparts,
                      fcs_float* pos, fcs_float* qs, fcs_float* es, fcs_float* pot, fcs_float boxa[3], fcs_float boxb[3], fcs_float boxc[3],
                      MPI_Comm* comm)
{
  fcs_float cell_sz [3]   = { 4.9134, 8.51025844, 5.4052 };
  fcs_float cell_chg[NUM_SIO2] = { 2.4, 2.4, 2.4, 2.4, 2.4, 2.4,-1.2,-1.2,-1.2,
                             -1.2,-1.2,-1.2,-1.2,-1.2,-1.2,-1.2,-1.2,-1.2 };
  fcs_float cell_pos[NUM_SIO2][3] = {
    { 0.677893, 5.145130, 0.900000 }, { 3.134590, 0.890000, 0.900000 },
    { 1.684400, 2.889490, 2.701730 }, { 4.141100, 7.144610, 2.701730 },
    { 1.684400, 7.400770, 4.503470 }, { 4.141100, 3.145640, 4.503470 },
    { 4.067400, 8.259460, 1.541777 }, { 1.610700, 4.004330, 1.541777 },
    { 2.205960, 1.511250, 2.059960 }, { 4.662660, 5.766380, 2.059960 },
    { 0.230040, 2.652050, 3.343510 }, { 2.686740, 6.907180, 3.343510 },
    { 2.686740, 3.383080, 3.861690 }, { 0.230040, 7.638210, 3.861690 },
    { 2.205960, 0.268752, 5.145240 }, { 4.662660, 4.523880, 5.145240 },
    { 1.610700, 6.285930, 0.258220 }, { 4.067400, 2.030800, 0.258220 }  };

  setup_crystal(my_rank, num_ranks, ntotal, nparts, max_nparts, pos, qs, es, pot, boxa, boxb, boxc, comm, cell_sz, cell_chg, cell_pos, NUM_SIO2);
}
