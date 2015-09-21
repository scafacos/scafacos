
#include <stdio.h>
#include <stdlib.h>

#include <p2nfft.h>
#include <fcs.h>
#include <FCSCommon.h>

#define TEST_BOX_SIZE 2 
#define TEST_N_PARTICLES (TEST_BOX_SIZE * TEST_BOX_SIZE * TEST_BOX_SIZE)

static void assert_fcs(FCSResult r)
{
  if (r) {
    fcs_result_print_result(r);
    MPI_Finalize();
    exit(-1);
  }
}

int main (int argc, char** argv) {
  fcs_int num_particles = TEST_N_PARTICLES;
  fcs_float box_size = TEST_BOX_SIZE;
  fcs_int i, px, py, pz;
  fcs_float positions[3*TEST_N_PARTICLES];
  fcs_float charges[TEST_N_PARTICLES];
  fcs_float field[3*TEST_N_PARTICLES];
  fcs_float potentials[TEST_N_PARTICLES];
  MPI_Init(&argc, &argv);
  MPI_Comm comm = MPI_COMM_WORLD;

  i = 0;
  for (px = 0; px < TEST_BOX_SIZE; ++px) {
    for (py = 0; py < TEST_BOX_SIZE; ++py) {
      for (pz = 0; pz < TEST_BOX_SIZE; ++pz) {
        positions[3*i] = px + 0.5;
        positions[3*i + 1] = py + 0.5;
        positions[3*i + 2] = pz + 0.5;
        charges[i] = 1.0-((px + py + pz) % 2)*2;
        ++i;
      }
    }
  }

  /* Calculate this system via FCS direct solver */
  fcs_int near_field_flag = 1;
  fcs_float box_a[] = { box_size, 0.0, 0.0 };
  fcs_float box_b[] = { 0.0, box_size, 0.0 };
  fcs_float box_c[] = { 0.0, 0.0, box_size };
  fcs_float offset[] = {0.0, 0.0, 0.0};
  fcs_int periodicity[] = {1, 1, 1};

  FCS fcs_handle;
  FCSResult fcs_result;

  fcs_result = fcs_init(&fcs_handle, "DIRECT", comm);
  assert_fcs(fcs_result);
  fcs_result = fcs_set_common(fcs_handle, near_field_flag, box_a, box_b, box_c, offset, periodicity, num_particles);
  assert_fcs(fcs_result);
  fcs_result = fcs_tune(fcs_handle, num_particles, positions, charges);
  assert_fcs(fcs_result);
  fcs_result = fcs_run(fcs_handle, num_particles, positions, charges, field, potentials);
  assert_fcs(fcs_result);

  printf("Potentials via FCS DIRECT:\n");
  printf("%" FCS_LMOD_FLOAT "f\n", potentials[0]);
//   for (i = 0; i < num_particles; ++i)
//     printf("%" FCS_LMOD_FLOAT "f\n", potentials[i]);
  printf("\n");
  fcs_destroy(fcs_handle);

  /* Calculate this system via P2NFFT */
  fcs_int k;
  fcs_float tolerance = 1.0;
  void* rd = NULL;
  ifcs_p2nfft_init(&rd, comm);
  for (k = 0; k < 8; ++k) {
    printf("===================================\n");
    printf("Trying tolerance %" FCS_LMOD_FLOAT "f.\n", tolerance);
    fcs_p2nfft_set_tolerance(rd, FCS_TOLERANCE_TYPE_FIELD, tolerance);
    ifcs_p2nfft_tune(rd, periodicity, num_particles, positions, charges, box_a, box_b, box_c, offset, near_field_flag);
    ifcs_p2nfft_run(rd, num_particles, num_particles, positions, charges, potentials, field, 0, 0, NULL, NULL, NULL, NULL);

    printf("Potentials via P2NFFT\n");
//    for(i = 0; i<num_particles; ++i) {
      printf("%.10" FCS_LMOD_FLOAT "f\n", potentials[0]);
//    }
    tolerance /= 10.0;
  }
  ifcs_p2nfft_destroy(rd);
  MPI_Finalize();
  return 0;
}
