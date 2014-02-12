#include "types.hpp"
#include "Communication.hpp"
#include "ErrorEstimate.hpp"
#include <cstdlib>
#include <cmath>

using namespace P3M;

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);

  if (argc <= 1) {
    printf("Usage: compute_error <grid>\n");
    return 2;
  }
  p3m_int grid = atoi(argv[1]);

  Communication comm(MPI_COMM_WORLD);
  ErrorEstimate *error = ErrorEstimate::create(comm);
  p3m_float L[3] = { 10.0, 10.0, 10.0 };
  p3m_float num_charges = 300.0;
  p3m_float sum_q2 = 300.0;

  Parameters p;
  p.cao = 4;
  p.grid[0] = grid;
  p.grid[1] = grid;
  p.grid[2] = grid;

  for (p.alpha = 0.1; p.alpha < 10.0; p.alpha += 0.1) {
	  p3m_float ks_error = error->compute_ks_error(p, num_charges, sum_q2, L);

	  p3m_float alpha_L = p.alpha*L[0];
	  p3m_float Q = sqrt(ks_error/(L[0]*L[0]));
	  printf("%le %le %le %le\n", p.alpha, ks_error, alpha_L, Q);
  }

  MPI_Finalize();
}
