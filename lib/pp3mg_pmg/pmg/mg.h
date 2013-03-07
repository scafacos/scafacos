#ifndef _MG__H_
#define _MG__H_

#include "mpi.h"

typedef struct {
  double ***v;
  double ***f;
  double ***r;
  double ***e;
  double ***tmp;

  /* buffers for send */
  double ***sbufxy;
  double ***sbufxz;
  double ***sbufyz;

  /* buffers for receive */
  double ***rbufxy;
  double ***rbufxz;
  double ***rbufyz;

  /* global dimensions */
  int m;
  int n;
  int o;

  /* local dimensions */
  int m_l;
  int n_l;
  int o_l;

  /* periodic? */
  int periodic;

  /* cartesian communicator */
  MPI_Comm cart_comm;

  int x_off;
  int y_off;
  int z_off;

  /* neighbors */
  int left;
  int right;
  int lower;
  int upper;
  int back;
  int front;

  /* ghosts */
  int x_ghosts;
  int y_ghosts;
  int z_ghosts;

  /* number of pre and post smoothing steps */
  int nu1;
  int nu2;

  /* position of zero */
  _Bool zero_at_pi3;

  /* relaxation coefficient */
  double omega;

  /* stencil data */
  int size;
  double* values;
  int* x_offsets;
  int* y_offsets;
  int* z_offsets;
} mg_data;

void mg_setup( mg_data **outdata, int maxlevel, int m, int n, int o,
	       int xstart, int xend, int ystart, int yend, int zstart, int zend,
	       int p, int nu1, int nu2, double omega, int size, double* values, 
	       int* xoff, int* yoff, int* zoff, MPI_Comm cart_comm);

void mg_init(double ***v, double ***f, mg_data *data );

void mg_free(mg_data *data, int maxlevel);

double mg_vcycle( mg_data *data, int level, int maxlevel );

double mg(double ***u, double ***f, int maxiter, double tol, int m, int n, int o,
	  int xstart, int xend, int ystart, int yend, int zstart, int zend,
	  int p, int nu1, int nu2, double omega, int size, double* values,
	  int* xoff, int* yoff, int* zoff, MPI_Comm cart_comm, int verbose);

#endif /* ifndef _MG__H_ */
