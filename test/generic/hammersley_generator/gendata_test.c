/*
 * Author: Toni Volkmer
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "gendata_hammersley.h"
#include "gendata_grid.h"
#include "gendata_plummer.h"
#include "gendata_util.h"

/* Prints positions and charges to stdout */
void print_data(int d, int n, double *x, double *q)
{
  int j;
  int t;

  for (j = 0; j < n; j++)
  {
    for (t = 0; t < d; t++)
      printf("%g ", x[d*j+t]);
    printf("%g\n", q[j]);
  }
}

/** Prints positions and charges to stdout in new libfcs input format.
  * First line contains particle numbers, second line reserved for comments.
  * Then each line consists of positions x y z and charge.
  */
void print_data_fcs(int d, int n, double *x, double *q)
{
  int j;
  int t;

  printf("%d\n", n);
  printf("\n");
  for (j = 0; j < n; j++)
  {
    for (t = 0; t < d; t++)
      printf("%g ", x[d*j+t]);
    printf("%g\n", q[j]);
  }
}

/* Prints positions and charges to file */
void print_data_to_file(char *filename, int d, int n, double *x, double *q)
{
  int j;
  int t;
  FILE *outfile = NULL;

  outfile = fopen(filename, "w");
  if (outfile <= 0)
    return;
  for (j = 0; j < n; j++)
  {
    for (t = 0; t < d; t++)
      fprintf(outfile, "%g ", x[d*j+t]);
    fprintf(outfile, "%g\n", q[j]);
  }
  fclose(outfile);

}

void show_usage_exit(char *pname)
{
  fprintf(stderr, "n means number of particles\n");
  fprintf(stderr, "charge_type: 0 for alternating charges, 1 for positive charges\n"); 
  fprintf(stderr, "usage:\n");
  fprintf(stderr, "%s hammersley_cube n charge_type\n", pname);
  fprintf(stderr, "%s hammersley_ball n charge_type\n", pname);
  fprintf(stderr, "%s hammersley_two_balls n charge_type n2_rel_to_n1 dist_rel_to_r1\n", pname);
  fprintf(stderr, "%s hammersley_sphere n charge_type\n", pname);
  fprintf(stderr, "%s hammersley_circle n charge_type\n", pname);
  fprintf(stderr, "%s hammersley_square n charge_type\n", pname);
  fprintf(stderr, "%s halton_ellipsoid n charge_type a b c\n", pname);
  fprintf(stderr, "%s halton_cylinder n charge_type ratio_radius_length\n", pname);
  fprintf(stderr, "%s grid_face_centered_cube n charge_type\n", pname);
  fprintf(stderr, "%s grid_body_centered_cube n charge_type\n", pname);
  fprintf(stderr, "%s grid_nacl_cube n\n", pname);
  fprintf(stderr, "%s grid_nacl_cube_periodic n\n", pname);
  fprintf(stderr, "%s plummer_ball n charge_type\n", pname);
  fprintf(stderr, "\n");
  exit(1);
}

/** Sets charges to *q.
  * Use 0 for alternating charges (-1 and 1), where sum of all charges is
  * in {-1, 0, 1}.
  * Use 1 for positive charges 1.
  * Use -1 for negative charges -1.
  */
void set_charges(char *argv0, int type, int n, double **q)
{
  *q = (double *) malloc(n*sizeof(double));
  if (type == 0)
  {
    set_random_charge(n, *q);
  }
  else if (type == 1)
  {
    set_positive_charge(n, *q);
  }
  else if (type == -1)
  {
    set_negative_charge(n, *q);
  }
  else
  {
    printf("Charge type must be one of {-1, 0, 1}\n");
    show_usage_exit(argv0);
  }
}

/** IMPORTANT: The total number of particles can be different from the
  * suggested number "n" specified.
  * (e.g. in grid_face_centered_cube, grid_body_centered_cube, grid_nacl_cube
  *  AT LEAST "n" particles are generated)
  */
int main(int argc, char **argv)
{
  double *x, *q;
  int d, n, n_component, j;
  double temp;
  int charge_default_type = 0; /* alternating charges, sum = 0 */
  double n2_rel_to_n1 = 0.015625; /* for hammersley_two_balls */
  double dist_rel_to_r1 = 10; /* for hammersley_two_balls */
  double a = 3.0; /* for halton_ellipsoid */
  double b = 2.0; /* for halton_ellipsoid */
  double c = 1.0; /* for halton_ellipsoid */
  double ratio_radius_length = 0.1; /* for halton_cylinder */
  
  if (argc < 3)
  {
    show_usage_exit(argv[0]);
  }

  n = atoi(argv[2]);
  if (n < 1)
  {
    fprintf(stderr, "n must be >= 1\n");
    show_usage_exit(argv[0]);
  }

  if(strcmp(argv[1], "hammersley_cube")==0)
  {
    d = 3;
    gendata_hammersley_cube(n, &x, &n);

    if (argc < 4)
      set_charges(argv[0], charge_default_type, n, &q);
    else
      set_charges(argv[0], atoi(argv[3]), n, &q);
  }
  else if (strcmp(argv[1], "hammersley_ball")==0)
  {
    d = 3;
    gendata_hammersley_ball(n, &x, &n);

    if (argc < 4)
      set_charges(argv[0], charge_default_type, n, &q);
    else
      set_charges(argv[0], atoi(argv[3]), n, &q);
  }
  else if (strcmp(argv[1], "hammersley_two_balls")==0)
  {
    d = 3;

    if (argc > 4)
      n2_rel_to_n1 = atof(argv[4]);

    if (argc > 5)
      dist_rel_to_r1 = atof(argv[5]);


    gendata_hammersley_two_balls(n, &x, n2_rel_to_n1, dist_rel_to_r1, &n);
    auto_shift_scale(d, n, x);

    if (argc < 4)
      set_charges(argv[0], charge_default_type, n, &q);
    else
      set_charges(argv[0], atoi(argv[3]), n, &q);
  }
  else if (strcmp(argv[1], "hammersley_sphere")==0)
  {
    d = 3;
    gendata_hammersley_sphere(n, &x, &n);

    if (argc < 4)
      set_charges(argv[0], charge_default_type, n, &q);
    else
      set_charges(argv[0], atoi(argv[3]), n, &q);
  }
  else if (strcmp(argv[1], "hammersley_circle")==0)
  {
    d = 3;
    gendata_hammersley_circle(n, &x, &n);

    if (argc < 4)
      set_charges(argv[0], charge_default_type, n, &q);
    else
      set_charges(argv[0], atoi(argv[3]), n, &q);
  }
  else if (strcmp(argv[1], "hammersley_square")==0)
  {
    d = 3;
    gendata_hammersley_square(n, &x, &n);

    if (argc < 4)
      set_charges(argv[0], charge_default_type, n, &q);
    else
      set_charges(argv[0], atoi(argv[3]), n, &q);
  }
  else if (strcmp(argv[1], "grid_face_centered_cube")==0)
  {
    d = 3;
    temp = pow(n / 4.0, 1.0/3.0);
    n_component = ceil(temp);

    gendata_grid_face_centered_cube(n_component, &x, &n);

    if (argc < 4)
      set_charges(argv[0], charge_default_type, n, &q);
    else
      set_charges(argv[0], atoi(argv[3]), n, &q);
  }
  else if (strcmp(argv[1], "grid_body_centered_cube")==0)
  {
    d = 3;
    temp = pow(n / 2.0, 1.0/3.0);
    n_component = ceil(temp);

    gendata_grid_body_centered_cube(n_component, &x, &n);

    if (argc < 4)
      set_charges(argv[0], charge_default_type, n, &q);
    else
      set_charges(argv[0], atoi(argv[3]), n, &q);
  }
  else if (strcmp(argv[1], "grid_nacl_cube")==0)
  {
    d = 3;
    temp = pow(n, 1.0/3.0);
    n_component = ceil(temp);

    gendata_nacl_cube(n_component, &x, &q, &n);
  }
  else if (strcmp(argv[1], "grid_nacl_cube_periodic")==0)
  {
    d = 3;
    temp = pow(n, 1.0/3.0);
    n_component = ceil(temp);

    gendata_nacl_cube_periodic(n_component, &x, &q, &n);
  }
  else if (strcmp(argv[1], "halton_ellipsoid")==0)
  {
    d = 3;

    if (argc > 4)
      a = atof(argv[4]);

    if (argc > 5)
      b = atof(argv[5]);

    if (argc > 6)
      c = atof(argv[6]);

    gendata_halton_ellipsoid(n, &x, a, b, c, &n);

    if (argc < 4)
      set_charges(argv[0], charge_default_type, n, &q);
    else
      set_charges(argv[0], atoi(argv[3]), n, &q);
  }
  else if (strcmp(argv[1], "halton_cylinder")==0)
  {
    d = 3;

    if (argc > 4)
      ratio_radius_length = atof(argv[4]);

    gendata_halton_cylinder(n, &x, ratio_radius_length, &n);

    if (argc < 4)
      set_charges(argv[0], charge_default_type, n, &q);
    else
      set_charges(argv[0], atoi(argv[3]), n, &q);
  }
  else if (strcmp(argv[1], "plummer_ball")==0)
  {
    d = 3;

    gendata_plummer_ball(n, &x, &n);

    if (argc < 4)
      set_charges(argv[0], charge_default_type, n, &q);
    else
      set_charges(argv[0], atoi(argv[3]), n, &q);
  }
  else
  {
    fprintf(stderr, "Unknown distribution\n");
    show_usage_exit(argv[0]);
  }

  print_data_fcs(d, n, x, q);

  free(q);
  free(x);

  return 0;
}
