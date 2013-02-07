/*
 * Author: Toni Volkmer
 *
 */


#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

extern "C" {
  #include "gendata_util.h"
}

#include "hammersley.H"

extern "C" {
  int seed[3];
  int leap[3];
  int base[3];

  void set_hammersley_default_params(int d, int n)
  {
    int j;

    assert(d > 0 && d <= 3);

    seed[0] = 0;
    leap[0] = 1;
    base[0] = -n;
    for (j = 1; j < d; j++)
    {
      seed[j] = 0;
      leap[j] = 1;
      base[j] = get_prime(j);
    }

    hammersley_dim_num_set(d);
    hammersley_seed_set(seed);
    hammersley_leap_set(leap);
    hammersley_base_set(base); 
  }

  void set_halton_default_params(int d)
  {
    int j;

    assert(d > 0 && d <= 3);

    for (j = 0; j < d; j++)
    {
      seed[j] = 0;
      leap[j] = 1;
      base[j] = get_prime(j+1);
    }

    hammersley_dim_num_set(d);
    hammersley_seed_set(seed);
    hammersley_leap_set(leap);
    hammersley_base_set(base);
  }

  void gendata_hammersley_cube(int n, double **x, int *n_total)
  {
    *x = (double *) malloc(3*n*sizeof(double));
    set_hammersley_default_params(3, n);
    hammersley_sequence(n, *x);
    *n_total = n;
  }

  void gendata_hammersley_ball(int n, double **x, int *n_total)
  {
    int j;

    *x = (double *) malloc(3*n*sizeof(double));

    set_hammersley_default_params(3, n);
    hammersley_sequence(n, *x);

    for (j = 0; j < n; j++)
      u3_to_ball_unit_3d((*x)+3*j, (*x)+3*j);

    shift_scale(3, n, *x, 1, 0.5);

    *n_total = n;
  }

  void gendata_hammersley_two_balls(int n, double **x, double n_2_rel_n_1,
                                    double distance_rel_r_1, int *n_total)
  {
    double dist; 
    double *x_2;
    int j;
    double r_2;
    double r_1 = 1.0;
    int n_1;
    int n_2;
    int diff;

    dist = r_1 * distance_rel_r_1;
    
    r_2 = pow(n_2_rel_n_1, 1.0/3.0) * r_1;
    n_1 = (int) ceil(n * (1.0-n_2_rel_n_1));
    n_2 = (int) ceil(n * n_2_rel_n_1);

    if (n_1 + n_2 > n)
    {
      diff = (n_1 + n_2) - n;
      if (diff <= n_1)
        n_1 -= diff;
    }

    *n_total = n_1 + n_2;

    *x = (double *) malloc(3 * (*n_total) * sizeof(double));

    set_hammersley_default_params(3, n_1);
    hammersley_sequence(n_1, *x);

    for (j = 0; j < n_1; j++)
      u3_to_ball_unit_3d((*x)+3*j, (*x)+3*j);

    if (n_2 > 0)
    {
      x_2 = (double *) malloc(3*n_2*sizeof(double));

      set_hammersley_default_params(3, n_2);
      hammersley_sequence(n_2, x_2);

      for (j = 0; j < n_2; j++)
      {
        u3_to_ball_unit_3d(x_2+3*j, x_2+3*j);
        (*x)[3*(n_1+j)] = x_2[3*j] * r_2 + dist;
        (*x)[3*(n_1+j)+1] = x_2[3*j+1] * r_2;
        (*x)[3*(n_1+j)+2] = x_2[3*j+2] * r_2;
      }

      free(x_2);
    }

  }

  void gendata_halton_ellipsoid(int n, double **x, double a, double b,
                                double c, int *n_total)
  {
    double max = a;
    int j = 0;
    double u[3];

    if (max < b)
      max = b;

    if (max < c)
      max = c;

    a = a / max;
    b = b / max;
    c = c / max;

    *x = (double *) malloc(3*n*sizeof(double));

    set_halton_default_params(3);

    while (j < n)
    {
      hammersley(u);
      u[0] -= 0.5;
      u[1] -= 0.5;
      u[2] -= 0.5;
      if (u[0]*u[0]/a/a + u[1]*u[1]/b/b + u[2]*u[2]/c/c > 0.25)
        continue;

      (*x)[3*j] = u[0]+0.5;
      (*x)[3*j+1] = u[1]+0.5;
      (*x)[3*j+2] = u[2]+0.5;
      j++;
    }

    *n_total = n;
  }

  void gendata_halton_cylinder(int n, double **x, double r_div_len,
                           int *n_total)
  {
    int j = 0;
    double u[3];
    double r, len;

    *x = (double *) malloc(3*n*sizeof(double));

    if (r_div_len >= 0.5)
    {
      r = 0.5;
      len = r / r_div_len;
    }
    else
    {
      len = 1.0;
      r = len * r_div_len;
    }


    set_halton_default_params(3);

    while (j < n)
    {
      hammersley(u);
      u[0] -= 0.5;
      u[1] -= 0.5;
      u[2] -= 0.5;

      if (2*fabs(u[0]) > len || u[1]*u[1] + u[2]*u[2] > r*r)
        continue;

      (*x)[3*j] = u[0]+0.5;
      (*x)[3*j+1] = u[1]+0.5;
      (*x)[3*j+2] = u[2]+0.5;
      j++;
    }

    *n_total = n;
  }

  void gendata_hammersley_ball_neg_charge(int n, double **x, int *n_total)
  {
    int j;

    *x = (double *) malloc(3*n*sizeof(double));

    set_hammersley_default_params(3, n);
    hammersley_sequence(n, *x);

    for (j = 0; j < n; j++)
      u3_to_ball_unit_3d((*x)+3*j, (*x)+3*j);

    shift_scale(3, n, *x, 1, 0.5);

    *n_total = n;
  }

  void gendata_hammersley_sphere(int n, double **x, int *n_total)
  {
    double u[2*n];
    int j;

    *x = (double *) malloc(3*n*sizeof(double));

    set_hammersley_default_params(2, n);
    hammersley_sequence(n, u);

    for (j = 0; j < n; j++)
      u2_to_sphere_unit_3d(u+2*j, (*x)+3*j);

    shift_scale(3, n, *x, 1, 0.5);

    *n_total = n;
  }

  void gendata_hammersley_circle(int n, double **x, int *n_total)
  {
    double u[2*n];
    int j;

    *x = (double *) malloc(3*n*sizeof(double));

    set_hammersley_default_params(2, n);
    hammersley_sequence(n, u);

    for (j = 0; j < n; j++)
      u2_to_ball_unit_2d(u+2*j, u+2*j);

    for (j = 0; j < n; j++)
    {
      (*x)[3*j] = u[2*j] * 0.5 + 0.5;
      (*x)[3*j+1] = u[2*j+1] * 0.5 + 0.5;
      (*x)[3*j+2] = 0.5;
    }

    *n_total = n;
  }

  void gendata_hammersley_square(int n, double **x, int *n_total)
  {
    double u[2*n];
    int j;

    *x = (double *) malloc(3*n*sizeof(double));

    set_hammersley_default_params(2, n);
    hammersley_sequence(n, u);

    for (j = 0; j < n; j++)
    {
      (*x)[3*j] = u[2*j];
      (*x)[3*j+1] = u[2*j+1];
      (*x)[3*j+2] = 0.5;
    }

    *n_total = n;
  }
}
