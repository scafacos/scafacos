/*
 * Author: Toni Volkmer
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <float.h>

#include "gendata_util.h"
#include "gendata_plummer.h"

#define PI 4*atan(1)
#define PLUM_NORM 0.033510321638291128 /*See Maple worksheet */

/* From libfcs, generateInput, hilfsfunktionen.c */
/*Int 4*PI*p(x)*x*x, ausgerechnet mit MAPLE*/
double int_plummer(double r) {
  return 4*PI*pow(r,3) / (3* pow(1.0+25*r*r,1.5)) / PLUM_NORM;
}

void gendata_plummer_ball(int n, double **x, int *n_total)
{
  int j, fib_i, fib_j, i_ctrl;
  double fib_x[17];
  double rnd_seed;
  double step = 1.0 / n;
  double b_left, b_right, integral_left, integral_right, integral_val, b_m;
  double radius, phi, theta;
  double r_max = 1.0e8;
  double b_eps = 1.0 / n * 1.0e-3;

  if (b_eps < DBL_EPSILON)
    b_eps = DBL_EPSILON;

  *x = (double *) malloc(3*n*sizeof(double));

  init_fib_rnd(&rnd_seed, &fib_i, &fib_j, &i_ctrl);

  b_left = 0.0;
  b_right = r_max;
  integral_right = 1.0 - step;
  while (b_right - b_left > DBL_EPSILON)
  {
    integral_val = int_plummer(b_right);
    if (integral_val < integral_right)
      break;
    b_right = 0.5 * (b_right - b_left);
  }
  /* Get upper boundary of search radius */
  r_max = b_right;

  integral_val = 0.0;

  for (j = 0; j < n; j++)
  {
    integral_left = step * j;
    integral_right = step * (j + 1);

    b_left = 0.0;
    b_right = r_max;

    /* Get radius by binary search */    
    while (b_right - b_left > b_eps)
    {
      b_m = 0.5 * (b_left + b_right);
      integral_val = int_plummer(b_m);
      if (integral_val < integral_left)
        b_left = b_m;
      else if (integral_val > integral_right)
        b_right = b_m;
      else
      {
        b_left = b_m - 0.5 * (b_m - b_left);
        b_right = b_m + 0.5 * (b_m - b_left);
      }
    }
    radius = b_m;

    /* Random phi in [0, 2*PI] */
    phi = 2 * PI * fib_rnd(fib_x, &rnd_seed, &fib_i, &fib_j, &i_ctrl);

    /* Random theta in [0, PI] */
    theta = acos(1.0 - 2*fib_rnd(fib_x, &rnd_seed, &fib_i, &fib_j, &i_ctrl));

    (*x)[3*j]=radius*sin(theta)*cos(phi);
    (*x)[3*j+1]=radius*sin(theta)*sin(phi);
    (*x)[3*j+2]=radius*cos(theta);
  }

  auto_shift_scale(3, n, *x);

  *n_total = n;
}
