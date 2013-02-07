/*
 * Author: Toni Volkmer
 *
 */

#include <stdlib.h>

#include "gendata_grid.h"
#include "gendata_util.h"

void gendata_grid_face_centered_cube(int n_component, double **x, int *n)
{
  int j = 0;
  int u, v, w;

  double a = 1.0 / (n_component - 0.5);
  *n = 4 * n_component * n_component * n_component;

  *x = (double *) malloc(3 * (*n) * sizeof(double));

  for (u = 0; u < n_component; u++)
    for (v = 0; v < n_component; v++)
      for (w = 0; w < n_component; w++)
      {
        (*x)[j++] = u * a;
        (*x)[j++] = v * a;
        (*x)[j++] = w * a;

        (*x)[j++] = (u+0.5) * a;
        (*x)[j++] = (v+0.5) * a;
        (*x)[j++] = w * a;

        (*x)[j++] = (u+0.5) * a;
        (*x)[j++] = v * a;
        (*x)[j++] = (w+0.5) * a;

        (*x)[j++] = u * a;
        (*x)[j++] = (v+0.5) * a;
        (*x)[j++] = (w+0.5) * a;
      }
}

void gendata_grid_body_centered_cube(int n_component, double **x, int *n)
{
  int j = 0;
  int u, v, w;

  double a = 1.0 / (n_component - 0.5);
  *n = 2 * n_component * n_component * n_component;

  *x = (double *) malloc(3 * (*n) * sizeof(double));

  for (u = 0; u < n_component; u++)
    for (v = 0; v < n_component; v++)
      for (w = 0; w < n_component; w++)
      {
        (*x)[j++] = u * a;
        (*x)[j++] = v * a;
        (*x)[j++] = w * a;

        (*x)[j++] = (u+0.5) * a;
        (*x)[j++] = (v+0.5) * a;
        (*x)[j++] = (w+0.5) * a;
      }
}

void gendata_nacl_cube(int n_component, double **x, double **q, int *n)
{
  int j = 0;
  int k = 0;
  int u, v, w;
  double a;

  *n = n_component * n_component * n_component;

  *x = (double *) malloc(3 * (*n) * sizeof(double));
  *q = (double *) malloc((*n) * sizeof(double));

  if (n_component > 1)
  {
    a = 1.0 / (n_component - 1);
  }
  else
  {
    a = 1.0;
  }

  for (u = 0; u < n_component; u++)
    for (v = 0; v < n_component; v++)
      for (w = 0; w < n_component; w++)
      {
        if ((u+v+w) % 2)
          (*q)[k++] = 1.0;
        else
          (*q)[k++] = -1.0;
        (*x)[j++] = u * a;
        (*x)[j++] = v * a;
        (*x)[j++] = w * a;
      }
}

void gendata_nacl_cube_periodic(int n_component, double **x, double **q, int *n)
{
  int j = 0;
  int k = 0;
  int u, v, w;
  double a;

  /* assure even number of charges in each direction */
  n_component = 2 * (n_component/2);

  *n = n_component * n_component * n_component;

  *x = (double *) malloc(3 * (*n) * sizeof(double));
  *q = (double *) malloc((*n) * sizeof(double));

  a = 1.0;

  for (u = 0; u < n_component; u++)
    for (v = 0; v < n_component; v++)
      for (w = 0; w < n_component; w++)
      {
        if ((u+v+w) % 2)
          (*q)[k++] = 1.0;
        else
          (*q)[k++] = -1.0;
        (*x)[j++] = u * a + 0.5;
        (*x)[j++] = v * a + 0.5;
        (*x)[j++] = w * a + 0.5;
      }
}
