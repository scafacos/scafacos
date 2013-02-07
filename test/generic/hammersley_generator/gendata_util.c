/*
 * Author: Toni Volkmer
 *
 */

#include <stdio.h>
#include <assert.h>
#include <math.h>

#include "gendata_util.h"

int prime_vec[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31};

int get_prime(int n)
{
  assert(n>=0);
  if (n == 0)
    return 1;

  if (n > 10)
    return prime_vec[10];

  return prime_vec[n - 1];
}

void init_fib_rnd(double *rnd_seed, int *fib_i, int *fib_j, int *i_ctrl)
{
  *rnd_seed = 1.0;
  *fib_i = 16;
  *fib_j = 4;
  *i_ctrl = 0;
}

double rnd(double *rnd_seed)
{
  double rand;
  double d2p31m = 2147483647.0;
  double d2p31  = 2147483648.0;

  *rnd_seed = fmod(16807.0 * (*rnd_seed), d2p31m);

  rand     = (*rnd_seed / d2p31);
  return rand;
}

double fib_rnd(double *x, double *rnd_seed, int *fib_i, int *fib_j, int *i_ctrl)
{
  double rand;
  int k;

  if(*i_ctrl == 0)
  {
    for (k = 0; k < 17; k++)
    {
      x[k] = rnd(rnd_seed);
    }
    *i_ctrl = 1;
  }
  rand = x[*fib_i] - x[*fib_j];
  if (rand < 0.0)
    rand = rand + 1.0;
  x[*fib_i] = rand;

  *fib_i = (*fib_i)-1;
  if( *fib_i < 0 )
    *fib_i = 16;

  *fib_j = (*fib_j)-1;
  if( *fib_j < 0 )
    *fib_j = 16;

  return rand;
}

void set_random_charge(int n, double *q)
{
    int j;
    int count_pos = 0;
    int count_neg = 0;
    int count_neg_max = n/2;
    int count_pos_max = n - count_neg_max;
    double rand;
    double rnd_seed;
    int fib_i, fib_j, i_ctrl;
    double x[17];

    init_fib_rnd(&rnd_seed, &fib_i, &fib_j, &i_ctrl);

    for (j = 0; j < n; j++)
    {
      rand = fib_rnd(x, &rnd_seed, &fib_i, &fib_j, &i_ctrl);
      if (rand < 0.5)
      {
        q[j] = 1.0;
        count_pos++;
      }
      else
      {
        q[j] = -1.0;
        count_neg++;
      }
    }

    while (count_pos > count_pos_max)
    {
      j = fib_rnd(x, &rnd_seed, &fib_i, &fib_j, &i_ctrl) * (n-1);
      assert(j >= 0 && j < n);
      if (q[j] > 0)
      {
        q[j] = -1.0;
        count_pos--;
      }
    }

    while (count_neg > count_neg_max)
    {
      j = fib_rnd(x, &rnd_seed, &fib_i, &fib_j, &i_ctrl) * (n-1);
      assert(j >= 0 && j < n);
      if (q[j] < 0)
      {
        q[j] = 1.0;
        count_neg--;
      }
    }
 }

  void set_negative_charge(int n, double *q)
  {
    int j;
    for (j = 0; j < n; j++)
    {
      q[j] = -1.0;
    }
  }

  void set_positive_charge(int n, double *q)
  {
    int j;
    for (j = 0; j < n; j++)
    {
      q[j] = 1.0;
    }
  }

  void shift_scale(int d, int n, double *x, double shift, double scale)
  {
    int j;
    for (j = 0; j < d*n; j++)
    {
      x[j] = (x[j] + shift) * scale;
    }
  }

  void auto_shift_scale(int d, int n, double *x)
  {
    double min[d], max[d], shift[d], scale[d];
    double min_scale;
    int j, t;

    for (t = 0; t < d; t++)
    {
      min[t] = x[t];
      max[t] = x[t];
    }

    for (j = 1; j < n; j++)
    {
      for (t = 0; t < d; t++)
      {
        if (min[t] > x[3*j+t])
          min[t] = x[3*j+t];
        if (max[t] < x[3*j+t])
          max[t] = x[3*j+t];
      }
    }

    for (t = 0; t < d; t++)
    {
      scale[t] = 1.0 / (max[t] - min[t]);
    }

    min_scale = scale[0];
    for (t = 1; t < d; t++)
    {
      if (min_scale > scale[t])
        min_scale = scale[t];
    }

    for (t = 0; t < d; t++)
    {
      shift[t] = -min[t];
      if (scale[t] > min_scale)
        shift[t] += 0.5 * (1.0 - min_scale / scale[t]) / min_scale;
    }

    for (j = 0; j < n; j++)
      for (t = 0; t < d; t++)
        x[3*j+t] = (x[3*j+t] + shift[t]) * min_scale;
  }

