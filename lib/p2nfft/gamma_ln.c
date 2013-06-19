/*
 * Copyright (C) 2013 Michael Pippig
 *
 * This file is part of ScaFaCoS.
 * 
 * ScaFaCoS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ScaFaCoS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *	
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


/* Compute the function 
 *   EulerGamma + Gamma[0,x] + Log[x]
 * on the interval [0,1] with a Chebyshev approximation to
 * avoids numerical difficulties due to the opposite singularities
 * of Gamma[0,x] and Log[x]. */

#include <math.h>
#include <stdio.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_gamma.h>


static const fcs_float P[] =
{
  K(0.22856666525896971348219207753017),
  K(0.22174041832582477887921648885903),
  K(-0.0066414479047707522728438675631585),
  K(0.00018053841266890865803650876681953),
  K(-4.1763665459723693706279837312874e-6),
  K(8.2798615179621263672930538970763e-8),
  K(-1.4284823341958152541342579101148e-9),
  K(2.176166497489546124276791622805e-11),
  K(-2.9643222053483918846152412326154e-13),
  K(3.6489123850309009562831771890637e-15),
  K(-4.0951564610143524156809346199965e-17),
  K(4.2220511711581257368276994505152e-19),
  K(-4.0243338005862861421852849127205e-21)
}


static const INT N = sizeof(P)/sizeof(P[0]);

static inline R evaluate_chebyshev(const INT n, const fcs_float *c, const fcs_float x)
{
  R a = c[n-2], b = c[n-1], t;
  int j;
  
  A(n >= 2);
  
  for (j = n - 2; j > 0; j--)
  {
    t = c[j-1] - b;
    b = a + K(2.0) * x * b;
    a = t;
  }
  return a + x * b;
}


fcs_float ifcs_p2nfft_gamma_ln(
    fcs_float x
    )
{
  /* Use a Chebyshev approximation on the interval [0,1] to
   * avoid numerical difficulties due to the opposite singularities
   * of Gamma[0,x] and Log[x]. */
  if(x<1)
    return evaluate_chebyshev(N, P, x);

  return  






}
