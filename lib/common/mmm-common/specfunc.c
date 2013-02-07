/*
  Copyright (C) 2010,2011 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 
    Max-Planck-Institute for Polymer Research, Theory Group
  
  This file is part of ScaFaCoS.
  
  ScaFaCoS is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  ScaFaCoS is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>. 
*/

/** \file specfunc.c
    Special functions, see \ref specfunc.h "specfunc.h"
*/
#include <math.h>
#include <stdlib.h>
#include "specfunc.h"

static void mmm_preparePolygammaEven(fcs_int n, fcs_float binom, SizedList *series);

static fcs_float mmm_bk0_data[11] = {
  -.5 -0.03532739323390276872,
   0.3442898999246284869, 
   0.03597993651536150163,
   0.00126461541144692592,
   0.00002286212103119451,
   0.00000025347910790261,
   0.00000000190451637722,
   0.00000000001034969525,
   0.00000000000004259816,
   0.00000000000000013744,
   0.00000000000000000035
};
static SizedList mmm_bk0_cs = {mmm_bk0_data, 11, 11 };

static fcs_float mmm_ak0_data[17] = {
  2.5 -0.07643947903327941,
  -0.02235652605699819,
   0.00077341811546938,
  -0.00004281006688886,
   0.00000308170017386,
  -0.00000026393672220,
   0.00000002563713036,
  -0.00000000274270554,
   0.00000000031694296,
  -0.00000000003902353,
   0.00000000000506804,
  -0.00000000000068895,
   0.00000000000009744,
  -0.00000000000001427,
   0.00000000000000215,
  -0.00000000000000033,
   0.00000000000000005
};
static SizedList mmm_ak0_cs = { mmm_ak0_data, 16, 16 };

static fcs_float mmm_ak02_data[14] = {
  2.5 -0.01201869826307592,
  -0.00917485269102569,
   0.00014445509317750,
  -0.00000401361417543,
   0.00000015678318108,
  -0.00000000777011043,
   0.00000000046111825,
  -0.00000000003158592,
   0.00000000000243501,
  -0.00000000000020743,
   0.00000000000001925,
  -0.00000000000000192,
   0.00000000000000020,
  -0.00000000000000002
};
static SizedList mmm_ak02_cs = { mmm_ak02_data, 13, 13 };

static fcs_float mmm_ak1_data[17] = {
  2.5 +0.27443134069738830, 
   0.07571989953199368,
  -0.00144105155647540,
   0.00006650116955125,
  -0.00000436998470952,
   0.00000035402774997,
  -0.00000003311163779,
   0.00000000344597758,
  -0.00000000038989323,
   0.00000000004720819,
  -0.00000000000604783,
   0.00000000000081284,
  -0.00000000000011386,
   0.00000000000001654,
  -0.00000000000000248,
   0.00000000000000038,
  -0.00000000000000006
};
static SizedList mmm_ak1_cs = { mmm_ak1_data, 17, 17 };

static fcs_float mmm_ak12_data[14] = {
  2.5 +0.06379308343739001,
   0.02832887813049721,
  -0.00024753706739052,
   0.00000577197245160,
  -0.00000020689392195,
   0.00000000973998344,
  -0.00000000055853361,
   0.00000000003732996,
  -0.00000000000282505,
   0.00000000000023720,
  -0.00000000000002176,
   0.00000000000000215,
  -0.00000000000000022,
   0.00000000000000002
};
static SizedList mmm_ak12_cs = { mmm_ak12_data, 14, 14 };

/* based on SLATEC besi1(), besi1e() */

/* chebyshev expansions

 series for bi1        on the interval  0.      to  9.00000d+00
          with weighted error   2.40e-17
           log weighted error  16.62
             significant figures required  16.23
            decimal places required  17.14

 series for ai1        on the interval  1.25000d-01 to  3.33333d-01
          with weighted error   6.98e-17
           log weighted error  16.16
             significant figures required  14.53
            decimal places required  16.82

 series for ai12       on the interval  0.      to  1.25000d-01
               with weighted error   3.55e-17
          log weighted error  16.45
            significant figures required  14.69
           decimal places required  17.12
*/

static fcs_float mmm_bi0_data[12] = {
  5.5 -.07660547252839144951,
  1.92733795399380827000,
   .22826445869203013390, 
   .01304891466707290428,
   .00043442709008164874,
   .00000942265768600193,
   .00000014340062895106,
   .00000000161384906966,
   .00000000001396650044,
   .00000000000009579451,
   .00000000000000053339,
   .00000000000000000245
};
static SizedList mmm_bi0_cs = { mmm_bi0_data, 12, 12 };

static fcs_float mmm_bk1_data[11] = {
  1.5 +0.0253002273389477705,
  -0.3531559607765448760, 
  -0.1226111808226571480, 
  -0.0069757238596398643,
  -0.0001730288957513052,
  -0.0000024334061415659,
  -0.0000000221338763073,
  -0.0000000001411488392,
  -0.0000000000006666901,
  -0.0000000000000024274,
  -0.0000000000000000070
};
static SizedList mmm_bk1_cs = { mmm_bk1_data, 11, 11 };

static fcs_float mmm_bi1_data[11] = {
  1.75 -0.001971713261099859,
   0.407348876675464810,
   0.034838994299959456,
   0.001545394556300123,
   0.000041888521098377,
   0.000000764902676483,
   0.000000010042493924,
   0.000000000099322077,
   0.000000000000766380,
   0.000000000000004741,
   0.000000000000000024
};
static SizedList mmm_bi1_cs = { mmm_bi1_data, 11, 11 };

/************************************************
 * chebychev expansions
 ************************************************/

/* coefficients for Maclaurin summation in hzeta()
 * B_{2j}/(2j)!
 */
static fcs_float mmm_hzeta_c[15] = {
  1.00000000000000000000000000000,
  0.083333333333333333333333333333,
 -0.00138888888888888888888888888889,
  0.000033068783068783068783068783069,
 -8.2671957671957671957671957672e-07,
  2.0876756987868098979210090321e-08,
 -5.2841901386874931848476822022e-10,
  1.3382536530684678832826980975e-11,
 -3.3896802963225828668301953912e-13,
  8.5860620562778445641359054504e-15,
 -2.1748686985580618730415164239e-16,
  5.5090028283602295152026526089e-18,
 -1.3954464685812523340707686264e-19,
  3.5347070396294674716932299778e-21,
 -8.9535174270375468504026113181e-23
};

fcs_float mmm_hzeta(fcs_float s, fcs_float q)
{
  fcs_float max_bits = 54.0;
  fcs_int jmax = 12, kmax = 10;
  fcs_int j, k;
  fcs_float pmax, scp, pcp, ans;

  if((s > max_bits && q < 1.0) || (s > 0.5*max_bits && q < 0.25))
    return pow(q, -s);
  if(s > 0.5*max_bits && q < 1.0) {
    fcs_float p1 = pow(q, -s);
    fcs_float p2 = pow(q/(1.0+q), s);
    fcs_float p3 = pow(q/(2.0+q), s);
    return p1 * (1.0 + p2 + p3);
  }
  /* Euler-Maclaurin summation formula 
   * [Moshier, p. 400, with several typo corrections]
   */
  pmax  = pow(kmax + q, -s);
  scp = s;
  pcp = pmax / (kmax + q);
  ans = pmax*((kmax+q)/(s-1.0) + 0.5);

  for(k=0; k<kmax; k++)
    ans += pow(k + q, -s);


  for(j=0; j<=jmax; j++) {
    fcs_float delta = mmm_hzeta_c[j+1] * scp * pcp;
    ans += delta;
    scp *= (s+2*j+1)*(s+2*j+2);
    pcp /= (kmax + q)*(kmax + q);
  }

  return ans;
}

fcs_float mmm_K0(fcs_float x)
{
  fcs_float c, I0;
  if(x <= 2.0) {
    c  = mmm_evaluateAsChebychevSeriesAt(&mmm_bk0_cs, 0.5*x*x-1.0);
    I0 = mmm_evaluateAsChebychevSeriesAt(&mmm_bi0_cs, x*x/4.5-1.0);
    return (-log(x) + M_LN2)*I0 + c;
  }
  c = (x <= 8.0) ?
    mmm_evaluateAsChebychevSeriesAt(&mmm_ak0_cs, (16.0/x-5.0)/3.0) :
    mmm_evaluateAsChebychevSeriesAt(&mmm_ak02_cs, 16.0/x-1.0);
  return exp(-x)*c/sqrt(x); 
}

fcs_float mmm_K1(fcs_float x)
{
  fcs_float c, I1;
  if(x <= 2.0) {
    c = mmm_evaluateAsChebychevSeriesAt(&mmm_bk1_cs, 0.5*x*x-1.0);
    I1 = x * mmm_evaluateAsChebychevSeriesAt(&mmm_bi1_cs, x*x/4.5-1.0);
    return (log(x) - M_LN2)*I1 + c/x;
  }
  c = (x <= 8.0) ?
    mmm_evaluateAsChebychevSeriesAt(&mmm_ak1_cs, (16.0/x-5.0)/3.0) :
    mmm_evaluateAsChebychevSeriesAt(&mmm_ak12_cs, 16.0/x-1.0);
  return exp(-x)*c/sqrt(x);
}

/** evaluate the polynomial interpreted as a Chebychev series. Requires a series with at least
    three coefficients, i.e. no linear approximations! */
fcs_float mmm_evaluateAsChebychevSeriesAt(SizedList *series, fcs_float x)
{
  fcs_int j;
  fcs_float *c = series->e;
  fcs_float x2 = 2.0 * x;
  fcs_float dd = c[series->n - 1];
  fcs_float d  = x2*dd + c[series->n - 2];
  for(j = series->n - 3; j >= 1; j--) {
    fcs_float tmp = d;
    d = x2*d - dd + c[j];
    dd = tmp;
  }
  return x*d - dd + 0.5 * c[0];
}

/** evaluate the polynomial interpreted as a Taylor series via the Horner scheme */
fcs_float mmm_evaluateAsTaylorSeriesAt(SizedList *series, fcs_float x)
{
  fcs_int cnt   = (series->n) - 1;
  fcs_float *c = series->e;
  fcs_float r  = c[cnt];
  while (--cnt >= 0)
    r = r*x + c[cnt];
  return r;
}

/***********************************************************
 * optimized K0/1 implementations for 10^(-14) precision
 ***********************************************************/

/** necessary orders for K0/1 from 2 up to 22 for 10^-14 precision. Note that at 8
    the expansion changes. From 23 to 26 order 2 is used, above order 1. For the
    latter cases separate implementations are necessary. */
static int mmm_ak01_orders[] = {
  /* 2 - 8 */
     11, 11, 10,
  10, 9, 9,
  /* 8 - 26 */
           6, 6,
  5, 5, 5, 4, 4,
  4, 3, 3, 2, 2,
  2, 2, 2
};

void mmm_LPK01(fcs_float x, fcs_float *K0, fcs_float *K1)
{
  if (x >= 27.) {
    fcs_float tmp = .5*exp(-x)/sqrt(x);
    *K0 = tmp*mmm_ak0_data[0];
    *K1 = tmp*mmm_ak1_data[0];
    return;
  }
  if (x >= 23.) {
    fcs_float tmp = exp(-x)/sqrt(x), xx = (16./3.)/x - 5./3.;
    *K0 = tmp*(xx*mmm_ak0_data[1] + 0.5*mmm_ak0_data[0]);
    *K1 = tmp*(xx*mmm_ak1_data[1] + 0.5*mmm_ak1_data[0]);
    return;    
  }
  if (x > 2) {
    fcs_int j = mmm_ak01_orders[((int)x) - 2];
    fcs_float tmp, x2;
    fcs_float dd0, dd1, d0, d1;
    fcs_float *s0, *s1;
    if (x <= 8) {
      s0 = mmm_ak0_data; s1 = mmm_ak1_data;
      x2 = (2.*16./3.)/x - 2.*5./3.;
    } else {
      s0 = mmm_ak02_data; s1 = mmm_ak12_data;
      x2 = (2.*16.)/x - 2.;
    }
    dd0 = s0[j];
    dd1 = s1[j];
    d0  = x2*dd0 + s0[j - 1];
    d1  = x2*dd1 + s1[j - 1];
    for(j -= 2; j >= 1; j--) {
      fcs_float tmp0 = d0, tmp1 = d1;
      d0 = x2*d0 - dd0 + s0[j];
      d1 = x2*d1 - dd1 + s1[j];
      dd0 = tmp0;
      dd1 = tmp1;      
    }
    tmp = exp(-x)/sqrt(x);
    *K0 = tmp*(0.5*(s0[0] + x2*d0) - dd0);
    *K1 = tmp*(0.5*(s1[0] + x2*d1) - dd1);
    return;
  }
  /* x <= 2 */
  {
    /* I0/1 series */
    int j = 10;
    fcs_float tmp, x2 = (2./4.5)*x*x - 2.;
    fcs_float dd0, dd1, d0, d1;
    dd0 = mmm_bi0_data[j];
    dd1 = mmm_bi1_data[j];
    d0  = x2*dd0 + mmm_bi0_data[j - 1];
    d1  = x2*dd1 + mmm_bi1_data[j - 1];
    for(j -= 2; j >= 1; j--) {
      fcs_float tmp0 = d0, tmp1 = d1;
      d0 = x2*d0 - dd0 + mmm_bi0_data[j];
      d1 = x2*d1 - dd1 + mmm_bi1_data[j];
      dd0 = tmp0;
      dd1 = tmp1;      
    }
    tmp = log(x) - M_LN2;
    *K0 =  -tmp*(0.5*(mmm_bi0_data[0] + x2*d0)- dd0);
    *K1 = x*tmp*(0.5*(mmm_bi1_data[0] + x2*d1) - dd1);

    /* K0/K1 correction */
    j = 9;
    x2 = x*x - 2.;
    dd0 = mmm_bk0_data[j];
    dd1 = mmm_bk1_data[j];
    d0  = x2*dd0 + mmm_bk0_data[j - 1];
    d1  = x2*dd1 + mmm_bk1_data[j - 1];
    for(j -= 2; j >= 1; j--) {
      fcs_float tmp0 = d0, tmp1 = d1;
      d0 = x2*d0 - dd0 + mmm_bk0_data[j];
      d1 = x2*d1 - dd1 + mmm_bk1_data[j];
      dd0 = tmp0;
      dd1 = tmp1;      
    }
    *K0 += (0.5*(x2*d0 + mmm_bk0_data[0]) - dd0);
    *K1 += (0.5*(x2*d1 + mmm_bk1_data[0]) - dd1)/x;
    return;
  }
}

static void mmm_preparePolygammaEven(fcs_int n, fcs_float binom, SizedList *series)
{
  /* (-0.5 n) psi^2n/2n! (-0.5 n) and psi^(2n+1)/(2n)! series expansions
     note that BOTH carry 2n! */
  fcs_int order;
  fcs_float deriv;
  fcs_float maxx, x_order, coeff, pref;

  deriv = 2*n;
  if (n == 0) {
    // psi^0 has a slightly different series expansion
    maxx = 0.25;
    mmm_alloc_doublelist(series, 1);
    series->e[0] = 2*(1 - MMM_COMMON_C_GAMMA);
    for (order = 1;; order += 1) {
      x_order = 2*order;
      coeff = -2*mmm_hzeta(x_order + 1, 2);
      if (fcs_float_is_zero(fabs(maxx*coeff)*(4.0/3.0)))
  break;
      mmm_realloc_doublelist(series, order + 1);
      series->e[order] = coeff;
      maxx *= 0.25;
    }
    series->n = order;
  }
  else {
    // even, n > 0
    maxx = 1;
    pref = 2;
    mmm_init_doublelist(series);
    for (order = 0;; order++) {
      // only even exponents of x
      x_order = 2*order;
      coeff = pref*mmm_hzeta(1 + deriv + x_order, 2);
      if ((fcs_float_is_zero(fabs(maxx*coeff)*(4.0/3.0))) && (x_order > deriv))
  break;
      mmm_realloc_doublelist(series, order + 1);
      series->e[order] = -binom*coeff;
      maxx *= 0.25;
      pref *= (1.0 + deriv/(x_order + 1));
      pref *= (1.0 + deriv/(x_order + 2));
    }
    series->n = order;
  }
}

static void mmm_preparePolygammaOdd(fcs_int n, fcs_float binom, SizedList *series)
{
  fcs_int order;
  fcs_float deriv;
  fcs_float maxx, x_order, coeff, pref;
  deriv  = 2*n + 1;
  maxx = 0.5;
  // to get 1/(2n)! instead of 1/(2n+1)!
  pref = 2*deriv*(1 + deriv);
  mmm_init_doublelist(series);
  for (order = 0;; order++) {
    // only odd exponents of x
    x_order = 2*order + 1;
    coeff = pref*mmm_hzeta(1 + deriv + x_order, 2);
    if ((fcs_float_is_zero(fabs(maxx*coeff)*(4.0/3.0))) && (x_order > deriv))
      break;

    mmm_realloc_doublelist(series, order + 1);
    series->e[order] = -binom*coeff;
    maxx *= 0.25;
    pref *= (1.0 + deriv/(x_order + 1));
    pref *= (1.0 + deriv/(x_order + 2));
  }

  series->n = order;
}

void mmm_create_mod_psi_up_to(mmm_data_struct *polTaylor, fcs_int new_n)
{
  fcs_int n;
  fcs_float binom;

  if (new_n > polTaylor->n_modPsi) {
    fcs_int old = polTaylor->n_modPsi;
    polTaylor->n_modPsi = new_n;
    polTaylor->modPsi = realloc(polTaylor->modPsi, 2*(polTaylor->n_modPsi)*sizeof(SizedList));

    binom = 1.0;
    for (n = 0; n < old; n++)
      binom *= (-0.5 - n)/(fcs_float)(n+1);

    for (; n < polTaylor->n_modPsi; n++) {
      mmm_preparePolygammaEven(n, binom, &(polTaylor->modPsi)[2*n]);
      mmm_preparePolygammaOdd(n, binom, &(polTaylor->modPsi)[2*n + 1]);
      binom *= (-0.5 - n)/(fcs_float)(n+1);
    }
  }
}

/** modified polygamma for even order 2*n, n >= 0 */
fcs_float mmm_mod_psi_even(mmm_data_struct *polTaylor, fcs_int n, fcs_float x) {
   return mmm_evaluateAsTaylorSeriesAt(&(polTaylor->modPsi)[2*n], x*x);
}

/** modified polygamma for odd order 2*n+1, n>= 0 */
fcs_float mmm_mod_psi_odd(mmm_data_struct *polTaylor, fcs_int n, fcs_float x) {
   return x*mmm_evaluateAsTaylorSeriesAt(&(polTaylor->modPsi)[2*n+1], x*x);
}
