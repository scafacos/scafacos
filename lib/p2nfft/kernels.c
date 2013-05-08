/*
 * Copyright (C) 2011-2013 Michael Pippig
 * Copyright (c) 2002, 2009 Jens Keiner, Stefan Kunis, Daniel Potts
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

#include "config.h"

#include <stdio.h>
#include <math.h>
#include <float.h>
#ifdef HAVE_COMPLEX_H
#include <complex.h>
#endif

#include "kernels.h"

FCS_P2NFFT_KERNEL_TYPE ifcs_p2nfft_gaussian(fcs_float x, fcs_int der, const fcs_float *param)    /* K(x)=exp(-x^2/c^2) */
{
  fcs_float c=param[0];
  fcs_float value=0.0;

  switch (der)
  {
    case  0 : value=exp(-x*x/(c*c)); break;
    case  1 : value=-2.0*x/(c*c)*exp(-x*x/(c*c)); break;
    case  2 : value=2.0*exp(-x*x/(c*c))*(-c*c+2.0*x*x)/(c*c*c*c); break;
    case  3 : value=-4.0*x*exp(-x*x/(c*c))*(-3.0*c*c+2.0*x*x)/(c*c*c*c*c*c); break;
    case  4 : value=4.0*exp(-x*x/(c*c))*(3.0*c*c*c*c-12.0*c*c*x*x+4.0*x*x*x*x)/(c*c*c*c*c*c*c*c); break;
    case  5 : value=-8.0*x*exp(-x*x/(c*c))*(15.0*c*c*c*c-20.0*c*c*x*x+4.0*x*x*x*x)/pow(c,10.0); break;
    case  6 : value=8.0*exp(-x*x/(c*c))*(-15.0*c*c*c*c*c*c+90.0*x*x*c*c*c*c-60.0*x*x*x*x*c*c+8.0*x*x*x*x*x*x)/pow(c,12.0); break;
    case  7 : value=-16.0*x*exp(-x*x/(c*c))*(-105.0*c*c*c*c*c*c+210.0*x*x*c*c*c*c-84.0*x*x*x*x*c*c+8.0*x*x*x*x*x*x)/pow(c,14.0); break;
    case  8 : value=16.0*exp(-x*x/(c*c))*(105.0*c*c*c*c*c*c*c*c-840.0*x*x*c*c*c*c*c*c+840.0*x*x*x*x*c*c*c*c-224.0*x*x*x*x*x*x*c*c+16.0*x*x*x*x*x*x*x*x)/pow(c,16.0); break;
    case  9 : value=-32.0*x*exp(-x*x/(c*c))*(945.0*c*c*c*c*c*c*c*c-2520.0*x*x*c*c*c*c*c*c+1512.0*x*x*x*x*c*c*c*c-288.0*x*x*x*x*x*x*c*c+16.0*x*x*x*x*x*x*x*x)/pow(c,18.0); break;
    case 10 : value=32.0*exp(-x*x/(c*c))*(-945.0*pow(c,10.0)+9450.0*x*x*c*c*c*c*c*c*c*c-12600.0*x*x*x*x*c*c*c*c*c*c+5040.0*x*x*x*x*x*x*c*c*c*c-720.0*x*x*x*x*x*x*x*x*c*c+32.0*pow(x,10.0))/pow(c,20.0); break;
    case 11 : value=-64.0*x*exp(-x*x/(c*c))*(-10395.0*pow(c,10.0)+34650.0*x*x*c*c*c*c*c*c*c*c-27720.0*x*x*x*x*c*c*c*c*c*c+7920.0*x*x*x*x*x*x*c*c*c*c-880.0*x*x*x*x*x*x*x*x*c*c+32.0*pow(x,10.0))/pow(c,22.0); break;
    case 12 : value=64.0*exp(-x*x/(c*c))*(10395.0*pow(c,12.0)-124740.0*x*x*pow(c,10.0)+207900.0*x*x*x*x*c*c*c*c*c*c*c*c-110880.0*x*x*x*x*x*x*c*c*c*c*c*c+23760.0*x*x*x*x*x*x*x*x*c*c*c*c-2112.0*pow(x,10.0)*c*c+64.0*pow(x,12.0))/pow(c,24.0); break;
    default : value=0.0;
  }

  return value;
}

FCS_P2NFFT_KERNEL_TYPE ifcs_p2nfft_multiquadric(fcs_float x, fcs_int der, const fcs_float *param)    /* K(x)=sqrt(x^2+c^2) */
{
  fcs_float c=param[0];
  fcs_float value=0.0;

  switch (der)
  {
    case  0 : value=sqrt(x*x+c*c); break;
    case  1 : value=1.0/(sqrt(x*x+c*c))*x; break;
    case  2 : value=c*c/sqrt(pow(x*x+c*c,3.0)); break;
    case  3 : value=-3.0*x*c*c/sqrt(pow(x*x+c*c,5.0)); break;
    case  4 : value=3.0*c*c*(4.0*x*x-c*c)/sqrt(pow(x*x+c*c,7.0)); break;
    case  5 : value=-15.0*x*c*c*(4.0*x*x-3.0*c*c)/sqrt(pow(x*x+c*c,9.0)); break;
    case  6 : value=45.0*c*c*(8.0*x*x*x*x-12.0*x*x*c*c+c*c*c*c)/sqrt(pow(x*x+c*c,11.0)); break;
    case  7 : value=-315.0*x*c*c*(8.0*x*x*x*x-20.0*x*x*c*c+5.0*c*c*c*c)/sqrt(pow(x*x+c*c,13.0)); break;
    case  8 : value=315.0*c*c*(64.0*x*x*x*x*x*x-240.0*x*x*x*x*c*c+120.0*x*x*c*c*c*c-5.0*c*c*c*c*c*c)/sqrt(pow(x*x+c*c,15.0)); break;
    case  9 : value=-2835.0*x*c*c*(64.0*x*x*x*x*x*x-336.0*x*x*x*x*c*c+280.0*x*x*c*c*c*c-35.0*c*c*c*c*c*c)/sqrt(pow(x*x+c*c,17.0)); break;
    case 10 : value=14175.0*c*c*(128.0*x*x*x*x*x*x*x*x-896.0*x*x*x*x*x*x*c*c+1120.0*x*x*x*x*c*c*c*c-280.0*x*x*c*c*c*c*c*c+7.0*c*c*c*c*c*c*c*c)/sqrt(pow(x*x+c*c,19.0)); break;
    case 11 : value=-155925.0*x*c*c*(128.0*x*x*x*x*x*x*x*x-1152.0*x*x*x*x*x*x*c*c+2016.0*x*x*x*x*c*c*c*c-840.0*x*x*c*c*c*c*c*c+63.0*c*c*c*c*c*c*c*c)/sqrt(pow(x*x+c*c,21.0)); break;
    case 12 : value=467775.0*c*c*(1260.0*x*x*c*c*c*c*c*c*c*c-21.0*pow(c,10.0)+512.0*pow(x,10.0)-5760.0*x*x*x*x*x*x*x*x*c*c+13440.0*x*x*x*x*x*x*c*c*c*c-8400.0*x*x*x*x*c*c*c*c*c*c)/sqrt(pow(x*x+c*c,23.0)); break;
    default : value=0.0;
  }

  return value;
}

FCS_P2NFFT_KERNEL_TYPE ifcs_p2nfft_inverse_multiquadric(fcs_float x, fcs_int der, const fcs_float *param)    /* K(x)=1/sqrt(x^2+c^2) */
{
  fcs_float c=param[0];
  fcs_float value=0.0;

  switch (der)
  {
    case  0 : value=1.0/sqrt(x*x+c*c); break;
    case  1 : value=-1.0/(sqrt(pow(x*x+c*c,3.0)))*x; break;
    case  2 : value=(2.0*x*x-c*c)/sqrt(pow(x*x+c*c,5.0)); break;
    case  3 : value=-3.0*x*(2.0*x*x-3.0*c*c)/sqrt(pow(x*x+c*c,7.0)); break;
    case  4 : value=3.0*(8.0*x*x*x*x-24.0*x*x*c*c+3.0*c*c*c*c)/sqrt(pow(x*x+c*c,9.0)); break;
    case  5 : value=-15.0*x*(8.0*x*x*x*x-40.0*x*x*c*c+15.0*c*c*c*c)/sqrt(pow(x*x+c*c,11.0)); break;
    case  6 : value=45.0*(16.0*x*x*x*x*x*x-120.0*x*x*x*x*c*c+90.0*x*x*c*c*c*c-5.0*c*c*c*c*c*c)/sqrt(pow(x*x+c*c,13.0)); break;
    case  7 : value=-315.0*x*(16.0*x*x*x*x*x*x-168.0*x*x*x*x*c*c+210.0*x*x*c*c*c*c-35.0*c*c*c*c*c*c)/sqrt(pow(x*x+c*c,15.0)); break;
    case  8 : value=315.0*(128.0*x*x*x*x*x*x*x*x-1792.0*x*x*x*x*x*x*c*c+3360.0*x*x*x*x*c*c*c*c-1120.0*x*x*c*c*c*c*c*c+35.0*c*c*c*c*c*c*c*c)/sqrt(pow(x*x+c*c,17.0)); break;
    case  9 : value=-2835.0*x*(128.0*x*x*x*x*x*x*x*x-2304.0*x*x*x*x*x*x*c*c+6048.0*x*x*x*x*c*c*c*c-3360.0*x*x*c*c*c*c*c*c+315.0*c*c*c*c*c*c*c*c)/sqrt(pow(x*x+c*c,19.0)); break;
    case 10 : value=14175.0*(256.0*pow(x,10.0)-5760.0*x*x*x*x*x*x*x*x*c*c+20160.0*x*x*x*x*x*x*c*c*c*c-16800.0*x*x*x*x*c*c*c*c*c*c+3150.0*x*x*c*c*c*c*c*c*c*c-63.0*pow(c,10.0))/sqrt(pow(x*x+c*c,21.0)); break;
    case 11 : value=-155925.0*x*(256.0*pow(x,10.0)-7040.0*x*x*x*x*x*x*x*x*c*c+31680.0*x*x*x*x*x*x*c*c*c*c-36960.0*x*x*x*x*c*c*c*c*c*c+11550.0*x*x*c*c*c*c*c*c*c*c-693.0*pow(c,10.0))/sqrt(pow(x*x+c*c,23.0)); break;
    case 12 : value=467775.0*(231.0*pow(c,12.0)+190080.0*x*x*x*x*x*x*x*x*c*c*c*c-16632.0*x*x*pow(c,10.0)-295680.0*x*x*x*x*x*x*c*c*c*c*c*c+138600.0*x*x*x*x*c*c*c*c*c*c*c*c+1024.0*pow(x,12.0)-33792.0*pow(x,10.0)*c*c)/sqrt(pow(x*x+c*c,25.0)); break;
    default : value=0.0;
  }

  return value;
}

FCS_P2NFFT_KERNEL_TYPE ifcs_p2nfft_logarithm(fcs_float x, fcs_int der, const fcs_float *param)    /* K(x)=log |x| */
{
  fcs_float value=0.0;

  (void)param;

  if (fabs(x)<DBL_EPSILON) value=0.0;
  else switch (der)
  {
    case  0 : value=log(fabs(x)); break;
    case  1 : value=(x<0 ? -1 : 1)/fabs(x); break;
    case  2 : value=-1/(x*x); break;
    case  3 : value=2.0*(x<0 ? -1 : 1)/pow(fabs(x),3.0); break;
    case  4 : value=-6.0/(x*x*x*x); break;
    case  5 : value=24.0*(x<0 ? -1 : 1)/pow(fabs(x),5.0); break;
    case  6 : value=-120.0/(x*x*x*x*x*x); break;
    case  7 : value=720.0*(x<0 ? -1 : 1)/pow(fabs(x),7.0); break;
    case  8 : value=-5040.0/(x*x*x*x*x*x*x*x); break;
    case  9 : value=40320.0*(x<0 ? -1 : 1)/pow(fabs(x),9.0); break;
    case 10 : value=-362880.0/pow(x,10.0); break;
    case 11 : value=3628800.0*(x<0 ? -1 : 1)/pow(fabs(x),11.0); break;
    case 12 : value=-39916800.0/pow(x,12.0); break;
    case 13 : value=479001600.0/pow(x,13.0); break;
    case 14 : value=-6227020800.0/pow(x,14.0); break;
    case 15 : value=87178291200.0/pow(x,15.0); break;
    case 16 : value=-1307674368000.0/pow(x,16.0); break;
    case 17 : value=20922789888000.0/pow(x,17.0); break;
    default : value=0.0;
  }

  return value;
}

FCS_P2NFFT_KERNEL_TYPE ifcs_p2nfft_thinplate_spline(fcs_float x, fcs_int der, const fcs_float *param)    /* K(x) = x^2 log |x| */
{
  fcs_float value=0.0;

  (void)param;

  if (fabs(x)<DBL_EPSILON) value=0.0;
  else switch (der)
  {
    case  0 : value=x*x*log(fabs(x)); break;
    case  1 : value=2.0*x*log(fabs(x))+x; break;
    case  2 : value=2.0*log(fabs(x))+3.0; break;
    case  3 : value=2.0/x; break;
    case  4 : value=-2.0/(x*x); break;
    case  5 : value=4.0/(x*x*x); break;
    case  6 : value=-12.0/(x*x*x*x); break;
    case  7 : value=48.0/(x*x*x*x*x); break;
    case  8 : value=-240.0/(x*x*x*x*x*x); break;
    case  9 : value=1440.0/(x*x*x*x*x*x*x); break;
    case 10 : value=-10080.0/(x*x*x*x*x*x*x*x); break;
    case 11 : value=80640.0/(x*x*x*x*x*x*x*x*x); break;
    case 12 : value=-725760.0/pow(x,10.0); break;
    default : value=0.0;
  }

  return value;
}

FCS_P2NFFT_KERNEL_TYPE ifcs_p2nfft_one_over_square(fcs_float x, fcs_int der, const fcs_float *param)    /* K(x) = 1/x^2 */
{
  fcs_float value=0.0;

  (void)param;

  if (fabs(x)<DBL_EPSILON) value=0.0;
  else switch (der)
  {
    case  0 : value=1.0/(x*x); break;
    case  1 : value=-2.0/(x*x*x); break;
    case  2 : value=6.0/(x*x*x*x); break;
    case  3 : value=-24.0/(x*x*x*x*x); break;
    case  4 : value=120.0/(x*x*x*x*x*x); break;
    case  5 : value=-720.0/(x*x*x*x*x*x*x); break;
    case  6 : value=5040.0/(x*x*x*x*x*x*x*x); break;
    case  7 : value=-40320.0/(x*x*x*x*x*x*x*x*x); break;
    case  8 : value=362880.0/pow(x,10.0); break;
    case  9 : value=-3628800.0/pow(x,11.0); break;
    case 10 : value=39916800.0/pow(x,12.0); break;
    case 11 : value=-479001600.0/pow(x,13.0); break;
    case 12 : value=6227020800.0/pow(x,14.0); break;
    default : value=0.0;
  }

  return value;
}

FCS_P2NFFT_KERNEL_TYPE ifcs_p2nfft_one_over_modulus(fcs_float x, fcs_int der, const fcs_float *param)    /* K(x) = 1/|x| */
{
  fcs_float value=0.0;

  (void)param;

  if (fabs(x)<DBL_EPSILON) value=0.0;
  else switch (der)
  {
    case  0 : value=1.0/fabs(x); break;
    case  1 : value=-1/x/fabs(x); break;
    case  2 : value=2.0/pow(fabs(x),3.0); break;
    case  3 : value=-6.0/(x*x*x)/fabs(x); break;
    case  4 : value=24.0/pow(fabs(x),5.0); break;
    case  5 : value=-120.0/(x*x*x*x*x)/fabs(x); break;
    case  6 : value=720.0/pow(fabs(x),7.0); break;
    case  7 : value=-5040.0/(x*x*x*x*x*x*x)/fabs(x); break;
    case  8 : value=40320.0/pow(fabs(x),9.0); break;
    case  9 : value=-362880.0/(x*x*x*x*x*x*x*x*x)/fabs(x); break;
    case 10 : value=3628800.0/pow(fabs(x),11.0); break;
    case 11 : value=-39916800.0/pow(x,11.0)/fabs(x); break;
    case 12 : value=479001600.0/pow(fabs(x),13.0); break;
    default : value=0.0;
  }

  return value;
}

FCS_P2NFFT_KERNEL_TYPE ifcs_p2nfft_one_over_x(fcs_float x, fcs_int der, const fcs_float *param)    /* K(x) = 1/x */
{
  fcs_float value=0.0;

  (void)param;

  if (fabs(x)<DBL_EPSILON) value=0.0;
  else switch (der)
  {
    case  0 : value=1.0/x; break;
    case  1 : value=-1.0/(x*x); break;
    case  2 : value=2.0/(x*x*x); break;
    case  3 : value=-6.0/(x*x*x*x); break;
    case  4 : value=24.0/(x*x*x*x*x); break;
    case  5 : value=-120.0/(x*x*x*x*x*x); break;
    case  6 : value=720.0/(x*x*x*x*x*x*x); break;
    case  7 : value=-5040.0/(x*x*x*x*x*x*x*x); break;
    case  8 : value=40320.0/(x*x*x*x*x*x*x*x*x); break;
    case  9 : value=-362880.0/pow(x,10.0); break;
    case 10 : value=3628800.0/pow(x,11.0); break;
    case 11 : value=-39916800.0/pow(x,12.0); break;
    case 12 : value=479001600.0/pow(x,13.0); break;
    default : value=0.0;
  }

  return value;
}

FCS_P2NFFT_KERNEL_TYPE ifcs_p2nfft_inverse_multiquadric3(fcs_float x, fcs_int der, const fcs_float *param)    /* K(x) = 1/sqrt(x^2+c^2)^3 */
{
  fcs_float c=param[0];
  fcs_float value=0.0;

  switch (der)
  {
    case  0 : value=1.0/(sqrt(pow(x*x+c*c,3.0))); break;
    case  1 : value=-3.0/sqrt(pow(x*x+c*c,5.0))*x; break;
    case  2 : value=3.0*(4.0*x*x-c*c)/sqrt(pow(x*x+c*c,7.0)); break;
    case  3 : value=-15.0*x*(4.0*x*x-3.0*c*c)/sqrt(pow(x*x+c*c,9.0)); break;
    case  4 : value=45.0*(8.0*x*x*x*x-12.0*x*x*c*c+c*c*c*c)/sqrt(pow(x*x+c*c,11.0)); break;
    case  5 : value=-315.0*x*(8.0*x*x*x*x-20.0*x*x*c*c+5.0*c*c*c*c)/sqrt(pow(x*x+c*c,13.0)); break;
    case  6 : value=315.0*(64.0*x*x*x*x*x*x-240.0*x*x*x*x*c*c+120.0*x*x*c*c*c*c-5.0*c*c*c*c*c*c)/sqrt(pow(x*x+c*c,15.0)); break;
    case  7 : value=-2835.0*x*(64.0*x*x*x*x*x*x-336.0*x*x*x*x*c*c+280.0*x*x*c*c*c*c-35.0*c*c*c*c*c*c)/sqrt(pow(x*x+c*c,17.0)); break;
    case  8 : value=14175.0*(128.0*x*x*x*x*x*x*x*x-896.0*x*x*x*x*x*x*c*c+1120.0*x*x*x*x*c*c*c*c-280.0*x*x*c*c*c*c*c*c+7.0*c*c*c*c*c*c*c*c)/sqrt(pow(x*x+c*c,19.0)); break;
    case  9 : value=-155925.0*x*(128.0*x*x*x*x*x*x*x*x-1152.0*x*x*x*x*x*x*c*c+2016.0*x*x*x*x*c*c*c*c-840.0*x*x*c*c*c*c*c*c+63.0*c*c*c*c*c*c*c*c)/sqrt(pow(x*x+c*c,21.0)); break;
    case 10 : value=467775.0*(512.0*pow(x,10.0)-5760.0*x*x*x*x*x*x*x*x*c*c+13440.0*x*x*x*x*x*x*c*c*c*c-8400.0*x*x*x*x*c*c*c*c*c*c+1260.0*x*x*c*c*c*c*c*c*c*c-21.0*pow(c,10.0))/sqrt(pow(x*x+c*c,23.0)); break;
    case 11 : value=-6081075.0*x*(512.0*pow(x,10.0)-7040.0*x*x*x*x*x*x*x*x*c*c+21120.0*x*x*x*x*x*x*c*c*c*c-18480.0*x*x*x*x*c*c*c*c*c*c+4620.0*x*x*c*c*c*c*c*c*c*c-231.0*pow(c,10.0))/sqrt(pow(x*x+c*c,25.0)); break;
    case 12 : value=42567525.0*(1024.0*pow(x,12.0)+27720.0*x*x*x*x*c*c*c*c*c*c*c*c+33.0*pow(c,12.0)-2772.0*x*x*pow(c,10.0)-73920.0*x*x*x*x*x*x*c*c*c*c*c*c+63360.0*x*x*x*x*x*x*x*x*c*c*c*c-16896.0*pow(x,10.0)*c*c)/sqrt(pow(x*x+c*c,27.0)); break;
    default : value=0.0;
  }

  return value;
}

FCS_P2NFFT_KERNEL_TYPE ifcs_p2nfft_sinc_kernel(fcs_float x, fcs_int der, const fcs_float *param)    /* K(x) = sin(cx)/x */
{
  fcs_float c=param[0];
  fcs_float value=0.0;

  if (fabs(x)<DBL_EPSILON) value=0.0;
  else switch (der)
  {
    case  0 : value=sin(c*x)/x; break;
    case  1 : value=(cos(c*x)*c*x-sin(c*x))/(x*x); break;
    case  2 : value=-(sin(c*x)*c*c*x*x+2.0*cos(c*x)*c*x-2.0*sin(c*x))/(x*x*x); break;
    case  3 : value=-(cos(c*x)*c*c*c*x*x*x-3.0*sin(c*x)*c*c*x*x-6.0*cos(c*x)*c*x+6.0*sin(c*x))/(x*x*x*x); break;
    case  4 : value=(sin(c*x)*c*c*c*c*x*x*x*x+4.0*cos(c*x)*c*c*c*x*x*x-12.0*sin(c*x)*c*c*x*x-24.0*cos(c*x)*c*x+24.0*sin(c*x))/(x*x*x*x*x); break;
    case  5 : value=(cos(c*x)*c*c*c*c*c*x*x*x*x*x-5.0*sin(c*x)*c*c*c*c*x*x*x*x-20.0*cos(c*x)*c*c*c*x*x*x+60.0*sin(c*x)*c*c*x*x+120.0*cos(c*x)*c*x-120.0*sin(c*x))/(x*x*x*x*x*x); break;
    case  6 : value=-(sin(c*x)*c*c*c*c*c*c*x*x*x*x*x*x+6.0*cos(c*x)*c*c*c*c*c*x*x*x*x*x-30.0*sin(c*x)*c*c*c*c*x*x*x*x-120.0*cos(c*x)*c*c*c*x*x*x+360.0*sin(c*x)*c*c*x*x+720.0*cos(c*x)*c*x-720.0*sin(c*x))/(x*x*x*x*x*x*x); break;
    case  7 : value=-(cos(c*x)*c*c*c*c*c*c*c*x*x*x*x*x*x*x-7.0*sin(c*x)*c*c*c*c*c*c*x*x*x*x*x*x-42.0*cos(c*x)*c*c*c*c*c*x*x*x*x*x+210.0*sin(c*x)*c*c*c*c*x*x*x*x+840.0*cos(c*x)*c*c*c*x*x*x-2520.0*sin(c*x)*c*c*x*x-5040.0*cos(c*x)*c*x+5040.0*sin(c*x))/(x*x*x*x*x*x*x*x); break;
    case  8 : value=(sin(c*x)*c*c*c*c*c*c*c*c*x*x*x*x*x*x*x*x+8.0*cos(c*x)*c*c*c*c*c*c*c*x*x*x*x*x*x*x-56.0*sin(c*x)*c*c*c*c*c*c*x*x*x*x*x*x-336.0*cos(c*x)*c*c*c*c*c*x*x*x*x*x+1680.0*sin(c*x)*c*c*c*c*x*x*x*x+6720.0*cos(c*x)*c*c*c*x*x*x-20160.0*sin(c*x)*c*c*x*x-40320.0*cos(c*x)*c*x+40320.0*sin(c*x))/(x*x*x*x*x*x*x*x*x); break;
    case  9 : value=(cos(c*x)*c*c*c*c*c*c*c*c*c*x*x*x*x*x*x*x*x*x-9.0*sin(c*x)*c*c*c*c*c*c*c*c*x*x*x*x*x*x*x*x-72.0*cos(c*x)*c*c*c*c*c*c*c*x*x*x*x*x*x*x+504.0*sin(c*x)*c*c*c*c*c*c*x*x*x*x*x*x+3024.0*cos(c*x)*c*c*c*c*c*x*x*x*x*x-15120.0*sin(c*x)*c*c*c*c*x*x*x*x-60480.0*cos(c*x)*c*c*c*x*x*x+181440.0*sin(c*x)*c*c*x*x+362880.0*cos(c*x)*c*x-362880.0*sin(c*x))/pow(x,10.0); break;
    case 10 : value=-(sin(c*x)*pow(c,10.0)*pow(x,10.0)+10.0*cos(c*x)*c*c*c*c*c*c*c*c*c*x*x*x*x*x*x*x*x*x-90.0*sin(c*x)*c*c*c*c*c*c*c*c*x*x*x*x*x*x*x*x-720.0*cos(c*x)*c*c*c*c*c*c*c*x*x*x*x*x*x*x+5040.0*sin(c*x)*c*c*c*c*c*c*x*x*x*x*x*x+30240.0*cos(c*x)*c*c*c*c*c*x*x*x*x*x-151200.0*sin(c*x)*c*c*c*c*x*x*x*x-604800.0*cos(c*x)*c*c*c*x*x*x+1814400.0*sin(c*x)*c*c*x*x+3628800.0*cos(c*x)*c*x-3628800.0*sin(c*x))/pow(x,11.0); break;
    case 11 : value=-(cos(c*x)*pow(c,11.0)*pow(x,11.0)-11.0*sin(c*x)*pow(c,10.0)*pow(x,10.0)-110.0*cos(c*x)*c*c*c*c*c*c*c*c*c*x*x*x*x*x*x*x*x*x+990.0*sin(c*x)*c*c*c*c*c*c*c*c*x*x*x*x*x*x*x*x+7920.0*cos(c*x)*c*c*c*c*c*c*c*x*x*x*x*x*x*x-55440.0*sin(c*x)*c*c*c*c*c*c*x*x*x*x*x*x-332640.0*cos(c*x)*c*c*c*c*c*x*x*x*x*x+1663200.0*sin(c*x)*c*c*c*c*x*x*x*x+6652800.0*cos(c*x)*c*c*c*x*x*x-19958400.0*sin(c*x)*c*c*x*x-39916800.0*cos(c*x)*c*x+39916800.0*sin(c*x))/pow(x,12.0); break;
    case 12 : value=(sin(c*x)*pow(c,12.0)*pow(x,12.0)+12.0*cos(c*x)*pow(c,11.0)*pow(x,11.0)-132.0*sin(c*x)*pow(c,10.0)*pow(x,10.0)-1320.0*cos(c*x)*c*c*c*c*c*c*c*c*c*x*x*x*x*x*x*x*x*x+11880.0*sin(c*x)*c*c*c*c*c*c*c*c*x*x*x*x*x*x*x*x+95040.0*cos(c*x)*c*c*c*c*c*c*c*x*x*x*x*x*x*x-665280.0*sin(c*x)*c*c*c*c*c*c*x*x*x*x*x*x-3991680.0*cos(c*x)*c*c*c*c*c*x*x*x*x*x+19958400.0*sin(c*x)*c*c*c*c*x*x*x*x+79833600.0*cos(c*x)*c*c*c*x*x*x-239500800.0*sin(c*x)*c*c*x*x-479001600.0*cos(c*x)*c*x+479001600.0*sin(c*x))/pow(x,13.0); break;
    default : value=0.0;
  }

  return value;
}

FCS_P2NFFT_KERNEL_TYPE ifcs_p2nfft_cosc(fcs_float x, fcs_int der, const fcs_float *param)    /* K(x) = cos(cx)/x */
{
  fcs_float c=param[0];
  fcs_float value=0.0;
  fcs_float sign;

  if (x<0) sign=-1.0; else sign=1.0;
  x=fabs(x);

  if (fabs(x)<DBL_EPSILON) value=0.0;
  else switch (der)
  {
    case  0 : value=cos(c*x)/x; break;
    case  1 : value=-(sin(c*x)*c*x+cos(c*x))/(x*x); break;
    case  2 : value=(-cos(c*x)*c*c*x*x+2.0*sin(c*x)*c*x+2.0*cos(c*x))/(x*x*x); break;
    case  3 : value=(sin(c*x)*c*c*c*x*x*x+3.0*cos(c*x)*c*c*x*x-6.0*sin(c*x)*c*x-6.0*cos(c*x))/(x*x*x*x); break;
    case  4 : value=(cos(c*x)*c*c*c*c*x*x*x*x-4.0*sin(c*x)*c*c*c*x*x*x-12.0*cos(c*x)*c*c*x*x+24.0*sin(c*x)*c*x+24.0*cos(c*x))/(x*x*x*x*x); break;
    case  5 : value=-(sin(c*x)*c*c*c*c*c*x*x*x*x*x+5.0*cos(c*x)*c*c*c*c*x*x*x*x-20.0*sin(c*x)*c*c*c*x*x*x-60.0*cos(c*x)*c*c*x*x+120.0*sin(c*x)*c*x+120.0*cos(c*x))/(x*x*x*x*x*x); break;
    case  6 : value=-(cos(c*x)*c*c*c*c*c*c*x*x*x*x*x*x-6.0*sin(c*x)*c*c*c*c*c*x*x*x*x*x-30.0*cos(c*x)*c*c*c*c*x*x*x*x+120.0*sin(c*x)*c*c*c*x*x*x+360.0*cos(c*x)*c*c*x*x-720.0*sin(c*x)*c*x-720.0*cos(c*x))/(x*x*x*x*x*x*x); break;
    case  7 : value=(sin(c*x)*c*c*c*c*c*c*c*x*x*x*x*x*x*x+7.0*cos(c*x)*c*c*c*c*c*c*x*x*x*x*x*x-42.0*sin(c*x)*c*c*c*c*c*x*x*x*x*x-210.0*cos(c*x)*c*c*c*c*x*x*x*x+840.0*sin(c*x)*c*c*c*x*x*x+2520.0*cos(c*x)*c*c*x*x-5040.0*sin(c*x)*c*x-5040.0*cos(c*x))/(x*x*x*x*x*x*x*x); break;
    case  8 : value=(cos(c*x)*c*c*c*c*c*c*c*c*x*x*x*x*x*x*x*x-8.0*sin(c*x)*c*c*c*c*c*c*c*x*x*x*x*x*x*x-56.0*cos(c*x)*c*c*c*c*c*c*x*x*x*x*x*x+336.0*sin(c*x)*c*c*c*c*c*x*x*x*x*x+1680.0*cos(c*x)*c*c*c*c*x*x*x*x-6720.0*sin(c*x)*c*c*c*x*x*x-20160.0*cos(c*x)*c*c*x*x+40320.0*sin(c*x)*c*x+40320.0*cos(c*x))/(x*x*x*x*x*x*x*x*x); break;
    case  9 : value=-(sin(c*x)*c*c*c*c*c*c*c*c*c*x*x*x*x*x*x*x*x*x+9.0*cos(c*x)*c*c*c*c*c*c*c*c*x*x*x*x*x*x*x*x-72.0*sin(c*x)*c*c*c*c*c*c*c*x*x*x*x*x*x*x-504.0*cos(c*x)*c*c*c*c*c*c*x*x*x*x*x*x+3024.0*sin(c*x)*c*c*c*c*c*x*x*x*x*x+15120.0*cos(c*x)*c*c*c*c*x*x*x*x-60480.0*sin(c*x)*c*c*c*x*x*x-181440.0*cos(c*x)*c*c*x*x+362880.0*sin(c*x)*c*x+362880.0*cos(c*x))/pow(x,10.0); break;
    case 10 : value=-(cos(c*x)*pow(c,10.0)*pow(x,10.0)-10.0*sin(c*x)*c*c*c*c*c*c*c*c*c*x*x*x*x*x*x*x*x*x-90.0*cos(c*x)*c*c*c*c*c*c*c*c*x*x*x*x*x*x*x*x+720.0*sin(c*x)*c*c*c*c*c*c*c*x*x*x*x*x*x*x+5040.0*cos(c*x)*c*c*c*c*c*c*x*x*x*x*x*x-30240.0*sin(c*x)*c*c*c*c*c*x*x*x*x*x-151200.0*cos(c*x)*c*c*c*c*x*x*x*x+604800.0*sin(c*x)*c*c*c*x*x*x+1814400.0*cos(c*x)*c*c*x*x-3628800.0*sin(c*x)*c*x-3628800.0*cos(c*x))/pow(x,11.0); break;
    case 11 : value=(sin(c*x)*pow(c,11.0)*pow(x,11.0)+11.0*cos(c*x)*pow(c,10.0)*pow(x,10.0)-110.0*sin(c*x)*c*c*c*c*c*c*c*c*c*x*x*x*x*x*x*x*x*x-990.0*cos(c*x)*c*c*c*c*c*c*c*c*x*x*x*x*x*x*x*x+7920.0*sin(c*x)*c*c*c*c*c*c*c*x*x*x*x*x*x*x+55440.0*cos(c*x)*c*c*c*c*c*c*x*x*x*x*x*x-332640.0*sin(c*x)*c*c*c*c*c*x*x*x*x*x-1663200.0*cos(c*x)*c*c*c*c*x*x*x*x+6652800.0*sin(c*x)*c*c*c*x*x*x+19958400.0*cos(c*x)*c*c*x*x-39916800.0*sin(c*x)*c*x-39916800.0*cos(c*x))/pow(x,12.0); break;
    case 12 : value=(cos(c*x)*pow(c,12.0)*pow(x,12.0)-12.0*sin(c*x)*pow(c,11.0)*pow(x,11.0)-132.0*cos(c*x)*pow(c,10.0)*pow(x,10.0)+1320.0*sin(c*x)*c*c*c*c*c*c*c*c*c*x*x*x*x*x*x*x*x*x+11880.0*cos(c*x)*c*c*c*c*c*c*c*c*x*x*x*x*x*x*x*x-95040.0*sin(c*x)*c*c*c*c*c*c*c*x*x*x*x*x*x*x-665280.0*cos(c*x)*c*c*c*c*c*c*x*x*x*x*x*x+3991680.0*sin(c*x)*c*c*c*c*c*x*x*x*x*x+19958400.0*cos(c*x)*c*c*c*c*x*x*x*x-79833600.0*sin(c*x)*c*c*c*x*x*x-239500800.0*cos(c*x)*c*c*x*x+479001600.0*sin(c*x)*c*x+479001600.0*cos(c*x))/pow(x,13.0); break;
    default : value=0.0;
  }
  value*=pow(sign,der);

  return value;
}

FCS_P2NFFT_KERNEL_TYPE ifcs_p2nfft_kcot(fcs_float x, fcs_int der, const fcs_float *param)   /* K(x) = cot(cx) */
{
  fcs_float c=param[0];
  fcs_float value=0.0;

  if (fabs(x)<DBL_EPSILON) value=0.0;
  else switch (der)
  {
    case  0 : value = 1.0/tan(c * x); break;
    case  1 : value = -(1.0 + pow(1.0/tan(c * x), 2.0)) * c; break;
    case  2 : value = 2.0 * 1.0/tan(c * x) * (1.0 + pow(1.0/tan(c * x), 2.0)) * c * c; break;
    case  3 : value = -2.0 * (1.0 + pow(1.0/tan(c * x), 2.0)) * pow(c, 3.0) * (1.0 + 3.0 * pow(1.0/tan(c * x), 2.0)); break;
    case  4 : value = 8.0 * (1.0 + pow(1.0/tan(c * x), 2.0)) * pow(c, 4.0) * 1.0/tan(c * x) * (2.0 + 3.0 * pow(1.0/tan(c * x), 2.0)); break;
    case  5 : value = -0.8e1 * (0.1e1 + pow(1.0/tan(c * x), 0.2e1)) * pow(c, 0.5e1) * (0.15e2 * pow(1.0/tan(c * x), 0.2e1) + 0.15e2 * pow(1.0/tan(c * x), 0.4e1) + 0.2e1); break;
    case  6 : value = 0.16e2 * (0.1e1 + pow(1.0/tan(c * x), 0.2e1)) * pow(c, 0.6e1) * 1.0/tan(c * x) * (0.60e2 * pow(1.0/tan(c * x), 0.2e1) + 0.45e2 * pow(1.0/tan(c * x), 0.4e1) + 0.17e2); break;
    case  7 : value = -0.16e2 * (0.1e1 + pow(1.0/tan(c * x), 0.2e1)) * pow(c, 0.7e1) * (0.525e3 * pow(1.0/tan(c * x), 0.4e1) + 0.315e3 * pow(1.0/tan(c * x), 0.6e1) + 0.231e3 * pow(1.0/tan(c * x), 0.2e1) + 0.17e2); break;
    case  8 : value = 0.128e3 * (0.1e1 + pow(1.0/tan(c * x), 0.2e1)) * pow(c, 0.8e1) * 1.0/tan(c * x) * (0.630e3 * pow(1.0/tan(c * x), 0.4e1) + 0.315e3 * pow(1.0/tan(c * x), 0.6e1) + 0.378e3 * pow(1.0/tan(c * x), 0.2e1) + 0.62e2); break;
    case  9 : value = -0.128e3 * (0.1e1 + pow(1.0/tan(c * x), 0.2e1)) * pow(c, 0.9e1) * (0.6615e4 * pow(1.0/tan(c * x), 0.6e1) + 0.2835e4 * pow(1.0/tan(c * x), 0.8e1) + 0.5040e4 * pow(1.0/tan(c * x), 0.4e1) + 0.1320e4 * pow(1.0/tan(c * x), 0.2e1) + 0.62e2); break;
    case 10 : value = 0.256e3 * (0.1e1 + pow(1.0/tan(c * x), 0.2e1)) * pow(c, 0.10e2) * 1.0/tan(c * x) * (0.37800e5 * pow(1.0/tan(c * x), 0.6e1) + 0.14175e5 * pow(1.0/tan(c * x), 0.8e1) + 0.34965e5 * pow(1.0/tan(c * x), 0.4e1) + 0.12720e5 * pow(1.0/tan(c * x), 0.2e1) + 0.1382e4); break;
    case 11 : value = -0.256e3 * (0.1e1 + pow(1.0/tan(c * x), 0.2e1)) * pow(c, 0.11e2) * (0.467775e6 * pow(1.0/tan(c * x), 0.8e1) + 0.155925e6 * pow(1.0/tan(c * x), 0.10e2) + 0.509355e6 * pow(1.0/tan(c * x), 0.6e1) + 0.238425e6 * pow(1.0/tan(c * x), 0.4e1) + 0.42306e5 * pow(1.0/tan(c * x), 0.2e1) + 0.1382e4); break;
    case 12 : value = 0.1024e4 * (0.1e1 + pow(1.0/tan(c * x), 0.2e1)) * pow(c, 0.12e2) * 1.0/tan(c * x) * (0.1559250e7 * pow(1.0/tan(c * x), 0.8e1) + 0.467775e6 * pow(1.0/tan(c * x), 0.10e2) + 0.1954260e7 * pow(1.0/tan(c * x), 0.6e1) + 0.1121670e7 * pow(1.0/tan(c * x), 0.4e1) + 0.280731e6 * pow(1.0/tan(c * x), 0.2e1) + 0.21844e5); break;
    default : value=0.0;
  }

  return value;
}


FCS_P2NFFT_KERNEL_TYPE ifcs_p2nfft_one_over_cube(fcs_float x, fcs_int der, const fcs_float *param) /* K(x) = 1/x^3 */
{
  fcs_float value=0.0;

  if (fabs(x)<DBL_EPSILON) value=0.0;
  else switch (der)
  {
    case  0 : value = 1.0/(x*x*x); break;
    case  1 : value = -3.0/(x*x*x*x); break;
    case  2 : value = 12.0/(x*x*x*x*x); break;
    case  3 : value = -60.0/(x*x*x*x*x*x); break;
    case  4 : value = 360.0/(x*x*x*x*x*x*x); break;
    case  5 : value = -2520.0/(x*x*x*x*x*x*x*x); break;
    case  6 : value = 20160.0/pow(x, 9.0); break;
    case  7 : value = -181440.0/pow(x, 10.0); break;
    case  8 : value = 1814400.0/pow(x, 11.0); break;
    case  9 : value = -19958400.0/pow(x, 12.0); break;
    case  10 : value = 239500800.0/pow(x, 13.0); break;
    case  11 : value = -3113510400.0/pow(x, 14.0); break;
    case  12 : value = 43589145600.0/pow(x, 15.0); break;
    default : value=0.0;
  }

  return value;
}

static fcs_float theta(
    fcs_float x, fcs_float k, fcs_float alpha
    )
{
  return exp(2*FCS_PI*k*x) * (1-erf(FCS_PI*k/alpha + alpha*x)); /* use erf instead of erfc to fix ICC performance problems */
}

static fcs_float theta_p(
    fcs_float x, fcs_float k, fcs_float alpha
    )
{
  return theta(x,k,alpha) + theta(-x,k,alpha);
}

static fcs_float theta_m(
    fcs_float x, fcs_float k, fcs_float alpha
    )
{
  return theta(x,k,alpha) - theta(-x,k,alpha);
}

/* Fourier space part of 2d-periodic Ewald and k<>0 */
FCS_P2NFFT_KERNEL_TYPE ifcs_p2nfft_ewald_2dp_kneq0(fcs_float x, fcs_int der, const fcs_float *param) /* K(x) = exp(2*pi*k*x) * erf(pi*k/alpha + alpha*x) */
{
  fcs_float alpha = param[0]; /* Ewald splitting parameter alpha */
  fcs_float k     = param[1]; /* norm of (k_0/B_0,k_1/B_1) */
  fcs_float z     = param[2]; /* norm of (x_2/box_scale[2]) */
  fcs_float one_over_alpha = 1.0 / alpha;

  fcs_float b = FCS_PI*k/alpha;
  fcs_float c = 2*FCS_PI*k;

  fcs_float value=0.0;

  if(fcs_float_is_zero(k))
    return 0.0;

  switch (der)
  {
    case  0 : value = theta_p(x,k,alpha); break;
    case  1 : value = c*theta_m(x,k,alpha); break;
    default : value = c*c * ifcs_p2nfft_ewald_2dp_kneq0(x,der-2,param) - 8*alpha*FCS_SQRTPI*k* exp(-b*b) * ifcs_p2nfft_gaussian(x,der-2,&one_over_alpha);
  }

  return value;
}

/* Fourier space part of 2d-periodic Ewald and k==0 */
FCS_P2NFFT_KERNEL_TYPE ifcs_p2nfft_ewald_2dp_keq0(fcs_float x, fcs_int der, const fcs_float *param) /* K(x) = 1/alpha * exp(alpha*x) + sqrt(pi)*x*erf(alpha*x) */
{
  fcs_float alpha = param[0]; /* Ewald splitting parameter alpha */
  fcs_float one_over_alpha = 1.0 / alpha;

  return one_over_alpha * ifcs_p2nfft_gaussian(x, der, &one_over_alpha) + ifcs_p2nfft_x_times_erf(x, der, &alpha);
}

FCS_P2NFFT_KERNEL_TYPE ifcs_p2nfft_x_times_erf(fcs_float x, fcs_int der, const fcs_float *param) /* K(x) = sqrt(PI) * x * erf(alpha*x) */
{
  fcs_float alpha=param[0]; /* Ewald splitting parameter alpha */

  /* substitute x := alpha*x */
  x = alpha * x;

  fcs_float value=0.0;
  switch (der)
  {
    /* This code was generated with Maple */
    case  0 : value = sqrt(FCS_PI) * erf(x) * x; break;
    case  1 : value = 0.2e1 * exp(-x * x) * x + sqrt(FCS_PI) * erf(x); break;
    case  2 : value = -0.4e1 * exp(-x * x) * (x * x - 0.1e1); break;
    case  3 : value = 0.8e1 * exp(-x * x) * (x * x - 0.2e1) * x; break;
    case  4 : value = -0.8e1 * exp(-x * x) * ((0.2e1 * x * x - 0.7e1) * x * x + 0.2e1); break;
    case  5 : value = 0.16e2 * exp(-x * x) * ((0.2e1 * x * x - 0.11e2) * x * x + 0.9e1) * x; break;
    case  6 : value = -0.16e2 * exp(-x * x) * (((0.4e1 * x * x - 0.32e2) * x * x + 0.51e2) * x * x - 0.9e1); break;
    case  7 : value = 0.32e2 * exp(-x * x) * (((0.4e1 * x * x - 0.44e2) * x * x + 0.115e3) * x * x - 0.60e2) * x; break;
    case  8 : value = -0.32e2 * exp(-x * x) * ((((0.8e1 * x * x - 0.116e3) * x * x + 0.450e3) * x * x - 0.465e3) * x * x + 0.60e2); break;
    case  9 : value = 0.64e2 * exp(-x * x) * ((((0.8e1 * x * x - 0.148e3) * x * x + 0.798e3) * x * x - 0.1365e4) * x * x + 0.525e3) * x; break;
    case 10 : value = -0.64e2 * exp(-x * x) * (((((0.16e2 * x * x - 0.368e3) * x * x + 0.2632e4) * x * x - 0.6720e4) * x * x + 0.5145e4) * x * x - 0.525e3); break;
    case 11 : value = 0.128e3 * exp(-x * x) * (((((0.16e2 * x * x - 0.448e3) * x * x + 0.4104e4) * x * x - 0.14616e5) * x * x + 0.18585e5) * x * x - 0.5670e4) * x; break;
    case 12 : value = -0.128e3 * exp(-x * x) * ((((((0.32e2 * x * x - 0.1072e4) * x * x + 0.12240e5) * x * x - 0.57960e5) * x * x + 0.110250e6) * x * x - 0.67095e5) * x * x + 0.5670e4); break;
    case 13 : value = 0.256e3 * exp(-x * x) * ((((((0.32e2 * x * x - 0.1264e4) * x * x + 0.17600e5) * x * x - 0.106920e6) * x * x + 0.284130e6) * x * x - 0.287595e6) * x * x + 0.72765e5) * x; break;
    case 14 : value = -0.256e3 * exp(-x * x) * (((((((0.64e2 * x * x - 0.2944e4) * x * x + 0.49104e5) * x * x - 0.372240e6) * x * x + 0.1316700e7) * x * x - 0.1995840e7) * x * x + 0.1008315e7) * x * x - 0.72765e5); break;
    case 15 : value = 0.512e3 * exp(-x * x) * (((((((0.64e2 * x * x - 0.3392e4) * x * x + 0.66768e5) * x * x - 0.617760e6) * x * x + 0.2805660e7) * x * x - 0.5945940e7) * x * x + 0.4999995e7) * x * x - 0.1081080e7) * x; break;
    case 16 : value = -0.512e3 * exp(-x * x) * ((((((((0.128e3 * x * x - 0.7744e4) * x * x + 0.177632e6) * x * x - 0.1969968e7) * x * x + 0.11171160e8) * x * x - 0.31531500e8) * x * x + 0.39729690e8) * x * x - 0.17162145e8) * x * x + 0.1081080e7); break;
    case 17 : value = 0.1024e4 * exp(-x * x) * ((((((((0.128e3 * x * x - 0.8768e4) * x * x + 0.231840e6) * x * x - 0.3035760e7) * x * x + 0.21021000e8) * x * x - 0.76216140e8) * x * x + 0.134324190e9) * x * x - 0.96621525e8) * x * x + 0.18243225e8) * x; break;
    case 18 : value = -0.1024e4 * exp(-x * x) * (((((((((0.256e3 * x * x - 0.19712e5) * x * x + 0.595200e6) * x * x - 0.9085440e7) * x * x + 0.75435360e8) * x * x - 0.341621280e9) * x * x + 0.802161360e9) * x * x - 0.864864000e9) * x * x + 0.326351025e9) * x * x - 0.18243225e8); break;
    case 19 : value = 0.2048e4 * exp(-x * x) * (((((((((0.256e3 * x * x - 0.22016e5) * x * x + 0.752896e6) * x * x - 0.13251840e8) * x * x + 0.129948000e9) * x * x - 0.718798080e9) * x * x + 0.2168646480e10) * x * x - 0.3271348080e10) * x * x + 0.2056079025e10) * x * x - 0.344594250e9) * x; break;
    case 20 : value = -0.2048e4 * exp(-x * x) * ((((((((((0.512e3 * x * x - 0.48896e5) * x * x + 0.1880064e7) * x * x - 0.37797120e8) * x * x + 0.432169920e9) * x * x - 0.2867024160e10) * x * x + 0.10806475680e11) * x * x - 0.21723221520e11) * x * x + 0.20468898450e11) * x * x - 0.6857425575e10) * x * x + 0.344594250e9); break;
    case 21 : value = 0.4096e4 * exp(-x * x) * ((((((((((0.512e3 * x * x - 0.54016e5) * x * x + 0.2320128e7) * x * x - 0.52837632e8) * x * x + 0.696749760e9) * x * x - 0.5460043680e10) * x * x + 0.25141596480e11) * x * x - 0.64949124240e11) * x * x + 0.85638563010e11) * x * x - 0.47795222475e11) * x * x + 0.7202019825e10) * x; break;
    default : value = 0.0;
  }

  return pow(alpha, der) * value;
}

