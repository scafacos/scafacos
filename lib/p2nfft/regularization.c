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

#include <stdlib.h>
#include <math.h>
#ifdef HAVE_COMPLEX_H
#include <complex.h>
#endif

#include "regularization.h"

static fcs_float intpol_even_const(
    const fcs_float x, const fcs_float *table,
    const fcs_int num_nodes, const fcs_float one_over_eps_I);
static fcs_float intpol_even_lin(
    const fcs_float x, const fcs_float *table,
    const fcs_int num_nodes, const fcs_float one_over_eps_I);
static fcs_float intpol_even_quad(
    const fcs_float x, const fcs_float *table,
    const fcs_int num_nodes, const fcs_float one_over_eps_I);
static fcs_float intpol_even_cub(
    const fcs_float x, const fcs_float *table,
    const fcs_int num_nodes, const fcs_float one_over_eps_I);




/** factorial */
static fcs_float fak(fcs_int n)
{
  if (n<=1) return 1.0;
  else return (fcs_float)n*fak(n-1);
}

/** binomial coefficient */
static fcs_float binom(fcs_int n, fcs_int m)
{
  return fak(n)/fak(m)/fak(n-m);
}

/** basis polynomial for regularized kernel */
static fcs_float BasisPoly(fcs_int m, fcs_int r, fcs_float xx)
{
  fcs_int k;
  fcs_float sum=0.0;

  for (k=0; k<=m-r; k++) {
    sum+=binom(m+k,k)*pow((xx+1.0)/2.0,(fcs_float)k);
  }
  return sum*pow((xx+1.0),(fcs_float)r)*pow(1.0-xx,(fcs_float)(m+1))/(1<<(m+1))/fak(r); /* 1<<(m+1) = 2^(m+1) */
}

/** integrated basis polynomial for regularized kernel */
static fcs_float IntBasisPoly(fcs_int p, fcs_int j, fcs_float y)
{
  fcs_int k,l;
  fcs_float sum1=0.0, sum2=0.0;
  
  for (l=0; l<=p; l++) {
    sum1 = 0.0;
    for (k=0; k<=p-j-1; k++) {
      sum1 += binom(p+k-1,k)*fak(j+k)/fak(j+k+1+l)*pow((1.0+y)/2.0,(fcs_float)k);
    }
    sum2 += pow(1.0+y,(fcs_float)l)*pow(1.0-y,(fcs_float)(p-l))/fak(p-l)*sum1;
  }
  return sum2 * fak(p)/fak(j)/(1<<p)*pow(1.0+y,(fcs_float)(j+1)); /* 1<<p = 2^p */
}


/** regularized kernel for even kernels with K_I even
 *  and K_B mirrored smooth to K(1/2) (used in dD, d>1)
 */
// static double _Complex regkern3(ifcs_p2nfft_kernel k, fcs_float xx, fcs_int p, const fcs_float *param, fcs_float a, fcs_float b)
// {
//   fcs_int r;
//   double _Complex sum=0.0;
// 
//   xx=fabs(xx);
// 
//   /* constant continuation for radii > 0.5 */
//   if (xx>=0.5)
//     xx=0.5;
// 
//   /* regularization at farfield border */
//   if ( 0.5-b < xx ) {
//     sum=k(0.5,0,param)*BasisPoly(p-1,0,-2.0*xx/b+(1.0-b)/b);
//     /* sum=regkern2(typ,c,p,a,b, 0.5)*BasisPoly(p-1,0,-2.0*xx/b+(1.0-b)/b); */
//     for (r=0; r<p; r++) {
//       sum+=pow(b/2.0,(fcs_float)r)
//           *k(0.5-b,r,param)
//           *BasisPoly(p-1,r,2.0*xx/b-(1.0-b)/b);
//     }
//     return sum;
//   }
//   
//   /* farfield: original kernel function */ 
//   if ( a <= xx )
//     return k(xx,0,param);
//   
//   /* nearfield regularization */
//   for (r=0; r<p; r++) {
//     sum+=pow(-a,(fcs_float)r)*k(a,r,param)
//         *(BasisPoly(p-1,r,xx/a)+BasisPoly(p-1,r,-xx/a));
//   }
//   return sum;
// }


/*
 * Calculates regularized kernel 1/x at x. Uses 2p-taylor polynomial for
 * near field part.
 * Implemented for p=3,4,5,6.
 */
fcs_float ifcs_regkern3_2ptaylor_x_inv(fcs_float *x, fcs_int p, fcs_float a, fcs_float b, fcs_float *n)
{
  fcs_int r;
  fcs_int t;
  fcs_float sum=0.0;
  fcs_float xx_2=0.0;
  fcs_float xx_e=0.0;
  fcs_float xx_outer=0.0;
  fcs_float param[] = {0.0};

  for (t=0; t<3; t++) {
    xx_2 += x[t] * x[t];
    xx_e += x[t] * x[t] / n[t] / n[t];
    if (n[t] / 2.0 > xx_outer)
      xx_outer = n[t] / 2.0;
  }

  xx_2 = fcs_sqrt(xx_2);
  xx_e = fcs_sqrt(xx_e);

  if (xx_e>=0.5) {
    return ifcs_p2nfft_one_over_modulus(xx_outer,0,param);
  }
  /* else */
  if ((a<=xx_2) && (xx_e<=0.5-b)) {
    return ifcs_p2nfft_one_over_modulus(xx_2,0,param);
  }
  else if (xx_2<a) {
    for (r=0; r<p; r++) {
      sum += fcs_pow(-a,(fcs_float)r) * ifcs_p2nfft_one_over_modulus(a,r,param)
          * (BasisPoly(p-1,r,xx_2/a)+BasisPoly(p-1,r,-xx_2/a));
    }
    return sum;
  }
  else if ((0.5-b<xx_e) && (xx_e<=0.5)) {
    fcs_float beta;
    fcs_float x_outer;

    x_outer = 0.5 / xx_e * xx_2;
    beta = x_outer - (0.5 - b) / xx_e * xx_2;

    sum = ifcs_p2nfft_one_over_modulus(xx_outer,0,param)*BasisPoly(p-1,0,-2.0*xx_2/beta+(2*x_outer-beta)/beta);
    for (r=0; r<p; r++) {
      sum += fcs_pow(beta/2.0,(fcs_float)r)
          * ifcs_p2nfft_one_over_modulus(x_outer-beta,r,param)
          * BasisPoly(p-1,r,2.0*xx_2/beta-(2*x_outer-beta)/beta);
    }
    return sum;
  }
  return 0.0;
}





/** regularized kernel for even kernels with K_I even
 *  and K_B mirrored smooth into x=1/2 (used in dD, d>1)
 */
// static fcs_float const_regkern4(ifcs_p2nfft_kernel k, fcs_int p, const fcs_float *param, fcs_float b)
// {
//   fcs_int j;
//   fcs_float xx, sum=0.0;
// 
//   xx=0.5-b;
//   for (j=0; j<=p-2; j++) {
//     sum+=pow(b/2.0,(fcs_float)j+1)
//         *k(0.5-b,j+1,param)
//         *IntBasisPoly(p-1,j,2.0*xx/b-(1.0-b)/b);
//   }
//   return k(xx,0,param)-sum;
// }


/** regularized kernel for even kernels with K_I even
 *  and K_B mirrored smooth into x=1/2 (used in dD, d>1)
 */
static fcs_float regkern4(ifcs_p2nfft_kernel k, fcs_float xx, fcs_int p, const fcs_float *param, fcs_float a, fcs_float b)
{
  fcs_int j;
  fcs_float sum=0.0;

  xx=fabs(xx);

  /* constant continuation for radii > 0.5 */
  if (xx>=0.5)
    xx=0.5;

  /* regularization at farfield border */
  if ( 0.5-b<xx ) {
    for (j=0; j<=p-2; j++) {
      sum += 
        creal(pow(b/2.0,(fcs_float)j+1)
          *k(0.5-b,j+1,param)
          *(IntBasisPoly(p-1,j,2.0*xx/b-(1.0-b)/b)
              -IntBasisPoly(p-1,j,-1))
//       *IntBasisPoly(p-1,j,2.0*xx/b-(1.0-b)/b)
       );
    }
//     return sum+const_regkern4(k, p, param, b);
    return sum + k(0.5-b,0,param);
  }
 
  /* farfield: original kernel function */ 
  if ( a<=xx )
    return k(xx,0,param);
  
  /* nearfield regularization */
  for (j=0; j<p; j++) {
    sum+=pow(-a,(fcs_float)j)*k(a,j,param)
        *(BasisPoly(p-1,j,xx/a)+BasisPoly(p-1,j,-xx/a));
  }
  return sum;
} 


fcs_float ifcs_p2nfft_regkernel(
    ifcs_p2nfft_kernel k, fcs_float xx, fcs_int p, const fcs_float *param, fcs_float epsI, fcs_float epsB
    )
{
//  return creal(regkern3(k, xx, p, param, a, b));
  return regkern4(k, xx, p, param, epsI, epsB);
}


fcs_float ifcs_p2nfft_interpolation(
    fcs_float x, fcs_float one_over_epsI,
    fcs_int order, fcs_int num_nodes,
    const fcs_float *table
    )
{
  switch(order){
    case 0: return intpol_even_const(x, table, num_nodes, one_over_epsI);
    case 1: return intpol_even_lin(x, table, num_nodes, one_over_epsI);
    case 2: return intpol_even_quad(x, table, num_nodes, one_over_epsI);
    default: return intpol_even_cub(x, table, num_nodes, one_over_epsI);
  }
}




/** regularized kernel for even kernels with K_I even
 *  and K_B mirrored smooth into x=1/2 (used in dD, d>1)
 */
fcs_float ifcs_p2nfft_nearfield_correction_taylor2p(
    fcs_float xx, fcs_int p, const fcs_float *param
    )
{
  fcs_float horner = param[0];

  xx=fabs(xx);

  /* use horner scheme for next p-2 coefficients */
  for(fcs_int t=1; t<p; t++)
    horner = xx * xx * horner + param[t];

  return horner;
}

/** regularized kernel for even kernels with K_I even
 *  and K_B mirrored smooth into x=1/2 (used in dD, d>1)
 */
fcs_float ifcs_p2nfft_nearfield_correction_taylor2p_derive(
    fcs_float xx, fcs_int p, const fcs_float *param
    )
{
  fcs_float horner = param[0];

  xx=fabs(xx);

  /* use horner schema for next p-2 coefficients */
  for(fcs_int t=1; t<p-1; t++)
    horner = xx * xx * horner + param[t];

  /* multiply last time by xx, since we have an odd polynomial */
  return horner * xx;
}


fcs_int ifcs_p2nfft_load_taylor2p_coefficients(
   fcs_float epsI, fcs_int p,
   fcs_float *param
   )
{
  fcs_float a = epsI;

  /* Case p<2 is senseless, for p>16 we need to add the coefficients. */
  /* These coefficients are calculated analytically with a Maple script. */ 
  switch(p){
    case 2:
      param[0] = -1/(a*a*a)/2.0;
      param[1] = 3.0/2.0/a;
      break;
    case 3:
      param[0] = 3.0/8.0/(a*a*a*a*a);
      param[1] = -5.0/4.0/(a*a*a);
      param[2] = 15.0/8.0/a;
      break;
    case 4:
      param[0] = -5.0/16.0/(a*a*a*a*a*a*a);
      param[1] = 21.0/16.0/(a*a*a*a*a);
      param[2] = -35.0/16.0/(a*a*a);
      param[3] = 35.0/16.0/a;
      break;
    case 5:
      param[0] = 35.0/128.0/(a*a*a*a*a*a*a*a*a);
      param[1] = -45.0/32.0/(a*a*a*a*a*a*a);
      param[2] = 189.0/64.0/(a*a*a*a*a);
      param[3] = -105.0/32.0/(a*a*a);
      param[4] = 315.0/128.0/a;
      break;
    case 6:
      param[0] = -63.0/256.0/pow(a,11.0);
      param[1] = 385.0/256.0/(a*a*a*a*a*a*a*a*a);
      param[2] = -495.0/128.0/(a*a*a*a*a*a*a);
      param[3] = 693.0/128.0/(a*a*a*a*a);
      param[4] = -1155.0/256.0/(a*a*a);
      param[5] = 693.0/256.0/a;
      break;
    case 7:
      param[0] = 231.0/1024.0/pow(a,13.0);
      param[1] = -819.0/512.0/pow(a,11.0);
      param[2] = 5005.0/1024.0/(a*a*a*a*a*a*a*a*a);
      param[3] = -2145.0/256.0/(a*a*a*a*a*a*a);
      param[4] = 9009.0/1024.0/(a*a*a*a*a);
      param[5] = -3003.0/512.0/(a*a*a);
      param[6] = 3003.0/1024.0/a;
      break;
    case 8:
      param[0] = -429.0/2048.0/pow(a,15.0);
      param[1] = 3465.0/2048.0/pow(a,13.0);
      param[2] = -12285.0/2048.0/pow(a,11.0);
      param[3] = 25025.0/2048.0/(a*a*a*a*a*a*a*a*a);
      param[4] = -32175.0/2048.0/(a*a*a*a*a*a*a);
      param[5] = 27027.0/2048.0/(a*a*a*a*a);
      param[6] = -15015.0/2048.0/(a*a*a);
      param[7] = 6435.0/2048.0/a;
      break;
    case 9:
      param[0] = 6435.0/32768.0/pow(a,17.0);
      param[1] = -7293.0/4096.0/pow(a,15.0);
      param[2] = 58905.0/8192.0/pow(a,13.0);
      param[3] = -69615.0/4096.0/pow(a,11.0);
      param[4] = 425425.0/16384.0/(a*a*a*a*a*a*a*a*a);
      param[5] = -109395.0/4096.0/(a*a*a*a*a*a*a);
      param[6] = 153153.0/8192.0/(a*a*a*a*a);
      param[7] = -36465.0/4096.0/(a*a*a);
      param[8] = 109395.0/32768.0/a;
      break;
    case 10:
      param[0] = -12155.0/65536.0/pow(a,19.0);
      param[1] = 122265.0/65536.0/pow(a,17.0);
      param[2] = -138567.0/16384.0/pow(a,15.0);
      param[3] = 373065.0/16384.0/pow(a,13.0);
      param[4] = -1322685.0/32768.0/pow(a,11.0);
      param[5] = 1616615.0/32768.0/(a*a*a*a*a*a*a*a*a);
      param[6] = -692835.0/16384.0/(a*a*a*a*a*a*a);
      param[7] = 415701.0/16384.0/(a*a*a*a*a);
      param[8] = -692835.0/65536.0/(a*a*a);
      param[9] = 230945.0/65536.0/a;
      break;
    case 11:
      param[0] = 46189.0/262144.0/pow(a,21.0);
      param[1] = -255255.0/131072.0/pow(a,19.0);
      param[2] = 2567565.0/262144.0/pow(a,17.0);
      param[3] = -969969.0/32768.0/pow(a,15.0);
      param[4] = 7834365.0/131072.0/pow(a,13.0);
      param[5] = -5555277.0/65536.0/pow(a,11.0);
      param[6] = 11316305.0/131072.0/(a*a*a*a*a*a*a*a*a);
      param[7] = -2078505.0/32768.0/(a*a*a*a*a*a*a);
      param[8] = 8729721.0/262144.0/(a*a*a*a*a);
      param[9] = -1616615.0/131072.0/(a*a*a);
      param[10] = 969969.0/262144.0/a;
      break;
    case 12:
      param[0] = -88179.0/524288.0/pow(a,23.0);
      param[1] = 1062347.0/524288.0/pow(a,21.0);
      param[2] = -5870865.0/524288.0/pow(a,19.0);
      param[3] = 19684665.0/524288.0/pow(a,17.0);
      param[4] = -22309287.0/262144.0/pow(a,15.0);
      param[5] = 36038079.0/262144.0/pow(a,13.0);
      param[6] = -42590457.0/262144.0/pow(a,11.0);
      param[7] = 37182145.0/262144.0/(a*a*a*a*a*a*a*a*a);
      param[8] = -47805615.0/524288.0/(a*a*a*a*a*a*a);
      param[9] = 22309287.0/524288.0/(a*a*a*a*a);
      param[10] = -7436429.0/524288.0/(a*a*a);
      param[11] = 2028117.0/524288.0/a;
      break;
    case 13:
      param[0] = 676039.0/4194304.0/pow(a,25.0);
      param[1] = -2204475.0/1048576.0/pow(a,23.0);
      param[2] = 26558675.0/2097152.0/pow(a,21.0);
      param[3] = -48923875.0/1048576.0/pow(a,19.0);
      param[4] = 492116625.0/4194304.0/pow(a,17.0);
      param[5] = -111546435.0/524288.0/pow(a,15.0);
      param[6] = 300317325.0/1048576.0/pow(a,13.0);
      param[7] = -152108775.0/524288.0/pow(a,11.0);
      param[8] = 929553625.0/4194304.0/(a*a*a*a*a*a*a*a*a);
      param[9] = -132793375.0/1048576.0/(a*a*a*a*a*a*a);
      param[10] = 111546435.0/2097152.0/(a*a*a*a*a);
      param[11] = -16900975.0/1048576.0/(a*a*a);
      param[12] = 16900975.0/4194304.0/a;
      break;
    case 14:
      param[0] = -1300075.0/8388608.0/pow(a,27.0);
      param[1] = 18253053.0/8388608.0/pow(a,25.0);
      param[2] = -59520825.0/4194304.0/pow(a,23.0);
      param[3] = 239028075.0/4194304.0/pow(a,21.0);
      param[4] = -1320944625.0/8388608.0/pow(a,19.0);
      param[5] = 2657429775.0/8388608.0/pow(a,17.0);
      param[6] = -1003917915.0/2097152.0/pow(a,15.0);
      param[7] = 1158366825.0/2097152.0/pow(a,13.0);
      param[8] = -4106936925.0/8388608.0/pow(a,11.0);
      param[9] = 2788660875.0/8388608.0/(a*a*a*a*a*a*a*a*a);
      param[10] = -717084225.0/4194304.0/(a*a*a*a*a*a*a);
      param[11] = 273795795.0/4194304.0/(a*a*a*a*a);
      param[12] = -152108775.0/8388608.0/(a*a*a);
      param[13] = 35102025.0/8388608.0/a;
      break;
    case 15:
      param[0] = 5014575.0/33554432.0/pow(a,29.0);
      param[1] = -37702175.0/16777216.0/pow(a,27.0);
      param[2] = 529338537.0/33554432.0/pow(a,25.0);
      param[3] = -575367975.0/8388608.0/pow(a,23.0);
      param[4] = 6931814175.0/33554432.0/pow(a,21.0);
      param[5] = -7661478825.0/16777216.0/pow(a,19.0);
      param[6] = 25688487825.0/33554432.0/pow(a,17.0);
      param[7] = -4159088505.0/4194304.0/pow(a,15.0);
      param[8] = 33592637925.0/33554432.0/pow(a,13.0);
      param[9] = -13233463425.0/16777216.0/pow(a,11.0);
      param[10] = 16174233075.0/33554432.0/(a*a*a*a*a*a*a*a*a);
      param[11] = -1890494775.0/8388608.0/(a*a*a*a*a*a*a);
      param[12] = 2646692685.0/33554432.0/(a*a*a*a*a);
      param[13] = -339319575.0/16777216.0/(a*a*a);
      param[14] = 145422675.0/33554432.0/a;
      break;
    case 16:
      param[0] = -9694845.0/67108864.0/pow(a,31.0);
      param[1] = 155451825.0/67108864.0/pow(a,29.0);
      param[2] = -1168767425.0/67108864.0/pow(a,27.0);
      param[3] = 5469831549.0/67108864.0/pow(a,25.0);
      param[4] = -17836407225.0/67108864.0/pow(a,23.0);
      param[5] = 42977247885.0/67108864.0/pow(a,21.0);
      param[6] = -79168614525.0/67108864.0/pow(a,19.0);
      param[7] = 113763303225.0/67108864.0/pow(a,17.0);
      param[8] = -128931743655.0/67108864.0/pow(a,15.0);
      param[9] = 115707975075.0/67108864.0/pow(a,13.0);
      param[10] = -82047473235.0/67108864.0/pow(a,11.0);
      param[11] = 45581929575.0/67108864.0/(a*a*a*a*a*a*a*a*a);
      param[12] = -19535112675.0/67108864.0/(a*a*a*a*a*a*a);
      param[13] = 6311344095.0/67108864.0/(a*a*a*a*a);
      param[14] = -1502700975.0/67108864.0/(a*a*a);
      param[15] = 300540195.0/67108864.0/a;
      break;
    default:
      return 1;
  }

  return 0;
}


fcs_int ifcs_p2nfft_load_taylor2p_derive_coefficients(
   fcs_float epsI, fcs_int p,
   fcs_float *param
   )
{
  fcs_float a = epsI;

  /* Case p<2 is senseless, for p>16 we need to add the coefficients. */
  /* These coefficients are calculated analytically with a Maple script. */ 
  switch(p){
    case 2:
      param[0] = -1/(a*a*a);
      break;
    case 3:
      param[0] = 3.0/2.0/(a*a*a*a*a);
      param[1] = -5.0/2.0/(a*a*a);
      break;
    case 4:
      param[0] = -15.0/8.0/(a*a*a*a*a*a*a);
      param[1] = 21.0/4.0/(a*a*a*a*a);
      param[2] = -35.0/8.0/(a*a*a);
      break;
    case 5:
      param[0] = 35.0/16.0/(a*a*a*a*a*a*a*a*a);
      param[1] = -135.0/16.0/(a*a*a*a*a*a*a);
      param[2] = 189.0/16.0/(a*a*a*a*a);
      param[3] = -105.0/16.0/(a*a*a);
      break;
    case 6:
      param[0] = -315.0/128.0/pow(a,11.0);
      param[1] = 385.0/32.0/(a*a*a*a*a*a*a*a*a);
      param[2] = -1485.0/64.0/(a*a*a*a*a*a*a);
      param[3] = 693.0/32.0/(a*a*a*a*a);
      param[4] = -1155.0/128.0/(a*a*a);
      break;
    case 7:
      param[0] = 693.0/256.0/pow(a,13.0);
      param[1] = -4095.0/256.0/pow(a,11.0);
      param[2] = 5005.0/128.0/(a*a*a*a*a*a*a*a*a);
      param[3] = -6435.0/128.0/(a*a*a*a*a*a*a);
      param[4] = 9009.0/256.0/(a*a*a*a*a);
      param[5] = -3003.0/256.0/(a*a*a);
      break;
    case 8:
      param[0] = -3003.0/1024.0/pow(a,15.0);
      param[1] = 10395.0/512.0/pow(a,13.0);
      param[2] = -61425.0/1024.0/pow(a,11.0);
      param[3] = 25025.0/256.0/(a*a*a*a*a*a*a*a*a);
      param[4] = -96525.0/1024.0/(a*a*a*a*a*a*a);
      param[5] = 27027.0/512.0/(a*a*a*a*a);
      param[6] = -15015.0/1024.0/(a*a*a);
      break;
    case 9:
      param[0] = 6435.0/2048.0/pow(a,17.0);
      param[1] = -51051.0/2048.0/pow(a,15.0);
      param[2] = 176715.0/2048.0/pow(a,13.0);
      param[3] = -348075.0/2048.0/pow(a,11.0);
      param[4] = 425425.0/2048.0/(a*a*a*a*a*a*a*a*a);
      param[5] = -328185.0/2048.0/(a*a*a*a*a*a*a);
      param[6] = 153153.0/2048.0/(a*a*a*a*a);
      param[7] = -36465.0/2048.0/(a*a*a);
      break;
    case 10:
      param[0] = -109395.0/32768.0/pow(a,19.0);
      param[1] = 122265.0/4096.0/pow(a,17.0);
      param[2] = -969969.0/8192.0/pow(a,15.0);
      param[3] = 1119195.0/4096.0/pow(a,13.0);
      param[4] = -6613425.0/16384.0/pow(a,11.0);
      param[5] = 1616615.0/4096.0/(a*a*a*a*a*a*a*a*a);
      param[6] = -2078505.0/8192.0/(a*a*a*a*a*a*a);
      param[7] = 415701.0/4096.0/(a*a*a*a*a);
      param[8] = -692835.0/32768.0/(a*a*a);
      break;
    case 11:
      param[0] = 230945.0/65536.0/pow(a,21.0);
      param[1] = -2297295.0/65536.0/pow(a,19.0);
      param[2] = 2567565.0/16384.0/pow(a,17.0);
      param[3] = -6789783.0/16384.0/pow(a,15.0);
      param[4] = 23503095.0/32768.0/pow(a,13.0);
      param[5] = -27776385.0/32768.0/pow(a,11.0);
      param[6] = 11316305.0/16384.0/(a*a*a*a*a*a*a*a*a);
      param[7] = -6235515.0/16384.0/(a*a*a*a*a*a*a);
      param[8] = 8729721.0/65536.0/(a*a*a*a*a);
      param[9] = -1616615.0/65536.0/(a*a*a);
      break;
    case 12:
      param[0] = -969969.0/262144.0/pow(a,23.0);
      param[1] = 5311735.0/131072.0/pow(a,21.0);
      param[2] = -52837785.0/262144.0/pow(a,19.0);
      param[3] = 19684665.0/32768.0/pow(a,17.0);
      param[4] = -156165009.0/131072.0/pow(a,15.0);
      param[5] = 108114237.0/65536.0/pow(a,13.0);
      param[6] = -212952285.0/131072.0/pow(a,11.0);
      param[7] = 37182145.0/32768.0/(a*a*a*a*a*a*a*a*a);
      param[8] = -143416845.0/262144.0/(a*a*a*a*a*a*a);
      param[9] = 22309287.0/131072.0/(a*a*a*a*a);
      param[10] = -7436429.0/262144.0/(a*a*a);
      break;
    case 13:
      param[0] = 2028117.0/524288.0/pow(a,25.0);
      param[1] = -24249225.0/524288.0/pow(a,23.0);
      param[2] = 132793375.0/524288.0/pow(a,21.0);
      param[3] = -440314875.0/524288.0/pow(a,19.0);
      param[4] = 492116625.0/262144.0/pow(a,17.0);
      param[5] = -780825045.0/262144.0/pow(a,15.0);
      param[6] = 900951975.0/262144.0/pow(a,13.0);
      param[7] = -760543875.0/262144.0/pow(a,11.0);
      param[8] = 929553625.0/524288.0/(a*a*a*a*a*a*a*a*a);
      param[9] = -398380125.0/524288.0/(a*a*a*a*a*a*a);
      param[10] = 111546435.0/524288.0/(a*a*a*a*a);
      param[11] = -16900975.0/524288.0/(a*a*a);
      break;
    case 14:
      param[0] = -16900975.0/4194304.0/pow(a,27.0);
      param[1] = 54759159.0/1048576.0/pow(a,25.0);
      param[2] = -654729075.0/2097152.0/pow(a,23.0);
      param[3] = 1195140375.0/1048576.0/pow(a,21.0);
      param[4] = -11888501625.0/4194304.0/pow(a,19.0);
      param[5] = 2657429775.0/524288.0/pow(a,17.0);
      param[6] = -7027425405.0/1048576.0/pow(a,15.0);
      param[7] = 3475100475.0/524288.0/pow(a,13.0);
      param[8] = -20534684625.0/4194304.0/pow(a,11.0);
      param[9] = 2788660875.0/1048576.0/(a*a*a*a*a*a*a*a*a);
      param[10] = -2151252675.0/2097152.0/(a*a*a*a*a*a*a);
      param[11] = 273795795.0/1048576.0/(a*a*a*a*a);
      param[12] = -152108775.0/4194304.0/(a*a*a);
      break;
    case 15:
      param[0] = 35102025.0/8388608.0/pow(a,29.0);
      param[1] = -490128275.0/8388608.0/pow(a,27.0);
      param[2] = 1588015611.0/4194304.0/pow(a,25.0);
      param[3] = -6329047725.0/4194304.0/pow(a,23.0);
      param[4] = 34659070875.0/8388608.0/pow(a,21.0);
      param[5] = -68953309425.0/8388608.0/pow(a,19.0);
      param[6] = 25688487825.0/2097152.0/pow(a,17.0);
      param[7] = -29113619535.0/2097152.0/pow(a,15.0);
      param[8] = 100777913775.0/8388608.0/pow(a,13.0);
      param[9] = -66167317125.0/8388608.0/pow(a,11.0);
      param[10] = 16174233075.0/4194304.0/(a*a*a*a*a*a*a*a*a);
      param[11] = -5671484325.0/4194304.0/(a*a*a*a*a*a*a);
      param[12] = 2646692685.0/8388608.0/(a*a*a*a*a);
      param[13] = -339319575.0/8388608.0/(a*a*a);
      break;
    case 16:
      param[0] = -145422675.0/33554432.0/pow(a,31.0);
      param[1] = 1088162775.0/16777216.0/pow(a,29.0);
      param[2] = -15193976525.0/33554432.0/pow(a,27.0);
      param[3] = 16409494647.0/8388608.0/pow(a,25.0);
      param[4] = -196200479475.0/33554432.0/pow(a,23.0);
      param[5] = 214886239425.0/16777216.0/pow(a,21.0);
      param[6] = -712517530725.0/33554432.0/pow(a,19.0);
      param[7] = 113763303225.0/4194304.0/pow(a,17.0);
      param[8] = -902522205585.0/33554432.0/pow(a,15.0);
      param[9] = 347123925225.0/16777216.0/pow(a,13.0);
      param[10] = -410237366175.0/33554432.0/pow(a,11.0);
      param[11] = 45581929575.0/8388608.0/(a*a*a*a*a*a*a*a*a);
      param[12] = -58605338025.0/33554432.0/(a*a*a*a*a*a*a);
      param[13] = 6311344095.0/16777216.0/(a*a*a*a*a);
      param[14] = -1502700975.0/33554432.0/(a*a*a);
      break;
    default: 
      return 1;
  }

  return 0;
}



/** linear spline interpolation in near field with even kernels */
static fcs_float intpol_even_const(
    const fcs_float x, const fcs_float *table,
    const fcs_int num_nodes, const fcs_float one_over_eps_I
    )
{
  fcs_float c;
  fcs_int r;
  fcs_float f1;

  c=fcs_fabs(x*num_nodes*one_over_eps_I);
  r=(fcs_int)c;
  f1=table[r];
  return f1;
}

/** linear spline interpolation in near field with even kernels */
static fcs_float intpol_even_lin(
    const fcs_float x, const fcs_float *table,
    const fcs_int num_nodes, const fcs_float one_over_eps_I
    )
{
  fcs_float c,c1,c3;
  fcs_int r;
  fcs_float f1,f2;

  c=fcs_fabs(x*num_nodes*one_over_eps_I);
  r=(fcs_int)c;
  f1=table[r];f2=table[r+1];
  c1=c-r;
  c3=c1-1.0;
  return (-f1*c3+f2*c1);
}

/** quadratic spline interpolation in near field with even kernels */
static fcs_float intpol_even_quad(
    const fcs_float x, const fcs_float *table,
    const fcs_int num_nodes, const fcs_float one_over_eps_I
    )
{
  fcs_float c,c1,c2,c3;
  fcs_int r;
  fcs_float f0,f1,f2;

  c=fcs_fabs(x*num_nodes*one_over_eps_I);
  r=(fcs_int)c;
  if (r==0) {f0=table[r+1];f1=table[r];f2=table[r+1];}
  else { f0=table[r-1];f1=table[r];f2=table[r+1];}
  c1=c-r;
  c2=c1+1.0;
  c3=c1-1.0;
  return (0.5*f0*c1*c3-f1*c2*c3+0.5*f2*c2*c1);
}

/** cubic spline interpolation in near field with even kernels */
static fcs_float intpol_even_cub(
    const fcs_float x, const fcs_float *table,
    const fcs_int num_nodes, const fcs_float one_over_eps_I
    )
{
  fcs_float c,c1,c2,c3,c4;
  fcs_int r;
  fcs_float f0,f1,f2,f3;

  c=fcs_fabs(x*num_nodes*one_over_eps_I);
  r=(fcs_int)c;
  if (r==0) {f0=table[r+1];f1=table[r];f2=table[r+1];f3=table[r+2];}
  else { f0=table[r-1];f1=table[r];f2=table[r+1];f3=table[r+2];}
  c1=c-r;
  c2=c1+1.0;
  c3=c1-1.0;
  c4=c1-2.0;
  return(-f0*c1*c3*c4+3.0*f1*c2*c3*c4-3.0*f2*c2*c1*c4+f3*c2*c1*c3)/6.0;
}



