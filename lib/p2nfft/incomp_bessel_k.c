

/* Compute the incomplete Bessel-K function of order nu according to the paper
   Richard M. Slevinsky and Hassan Safouhi. 2010.
   A recursive algorithm for the G transformation and accurate computation of incomplete Bessel functions.
   Appl. Numer. Math. 60, 12 (December 2010), 1411-1417.
   DOI=10.1016/j.apnum.2010.04.005 http://dx.doi.org/10.1016/j.apnum.2010.04.005 
*/

#include <math.h>
#include <stdio.h>



/* General recursion formula for coefficients of the G transform */
static fcs_int A(
    fcs_int i, fcs_int k
    )
{
  /* Set mu and m suitable to the Bessel-K approximation */
  fcs_int mu = -2;
  fcs_int m = 0;

  if(k<i)
    fprintf(stderr, "Error in computation of G transform coefficients A(i,k): k < i is not valid\n");

  if(i==k)
    return 1;

  if(i==0)
    return (n - nu - (k-1)*(mu+1)) * A(0,k-1);
  else
    return (n - nu + i*(m+1)-(k-1)*(mu+1)) * A(i,k-1) + A(i-1,k-1);
}

static fcs_float D(
    fcs_int n, fcs_float x, fcs_float y, fcs_int nu
    )
{




}


static fcs_float D_tilde(
    fcs_int n, fcs_float x, fcs_float y, fcs_int nu
    )
{
  fcs_float sum1, sum2;

  sum1 = 0.0;
  for(fcs_int r=0; r<=n; r++){
    sum2 = 0.0;
    for(fcs_int i=0; i<=r; i++)
      sum2 += A(i,r) * pow(x,i);
    sum1 += binom(n,r) * pow(-y,-r) * sum2;
  }

  return pow(-x*y,n) * pow(x,nu+1) * exp(x+y) * sum1;
}

static fcs_float N_tilde(
    fcs_int n, fcs_float x, fcs_float y, fcs_int nu
    )
{
  fcs_float sum1, sum2, sum3;

  sum1 = 0.0;
  for(fcs_int r=1; r<=n; r++){
    sum2 = 0.0;
    for(fcs_int s=0; s<=r-1; s++){
      sum3 = 0.0;
      for(fcs_int i=0; i<=s; i++)
        sum3 += A(i,s)*pow(-x,i);
      sum2 += binom(r-1,s) * pow(y,-s) * sum3;
    }
    sum1 += binom(n,r) * D(n-r,x,y,nu) * pow(x*y,r) * sum2;
  }

  return exp(-x-y)/(pow(x,nu)*y) * sum1;
}

static fcs_float G_tilde(
    fcs_int n, fcs_float x, fcs_float y, fcs_int nu
    )
{
  return pow(x,nu) * N_tilde(n,x,y,nu) / D_tilde(n,x,y,nu);
} 

fcs_float inc_bessel_k(
    fcs_int nu, fcs_float x, fcs_float y
    )
{
  fcs_int n = 10;
  return G_tilde(n,x,y,nu);
}



