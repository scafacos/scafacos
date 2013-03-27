/*
 * interpolation.c
 *
 * This file contains the function "interp_poly" to interpolate from
 * gridded potential using polynomials.
 * 
 * Based on Fortran code by Matthias Bolten.
 *
 * Authors: Matthias Bolten, Stephanie Friedhoff
 * Created: 2009/06/02
 *
 * Copyright 2009, 2010, 2011, 2012 Matthias Bolten, Stephanie Friedhoff
 * All rights reserved.
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

#include "pp3mg.h"

#include "interpolation.h"

void interp_poly( int degree, pp3mg_particle* particles, int*** hoc, int* ll, double*** u, 
		  int n_particles, int m_start, int m_end, int n_start, int n_end, 
		  int o_start, int o_end, int boundary, double hx, double hy, double hz )
{

  /* Local variables */
  /* Loop variables */
  int i, j, k, ii, jj, kk, iii, jjj, kkk, iiii, jjjj, kkkk;

  /* Products */
  double prodiiii, prodjjjj, prodkkkk;

  /* Position in particle structure */
  int p;

  /* Coordinates */
  double x, y, z;

  /* Polynomials' coefficients */
  double* xi, *yj, *zk;

  /* Polynomials' values */
  double*** fijk;

  /* Lambdas */
  double* lambdax, *lambday, *lambdaz;

  /* Sums */
  double sumii, sumjj, sumkk;

  /* Value */
  double val;

/*   printf("m_start = %d, m_end = %d, n_start = %d, n_end = %d, o_start = %d, o_end = %d\n",m_start,m_end,n_start,n_end,o_start,o_end); */
	
  /* Checking sizes */
  if( boundary < ceil( 0.5*degree) )
    printf( "Error: Field not large enough for interpolation.\n" );

  /* Allocating memory for polynomial coefficients and values */
  xi = (double*) malloc( (degree+1) * sizeof( double ) );
  yj = (double*) malloc( (degree+1) * sizeof( double ) );
  zk = (double*) malloc( (degree+1) * sizeof( double ) );

  fijk = (double***) malloc( (degree+1) * sizeof( double** ) );
  for( i = 0; i <= degree; i++ ){
    fijk[i] = (double**) malloc( (degree+1) * sizeof( double* ) );
    for( j = 0; j <= degree; j++ )
      fijk[i][j] = (double*) malloc( (degree+1)*sizeof(double) );
  }

  lambdax = (double*) malloc( (degree+1) * sizeof( double ) );
  lambday = (double*) malloc( (degree+1) * sizeof( double ) );
  lambdaz = (double*) malloc( (degree+1) * sizeof( double ) );
	
  /* Interpolating data from grid */
  for( i = boundary; i <= m_end-m_start+boundary; i++ ){
    for( j = boundary; j <= n_end-n_start+boundary; j++ ){
      for( k = boundary; k <= o_end-o_start+boundary; k++){
	p = hoc[i][j][k];
	while( p >= 0 ){
	  /* Polynomial interpolation */
	  x = particles[p].x;
	  y = particles[p].y;
	  z = particles[p].z;

	  if( (degree+1-2*((degree+1)/2)) == 0 )
	    for( ii = 0; ii <= degree; ii++ )
	      xi[ii] = (floor(x/hx)-(degree+1)/2+ii+1) * hx;
	  else
	    for( ii = 0; ii <= degree; ii++ )
	      xi[ii] = ((int) round(x/hx)-(degree+1)/2+ii) * hx;

	  if( (degree+1-2*((degree+1)/2)) == 0 )                
	    for( jj = 0; jj <= degree; jj++ )
	      yj[jj] = (floor(y/hy)-(degree+1)/2+jj+1) * hy;
	  else
	    for( jj = 0; jj <= degree; jj++ )
	      yj[jj] = ((int) round(y/hy)-(degree+1)/2+jj) * hy;

	  if( (degree+1-2*((degree+1)/2)) == 0 )
	    for( kk = 0; kk <= degree; kk++ )
	      zk[kk] = (floor(z/hz)-(degree+1)/2+kk+1) * hz;
	  else
	    for( kk = 0; kk <= degree; kk++ )
	      zk[kk] = ((int) round(z/hz)-(degree+1)/2+kk) * hz;
					
	  for( ii = 0; ii <= degree; ii++ )
	    for( jj = 0; jj <= degree; jj++ )
	      for( kk = 0; kk <= degree; kk++ ) {
		/* printf("ii = %d, jj = %d, kk = %d, degree = %d\n",ii,jj,kk,degree); */
		/* printf("xi[ii] = %f, yj[jj] = %f, zk[kk] = %f\n",xi[ii],yj[jj],zk[kk]); */
		/* printf("x = %f, y = %f, z = %f, i = %d, j = %d, k = %d, p = %d\n",x,y,z,i,j,k,p); */
		/* printf("u[%d][%d][%d]\n",(int) round(xi[ii]/hx)-m_start+boundary,(int) round(yj[jj]/hy)-n_start+boundary,(int) round(zk[kk]/hz)-o_start+boundary); */
		fijk[ii][jj][kk] = u[(int) round(xi[ii]/hx)-m_start+boundary][(int) round(yj[jj]/hy)-n_start+boundary][(int) round(zk[kk]/hz)-o_start+boundary];
	      }
						
	  /* Compute lambdas (see Schwarz, "Numerische Mathematik", Teubner) */
	  /* For x */
	  lambdax[0] = 1.0;
	  for( iii = 1; iii <= degree; iii++ ){
	    lambdax[iii] = 0.0;
	    for( ii = 0; ii <= (iii-1); ii++ ){
	      lambdax[ii]  = lambdax[ii]/(xi[ii] - xi[iii]);
	      lambdax[iii] = lambdax[iii] - lambdax[ii];
	    }
	  }
										
	  /* For y */
	  lambday[0] = 1.0;
	  for( jjj = 1; jjj <= degree; jjj++ ){
	    lambday[jjj] = 0.0;
	    for( jj = 0; jj <= (jjj-1); jj++ ){
	      lambday[jj]  = lambday[jj]/(yj[jj] - yj[jjj]);
	      lambday[jjj] = lambday[jjj] - lambday[jj];
	    }
	  }	
									
	  /* For z */
	  lambdaz[0] = 1.0;
	  for( kkk = 1; kkk <= degree; kkk++ ){
	    lambdaz[kkk] = 0.0;
	    for( kk = 0; kk <= (kkk-1); kk++ ){
	      lambdaz[kk]  = lambdaz[kk]/(zk[kk] - zk[kkk]);
	      lambdaz[kkk] = lambdaz[kkk] - lambdaz[kk];
	    }
	  }
					
	  /* Compute interpolated value */
	  sumii = 0.0;
	  for( ii = 0; ii <= degree; ii++ ){
	    sumjj = 0.0;
	    for( jj = 0; jj <= degree; jj++ ){
	      sumkk = 0.0;
	      for( kk = 0; kk <= degree; kk++ ){
		val = 1.0;
		for( iii = 0; iii <= (ii-1); iii++ )
		  val = val * ( x - xi[iii] );
		for( iii = ii+1; iii <= degree; iii++ )
		  val = val * ( x - xi[iii] );
		for( jjj = 0; jjj <= (jj-1); jjj++ )
		  val = val * ( y - yj[jjj] );
		for( jjj = jj+1; jjj <= degree; jjj++ )
		  val = val * ( y - yj[jjj] );
		for( kkk = 0; kkk <= (kk-1); kkk++ )
		  val = val * ( z - zk[kkk] );
		for( kkk = kk+1; kkk <= degree; kkk++ )
		  val = val * ( z - zk[kkk] );
		sumkk += val * lambdaz[kk] * fijk[ii][jj][kk];									
	      }
	      sumjj += lambday[jj] * sumkk;
	    }
	    sumii += lambdax[ii] * sumjj;
	  }
	  particles[p].e += particles[p].q * sumii;
/* 	  printf("particles[%d].e = %f\n",p,particles[p].e); */
					
	  /* Compute derivative in x-direction */
	  sumii = 0.0;
	  for( ii = 0; ii <= degree; ii++ ){
	    sumjj = 0.0;
	    for( jj = 0; jj <= degree; jj++ ){
	      sumkk = 0.0;
	      for( kk = 0; kk <= degree; kk++ ){
		val = 0.0;
		for( iii = 0; iii <= (ii-1); iii++ ){
		  prodiiii = 1.0;
		  for( iiii = 0; iiii <= degree; iiii++ )
		    if( (iiii != ii) && (iiii != iii) )
		      prodiiii *= ( x - xi[iiii] );
		  val += prodiiii;
		}
		for( iii = ii+1; iii <= degree; iii++ ){
		  prodiiii = 1.0;
		  for( iiii = 0; iiii <= degree; iiii++ )
		    if( (iiii != ii) && (iiii != iii) )
		      prodiiii *= ( x - xi[iiii] );
		  val += prodiiii;
		}
								
		for( jjj = 0; jjj <= (jj-1); jjj++ )
		  val *= ( y - yj[jjj] );
		for( jjj = jj+1; jjj <= degree; jjj++ )
		  val *= ( y - yj[jjj] );

		for( kkk = 0; kkk <= (kk-1); kkk++ )
		  val *= ( z - zk[kkk] );
		for( kkk = kk+1; kkk <= degree; kkk++ )
		  val *= ( z - zk[kkk] );	
									
		sumkk += val * lambdaz[kk] * fijk[ii][jj][kk];						
	      }
	      sumjj += lambday[jj] * sumkk;
	    }
	    sumii += lambdax[ii] * sumjj;
	  }									

	  particles[p].fx -= particles[p].q * sumii;
					
	  /* Compute derivative in y-direction */
	  sumii = 0.0;
	  for( ii = 0; ii <= degree; ii++ ){
	    sumjj = 0.0;
	    for( jj = 0; jj <= degree; jj++ ){
	      sumkk = 0.0;
	      for( kk = 0; kk <= degree; kk++ ){
		val = 0.0;
		for( jjj = 0; jjj <= (jj-1); jjj++ ){
		  prodjjjj = 1.0;
		  for( jjjj = 0; jjjj <= degree; jjjj++ )
		    if( (jjjj != jj) && (jjjj != jjj) )
		      prodjjjj *= ( y - yj[jjjj] );
		  val += prodjjjj;
		}
		for( jjj = jj+1; jjj <= degree; jjj++ ){
		  prodjjjj = 1.0;
		  for( jjjj = 0; jjjj <= degree; jjjj++ )
		    if( (jjjj != jj) && (jjjj != jjj) )
		      prodjjjj *= ( y - yj[jjjj] );
		  val += prodjjjj;
		}
								
		for( iii = 0; iii <= (ii-1); iii++ )
		  val *= ( x - xi[iii] );
		for( iii = ii+1; iii <= degree; iii++ )
		  val *= ( x - xi[iii] );

		for( kkk = 0; kkk <= (kk-1); kkk++ )
		  val *= ( z - zk[kkk] );
		for( kkk = kk+1; kkk <= degree; kkk++ )
		  val *= ( z - zk[kkk] );
								
		sumkk += val * lambdaz[kk] * fijk[ii][jj][kk];
	      }
	      sumjj += lambday[jj] * sumkk;
	    }
	    sumii += lambdax[ii] * sumjj;
	  }

	  particles[p].fy -= particles[p].q * sumii;

	  /* Compute derivative in z-direction */
	  sumii = 0.0;
	  for( ii = 0; ii <= degree; ii++ ){
	    sumjj = 0.0;
	    for( jj = 0; jj <= degree; jj++ ){
	      sumkk = 0.0;
	      for( kk = 0; kk <= degree; kk++ ){
		val = 0.0;	
		for( kkk = 0; kkk <= (kk-1); kkk++ ){
		  prodkkkk = 1.0;
		  for( kkkk = 0; kkkk <= degree; kkkk++ )
		    if( (kkkk != kk) && (kkkk != kkk) )
		      prodkkkk *= ( z - zk[kkkk] );
		  val += prodkkkk;
		}
		for( kkk = kk+1; kkk <= degree; kkk++ ){
		  prodkkkk = 1.0;
		  for( kkkk = 0; kkkk <= degree; kkkk++ )
		    if( (kkkk != kk) && (kkkk != kkk) )
		      prodkkkk *= ( z - zk[kkkk] );
		  val += prodkkkk;
		}
								
		for( iii = 0; iii <= (ii-1); iii++ )
		  val *= ( x - xi[iii] );
		for( iii = ii+1; iii <= degree; iii++ )
		  val *= ( x - xi[iii] );

		for( jjj = 0; jjj <= (jj-1); jjj++ )
		  val *= ( y - yj[jjj] );
		for( jjj = jj+1; jjj <= degree; jjj++ )
		  val *= ( y - yj[jjj] );
								
		sumkk += val * lambdaz[kk] * fijk[ii][jj][kk];							
	      }
	      sumjj += lambday[jj] * sumkk;
	    }
	    sumii += lambdax[ii] * sumjj; 
	  }

	  particles[p].fz -= particles[p].q * sumii;
					
	  /* Getting next particle */
	  p = ll[p];
	}
      }
    }
  }

  for( i = 0; i <= degree; i++ ){
    for( j = 0; j <= degree; j++ )
      free(fijk[i][j]);
    free(fijk[i]);
  }
  free(fijk);
  free( xi );
  free( yj );
  free( zk );
  free( lambdax );
  free( lambday );
  free( lambdaz );	               
}
