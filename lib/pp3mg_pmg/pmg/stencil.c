#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include "stencil.h"
#include "cuboid.h"

inline int min(int a, int b) {
  return (a<b)?a:b;
}

inline int max(int a, int b) {
  return (a>b)?a:b;
}

/*
 *
 * conv3 - Computes the convolution of A and B and stores it in C.
 *
 */
void conv3(double ***A, double ***B, double ***C, int m_A, int n_A, int o_A, int m_B, int n_B, int o_B, int m_C, int n_C, int o_C) {
  int Aentriesx, Aentriesy, Aentriesz;
  int Bentriesx, Bentriesy, Bentriesz;
  int Centriesx, Centriesy, Centriesz;
  int i, j, k, ii, jj, kk;
  
  Aentriesx = (m_A-1)/2;
  Aentriesy = (n_A-1)/2;
  Aentriesz = (o_A-1)/2;
  Bentriesx = (m_B-1)/2;
  Bentriesy = (n_B-1)/2;
  Bentriesz = (o_B-1)/2;

  Centriesx = Aentriesx + Bentriesx;
  Centriesy = Aentriesy + Bentriesy;
  Centriesz = Aentriesz + Bentriesz;

  if (m_C != 2*Centriesx + 1) printf("Error! Size of array C does not match!");
  if (n_C != 2*Centriesy + 1) printf("Error! Size of array C does not match!");
  if (o_C != 2*Centriesz + 1) printf("Error! Size of array C does not match!");
  
  for (i=0;i<m_C;i++)
    for (j=0;j<n_C;j++)
      for (k=0;k<o_C;k++)
	C[i][j][k]=0.0;

  for (i=-Centriesx; i<=Centriesx; i++)
    for (j=-Centriesy; j<=Centriesy; j++)
      for (k=-Centriesz; k<=Centriesz; k++)
	for (ii=max(-Aentriesx,i-Bentriesx); ii<=min(Aentriesx,i+Bentriesx); ii++)
	  for (jj=max(-Aentriesy,j-Bentriesy); jj<=min(Aentriesy,j+Bentriesy); jj++)
	    for (kk=max(-Aentriesz,k-Bentriesz); kk<=min(Aentriesz,k+Bentriesz); kk++)
	      C[i+Centriesx][j+Centriesy][k+Centriesz] = 
		C[i+Centriesx][j+Centriesy][k+Centriesz] +
		B[ii-i+Bentriesx][jj-j+Bentriesy][kk-k+Bentriesz] * 
		A[ii+Aentriesx][jj+Aentriesy][kk+Aentriesz];

  return;
} 

/*
 *
 * set_stencil - calculates Galerkin stencil for next level
 * 
 */
void set_stencil( mg_data* data, int level, int dim )
{
  /* stencils */
  double ***A, ***P, ***Ac, ***Tmp, ***Ac_full;
  /* stencil sizes */
  int m_A, n_A, o_A;
  int A_entries_x, A_entries_y, A_entries_z;
  int m_P, n_P, o_P;
  int P_entries_x, P_entries_y, P_entries_z;
  int m_Ac, n_Ac, o_Ac;
  int Ac_entries_x, Ac_entries_y, Ac_entries_z;
  int m_Tmp, n_Tmp, o_Tmp;
  int Tmp_entries_x, Tmp_entries_y, Tmp_entries_z;
  int m_Ac_full, n_Ac_full, o_Ac_full;
  int Ac_full_entries_x, Ac_full_entries_y, Ac_full_entries_z;
  /* loop counter */
  int i, j, k;
  /* counter */
  int cnt;

  /* set omega (why do we do this here?) */
  data[level+1].omega = 0.7;

  /* create stencil for prolongation */
  /* WARNING: currently only zeros at (0.0,0.0,0.0) or (M_PI,M_PI,M_PI) are supported! */
  P_entries_x = P_entries_y = P_entries_z = 1;
  m_P = 2*P_entries_x + 1; n_P = 2*P_entries_y + 1; o_P = 2*P_entries_z + 1;
  P = cuboid_alloc(m_P,n_P,o_P);
  for (i=-1;i<=1;i++) {
    for (j=-1;j<=1;j++) {
      for (k=-1;k<=1;k++) {
	if (data[level].zero_at_pi3==false) {
	  P[i+1][j+1][k+1] = 0.015625 * (2 - abs(i)) * (2 - abs(j)) * (2 - abs(k));
	} else {
	  P[i+1][j+1][k+1] = 0.015625 * (2 - 3*abs(i)) * (2 - 3*abs(j)) * (2 - 3*abs(k));
	}
      }
    }
  }

  /* create stencil for fine grid operator */
  A_entries_x = A_entries_y = A_entries_z = 0;
  for (i=0;i<data[level].size;i++) {
    if (abs(data[level].x_offsets[i])>A_entries_x) A_entries_x = abs(data[level].x_offsets[i]);
    if (abs(data[level].y_offsets[i])>A_entries_y) A_entries_y = abs(data[level].y_offsets[i]);
    if (abs(data[level].z_offsets[i])>A_entries_z) A_entries_z = abs(data[level].z_offsets[i]);
  }
  m_A = 2*A_entries_x + 1; n_A = 2*A_entries_y + 1; o_A = 2*A_entries_z + 1;
  A = cuboid_alloc(m_A,n_A,o_A);
  for (i=0;i<m_A;i++) for (j=0;j<n_A;j++) for (k=0;k<o_A;k++) A[i][j][k] = 0.0;
  for (i=0;i<data[level].size;i++) {
    A[data[level].x_offsets[i]+A_entries_x][data[level].y_offsets[i]+A_entries_y][data[level].z_offsets[i]+A_entries_z] = data[level].values[i];
  }
  

  /* allocate memory for temporary stencil */
  Tmp_entries_x = A_entries_x + P_entries_x;
  Tmp_entries_y = A_entries_y + P_entries_y;
  Tmp_entries_z = A_entries_z + P_entries_z;
  m_Tmp = 2*Tmp_entries_x + 1;
  n_Tmp = 2*Tmp_entries_y + 1;
  o_Tmp = 2*Tmp_entries_z + 1;
  Tmp = cuboid_alloc(m_Tmp, n_Tmp, o_Tmp);

  /* calculate temporary stencil */
  conv3(A,P,Tmp,m_A,n_A,o_A,m_P,n_P,o_P,m_Tmp,n_Tmp,o_Tmp);

  /* allocate memory for full coarse grid stencil */
  Ac_full_entries_x = P_entries_x + Tmp_entries_x;
  Ac_full_entries_y = P_entries_y + Tmp_entries_y;
  Ac_full_entries_z = P_entries_z + Tmp_entries_z;
  m_Ac_full = 2*Ac_full_entries_x + 1;
  n_Ac_full = 2*Ac_full_entries_y + 1;
  o_Ac_full = 2*Ac_full_entries_z + 1;
  Ac_full = cuboid_alloc(m_Ac_full,n_Ac_full,o_Ac_full);

  /* calculate full coarse grid stencil */
  conv3(P,Tmp,Ac_full,m_P,n_P,o_P,m_Tmp,n_Tmp,o_Tmp,m_Ac_full,n_Ac_full,o_Ac_full);

  /* allocate memory for cut coarse grid stencil */
  Ac_entries_x = Ac_full_entries_x/2;
  Ac_entries_y = Ac_full_entries_y/2;
  Ac_entries_z = Ac_full_entries_z/2;
  m_Ac = 2*Ac_entries_x + 1;
  n_Ac = 2*Ac_entries_y + 1;
  o_Ac = 2*Ac_entries_z + 1;
  Ac = cuboid_alloc(m_Ac,n_Ac,o_Ac);

  /* cut coarse grid stencil */
  for (i=-Ac_entries_x;i<=Ac_entries_y;i++) {
    for (j=-Ac_entries_y;j<=Ac_entries_y;j++) {
      for (k=-Ac_entries_z;k<=Ac_entries_z;k++) {
	Ac[i+Ac_entries_x][j+Ac_entries_y][k+Ac_entries_z] =
	  Ac_full[2*i+Ac_full_entries_x][2*j+Ac_full_entries_y][2*k+Ac_full_entries_z];
      }
    }
  }

  /* /\* output coarse grid stencil *\/ */
  /* for (i=0;i<m_Ac;i++) { */
  /*   for (j=0;j<n_Ac;j++) { */
  /*     for (k=0;k<o_Ac;k++) */
  /* 	printf("%e\t",Ac[i][j][k]); */
  /*     printf("\n"); */
  /*   } */
  /*   printf("\n"); */
  /* } */

  /* Count non-zeros in stencil */
  data[level+1].size = 0;
  for (i=0;i<m_Ac;i++)
    for (j=0;j<n_Ac;j++)
      for (k=0;k<o_Ac;k++)
	if (Ac[i][j][k] != 0.0) data[level+1].size++;

  /* allocate space for stencil data */
  data[level+1].values    = (double*) realloc(data[level+1].values, data[level+1].size*sizeof(double) );
  data[level+1].x_offsets = (int*) realloc(data[level+1].x_offsets, data[level+1].size*sizeof(int) );
  data[level+1].y_offsets = (int*) realloc(data[level+1].y_offsets, data[level+1].size*sizeof(int) );
  data[level+1].z_offsets = (int*) realloc(data[level+1].z_offsets, data[level+1].size*sizeof(int) );
  
  /* save cut coarse grid stencil to data structure */
  cnt = 0;
  for (i=0;i<m_Ac;i++)
    for (j=0;j<n_Ac;j++)
      for (k=0;k<o_Ac;k++)
	if (Ac[i][j][k] != 0.0) {
	  data[level+1].values[cnt] = Ac[i][j][k];
	  data[level+1].x_offsets[cnt] = i - Ac_entries_x;
	  data[level+1].y_offsets[cnt] = j - Ac_entries_y;
	  data[level+1].z_offsets[cnt] = k - Ac_entries_z;
	  cnt++;
	}

  /* free all memory allocated for stencils */
  cuboid_free(P,m_P,n_P,o_P);
  cuboid_free(A,m_A,n_A,o_A);
  cuboid_free(Tmp,m_Tmp,n_Tmp,o_Tmp);
  cuboid_free(Ac,m_Ac,n_Ac,o_Ac);
  cuboid_free(Ac_full,m_Ac_full,n_Ac_full,o_Ac_full);

  return;
}


