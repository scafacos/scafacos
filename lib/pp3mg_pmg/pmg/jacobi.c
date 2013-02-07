/*
 *  jacobi.c
 *
 *  Copyright 2006 Matthias Bolten. All rights reserved.
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <mpi.h>

#include "jacobi.h"
#include "lueqf.h"
#include "ghosts.h"

double jacobi( double*** v, double*** f, double*** r, mg_data* data, int level, 
	       int maxiter )
{

  double res;
  int iter;
  int i, j, k;

  for (iter=0;iter<maxiter;iter++) {
    /* r = f - A*v */
    update_ghosts( v, data, level );
    lueqf_res( v, f, r, data, level );
    /* r = D^(-1)*r */
    lueqf_invd( r, data, level );
    /* v = v + omega * r */
    for (i=0;i<data[level].m_l;i++) {
      for (j=0;j<data[level].n_l;j++) {
	for (k=0;k<data[level].o_l;k++) {
	  v[i][j][k] = v[i][j][k] + data[level].omega * r[i][j][k];
	}
      }
    }
  }

  update_ghosts( v, data, level );
  res = lueqf_res( v, f, r, data, level );
  update_ghosts( v, data, level );

  return(res);
}
