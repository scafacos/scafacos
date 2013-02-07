/*
 *  ghosts.c
 *
 *  Copyright 2006 Matthias Bolten. All rights reserved.
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <stdio.h>

#include <mpi.h>

#include "ghosts.h"

void update_ghosts( double ***v, mg_data* data, int level )
{
  /* MPI variables */
  int myid;
  MPI_Status stats[2];
  MPI_Request reqs[2];

  /* Loop variables */
  int i, j, k;
  int m = data[level].m_l, n = data[level].n_l, o = data[level].o_l;

  /* Getting local rank */
  MPI_Comm_rank(data[level].cart_comm, &myid);

  /* Sending data to left */
  for (j=0;j<n;j++) {
    for (k=0;k<o;k++) {
      data[level].sbufyz[j][k] = v[1][j][k];
    }
  }
  reqs[0] = MPI_REQUEST_NULL; reqs[1] = MPI_REQUEST_NULL;
  if (myid!=data[level].left)
    MPI_Isend((void *) data[level].sbufyz[0],n*o,MPI_DOUBLE,data[level].left,0,data[level].cart_comm,&reqs[0]);
  if (myid!=data[level].right)
    MPI_Irecv((void *) data[level].rbufyz[0],n*o,MPI_DOUBLE,data[level].right,0,data[level].cart_comm,&reqs[1]);
  MPI_Waitall(2,reqs,stats);
  if( MPI_PROC_NULL != data[level].right )
    for (j=0;j<n;j++) {
      for (k=0;k<o;k++) {
        v[m-1][j][k] = data[level].rbufyz[j][k];
      }
  }

  /* Sending data to right */
  for (j=0;j<n;j++) {
    for (k=0;k<o;k++) {
      data[level].sbufyz[j][k] = v[m-2][j][k];
    }
  }
  reqs[0] = MPI_REQUEST_NULL; reqs[1] = MPI_REQUEST_NULL;
  if (myid!=data[level].right)
    MPI_Isend((void *) data[level].sbufyz[0],n*o,MPI_DOUBLE,data[level].right,0,data[level].cart_comm,&reqs[0]);
  if (myid!=data[level].left)
    MPI_Irecv((void *) data[level].rbufyz[0],n*o,MPI_DOUBLE,data[level].left,0,data[level].cart_comm,&reqs[1]);
  MPI_Waitall(2,reqs,stats);
  if( MPI_PROC_NULL != data[level].left )
    for (j=0;j<n;j++) {
      for (k=0;k<o;k++) {
        v[0][j][k] = data[level].rbufyz[j][k];
      }
  }

  /* Sending data downwards */
  for (i=0;i<m;i++) {
    for (k=0;k<o;k++) {
      data[level].sbufxz[i][k] = v[i][1][k];
    }
  }
  reqs[0] = MPI_REQUEST_NULL; reqs[1] = MPI_REQUEST_NULL;
  if (myid!=data[level].lower)
    MPI_Isend((void *) data[level].sbufxz[0],m*o,MPI_DOUBLE,data[level].lower,0,data[level].cart_comm,&reqs[0]);
  if (myid!=data[level].upper)
    MPI_Irecv((void *) data[level].rbufxz[0],m*o,MPI_DOUBLE,data[level].upper,0,data[level].cart_comm,&reqs[1]);
  MPI_Waitall(2,reqs,stats);
  if( MPI_PROC_NULL != data[level].upper )
    for (i=0;i<m;i++) {
      for (k=0;k<o;k++) {
        v[i][n-1][k] = data[level].rbufxz[i][k];
      }
  }
  
  /* Sending data upwards */
  for (i=0;i<m;i++) {
    for (k=0;k<o;k++) {
      data[level].sbufxz[i][k] = v[i][n-2][k];
    }
  }
  reqs[0] = MPI_REQUEST_NULL; reqs[1] = MPI_REQUEST_NULL;
  if (myid!=data[level].upper)
    MPI_Isend((void *) data[level].sbufxz[0],m*o,MPI_DOUBLE,data[level].upper,0,data[level].cart_comm,&reqs[0]);
  if (myid!=data[level].lower)
    MPI_Irecv((void *) data[level].rbufxz[0],m*o,MPI_DOUBLE,data[level].lower,0,data[level].cart_comm,&reqs[1]);
  MPI_Waitall(2,reqs,stats);
  if( MPI_PROC_NULL != data[level].lower )
    for (i=0;i<m;i++) {
      for (k=0;k<o;k++) {
        v[i][0][k] = data[level].rbufxz[i][k];
      }
  }

  /* Sending data backwards */
  for (i=0;i<m;i++) {
    for (j=0;j<n;j++) {
      data[level].sbufxy[i][j] = v[i][j][1];
    }
  }
  reqs[0] = MPI_REQUEST_NULL; reqs[1] = MPI_REQUEST_NULL;
  if (myid!=data[level].back)
    MPI_Isend((void *) data[level].sbufxy[0],m*n,MPI_DOUBLE,data[level].back,0,data[level].cart_comm,&reqs[0]);
  if (myid!=data[level].front)
    MPI_Irecv((void *) data[level].rbufxy[0],m*n,MPI_DOUBLE,data[level].front,0,data[level].cart_comm,&reqs[1]);
  MPI_Waitall(2,reqs,stats);
  if( MPI_PROC_NULL != data[level].front )
    for (i=0;i<m;i++) {
      for (j=0;j<n;j++) {
        v[i][j][o-1] = data[level].rbufxy[i][j];
      }
  }
  
  /* Sending data to front */
  for (i=0;i<m;i++) {
    for (j=0;j<n;j++) {
      data[level].sbufxy[i][j] = v[i][j][o-2];
    }
  }
  reqs[0] = MPI_REQUEST_NULL; reqs[1] = MPI_REQUEST_NULL;
  if (myid!=data[level].front)
    MPI_Isend((void *) data[level].sbufxy[0],m*n,MPI_DOUBLE,data[level].front,0,data[level].cart_comm,&reqs[0]);
  if (myid!=data[level].back)
    MPI_Irecv((void *) data[level].rbufxy[0],m*n,MPI_DOUBLE,data[level].back,0,data[level].cart_comm,&reqs[1]);
  MPI_Waitall(2,reqs,stats);
  if( MPI_PROC_NULL != data[level].back )
    for (i=0;i<m;i++) {
      for (j=0;j<n;j++) {
        v[i][j][0] = data[level].rbufxy[i][j];
      }
  }
  
  if (data[level].periodic) {
    if (myid==data[level].back && myid==data[level].front) {
      for (i=0;i<m;i++) {
	for (j=0;j<n;j++) {
	  v[i][j][0  ] = v[i][j][o-2];
	  v[i][j][o-1] = v[i][j][1  ];
	}
      }
    }

    if (myid==data[level].lower && myid==data[level].upper) {
      for (i=0;i<m;i++) {
	for (k=0;k<o;k++) {
	  v[i][0  ][k] = v[i][n-2][k];
	  v[i][n-1][k] = v[i][1  ][k];
	}
      }
    }

    if (myid==data[level].left && myid==data[level].right) {
      for (j=0;j<n;j++) {
	for (k=0;k<o;k++) {
	  v[0  ][j][k] = v[m-2][j][k];
	  v[m-1][j][k] = v[1  ][j][k];
	}
      }
    }
  }
    
  return;
}
