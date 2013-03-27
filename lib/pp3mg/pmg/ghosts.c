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
  for (i=0;i<data[level].x_ghosts;i++)
    for (j=0;j<n;j++)
      for (k=0;k<o;k++)
	data[level].sbufyz[i][j][k] = v[data[level].x_ghosts+i][j][k];

  reqs[0] = MPI_REQUEST_NULL; reqs[1] = MPI_REQUEST_NULL;
  if (myid!=data[level].left)
    MPI_Isend((void *) data[level].sbufyz[0][0],data[level].x_ghosts*n*o,MPI_DOUBLE,data[level].left,0,data[level].cart_comm,&reqs[0]);
  if (myid!=data[level].right)
    MPI_Irecv((void *) data[level].rbufyz[0][0],data[level].x_ghosts*n*o,MPI_DOUBLE,data[level].right,0,data[level].cart_comm,&reqs[1]);
  MPI_Waitall(2,reqs,stats);
  if( MPI_PROC_NULL != data[level].right )
    for (i=0;i<data[level].x_ghosts;i++)
      for (j=0;j<n;j++)
	for (k=0;k<o;k++)
	  v[m-data[level].x_ghosts+i][j][k] = data[level].rbufyz[i][j][k];

  /* Sending data to right */
  for (i=0;i<data[level].x_ghosts;i++)
    for (j=0;j<n;j++)
      for (k=0;k<o;k++)
	data[level].sbufyz[i][j][k] = v[m-2*data[level].x_ghosts+i][j][k];
  
  reqs[0] = MPI_REQUEST_NULL; reqs[1] = MPI_REQUEST_NULL;
  if (myid!=data[level].right)
    MPI_Isend((void *) data[level].sbufyz[0][0],data[level].x_ghosts*n*o,MPI_DOUBLE,data[level].right,0,data[level].cart_comm,&reqs[0]);
  if (myid!=data[level].left)
    MPI_Irecv((void *) data[level].rbufyz[0][0],data[level].x_ghosts*n*o,MPI_DOUBLE,data[level].left,0,data[level].cart_comm,&reqs[1]);
  MPI_Waitall(2,reqs,stats);
  if( MPI_PROC_NULL != data[level].left )
    for (i=0;i<data[level].x_ghosts;i++)
      for (j=0;j<n;j++)
	for (k=0;k<o;k++)
	  v[i][j][k] = data[level].rbufyz[i][j][k];

  /* Sending data downwards */
  for (i=0;i<m;i++)
    for (j=0;j<data[level].y_ghosts;j++)
      for (k=0;k<o;k++)
	data[level].sbufxz[i][j][k] = v[i][data[level].y_ghosts+j][k];

  reqs[0] = MPI_REQUEST_NULL; reqs[1] = MPI_REQUEST_NULL;
  if (myid!=data[level].lower)
    MPI_Isend((void *) data[level].sbufxz[0][0],m*data[level].y_ghosts*o,MPI_DOUBLE,data[level].lower,0,data[level].cart_comm,&reqs[0]);
  if (myid!=data[level].upper)
    MPI_Irecv((void *) data[level].rbufxz[0][0],m*data[level].y_ghosts*o,MPI_DOUBLE,data[level].upper,0,data[level].cart_comm,&reqs[1]);
  MPI_Waitall(2,reqs,stats);
  if( MPI_PROC_NULL != data[level].upper )
    for (i=0;i<m;i++)
      for (j=0;j<data[level].y_ghosts;j++)
	for (k=0;k<o;k++)
	  v[i][n-data[level].y_ghosts+j][k] = data[level].rbufxz[i][j][k];
  
  /* Sending data upwards */
  for (i=0;i<m;i++)
    for (j=0;j<data[level].y_ghosts;j++)
      for (k=0;k<o;k++)
	data[level].sbufxz[i][j][k] = v[i][n-2*data[level].y_ghosts+j][k];

  reqs[0] = MPI_REQUEST_NULL; reqs[1] = MPI_REQUEST_NULL;
  if (myid!=data[level].upper)
    MPI_Isend((void *) data[level].sbufxz[0][0],m*data[level].y_ghosts*o,MPI_DOUBLE,data[level].upper,0,data[level].cart_comm,&reqs[0]);
  if (myid!=data[level].lower)
    MPI_Irecv((void *) data[level].rbufxz[0][0],m*data[level].y_ghosts*o,MPI_DOUBLE,data[level].lower,0,data[level].cart_comm,&reqs[1]);
  MPI_Waitall(2,reqs,stats);
  if( MPI_PROC_NULL != data[level].lower )
    for (i=0;i<m;i++)
      for (j=0;j<data[level].y_ghosts;j++)
	for (k=0;k<o;k++)
	  v[i][j][k] = data[level].rbufxz[i][j][k];

  /* Sending data backwards */
  for (i=0;i<m;i++)
    for (j=0;j<n;j++)
      for (k=0;k<data[level].z_ghosts;k++)
	data[level].sbufxy[i][j][k] = v[i][j][data[level].z_ghosts+k];

  reqs[0] = MPI_REQUEST_NULL; reqs[1] = MPI_REQUEST_NULL;
  if (myid!=data[level].back)
    MPI_Isend((void *) data[level].sbufxy[0][0],m*n*data[level].z_ghosts,MPI_DOUBLE,data[level].back,0,data[level].cart_comm,&reqs[0]);
  if (myid!=data[level].front)
    MPI_Irecv((void *) data[level].rbufxy[0][0],m*n*data[level].z_ghosts,MPI_DOUBLE,data[level].front,0,data[level].cart_comm,&reqs[1]);
  MPI_Waitall(2,reqs,stats);
  if( MPI_PROC_NULL != data[level].front )
    for (i=0;i<m;i++)
      for (j=0;j<n;j++)
	for (k=0;k<data[level].z_ghosts;k++)
	  v[i][j][o-data[level].z_ghosts+k] = data[level].rbufxy[i][j][k];
  
  /* Sending data to front */
  for (i=0;i<m;i++)
    for (j=0;j<n;j++)
      for (k=0;k<data[level].z_ghosts;k++)
	data[level].sbufxy[i][j][k] = v[i][j][o-2*data[level].z_ghosts+k];

  reqs[0] = MPI_REQUEST_NULL; reqs[1] = MPI_REQUEST_NULL;
  if (myid!=data[level].front)
    MPI_Isend((void *) data[level].sbufxy[0][0],m*n*data[level].z_ghosts,MPI_DOUBLE,data[level].front,0,data[level].cart_comm,&reqs[0]);
  if (myid!=data[level].back)
    MPI_Irecv((void *) data[level].rbufxy[0][0],m*n*data[level].z_ghosts,MPI_DOUBLE,data[level].back,0,data[level].cart_comm,&reqs[1]);
  MPI_Waitall(2,reqs,stats);
  if( MPI_PROC_NULL != data[level].back )
    for (i=0;i<m;i++)
      for (j=0;j<n;j++)
	for (k=0;k<data[level].z_ghosts;k++)
	  v[i][j][k] = data[level].rbufxy[i][j][k];
  
  if (data[level].periodic) {
    if (myid==data[level].back && myid==data[level].front) {
      for (i=0;i<m;i++)
	for (j=0;j<n;j++) 
	  for (k=0;k<data[level].z_ghosts;k++) {
	    v[i][j][k] = v[i][j][o-2*data[level].z_ghosts+k];
	    v[i][j][o-data[level].z_ghosts+k] = v[i][j][data[level].z_ghosts+k];
	  }
    }
    
    if (myid==data[level].lower && myid==data[level].upper) {
      for (i=0;i<m;i++)
	for (j=0;j<data[level].y_ghosts;j++)
	  for (k=0;k<o;k++) {
	    v[i][j][k] = v[i][n-2*data[level].y_ghosts+j][k];
	    v[i][n-data[level].y_ghosts+j][k] = v[i][data[level].y_ghosts+j][k];
	  }
    }

    if (myid==data[level].left && myid==data[level].right) {
      for (i=0;i<data[level].x_ghosts;i++)
	for (j=0;j<n;j++)
	  for (k=0;k<o;k++) {
	    v[i][j][k] = v[m-2*data[level].x_ghosts+i][j][k];
	    v[m-data[level].x_ghosts+i][j][k] = v[data[level].x_ghosts+i][j][k];
	  }
    }
  }
    
  return;
}
