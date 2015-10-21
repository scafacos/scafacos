#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "mp3.h"
#include <armci.h>

int me, nprocs;
int LOOP=10;

int main(int argc, char **argv) {
  int i;
  double **myptrs;
  double t0, t1, tnbget=0, tnbwait=0, t2=0;

  MP_INIT(argc,argv);
  ARMCI_Init();

  MP_PROCS(&nprocs);
  MP_MYID(&me);

  if (nprocs < 2)
    ARMCI_Error("This program requires at least to processes", 1);

  myptrs = (double **)malloc(sizeof(double *)*nprocs);
  ARMCI_Malloc((void **)myptrs, LOOP*sizeof(double)); 
  
  MP_BARRIER();
  
  if(me == 0) {
    for(i = 0; i < 10; i++) {
      // This is a bug:
      // ARMCI_Get(myptrs[me]+i,myptrs[me+1]+i,sizeof(double),me+1);
      ARMCI_Get(myptrs[me+1]+i, myptrs[me]+i, sizeof(double), me+1);
    }

    t0 = MP_TIMER(); 
    for(i = 0; i < LOOP; i++) {
      // This is a bug:
      // ARMCI_Get(myptrs[me]+i,myptrs[me+1]+i,sizeof(double),me+1);
      ARMCI_Get(myptrs[me+1]+1, myptrs[me]+i, sizeof(double), me+1);
    }
    t1 = MP_TIMER(); 

    printf("\nGet Latency=%lf\n", 1e6*(t1-t0)/LOOP);
    fflush(stdout);

    t1 = t0 = 0;

    for(i = 0; i < LOOP; i++) {
      armci_hdl_t nbh;
      ARMCI_INIT_HANDLE(&nbh);

      t0 = MP_TIMER(); 
      //ARMCI_NbGet(myptrs[me]+i, myptrs[me+1]+i, sizeof(double), me+1, &nbh);
      ARMCI_NbGet(myptrs[me+1]+i, myptrs[me]+i, sizeof(double), me+1, &nbh);
      t1 = MP_TIMER(); 
      ARMCI_Wait(&nbh);
      t2 = MP_TIMER();

      tnbget  += (t1-t0);
      tnbwait += (t2-t1);
    }

    printf("\nNb Get Latency=%lf Nb Wait=%lf\n",1e6*tnbget/LOOP,1e6*tnbwait/LOOP);fflush(stdout);
  }

  else
    sleep(1);

  MP_BARRIER();

  ARMCI_Finalize();
  MP_FINALIZE();

  return 0;
}
