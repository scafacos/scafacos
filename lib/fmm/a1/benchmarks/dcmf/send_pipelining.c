/* 
* The following is a notice of limited availability of the code, and disclaimer
* which must be included in the prologue of the code and in all source listings
* of the code.
* 
* Copyright (c) 2010  Argonne Leadership Computing Facility, Argonne National 
* Laboratory
* 
* Permission is hereby granted to use, reproduce, prepare derivative works, and
* to redistribute to others.
* 
* 
*                          LICENSE
* 
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions are
* met:
* 
* - Redistributions of source code must retain the above copyright
*   notice, this list of conditions and the following disclaimer.
* 
* - Redistributions in binary form must reproduce the above copyright
*   notice, this list of conditions and the following disclaimer listed
*   in this license in the documentation and/or other materials
*   provided with the distribution.
* 
* - Neither the name of the copyright holders nor the names of its
*   contributors may be used to endorse or promote products derived from
*   this software without specific prior written permission.
* 
* The copyright holders provide no reassurances that the source code
* provided does not infringe any patent, copyright, or any other
* intellectual property rights of third parties.  The copyright holders
* disclaim any liability to any recipient for claims brought against
* recipient by any third party for infringement of that parties
* intellectual property rights.
* 
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
* "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
* LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
* A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
* OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
* SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
* LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
* DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
* THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
* (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
* OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. 
*/

#include "bench.h"

#define ITERATIONS_LOCAL 100
#define MAX_MSG_SIZE_LOCAL 1*1024*1024

void accumulate_pipelining() {

       DCMF_Request_t *snd_req;
       DCMF_Callback_t snd_done;
       DCQuad msginfo;
       volatile int done_count;
       int i, j, dst, msgsize;
       DCMF_NetworkCoord_t myaddr, dstaddr;
       DCMF_Network ntwk;
       int csize;

       csize = CHUNK_SIZE;

       snd_done.function = done;
       snd_done.clientdata = (void *) &done_count;

       snd_req = (DCMF_Request_t *) malloc(sizeof(DCMF_Request_t)*(MAX_MSG_SIZE_LOCAL/csize));

       if(snd_req == NULL) {
           printf("request memory allocation failed \n");
           fflush(stdout);
           exit(-1);
       }

       DCMF_Messager_rank2network(myrank, DCMF_TORUS_NETWORK, &myaddr);

       dstaddr.network = myaddr.network;
       dstaddr.torus.x = (myaddr.torus.x+3)%8;
       dstaddr.torus.y = (myaddr.torus.y+3)%8;
       dstaddr.torus.z = (myaddr.torus.z+3)%8;
       dstaddr.torus.t = myaddr.torus.t;

       DCMF_Messager_network2rank(&dstaddr, &dst, &ntwk);

       if(myrank == 0) {
          char buffer[100];
          printf("Accumulate Latency in usec, Chunk Size : %d \n", csize);
          sprintf(buffer,"%20s  %22s  %22s  %22s  %22s", "Msg Size", "(3,3,3)-Acc",\
               "(3,3,3) Pipelined-Acc", "(0,0,1)-Acc",
               "(0,0,1) Pipelined-Acc");
          printf("%s \n", buffer);
          fflush(stdout);
       }      
  
       for(msgsize=2*sizeof(double); msgsize<=MAX_MSG_SIZE_LOCAL; msgsize*=2) {
 
           t_avg = 0; t_avg1 = 0; done_count = 0; snd_rcv_active = 0;
           datasize = msgsize;
           ispipelined = 0;
           barrier();

           t_start = DCMF_Timebase();

           for(j=0; j<ITERATIONS_LOCAL; j++) {

                 done_count = 1;
                 snd_rcv_active = snd_rcv_active + 1;

                 DCMF_Send(&accumulate_snd_reg,
                    &snd_req[0],
                    snd_done,
                    DCMF_SEQUENTIAL_CONSISTENCY,
                    dst,
                    msgsize,
                    source + j*msgsize,
                    &msginfo,
                    1);

                 while(done_count > 0 || snd_rcv_active > 0) DCMF_Messager_advance();

           }

           t_stop = DCMF_Timebase();
           t_usec = (t_stop-t_start)/(clockMHz*ITERATIONS_LOCAL);

           barrier();
           allreduce(-1, (char *) &t_usec, (char *) &t_avg, 1, DCMF_DOUBLE, DCMF_SUM);
           target_index = 0;
           stage_index = 0;
           if(msgsize > csize)
               ispipelined = 1;
           barrier();

           t_start = DCMF_Timebase();

           for(j=0; j<ITERATIONS_LOCAL; j++) {

                 if(msgsize > csize) { 
                   done_count = msgsize/csize;
                   snd_rcv_active = snd_rcv_active + msgsize/csize;
                   for(i=0; i<done_count; i++) {
                       DCMF_Send(&accumulate_snd_reg,
                          &snd_req[i],
                          snd_done,
                          DCMF_SEQUENTIAL_CONSISTENCY,
                          dst,
                          csize,
                          source + j*msgsize + i*csize,
                          &msginfo,
                          1);
                   }
                 } else {
                   done_count = 1;
                   snd_rcv_active = snd_rcv_active + 1;
                   DCMF_Send(&accumulate_snd_reg,
                          &snd_req[0],
                          snd_done,
                          DCMF_SEQUENTIAL_CONSISTENCY,
                          dst,
                          msgsize,
                          source + j*msgsize,
                          &msginfo,
                          1);
                 }
                 while((done_count > 0) || (snd_rcv_active > 0)) DCMF_Messager_advance();

           }

           t_stop = DCMF_Timebase();
           t_usec1 = (t_stop-t_start)/(clockMHz*ITERATIONS_LOCAL);

           barrier();
           allreduce(-1, (char *) &t_usec1, (char *) &t_avg1, 1, DCMF_DOUBLE, DCMF_SUM);
           target_index = 0;
           stage_index = 0;
           ispipelined = 0;
           barrier();

           if(myrank == 0) {
              t_avg = t_avg/nranks;
              t_avg1 = t_avg1/nranks;
              printf("%20d %20.2f %20.2f ", msgsize, t_avg, t_avg1);
              fflush(stdout);
           } 

           barrier(); 

           t_start = DCMF_Timebase();

           for(j=0; j<ITERATIONS_LOCAL; j++) {

                 done_count = 1;
                 snd_rcv_active = snd_rcv_active + 1;

                 DCMF_Send(&accumulate_snd_reg,
                    &snd_req[0],
                    snd_done,
                    DCMF_SEQUENTIAL_CONSISTENCY,
                    (myrank+1)%nranks,
                    msgsize,
                    source + j*msgsize,
                    &msginfo,
                    1);

                 while(done_count > 0 || snd_rcv_active > 0) DCMF_Messager_advance();

           }

           t_stop = DCMF_Timebase();
           t_usec = (t_stop-t_start)/(clockMHz*ITERATIONS_LOCAL);

           barrier();
           allreduce(-1, (char *) &t_usec, (char *) &t_avg, 1, DCMF_DOUBLE, DCMF_SUM);
           target_index = 0;
           stage_index = 0;
           if(msgsize > csize)
               ispipelined = 1;
           barrier();

           t_start = DCMF_Timebase();

           for(j=0; j<ITERATIONS_LOCAL; j++) {

                 if(msgsize > csize) { 
                   done_count = msgsize/csize;
                   snd_rcv_active = snd_rcv_active + msgsize/csize;
                   for(i=0; i<done_count; i++) {
                       DCMF_Send(&accumulate_snd_reg,
                          &snd_req[i],
                          snd_done,
                          DCMF_SEQUENTIAL_CONSISTENCY,
                          (myrank+1)%nranks,
                          csize,
                          source + j*msgsize + i*csize,
                          &msginfo,
                          1);
                   }
                 } else {
                   done_count = 1;
                   snd_rcv_active = snd_rcv_active + 1;
                   DCMF_Send(&accumulate_snd_reg,
                          &snd_req[0],
                          snd_done,
                          DCMF_SEQUENTIAL_CONSISTENCY,
                          (myrank+1)%nranks,
                          msgsize,
                          source + j*msgsize,
                          &msginfo,
                          1);
                 }
                 while((done_count > 0) || (snd_rcv_active > 0)) DCMF_Messager_advance();

           }

           t_stop = DCMF_Timebase();
           t_usec1 = (t_stop-t_start)/(clockMHz*ITERATIONS_LOCAL);

           barrier();
           allreduce(-1, (char *) &t_usec1, (char *) &t_avg1, 1, DCMF_DOUBLE, DCMF_SUM);
           target_index = 0;
           stage_index = 0;
           ispipelined = 0;
           barrier();

           if(myrank == 0) {
              t_avg = t_avg/nranks;
              t_avg1 = t_avg1/nranks;
              printf("%20d %20.2f %20.2f \n", msgsize, t_avg, t_avg1);
              fflush(stdout);
           }

         }
}

int main ()
{

  DCMF_Messager_initialize();

  init();

  posix_memalign((void **) &source, 16, MAX_MSG_SIZE_LOCAL*ITERATIONS_LOCAL);
  posix_memalign((void **) &target, 16, MAX_MSG_SIZE_LOCAL*ITERATIONS_LOCAL);
  posix_memalign((void **) &stagebuf, 16, MAX_MSG_SIZE_LOCAL*ITERATIONS_LOCAL);

  DCMF_CriticalSection_enter(0);

  barrier_init(DCMF_DEFAULT_GLOBALBARRIER_PROTOCOL);

  allreduce_init(DCMF_DEFAULT_GLOBALALLREDUCE_PROTOCOL);

  accumulate_send_init(DCMF_DEFAULT_SEND_PROTOCOL, DCMF_TORUS_NETWORK);

  accumulate_pipelining();

  if(myrank == 0) {
    printf("[%d] Benchmark complete\n", myrank);
    fflush(stdout);
  }

  DCMF_CriticalSection_exit(0);

  DCMF_Messager_finalize ();

  free(source);
  free(target);
  free(stagebuf);  

  return 0;
}
