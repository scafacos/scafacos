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

#include <stdio.h>
#include <stdlib.h>
#include <a1.h>

#define MAX_DIM 1024 
#define ITERATIONS 100
#define SKIP 10
#define WINDOW 8 

int main() {

   int i, j, k, rank, nranks, msgsize;
   int dim;
   long bufsize;
   double **buffer;
   unsigned long long t_start, t_stop, t_latency, t_overlap;
   unsigned long long wait_start, wait_stop;
   int count[2], src_stride, trg_stride, stride_level, peer;
   double A[1024][1024], B[1024][1024], C[1024][1024];
   int m1,m2,m3;
   double expected, actual;
   A1_handle_t a1_handle;
   
   A1_Initialize(A1_THREAD_SINGLE); 

   rank = A1_Process_id(A1_GROUP_WORLD);
   nranks = A1_Process_total(A1_GROUP_WORLD);

   buffer = (double **) malloc (sizeof(double *) * nranks); 

   A1_Allocate_handle(&a1_handle);

   A1_Barrier_group(A1_GROUP_WORLD);

   bufsize = MAX_DIM * MAX_DIM * sizeof(double);
   A1_Alloc_segment((void **) &(buffer[rank]), bufsize);
   A1_Exchange_segments(A1_GROUP_WORLD, (void **) buffer);

   for(i=0; i< bufsize/sizeof(double); i++) {
       *(buffer[rank] + i) = 1.0 + rank;
   }

   if(rank == 0) {

      printf("A1_PutS Overlap - NbPutS + DGEMM + Wait. Time in cycles\n");
      printf("%30s %30s %22s %22s\n", "Msg Size", "Dimensions(array of doubles)", "Base Latency", "Overlaped Latency");
      fflush(stdout);

      src_stride = MAX_DIM*sizeof(double);
      trg_stride = MAX_DIM*sizeof(double);
      stride_level = 1;
 
      for(dim=1; dim<=MAX_DIM; dim*=2) {
 
         count[0] = dim*sizeof(double);
         count[1] = 512;
 
         peer = 1;          
  
         for(i=0; i<ITERATIONS+SKIP; i++) { 
 
            if(i == SKIP)
                t_start = A1_Time_cycles();              

             for(k=0; k<WINDOW; k++)      
             { 
                A1_NbPutS(peer, stride_level, count, (void *) buffer[rank],
                      &src_stride, (void *) buffer[peer], &trg_stride, a1_handle);
             }

             A1_Wait_handle(a1_handle);
 
         }
         t_stop = A1_Time_cycles();
         A1_Flush(peer);
         
         t_latency = (t_stop-t_start)/ITERATIONS;
  
         char temp[10];
         sprintf(temp,"%dX%d", count[1], dim);
         printf("%30d %30s %20lld", count[1]*count[0], temp, t_latency);
         fflush(stdout);
 
         t_start = A1_Time_cycles();
         for(i=0; i<ITERATIONS; i++) {

            for(k=0; k<WINDOW; k++)      
            {
               A1_NbPutS(peer, stride_level, count, (void *) buffer[rank],
                      &src_stride, (void *) buffer[peer], &trg_stride, a1_handle);
            } 
  
            wait_start = A1_Time_cycles();
            for(m1=0; m1<1024; m1++)
            {
               for(m2=0; m2<1024; m2++)
               {
                 for(m3=0; m3<1024; m3++)
                 {
                    C[m1][m2] +=  A[m1][m3] * B[m3][m2];
                    wait_stop = A1_Time_cycles();
                    if((wait_stop - wait_start) > t_latency)
                         break;
                 }
                 if((wait_stop - wait_start) > t_latency)
                      break;
               }
               if((wait_stop - wait_start) > t_latency)
                    break;
            } 
 
            A1_Wait_handle(a1_handle);
 
         }
         t_stop = A1_Time_cycles();
         A1_Flush(peer);
         t_overlap = (t_stop - t_start)/ITERATIONS;
 
         printf("%20lld \n", t_overlap);          
 
      }

   }

   A1_Barrier_group(A1_GROUP_WORLD);

   A1_Release_handle(a1_handle);

   A1_Release_segments(A1_GROUP_WORLD, (void *) buffer[rank]);

   A1_Free_segment((void *) buffer[rank]);

   A1_Finalize();

   return 0;
}
