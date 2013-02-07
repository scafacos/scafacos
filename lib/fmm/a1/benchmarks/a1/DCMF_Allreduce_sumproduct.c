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
#include <dcmf.h>
#include "dcmf_collectives.h"

#define MAX_MSG_SIZE 1024*1024

DCMF_Geometry_t geometry;

DCMF_Geometry_t *getGeometry (int comm)
{
  return &geometry;
}

void done(void *clientdata, DCMF_Error_t *error)
{
    --(*((uint32_t *) clientdata));
}

int main()
{

    int i, rank, nranks, msgsize, status, expected;
    long bufsize;
    int *src_buffer;
    int *trg_buffer;
    unsigned *ranks;
    DCMF_CollectiveProtocol_t barrier_protocol, lbarrier_protocol;
    DCMF_CollectiveProtocol_t allreduce_protocol, allreduce_notree_protocol;
    DCMF_Barrier_Configuration_t barrier_conf;
    DCMF_Allreduce_Configuration_t allreduce_conf;
    DCMF_CollectiveRequest_t crequest, crequest1, crequest2;
    DCMF_Callback_t done_callback;
    volatile unsigned allreduce_active = 0;

    DCMF_Messager_initialize();

    DCMF_Collective_initialize();

    rank = DCMF_Messager_rank();
    nranks = DCMF_Messager_size();

    ranks = (unsigned *) malloc(nranks * sizeof(int));
    for(i=0; i<nranks; i++)
         ranks[i] = i;

    bufsize = MAX_MSG_SIZE;
    src_buffer = (int *) malloc(bufsize);
    trg_buffer = (int *) malloc(bufsize);

    barrier_conf.protocol = DCMF_GI_BARRIER_PROTOCOL;
    barrier_conf.cb_geometry = getGeometry; 
    status = DCMF_Barrier_register(&barrier_protocol,
                                   &barrier_conf);

    barrier_conf.protocol = DCMF_LOCKBOX_BARRIER_PROTOCOL;
    barrier_conf.cb_geometry = getGeometry;
    status = DCMF_Barrier_register(&lbarrier_protocol,
                                   &barrier_conf);

    DCMF_CollectiveProtocol_t  *barrier_ptr, *lbarrier_ptr;
    barrier_ptr = &barrier_protocol;
    lbarrier_ptr  = &lbarrier_protocol;
    status = DCMF_Geometry_initialize(&geometry,
                                      0,
 				      ranks,
				      nranks,
       				      &barrier_ptr,
                                      1,
                                      &lbarrier_ptr,
                                      1,
       				      &crequest,
				      0, 
				      1);
       
    allreduce_conf.protocol = DCMF_TREE_ALLREDUCE_PROTOCOL;
    allreduce_conf.cb_geometry = getGeometry;
    allreduce_conf.reuse_storage = 1;
    status = DCMF_Allreduce_register(&allreduce_protocol,
                                      &allreduce_conf);
    if(status != DCMF_SUCCESS)
    { 
       printf("DCMF_Allreduce_register returned with error %d \n",
                 status);
       exit(-1);
    }

    allreduce_conf.protocol = DCMF_TORUS_BINOMIAL_ALLREDUCE_PROTOCOL;
    allreduce_conf.cb_geometry = getGeometry;
    allreduce_conf.reuse_storage = 1;
    status = DCMF_Allreduce_register(&allreduce_notree_protocol,
                                  &allreduce_conf);
    if(status != DCMF_SUCCESS)
    { 
       printf("DCMF_Allreduce_register returned with error %d \n",
                 status);
       exit(-1);
    }

    if(!DCMF_Geometry_analyze(&geometry, &allreduce_protocol))
    {
      printf("Not a supported geometry!! \n");
      fflush(stdout);
      return -1;
    }

    if(!DCMF_Geometry_analyze(&geometry, &allreduce_notree_protocol))
    {
      printf("Not a supported geometry!! \n");
      fflush(stdout);
      return -1;
    }

    done_callback.function = done;
    done_callback.clientdata = (void *) &allreduce_active;

    if (rank == 0)
    {
        printf("DCMF_Allreduce Test\n");
        fflush(stdout);
    }

    for (msgsize = sizeof(int); msgsize < MAX_MSG_SIZE; msgsize *= 2)
    {
            /*initializing buffer*/
            for (i = 0; i < bufsize/sizeof(int); i++)
            {
                 src_buffer[i] = rank;
                 trg_buffer[i] = 0;
            }

            allreduce_active += 1;

            /*sum reduce operation*/
            status = DCMF_Allreduce(&allreduce_protocol,
                                    &crequest1,
                                    done_callback,
                                    DCMF_SEQUENTIAL_CONSISTENCY,
                                    &geometry,
                                    (char *) src_buffer,
                                    (char *) trg_buffer,
                                    msgsize/sizeof(int),
                                    DCMF_SIGNED_INT,
                                    DCMF_SUM);

             while(allreduce_active > 0) DCMF_Messager_advance();

             expected = (nranks-1)*(nranks)/2;
             for (i = 0; i < msgsize/sizeof(int); i++)
             {
                if(trg_buffer[i] - expected != 0)
                {
                   printf("[%d] Validation has failed Expected: %d, Actual: %d, i: %d \n",
                               rank, expected, trg_buffer[i], i);
                   fflush(stdout);
                   exit(-1);
                }
             }

             printf("[%d] %d message sum allreduce successful \n", rank, msgsize);
             fflush(stdout);

             for (i = 0; i < bufsize/sizeof(int); i++)
             {
                   src_buffer[i] = 1;
                   trg_buffer[i] = 0;
             }

            allreduce_active += 1;

            /*sum reduce operation*/
            status = DCMF_Allreduce(&allreduce_notree_protocol,
                                    &crequest2,
                                    done_callback,
                                    DCMF_SEQUENTIAL_CONSISTENCY,
       				    &geometry,
                                    (char *) src_buffer,
                                    (char *) trg_buffer,
                                    msgsize/sizeof(int),
                                    DCMF_SIGNED_INT,
                                    DCMF_PROD);

             while(allreduce_active > 0) DCMF_Messager_advance();

             expected = 1;
             for (i = 0; i < msgsize/sizeof(int); i++)
             {
                if(trg_buffer[i] - expected != 0)
                {
                    printf("[%d] Validation has failed Expected: %d, Actual: %d, i: %d \n",
                                rank, expected, trg_buffer[i], i);
                    fflush(stdout);
                    exit(-1);
                }
             }

             printf("[%d] %d message product allreduce successful\n", rank, msgsize);
             fflush(stdout);

    }

    free(src_buffer);
    free(trg_buffer);
    DCMF_Messager_finalize();

    return 0;
}
