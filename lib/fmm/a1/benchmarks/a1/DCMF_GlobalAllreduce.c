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
#include <dcmf_globalcollectives.h>

#define MAX_MSG_SIZE 1024*1024

void done(void *clientdata, DCMF_Error_t *error)
{
    --(*((uint32_t *) clientdata));
}

int main()
{

    int i, rank, nranks, msgsize, status, expected;
    long bufsize;
    int *buffer;
    DCMF_Protocol_t ga_protocol;
    DCMF_GlobalAllreduce_Configuration_t ga_conf;
    DCMF_Request_t request;
    DCMF_Callback_t done_callback;
    volatile unsigned ga_active = 0;

    DCMF_Messager_initialize();

    rank = DCMF_Messager_rank();
    nranks = DCMF_Messager_size();

    bufsize = MAX_MSG_SIZE;
    buffer = (int *) malloc(bufsize);

    ga_conf.protocol = DCMF_DEFAULT_GLOBALALLREDUCE_PROTOCOL;
    status = DCMF_GlobalAllreduce_register(&ga_protocol,
                                           &ga_conf);
    if(status != DCMF_SUCCESS)
    { 
       printf("DCMF_GlobalAllreduce_register returned with error %d \n",
                 status);
       exit(-1);
    }

    done_callback.function = done;
    done_callback.clientdata = (void *) &ga_active;

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
                 buffer[i] = rank;
            }

            ga_active += 1;

            /*sum reduce operation*/
            status = DCMF_GlobalAllreduce(&ga_protocol,
                                          &request,
                                          done_callback,
                                          DCMF_SEQUENTIAL_CONSISTENCY,
                                          -1,
                                          (char *) buffer,
                                          (char *) buffer,
                                          msgsize/sizeof(int),
                                          DCMF_SIGNED_INT,
                                          DCMF_SUM);

             while(ga_active > 0) DCMF_Messager_advance();

             expected = (nranks-1)*(nranks)/2;
             for (i = 0; i < msgsize/sizeof(int); i++)
             {
                if(buffer[i] - expected != 0)
                {
                   printf("[%d] Validation has failed Expected: %d, Actual: %d, i: %d \n",
                               rank, expected, buffer[i], i);
                   fflush(stdout);
                   exit(-1);
                }
             }

             printf("[%d] %d message sum reduce successful \n", rank, msgsize);
             fflush(stdout);

             for (i = 0; i < bufsize/sizeof(int); i++)
             {
                   buffer[i] = 1;
             }

            ga_active += 1;

            status = DCMF_GlobalAllreduce(&ga_protocol,
                                          &request,
                                          done_callback,
                                          DCMF_SEQUENTIAL_CONSISTENCY,
                                          -1,
                                          (char *) buffer,
                                          (char *) buffer,
                                          msgsize/sizeof(int),
                                          DCMF_SIGNED_INT,
                                          DCMF_PROD);

             while(ga_active > 0) DCMF_Messager_advance();

             expected = 1;
             for (i = 0; i < msgsize/sizeof(int); i++)
             {
                if(buffer[i] - expected != 0)
                {
                    printf("[%d] Validation has failed Expected: %d, Actual: %d, i: %d \n",
                                rank, expected, buffer[i], i);
                    fflush(stdout);
                    exit(-1);
                }
             }

             printf("[%d] %d message product reduce successful\n", rank, msgsize);
             fflush(stdout);

    }

    free(buffer);
    DCMF_Messager_finalize();

    return 0;
}
