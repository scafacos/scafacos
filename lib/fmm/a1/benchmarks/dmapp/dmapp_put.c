/* -*- Mode: C; c-basic-offset:4 ; -*- */
/*
 *  Copyright (C) 2010 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 *
 * dmapp_put.c
 *
 *  Created on: Sep 18, 2010
 *      Author: jeff
 */

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include "pmi.h"
#include "dmapp.h"
#define MAX_NELEMS (128L*1024L)
/* If necessary, run the job with fewer than the maximum number of cores
 * per node so that enough memory is available for each PE. */
int main(int argc, char **argv)
{
    int pe = -1;
    int npes = -1;
    int target_pe;
    int fail_count = 0;
    long nelems = MAX_NELEMS;
    long *source;
    long *target;
    long i;
    dmapp_return_t status;
    dmapp_rma_attrs_t actual_args = { 0 }, rma_args = { 0 };
    dmapp_jobinfo_t job;
    dmapp_seg_desc_t *seg = NULL;
    /* Set the RMA parameters. */
    rma_args.put_relaxed_ordering = DMAPP_ROUTING_ADAPTIVE;
    rma_args.max_outstanding_nb = DMAPP_DEF_OUTSTANDING_NB;
    rma_args.offload_threshold = DMAPP_OFFLOAD_THRESHOLD;
    rma_args.max_concurrency = 1;
    /* Initialize DMAPP. */
    status = dmapp_init(&rma_args, &actual_args);
    if (status != DMAPP_RC_SUCCESS)
    {
        fprintf(stderr, " dmapp_init FAILED: %d\n", status);
        exit(1);
    }
    /* Allocate and initialize the source and target arrays. */
    source = (long *) dmapp_sheap_malloc(nelems * sizeof(long));
    target = (long *) dmapp_sheap_malloc(nelems * sizeof(long));
    if ((source == NULL) || (target == NULL))
    {
        fprintf(stderr,
                " sheap_malloc'd failed src 0x%lx targ 0x%lx\n",
                (long) source,
                (long) target);
        exit(1);
    }
    for (i = 0; i < nelems; i++)
    {
        source[i] = i;
        target[i] = -9L;
    }
    /* Wait for all PEs to complete array initialization. */
    PMI_Barrier();

    /* Get job related information. */
    status = dmapp_get_jobinfo(&job);
    if (status != DMAPP_RC_SUCCESS)
    {
        fprintf(stderr, " dmapp_get_jobinfo FAILED: %d\n", status);
        exit(1);
    }
    pe = job.pe;
    npes = job.npes;
    seg = &(job.sheap_seg);

    /* Send my data to my target PE. */
    target_pe = npes - pe - 1;
    status = dmapp_put(target, seg, target_pe, source, nelems, DMAPP_QW);
    if (status != DMAPP_RC_SUCCESS)
    {
        fprintf(stderr, " dmapp_put FAILED: %d\n", status);
        exit(1);
    }
    /* Wait for all PEs to complete their PUT. */
    PMI_Barrier();

    /* Check the results. */
    for (i = 0; i < nelems; i++)
    {
        if ((target[i] != i) && (fail_count < 10))
        {
            fprintf(stderr,
                    " PE %d: target[%ld] is %ld, should be %ld\n",
                    pe,
                    i,
                    target[i],
                    (long) i);
            fail_count++;
        }
    }
    if (fail_count == 0)
    {
        fprintf(stderr, " dmapp_put PASSED for PE %04d\n", pe);
    }
    else
    {
        fprintf(stderr, " dmapp_put FAILED for PE %04d: "
            "%d or more wrong values\n", pe, fail_count);
    }
    /* Finalize. */
    status = dmapp_finalize();
    return (0);
}
