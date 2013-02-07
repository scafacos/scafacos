/*
 * Copyright (c) 2009 Cray Inc. All Rights Reserved.
 *
 * The contents of this file is proprietary information of Cray Inc.
 * and may not be disclosed without prior written consent.
 *
 * $HeadURL: http://svn.us.cray.com/svn/baker/packages/dmapp/trunk/tests0/dmapp_sample_get.c $
 * $LastChangedRevision: 1671 $
 *
 * Simple DMAPP Get test - blocking version
 * This test uses the external job launch package PMI.
 *
 */

#include <stdio.h>
#include <stdlib.h>

#ifdef __CRAYXE
#include "pmi.h"
#include "dmapp.h"
#endif

int main(int argc,char **argv)
{
#ifdef __CRAYXE
        int               nelems = 128;
        int               i;
        int               pe = -1;
        int               npes = -1;
        int               fail_count = 0;
        long              *source = NULL;
        long              *target = NULL;
        dmapp_return_t    status;
        dmapp_rma_attrs_t actual_args;
        dmapp_jobinfo_t   job;
        dmapp_seg_desc_t  *seg = NULL;

        /* Initialize DMAPP resources before executing any other DMAPP calls. */
        status = dmapp_init(NULL, &actual_args);
        if (status != DMAPP_RC_SUCCESS) {
                fprintf(stderr,"\n dmapp_init FAILED: %d\n", status);
                exit(1);
        }

        /* Allocate remotely accessible memory for source and target buffers.
           Only memory in the data segment or the sheap is remotely accessible.
           Here we allocate from the sheap. */
        source = (long *)dmapp_sheap_malloc(nelems*sizeof(long));
        target = (long *)dmapp_sheap_malloc(nelems*sizeof(long));
        if ((source == NULL) || (target == NULL)) {
                fprintf(stderr,"\n dmapp_sheap_malloc FAILED\n");
                exit(1);
        }

        for (i=0; i<nelems; i++) {
                source[i] = i;
                target[i] = -9L;
        }

        /* Synchronize to make sure everyone's buffers are initialized before
           data transfer is started. */
        PMI_Barrier();

        /* Retrieve information about job details, such as PE id and number of PEs. */
        status = dmapp_get_jobinfo(&job);
        if (status != DMAPP_RC_SUCCESS) {
                fprintf(stderr,"\n dmapp_get_jobinfo FAILED: %d\n", status);
                exit(1);
        }
        pe = job.pe;
        npes = job.npes;

        /* Retrieve information about RMA attributes, such as offload_threshold
           and routing modes. */
        status = dmapp_get_rma_attrs(&actual_args);
        if (status != DMAPP_RC_SUCCESS) {
                fprintf(stderr,"\n dmapp_get_rma_attrs FAILED: %d\n", status);
                exit(1);
        }

        /* Specify in which segment the remote memory region (the source) lies.
           In this case, it is the sheap (see above). */
        seg = &(job.sheap_seg);

        fprintf(stderr," Hello from PE %d of %d, using seg start %p, seg size 0x%lx, offload_threshold %d\n",
                pe, npes, seg->addr, (unsigned long)seg->len, actual_args.offload_threshold);

        fprintf(stderr,"\n PE %d getting %d nelems from addr %p on PE %d to local addr %p",
                pe, nelems, (void *)source, npes-pe-1, (void *)source);

        /* Execute GET operation from remote memory region source on PE Y 
           into local memory region target on PE X. */
        status = dmapp_get(target, source, seg, npes-pe-1, nelems, DMAPP_QW);
        if (status != DMAPP_RC_SUCCESS) {
                fprintf(stderr,"\n dmapp_get FAILED: %d\n", status);
                exit(1);
        }

        /* Synchronize before verifying the data. */
        PMI_Barrier();

        /* Verify data received in target buffer. */
        for (i=0; i<nelems; i++) {
                if (target[i] != i) {
                        fprintf(stderr,"\n PE %d: target[%d] is %ld, should be %ld",
                                pe, i, target[i], (long)i);
                        fail_count++;
                }
        }
        if (fail_count == 0)
                fprintf(stderr,"\n dmapp_sample_get PASSED\n");
        else
                fprintf(stderr,"\n dmapp_sample_get FAILED: %d wrong values\n",
                        fail_count);

        /* Free buffers allocated from sheap. */
        dmapp_sheap_free(target);
        dmapp_sheap_free(source);

        /* Release DMAPP resources. This is a mandatory call. */
        status = dmapp_finalize();
        if (status != DMAPP_RC_SUCCESS) {
                fprintf(stderr,"\n dmapp_finalize FAILED: %d\n", status);
                exit(1);
        }
#endif
        return(0);
}

