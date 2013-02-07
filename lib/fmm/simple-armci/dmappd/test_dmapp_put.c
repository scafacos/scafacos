/* Code from which this is derived is copyright 2010 Cray Inc. */

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <assert.h>
#include <math.h>

#ifdef __CRAYXE
#include "pmi.h"
#include "dmapp.h"
#endif

#define MAX_NELEMS (128L*1024L)

/* If necessary, run the job with fewer than the maximum number of cores
 * per node so that enough memory is available for each PE. */

int main(int argc, char **argv)
{
#ifdef __CRAYXE
    int pe = -1;
    int npes = -1;
    int target_pe;
    int fail_count = 0;

    long i, nelems = MAX_NELEMS;

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
    assert(status==DMAPP_RC_SUCCESS);

    /* Get job related information. */
    status = dmapp_get_jobinfo(&job);
    assert(status==DMAPP_RC_SUCCESS);
    pe = job.pe;
    npes = job.npes;
    seg = &(job.sheap_seg);

    fprintf(stderr," Hello from PE %d of %d, using seg start %p, seg size 0x%lx \n",
                                    pe,  npes, seg->addr, (unsigned long)seg->len );
    fflush(stderr);
    PMI_Barrier();

    /* Allocate and initialize the source and target arrays. */
    long * source = (long *) dmapp_sheap_malloc(nelems * sizeof(long));
    assert(source!=NULL);
    fprintf(stderr,"dmapp_sheap_malloc returned source = %p \n", source );
    fflush(stderr);
    PMI_Barrier();

    long * target = (long *) dmapp_sheap_malloc(nelems * sizeof(long));
    assert(target!=NULL);
    fprintf(stderr,"dmapp_sheap_malloc returned target = %p \n", target );
    fflush(stderr);
    PMI_Barrier();

    for (i = 0; i < nelems; i++) source[i] = i;
    for (i = 0; i < nelems; i++) target[i] = -9L;

    /* Wait for all PEs to complete array initialization. */
    PMI_Barrier();

    /* Send my data to my target PE. */
    target_pe = npes - pe - 1;
    status = dmapp_put(target, seg, target_pe, source, nelems, DMAPP_QW);
    assert(status==DMAPP_RC_SUCCESS);

    /* Wait for all PEs to complete their PUT. */
    PMI_Barrier();

    /* Check the results. */
    for (i = 0; i < nelems; i++)
        if ((target[i] != i) && (fail_count < 10))
        {
            fprintf(stderr," PE %d: target[%ld] is %ld, should be %ld\n", pe, i, target[i], (long) i);
            fail_count++;
        }

    if (fail_count == 0) fprintf(stderr, " dmapp_put PASSED for PE %04d\n", pe);
    else fprintf(stderr, " dmapp_put FAILED for PE %04d: %d or more wrong values\n", pe, fail_count);
    fflush(stderr);

    /* Finalize. */
    status = dmapp_finalize();
    assert(status==DMAPP_RC_SUCCESS);

#endif
    return(0);
}
