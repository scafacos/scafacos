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
#include <assert.h>
#include <string.h>
#include <mpi.h>

#ifdef __CRAYXE
#include "pmi.h"
#include "dmapp.h"
#endif

int main(int argc,char **argv)
{
#ifdef __CRAYXE
    int                max;
    int                i;
    int                pe = -1;
    int                npes = -1;
    char *             source = NULL;
    char *             target = NULL;
    dmapp_return_t     status;
    //dmapp_rma_attrs_t  dmapp_config_in, dmapp_config_out;
    dmapp_rma_attrs_ext_t dmapp_config_in, dmapp_config_out;
    dmapp_jobinfo_t    job;
    dmapp_seg_desc_t * seg = NULL;

    double t0, t1, dt;
    double bw;

    MPI_Init(&argc, &argv);

    /* Initialize DMAPP resources before executing any other DMAPP calls. */
    //status = dmapp_init(NULL, &dmapp_config_out);

    dmapp_config_in.max_outstanding_nb   = DMAPP_DEF_OUTSTANDING_NB; /*  512 */
    dmapp_config_in.offload_threshold    = DMAPP_OFFLOAD_THRESHOLD;  /* 4096 */

    //dmapp_config_in.put_relaxed_ordering = DMAPP_ROUTING_DETERMINISTIC;
    //dmapp_config_in.get_relaxed_ordering = DMAPP_ROUTING_DETERMINISTIC;
    dmapp_config_in.put_relaxed_ordering = DMAPP_ROUTING_ADAPTIVE;
    dmapp_config_in.get_relaxed_ordering = DMAPP_ROUTING_ADAPTIVE;

    dmapp_config_in.max_concurrency      = 1; /* not thread-safe */

    //dmapp_config_in.PI_ordering          = DMAPP_PI_ORDERING_STRICT;
    dmapp_config_in.PI_ordering          = DMAPP_PI_ORDERING_RELAXED;

    status = dmapp_init_ext( &dmapp_config_in, &dmapp_config_out );
    assert(status==DMAPP_RC_SUCCESS);

    max = (argc>1) ? atoi(argv[1]) : 1000000;
    max *= 16; /* max must be a multiple of 16 for the test to work */

    /* Allocate remotely accessible memory for source and target buffers.
           Only memory in the data segment or the sheap is remotely accessible.
           Here we allocate from the sheap. */
    source = (char *)dmapp_sheap_malloc( max*sizeof(char) );
    target = (char *)dmapp_sheap_malloc( max*sizeof(char) );
    assert( (source!=NULL) && (target!=NULL));

    memset (source,'S',max);
    memset (target,'T',max);

    /* Retrieve information about job details, such as PE id and number of PEs. */
    status = dmapp_get_jobinfo(&job);
    assert(status==DMAPP_RC_SUCCESS);
    pe = job.pe;
    npes = job.npes;

    /* Retrieve information about RMA attributes, such as offload_threshold
           and routing modes. */
    //status = dmapp_get_rma_attrs(&dmapp_config_out);
    status = dmapp_get_rma_attrs_ext(&dmapp_config_out);
    assert(status==DMAPP_RC_SUCCESS);

    /* Specify in which segment the remote memory region (the source) lies.
           In this case, it is the sheap (see above). */
    seg = &(job.sheap_seg);

    if (pe == 0) fprintf(stderr," Hello from PE %d of %d, using seg start %p, seg size 0x%lx, offload_threshold %d \n",
            pe, npes, seg->addr, (unsigned long)seg->len, dmapp_config_out.offload_threshold);
    fflush(stderr);
    PMI_Barrier();

    if (pe == 0)
    {
        fprintf(stderr,"%d: max = %d bytes, dmapp_get using DMAPP_DQW \n", pe, max);
        for (i=1; i<(max/16); i*=2)
        {
            t0 = MPI_Wtime();
            status = dmapp_get(target, source, seg, 1, i, DMAPP_DQW);
            t1 = MPI_Wtime();
            assert(status==DMAPP_RC_SUCCESS);
            dt = t1-t0;
            bw = 16 * 1e-6 * (double)i / dt;
            fprintf(stderr,"%d: %12d bytes %12lf seconds = %lf MB/s \n", pe, 16*i, dt, bw);
        }

    }
    fflush(stderr);
    PMI_Barrier();

    if (pe == 0)
    {
        fprintf(stderr,"%d: max = %d bytes, dmapp_get using DMAPP_QW \n", pe, max);
        for (i=1; i<(max/8); i*=2)
        {
            t0 = MPI_Wtime();
            status = dmapp_get(target, source, seg, 1, i, DMAPP_QW);
            t1 = MPI_Wtime();
            assert(status==DMAPP_RC_SUCCESS);
            dt = t1-t0;
            bw = 8 * 1e-6 * (double)i / dt;
            fprintf(stderr,"%d: %12d bytes %12lf seconds = %lf MB/s \n", pe, 8*i, dt, bw);
        }

    }
    fflush(stderr);
    PMI_Barrier();

    if (pe == 0)
    {
        fprintf(stderr,"%d: max = %d bytes, dmapp_get using DMAPP_DW \n", pe, max);
        for (i=1; i<(max/4); i*=2)
        {
            t0 = MPI_Wtime();
            status = dmapp_get(target, source, seg, 1, i, DMAPP_DW);
            t1 = MPI_Wtime();
            assert(status==DMAPP_RC_SUCCESS);
            dt = t1-t0;
            bw = 4 * 1e-6 * (double)i / dt;
            fprintf(stderr,"%d: %12d bytes %12lf seconds = %lf MB/s \n", pe, 4*i, dt, bw);
        }

    }
    fflush(stderr);
    PMI_Barrier();

    if (pe == 0)
    {
        fprintf(stderr,"%d: max = %d bytes, dmapp_get using DMAPP_BYTE \n", pe, max);
        for (i=1; i<max; i*=2)
        {
            t0 = MPI_Wtime();
            status = dmapp_get(target, source, seg, 1, i, DMAPP_BYTE);
            t1 = MPI_Wtime();
            assert(status==DMAPP_RC_SUCCESS);
            dt = t1-t0;
            bw = 1 * 1e-6 * (double)i / dt;
            fprintf(stderr,"%d: %12d bytes %12lf seconds = %lf MB/s \n", pe, 1*i, dt, bw);
        }

    }
    fflush(stderr);
    PMI_Barrier();

    /* Free buffers allocated from sheap. */
    dmapp_sheap_free(target);
    dmapp_sheap_free(source);

    /* Release DMAPP resources. This is a mandatory call. */
    status = dmapp_finalize();
    assert(status==DMAPP_RC_SUCCESS);

    MPI_Finalize();
#endif
    return(0);
}

