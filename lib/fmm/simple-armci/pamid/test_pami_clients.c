#ifdef __bgq__

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <assert.h>
#include <string.h>
#include <mpi.h>
#include <pami.h>

static size_t world_size, world_rank = -1;

#define PRINT_SUCCESS 1

#define TEST_ASSERT(c,m) \
        do { \
        if (!(c)) { \
                    printf(m" FAILED on rank %ld\n", world_rank); \
                    fflush(stdout); \
                  } \
        else if (PRINT_SUCCESS) { \
                    printf(m" SUCCEEDED on rank %ld\n", world_rank); \
                    fflush(stdout); \
                  } \
        sleep(1); \
        /*assert(c);*/ \
        } \
        while(0);

int main(int argc, char *argv[])
{
    int provided = MPI_THREAD_SINGLE;
    int mpi_status = MPI_SUCCESS;

    mpi_status = MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    TEST_ASSERT( mpi_status==MPI_SUCCESS && provided==MPI_THREAD_MULTIPLE , "MPI_Init_thread" );

    /**********************************************************************/

    pami_result_t pami_result = PAMI_SUCCESS;

    /* initialize the client */
    char * clientname = "A1D";
    pami_client_t a1d_client;
    pami_result = PAMI_Client_create(clientname, &a1d_client, NULL, 0);
    TEST_ASSERT(pami_result == PAMI_SUCCESS,"PAMI_Client_create");

    /* query properties of the client */
    pami_configuration_t config;

    config.name = PAMI_CLIENT_NUM_TASKS;
    pami_result = PAMI_Client_query( a1d_client, &config,1);
    TEST_ASSERT(pami_result == PAMI_SUCCESS,"PAMI_Client_query");
    world_size = config.value.intval;
    TEST_ASSERT( world_size > 1 , "world_size > 1" );

    config.name = PAMI_CLIENT_TASK_ID;
    pami_result = PAMI_Client_query( a1d_client, &config,1);
    TEST_ASSERT(pami_result == PAMI_SUCCESS,"PAMI_Client_query");
    world_rank = config.value.intval;
    printf("hello world from rank %ld of %ld \n", world_rank, world_size );
    fflush(stdout);

    size_t num_contexts = 0;
    config.name = PAMI_CLIENT_NUM_CONTEXTS;
    pami_result = PAMI_Client_query( a1d_client, &config, 1);
    TEST_ASSERT(pami_result == PAMI_SUCCESS,"PAMI_Client_query");
    num_contexts = config.value.intval;

	printf(" Number of Contexts are %d on Rank %d \n", num_contexts, world_rank);
	fflush(stdout);	

    /* initialize the contexts */
    pami_context_t * contexts;
    contexts = (pami_context_t *) malloc( num_contexts * sizeof(pami_context_t) );
    TEST_ASSERT( contexts!=NULL , "malloc( num_contexts * sizeof(pami_context_t) )" );

    pami_result = PAMI_Context_createv( a1d_client, NULL, 0, contexts, num_contexts );
    TEST_ASSERT( pami_result == PAMI_SUCCESS , "PAMI_Context_createv" );

    //mpi_status = MPI_Barrier(MPI_COMM_WORLD);
    //TEST_ASSERT( mpi_status==MPI_SUCCESS , "MPI_Barrier" );

    MPI_Finalize();

    return(0);
}
#else
#error WTF

int main()
{
    return(1);
}

#endif
