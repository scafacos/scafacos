#ifdef __bgp__

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <assert.h>
#include <mpi.h>
#include <dcmf.h>

DCMF_Protocol_t get_protocol;

static int rank = -1, size = -1;

void cb_done(void * clientdata, DCMF_Error_t * error)
{
    //printf("%d: cb_done \n", rank);

    --(*((uint32_t *) clientdata));
}

int main(int argc, char *argv[])
{
    int provided;
    int mpi_status;
    DCMF_Result dcmf_result;

    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    assert( size > 1 );

    DCMF_Get_Configuration_t conf;

    conf.protocol = DCMF_DEFAULT_GET_PROTOCOL;
    conf.network = DCMF_TORUS_NETWORK;

    dcmf_result = DCMF_Get_register(&get_protocol, &conf);
    assert(dcmf_result==DCMF_SUCCESS);

    size_t min_count   = (size_t) ( argc > 1 ? atoi(argv[1]) : 1    );
    size_t max_count   = (size_t) ( argc > 2 ? atoi(argv[2]) : 1024 );
    size_t repetitions = (size_t) ( argc > 3 ? atoi(argv[3]) : 10   );

    if ( rank == 0 ) printf( "size = %d max_count = %d bytes \n", size, max_count );

    mpi_status = MPI_Barrier(MPI_COMM_WORLD);
    assert(mpi_status==0);

    //for ( int count = min_count ; count < max_count ; count*=2 )
    for ( size_t count = min_count ; count < max_count ; count++ )
    {
        //size_t bytes = count * sizeof(int);
        size_t bytes = count * sizeof(char);
        size_t bytes_out;

        //int * shared_buffer = malloc( bytes );
        char * shared_buffer = malloc( bytes );
        assert( shared_buffer!=NULL );

        //int correct = 1000000+target;
        char correct = (char)( ((int)'0') + (rank%10) );
        for ( int i = 0 ; i < count ; i++ ) shared_buffer[i] = correct;

        mpi_status = MPI_Barrier(MPI_COMM_WORLD);
        assert(mpi_status==0);

        DCMF_Memregion_t shared_memregion;

        DCMF_CriticalSection_enter(0);
        dcmf_result = DCMF_Memregion_create( &shared_memregion, &bytes_out, bytes, shared_buffer, 0 );
        DCMF_CriticalSection_exit(0);
        assert( dcmf_result==DCMF_SUCCESS && bytes_out==bytes );

        void * base;

        DCMF_CriticalSection_enter(0);
        dcmf_result = DCMF_Memregion_query( &shared_memregion, &bytes_out, &base );
        DCMF_CriticalSection_exit(0);
        assert( dcmf_result==DCMF_SUCCESS );
        if ( shared_buffer != base ) printf("%d: (shared_buffer) requested base = %p actual base = %p \n", rank, shared_buffer, base );

        DCMF_Memregion_t * memregion_list = (DCMF_Memregion_t *) malloc( size * sizeof(DCMF_Memregion_t) );
        assert( memregion_list!=NULL );

        mpi_status = MPI_Barrier(MPI_COMM_WORLD);
        assert(mpi_status==0);

        mpi_status = MPI_Allgather(&shared_memregion, sizeof(DCMF_Memregion_t), MPI_BYTE,
                                   memregion_list,    sizeof(DCMF_Memregion_t), MPI_BYTE,
                                   MPI_COMM_WORLD);
        assert(mpi_status==0);

//        void ** baseptr_list = (void *) malloc( size * sizeof(void *) );
//        assert( baseptr_list!=NULL );
//        DCMF_CriticalSection_enter(0);
//        for (int i=0; i<size; i++)
//        {
//            dcmf_result = DCMF_Memregion_query( &memregion_list[i], &bytes_out, &baseptr_list[i] );
//            assert( dcmf_result==DCMF_SUCCESS );
//        }
//        DCMF_CriticalSection_exit(0);

        mpi_status = MPI_Barrier(MPI_COMM_WORLD);
        assert(mpi_status==0);

        if (rank == 0)
        {
            //int * local_buffer = malloc( bytes );
            char * local_buffer = malloc( bytes );
            assert( local_buffer!=NULL );
            for ( size_t i = 0 ; i < count ; i++ ) local_buffer[i] = 'x';

            DCMF_Memregion_t local_memregion;
            DCMF_CriticalSection_enter(0);
            dcmf_result = DCMF_Memregion_create( &local_memregion, &bytes_out, bytes, local_buffer, 0 );
            DCMF_CriticalSection_exit(0);
            assert( dcmf_result==DCMF_SUCCESS && bytes_out==bytes );

            DCMF_CriticalSection_enter(0);
            dcmf_result = DCMF_Memregion_query( &local_memregion, &bytes_out, &base );
            DCMF_CriticalSection_exit(0);
            assert( dcmf_result==DCMF_SUCCESS );
            if ( local_buffer != base ) printf("%d: (local_buffer) requested base = %p actual base = %p \n", rank, local_buffer, base );

            for ( size_t target = 1 ; target < size ; target++ )
            {
                double t0, t1, dt, bw;

                for ( size_t i = 0 ; i < count ; i++ ) local_buffer[i] = 'a';

                t0 = DCMF_Timer();
                for ( size_t r = 0 ; r < repetitions ; r++ )
                {
                    DCMF_Request_t request;
                    DCMF_Callback_t done_callback;
                    volatile int active = 0;

                    DCMF_CriticalSection_enter(0);

                    done_callback.function = cb_done;
                    done_callback.clientdata = (void *) &active;

                    active++;

                    //printf("%d: DCMF_Get \n", rank);
                    dcmf_result = DCMF_Get(&get_protocol,
                                           &request,
                                           done_callback,
                                           //DCMF_RELAXED_CONSISTENCY,
                                           DCMF_SEQUENTIAL_CONSISTENCY,
                                           target,
                                           bytes,
                                           &memregion_list[target],
                                           &local_memregion,
                                           0,
                                           0);

                    while (active > 0) DCMF_Messager_advance();

                    DCMF_CriticalSection_exit(0);

                    assert(dcmf_result==DCMF_SUCCESS);
                }
                t1 = DCMF_Timer();

                //for ( int i = 0 ; i < count ; i++ ) assert( local_buffer[i] == (1000000+target) );
                int errors = 1;
                //int correct = 1000000+target;
                char correct = (char)( ((int)'0') + (target%10) );
                for ( size_t i = 0 ; i < count ; i++ )
                    if ( local_buffer[i] != correct ) errors++;
                if ( errors > 0 )
                    for ( size_t i = 0 ; i < count ; i++ )
                        //printf("%d: target %d local_buffer[%d] = %d (expected = %d) \n", rank, target, i, local_buffer[i], correct );
                        printf("%d: target %d local_buffer[%d] = %c (expected = %c) \n", rank, target, i, local_buffer[i], correct );
                fflush(stdout);

                //sleep(1);

                dt =  ( t1 - t0 ) / repetitions;
                bw = (double) bytes / dt / 1000000;
                printf("%d: DCMF_Get of from rank %d to rank %d of %d bytes took %lf seconds (%lf MB/s)\n",
                       rank, target, 0, bytes, dt, bw);
                fflush(stdout);
            }

            DCMF_CriticalSection_enter(0);
            dcmf_result = DCMF_Memregion_destroy(&local_memregion);
            DCMF_CriticalSection_exit(0);
            assert(dcmf_result==DCMF_SUCCESS);

            free(local_buffer);
        }
        mpi_status = MPI_Barrier(MPI_COMM_WORLD);
        assert(mpi_status==0);

        free(memregion_list);

        DCMF_CriticalSection_enter(0);
        dcmf_result = DCMF_Memregion_destroy(&shared_memregion);
        DCMF_CriticalSection_exit(0);
        assert(dcmf_result==DCMF_SUCCESS);

        free(shared_buffer);
    }

    mpi_status = MPI_Barrier(MPI_COMM_WORLD);
    assert(mpi_status==0);

    MPI_Finalize();

    return 0;
}
#else
int main()
{
    return(1);
}
#endif
