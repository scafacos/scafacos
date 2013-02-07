#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <mpi.h>

int main(int argc, char **argv)
{
    int status    = MPI_SUCCESS;
    int provided  = MPI_THREAD_SINGLE;
    int requested = MPI_THREAD_MULTIPLE;
    status = MPI_Init_thread(&argc, &argv, requested, &provided);
    assert( status == MPI_SUCCESS && provided == requested );

    int rank = -1;
    int size = -1;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    if (rank==0) printf("hello from rank %d of %d \n", rank, size );
    fflush(stdout);

    int max = -1, reps = -1;
    max  = ( argc>1 ? atoi(argv[1]) : 1000000 );
    reps = ( argc>2 ? atoi(argv[2]) : 100 );
    if (rank==0) printf("max = %d doubles, repetitions = %d \n", max, reps );
    fflush(stdout);

    int namelen = 0;
    char procname[MPI_MAX_PROCESSOR_NAME];

    MPI_Get_processor_name( procname, &namelen );
    printf("%5d: MPI_Get_processor_name = %s \n", rank, procname );
    fflush(stdout);

    status = MPI_Barrier(MPI_COMM_WORLD);
    assert( status == MPI_SUCCESS );

    for ( int s=1 ; size<max ; s*=2 )
    {
        int bytes = s * sizeof(double);

        /* allocate RMA buffers */
        double * winbuf = NULL; 
        status = MPI_Alloc_mem( s * sizeof(double), MPI_INFO_NULL, &winbuf );
        assert( status == MPI_SUCCESS && winbuf != NULL );
        memset( winbuf, 'w', s );

        /* register remote pointers */
        MPI_Win win;
        status = MPI_Win_create( winbuf, s * sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &win );
        assert( status == MPI_SUCCESS );

        double * tmp = NULL;
        status = MPI_Alloc_mem( s * sizeof(double), MPI_INFO_NULL, &tmp );
        assert( status == MPI_SUCCESS && tmp != NULL );
        memset( winbuf, 't', s );

        status = MPI_Barrier(MPI_COMM_WORLD);
        assert( status == MPI_SUCCESS );

        if (rank==0)
        {    
            for ( int target=0 ; target<size ; target++ )
            {
                double dt_put = 0.0, dt_acc = 0.0;
                double t0 = 0.0, t1 = 0.0;
             
                int mpiassert = 0;
             
                t0 = MPI_Wtime();
                for ( int r=0 ; r<reps ; r++ )
                {
                    status = MPI_Win_lock( MPI_LOCK_EXCLUSIVE, 0, mpiassert, win );
                    assert( status == MPI_SUCCESS );
             
                    status = MPI_Put( tmp, s, MPI_DOUBLE, target, 0, s, MPI_DOUBLE, win );
                    assert( status == MPI_SUCCESS );
             
                    status = MPI_Win_unlock( 0, win );
                    assert( status == MPI_SUCCESS );
                }
                t1 = MPI_Wtime();
                dt_put = (t1-t0)/reps;
             
                t0 = MPI_Wtime();
                for ( int r=0 ; r<reps ; r++ )
                {
                    status = MPI_Win_lock( MPI_LOCK_EXCLUSIVE, 0, mpiassert, win );
                    assert( status == MPI_SUCCESS );
             
                    status = MPI_Accumulate( tmp, s, MPI_DOUBLE, target, 0, s, MPI_DOUBLE, MPI_SUM, win );
                    assert( status == MPI_SUCCESS );
             
                    status = MPI_Win_unlock( 0, win );
                    assert( status == MPI_SUCCESS );
                }
                t1 = MPI_Wtime();
                dt_acc += (t1-t0)/reps;
             
                printf( "target = %4d size = %10d put: %e sec %e MB/s acc: %e sec  %e MB/s\n", 
                         target, s, dt_put, 1e-6 * bytes/dt_put, dt_acc, 1e-6 * bytes/dt_acc );
                fflush(stdout);
            }
        }

        status = MPI_Free_mem(tmp);

        status = MPI_Win_free(&win);
        status = MPI_Free_mem(winbuf);
    }

    status = MPI_Barrier(MPI_COMM_WORLD);
    assert( status == MPI_SUCCESS );

    if (rank==0) printf("%d: all done\n",rank);
    MPI_Finalize();

    return(0);
}



