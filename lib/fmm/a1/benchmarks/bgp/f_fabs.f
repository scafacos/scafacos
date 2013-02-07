        program f_fabs
        implicit none
        include 'mpif.h'
        integer ierror
        integer i,j,k
        double precision x(1024), y1(1024), y2(1024)
        complex(8) z1,z2
        double precision t0,t1

        call MPI_INIT(ierror)

        k = 1

        write(6,*) 'FORTRAN VECTOR ABS TEST'

        write(6,*) 'MPI_WTICK = ',MPI_WTICK()

        do i=1,1024
            x(i) = (i/(i+1))*((-1)**i)
        enddo
        do i=1,1024
            y1(i) = 0.0
        enddo
        do i=1,1024
            y2(i) = 0.0
        enddo

c       WARM-UP
        do i=1,1024
            y1(i) = dabs(x(i))
        enddo

c       TIMING
        t0 = MPI_WTIME()
        do j=1,k
            do i=1,1024
                y1(i) = dabs(x(i))
            enddo
        enddo
        t1 = MPI_WTIME()
        write(6,*) 'fabs version=',t1-t0

        call MPI_FINALIZE()

        write(6,*) 'ALL DONE'
        end
