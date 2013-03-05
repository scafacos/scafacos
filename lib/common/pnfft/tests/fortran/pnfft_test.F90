
      program main

      implicit none
      include "pnfft.f"

      integer np(3), m, window, window_flag, ierror, l
      integer(ptrdiff_t_kind) :: N(3), local_N(3), local_N_start(3)
      integer(ptrdiff_t_kind) :: num_nodes, local_num_nodes
      integer(8) pnfft
      integer comm_cart_3d
      complex(8), allocatable ::  f1(:), f2(:)
      real(8) f_hat_sum, x_max(3)
      real(8) :: lower_border(3), upper_border(3)
      complex(8), allocatable :: f_hat(:), f(:)
      real(8), allocatable :: x(:)
      integer myrank


      N = (/ 16,16,16 /)
      np = (/ 2,2,2 /)
      local_num_nodes = N(1)*N(2)*N(3) / (np(1)*np(2)*np(3))

!     Initialize MPI and PFFT
      call MPI_Init(ierror)
      call dpnfft_init();
      call MPI_Comm_rank(MPI_COMM_WORLD, myrank, ierror)
      
!     Create three-dimensional process grid of
!     size np(1) x np(2) x np(3), if possible
      call dpnfft_create_procmesh(ierror, 3, MPI_COMM_WORLD, &
     &     np, comm_cart_3d)
      if (ierror .ne. 0) then
        if(myrank .eq. 0) then
          write(*,*) "Error: This test file only works with 8 processes"
        endif
        call MPI_Finalize(ierror)
        call exit(1)
      endif

!     Get parameters of data distribution
      call dpnfft_local_size_3d( N, comm_cart_3d, &
     &     local_N, local_N_start, lower_border, upper_border)

!     Allocate memory

!     Plan parallel NFFT
      call dpnfft_init_3d(pnfft, &
     &     N, local_num_nodes, comm_cart_3d)

!     Get data pointers
      call dpnfft_get_f_hat(f_hat, pnfft)
      call dpnfft_get_f(f, pnfft)
      call dpnfft_get_x(x, pnfft)

!     Initialize Fourier coefficients
      call dpnfft_init_f_hat_3d(N, local_N, local_N_start, &
     &     f_hat)

!     Initialize nonequispaced nodes
      call init_random_x(lower_border, upper_border, local_num_nodes, &
     &    x)

!     Print input Fourier coefficents
      call vpr_complex(comm_cart_3d, 8, f_hat, &
     &      "Input Fourier coefficients on process 1:")

!     Execute parallel NFFT
      call dpnfft_trafo(pnfft)

!     Print NFFT results
      call vpr_complex(comm_cart_3d, 8, f, &
     &      "PNFFT Results on process 1:")

!     Execute parallel adjoint NFFT
      call dpnfft_adj(pnfft)

!     Scale data
      do l=1,local_N(1) * local_N(2) * local_N(3)
        f_hat(l) = f_hat(l) / (N(1)*N(2)*N(3))
      enddo

!     Print output Fourier coefficents
      call vpr_complex(comm_cart_3d, 8, f_hat, &
     &      "Fourier coefficients after one forward &
     & and backward PNFFT on process 1:")

!     Free mem and finalize      
      call dpnfft_finalize(pnfft, &
     &      PNFFT_FREE_X + PNFFT_FREE_F_HAT + PNFFT_FREE_F)
      call MPI_Comm_free(comm_cart_3d, ierror)
      call dpnfft_cleanup()
      call MPI_Finalize(ierror)
      end


      subroutine vpr_complex(comm, num, vec, info)
        implicit none
        include "pnfft.f"

        integer comm,k
        integer(8) num
        complex(8) vec(:)
        character info(:)
        integer myrank, ierror

        call MPI_Comm_rank(comm, myrank, ierror)

        if (myrank .eq. 0) then
          write (*,*) info
          do k=0,num
            if (MOD(k,4) .eq. 0) then
              write (*,*) k, "."
            endif
            write (*,*) vec(k)
            if (MOD(k,4) .eq. 3) then
              write (*,*) "\n"
            endif
          enddo
        endif

      end subroutine vpr_complex

      subroutine init_random_x(lo, up, local_num_nodes, x)
        implicit none
        include "pnfft.f"

        integer t
        integer(ptrdiff_t_kind) local_num_nodes
        integer(ptrdiff_t_kind) j
        real(8) lo(:), up(:), x(:)

        do j=1,local_num_nodes
          do t=1,3
            x(3*j+t) = (up(t) - lo(t)) * rand(0) + lo(t)
          enddo
        enddo
      end subroutine init_random_x




