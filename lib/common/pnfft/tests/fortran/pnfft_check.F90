
      program main

      implicit none

      include "mpif.h"
#include "fftw3.f"
      include "pfft.f"
      include "pnfft.f"

      integer np(3), m, window, window_flag, ierror
      integer(ptrdiff_t_kind) :: N(3), n(3), M, local_M
      integer(8) plan_forw, plan_back
      integer comm_cart_3d
      complex(8), allocatable ::  f1(:), f2(:)
      real(8) f_hat_sum, x_max(3)

      N = (/ 16,16,16 /)
      n = (/ 0,0,0 /)
      local_M = 0
      m = 6
      window = 0
      x_max = (/ 0.5,0.5,0.5 /)
      np = (/ 2,2,2 /)

!     Initialize MPI and PFFT
      call MPI_Init(ierror)
      call dpnfft_init();
      
!      call init_parameters()
!      if M or n are set to zero, we choose nice values
      local_M = (local_M==0) ? N[0]*N[1]*N[2]/(np[0]*np[1]*np[2]) : local_M;
      do t=0,3
        if (n(t) .eq. 0)
          n(t) = 2*N[t]
        endif
      enddo

      select case (window)
        case (1)
          window_flag = PNFFT_WINDOW_GAUSSIAN
        case (2)
          window_flag = PNFFT_WINDOW_BSPLINE
        case (3)
          window_flag = PNFFT_WINDOW_SINC_POWER
        case (4)
          window_flag = PNFFT_WINDOW_BESSEL_I0
        case default
          window_flag = PNFFT_WINDOW_KAISER_BESSEL
      end select





      call MPI_Comm_rank(MPI_COMM_WORLD, myrank, ierror)
      
!     Create two-dimensional process grid of
!     size np(1) x np(2), if possible
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
     &     N, local_M, comm_cart_3d)

!     Get data pointers
      call dpnfft_get_f_hat(f_hat, pnfft)
      call dpnfft_get_f(f, pnfft)
      call dpnfft_get_x(x, pnfft)

!     Initialize Fourier coefficients
      call dpnfft_init_f_hat_3d(N, local_N, local_N_start, &
     &     f_hat)

!     Initialize nonequispaced nodes
      call init_random_x(lower_border, upper_border, local_M, &
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

