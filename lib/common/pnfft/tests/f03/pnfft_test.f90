
program main

  use, intrinsic :: iso_c_binding
  implicit none
  include "mpif.h"
  include "fftw3-mpi.f03"
  include "pfft.f03"
  include "pnfft.f03"

  integer np(3), ierror
  integer(C_INTPTR_T) l1, l2, l3
  integer(C_INTPTR_T) :: N(3), local_N(3), local_N_start(3)
  integer(C_INTPTR_T) :: local_M, d=3
  real(C_DOUBLE) :: lower_border(3), upper_border(3)
  type(C_PTR) :: pnfft, cf_hat, cf, cx
  complex(C_DOUBLE_COMPLEX), pointer :: f_hat(:,:,:), f(:)
  real(C_DOUBLE), pointer :: x(:,:)

  integer comm_cart_3d, myrank

  N = (/ 16,16,16 /)
  np = (/ 2,2,2 /)
  local_M = N(1)*N(2)*N(3) / (np(1)*np(2)*np(3))

! Initialize MPI and PFFT
  call MPI_Init(ierror)
  call pnfft_init();
  call MPI_Comm_rank(MPI_COMM_WORLD, myrank, ierror)
  
! Create three-dimensional process grid of
! size np(1) x np(2) x np(3), if possible
  ierror =  pnfft_create_procmesh(3, MPI_COMM_WORLD, np, comm_cart_3d)
  if (ierror .ne. 0) then
    if(myrank .eq. 0) then
      write(*,*) "Error: This test file only works with ", np(1)*np(2)*np(3), " processes"
    endif
    call MPI_Finalize(ierror)
    call exit(1)
  endif

! Get parameters of data distribution
  call pnfft_local_size_3d(N, comm_cart_3d, PNFFT_TRANSPOSED_NONE, &
      local_N, local_N_start, lower_border, upper_border)

! Plan parallel NFFT
  pnfft = pnfft_init_3d(N, local_M, comm_cart_3d)

! Get data pointers in C format
  cf_hat = pnfft_get_f_hat(pnfft)
  cf     = pnfft_get_f(pnfft)
  cx     = pnfft_get_x(pnfft)

! Convert data pointers to Fortran format
  call c_f_pointer(cf_hat, f_hat, [local_N])
  call c_f_pointer(cf,     f,     [local_M])
  call c_f_pointer(cx,     x,     [d,local_M])

! Initialize Fourier coefficients
  call init_f_hat(N, local_N, local_N_start, &
      f_hat)

! Initialize nonequispaced nodes
  call init_random_x(lower_border, upper_border, local_M, &
     x)

! Print input Fourier coefficents
  if(myrank .eq. 0) then
    write (*,*) "Input Fourier coefficients on process 1:"
    write (*,*) f_hat(1:4,1,1)
  endif

! Execute parallel NFFT
  call pnfft_trafo(pnfft)

! Print NFFT results
  if(myrank .eq. 0) then
    write (*,*) ""
    write (*,*) "PNFFT Results on process 1:"
    write (*,*) f(1:4)
  endif

! Execute parallel adjoint NFFT
  call pnfft_adj(pnfft)

! Scale data
  do l3=1,local_N(3)
    do l2=1,local_N(2)
      do l1=1,local_N(1)
        f_hat(l1,l2,l3) = f_hat(l1,l2,l3) / (N(3)*N(2)*N(1))
      enddo
    enddo
  enddo

! Print output Fourier coefficents
  if(myrank .eq. 0) then
    write (*,*) ""
    write (*,*) "Fourier coefficients after one forward and backward PNFFT on process 1:"
    write (*,*) f_hat(1:4,1,1)
  endif

! Free mem and finalize      
  call pnfft_finalize(pnfft, &
      PNFFT_FREE_X + PNFFT_FREE_F_HAT + PNFFT_FREE_F)
  call MPI_Comm_free(comm_cart_3d, ierror)
  call pnfft_cleanup()
  call MPI_Finalize(ierror)
end program main


subroutine init_random_x(lo, up, local_M, x)
  use, intrinsic :: iso_c_binding
  implicit none

  integer(C_INTPTR_T) j, t
  integer(C_INTPTR_T), intent(in) :: local_M
  real(C_DOUBLE), intent(in) :: lo(3), up(3)
  real(C_DOUBLE), intent(out) :: x(3,local_M)

  do j=1,local_M
    do t=1,3
      x(t,j) = (up(t) - lo(t)) * rand(0) + lo(t)
    enddo
  enddo
end subroutine init_random_x

subroutine init_f_hat(N, local_N, local_N_start, f_hat)
  use, intrinsic :: iso_c_binding
  implicit none

  integer(C_INTPTR_T) :: l1, l2, l3
  integer(C_INTPTR_T), intent(in) :: N(3), local_N(3), local_N_start(3)
  complex(C_DOUBLE_COMPLEX), intent(out) :: f_hat(local_N(1),local_N(2),local_N(3))
  real :: func

! use C-like row-major order here  
  do l1 = 1,local_N(1)
    do l2 = 1,local_N(2)
      do l3 = 1,local_N(3)
        f_hat(l1,l2,l3) = &
            func(l1 + local_N_start(1), N(1)) & 
          * func(l2 + local_N_start(2), N(2)) &
          * func(l3 + local_N_start(3), N(3)) 
      enddo
    enddo
  enddo
end subroutine init_f_hat

real function func(l,N)
  use, intrinsic :: iso_c_binding
  implicit none

  integer(C_INTPTR_T), intent(in) :: l, N

  func = exp(-real(l*l)/real(N*N))
end function func

