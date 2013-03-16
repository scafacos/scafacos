
program main

  use, intrinsic :: iso_c_binding

  implicit none
  include "mpif.h"
  include "fftw3.f"
  include "pfft.f"
  include "pnfft.f03"

  integer np(3), m, window, window_flag, ierror, l1, l2, l3
  integer(C_INTPTR_T) :: N(3), local_N(3), local_N_start(3)
  integer(C_INTPTR_T) :: num_nodes, local_num_nodes
  real(C_DOUBLE) :: lower_border(3), upper_border(3)
  type(C_PTR) :: pnfft, cf_hat, cf, cx
  complex(C_DOUBLE_COMPLEX), pointer :: f_hat(:,:,:), f(:)
  real(C_DOUBLE), pointer :: x(:,:)

  real(C_DOUBLE) f_hat_sum, x_max(3)
  integer comm_cart_3d
  integer myrank

  N = (/ 16,16,16 /)
  np = (/ 2,2,2 /)
  local_num_nodes = N(1)*N(2)*N(3) / (np(1)*np(2)*np(3))

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
  call pnfft_local_size_3d(N, comm_cart_3d, &
      local_N, local_N_start, lower_border, upper_border)

! Plan parallel NFFT
  pnfft = pnfft_init_3d(N, local_num_nodes, comm_cart_3d)

! Get data pointers in C format
  cf_hat = pnfft_get_f_hat(pnfft)
  cf     = pnfft_get_f(pnfft)
  cx     = pnfft_get_x(pnfft)

! Convert data pointers to Fortran format
  call c_f_pointer(cf_hat, f_hat, [local_N])
  call c_f_pointer(cf,     f,     [local_num_nodes])
  call c_f_pointer(cx,     x,     [integer(C_INTPTR_T)::3,local_num_nodes])

! Initialize Fourier coefficients
  call init_f_hat(N, local_N, local_N_start, &
      f_hat)

! Initialize nonequispaced nodes
  call init_random_x(lower_border, upper_border, local_num_nodes, &
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


subroutine init_random_x(lo, up, local_num_nodes, x)
  use, intrinsic :: iso_c_binding
  implicit none

  integer(C_INTPTR_T) j, t
  integer(C_INTPTR_T), intent(in) :: local_num_nodes
  real(C_DOUBLE), intent(in) :: lo(3), up(3)
  real(C_DOUBLE), intent(out) :: x(3,local_num_nodes)

  do j=1,local_num_nodes
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
  complex(C_DOUBLE), intent(out) :: f_hat(local_N(1),local_N(2),local_N(3))
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


subroutine pnfft_perform_guru( &
    N, Nos, local_M, m, x_max, window_flag, &
    np, comm, &
    f, f_hat_sum &
    )

  use, intrinsic :: iso_c_binding

  implicit none
  include "mpif.h"
  include "fftw3.f"
  include "pfft.f"
  include "pnfft.f03"

  integer(C_INTPTR_T), intent(in)  :: N(3), Nos(3), local_M
  integer,             intent(in)  :: m, window_flag, np(3), comm
  real(C_DOUBLE),      intent(in)  :: x_max(3)
  complex(C_DOUBLE),   intent(out) :: f(local_M), f_hat_sum

  integer ierr, myrank, comm_cart_3d
  complex(C_DOUBLE_COMPLEX), pointer :: f_hat(:,:,:), f(:)
  real(C_DOUBLE), pointer :: x(:,:)
  type(C_PTR) :: pnfft, cf_hat, cf, cx

  call MPI_Comm_rank(comm, myrank, ierr)

! Create three-dimensional process grid of
! size np(1) x np(2) x np(3), if possible
  ierror =  pnfft_create_procmesh(3, comm, np, comm_cart_3d)
  if (ierror .ne. 0) then
    if(myrank .eq. 0) then
      write(*,*) "Error: This test file only works with ", np(1)*np(2)*np(3), " processes"
    endif
    call MPI_Finalize(ierror)
    call exit(1)
  endif

! Get parameters of data distribution
  call pnfft_local_size_guru(3, N, Nos, x_max, m, comm_cart_3d, &
      local_N, local_N_start, lower_border, upper_border)

! Plan parallel NFFT
  pnfft = pnfft_init_guru(3, N, Nos, x_max, local_M, m, &
      PNFFT_MALLOC_X + PNFFT_MALLOC_F_HAT + PNFFT_MALLOC_F + window_flag, &
      PFFT_ESTIMATE, comm_cart_3d)

! Get data pointers in C format
  cf_hat = pnfft_get_f_hat(pnfft)
  cf     = pnfft_get_f(pnfft)
  cx     = pnfft_get_x(pnfft)

! Convert data pointers to Fortran format
  call c_f_pointer(cf_hat, f_hat, [local_N])
  call c_f_pointer(cf,     f,     [local_M])
  call c_f_pointer(cx,     x,     [integer(C_INTPTR_T)::3,local_M])



end subroutine pnfft_perform_guru


