program main
  use, intrinsic :: iso_c_binding
  use mpi
  implicit none
  include "fftw3-mpi.f03"
  include "pfft.f03"
  include "pnfft.f03"
  interface
    subroutine compare_f(f1, f2, local_M, f_hat_sum, comm)
      use, intrinsic :: iso_c_binding
      integer(C_INTPTR_T),       intent(in) :: local_M
      complex(C_DOUBLE_COMPLEX), intent(in) :: f1(local_M), f2(local_M)
      real(C_DOUBLE),            intent(in) :: f_hat_sum
      integer, intent(in) :: comm
    end subroutine compare_f
    subroutine init_random_x(lo, up, local_M, x)
      use, intrinsic :: iso_c_binding
      integer(C_INTPTR_T), intent(in) :: local_M
      real(C_DOUBLE), intent(in) :: lo(3), up(3)
      real(C_DOUBLE), intent(out) :: x(3,local_M)
    end subroutine init_random_x
    subroutine init_f_hat(N, local_N, local_N_start, f_hat)
      use, intrinsic :: iso_c_binding
      integer(C_INTPTR_T),       intent(in) :: N(3), local_N(3), local_N_start(3)
      complex(C_DOUBLE_COMPLEX), intent(out) :: f_hat(local_N(1),local_N(3),local_N(2))
    end subroutine init_f_hat
    subroutine pnfft_perform_guru( &
        N, Nos, local_M, m, x_max, window_flag, &
        np, comm, &
        f, cf, f_hat_sum &
      )
      use, intrinsic :: iso_c_binding
      integer(C_INTPTR_T), intent(in)  :: N(3), Nos(3), local_M
      integer,             intent(in)  :: m, window_flag, np(3), comm
      real(C_DOUBLE),      intent(in)  :: x_max(3)
      complex(C_DOUBLE_COMPLEX), pointer, intent(inout) :: f(:)
      type(C_PTR),                        intent(inout) :: cf
      real(C_DOUBLE),      intent(out) :: f_hat_sum
    end subroutine pnfft_perform_guru
  end interface

  integer np(3), m, window_flag, ierror
  integer(C_INTPTR_T) :: N(3), Nos(3)
  integer(C_INTPTR_T) :: local_M
  complex(C_DOUBLE_COMPLEX), pointer :: f1(:), f2(:)
  type(C_PTR) :: cf1, cf2

  real(C_DOUBLE) f_hat_sum, x_max(3)
  integer comm_cart_3d
  integer myrank

  N   = (/ 16,16,16 /)
  Nos = (/ 32,32,32 /)
  np  = (/ 2,2,2 /)
  local_M = N(1)*N(2)*N(3) / (np(1)*np(2)*np(3))
  m = 6
  window_flag = PNFFT_WINDOW_KAISER_BESSEL
  x_max = (/ 0.5,0.5,0.5 /)

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

  if(myrank .eq. 0) then
    write(*,*) "******************************************************************************************************"
    write(*,*) "* Computation of parallel NFFT"
    write(*,*) "* for  N[0] x N[1] x N[2] = ", N(1), " x ", N(2), " x ", N(3), " Fourier coefficients"
    write(*,*) "* at   local_M = ", local_M, "nodes per process"
    write(*,*) "* with n[0] x n[1] x n[2] = ", Nos(1), " x ", Nos(2), " x ", Nos(3), " FFT grid size,"
    write(*,*) "*      m = ", m, " real space cutof,"
    write(*,*) "*      window_flag = PNFFT_WINDOW_KAISER_BESSEL,"
    write(*,*) "* on   np[0] x np[1] x np[2] = ", np(1), " x ", np(2), " x ", np(3), " processes"
    write(*,*) "*******************************************************************************************************"
  endif

  ! calculate parallel NFFT
  call srand(1)
  call pnfft_perform_guru(N, Nos, local_M, m, x_max, window_flag, np, MPI_COMM_WORLD, &
     f1, cf1, f_hat_sum)

  ! calculate parallel NFFT with higher accuracy
  call srand(1)
  call pnfft_perform_guru(N, Nos, local_M, m+2, x_max, window_flag, np, MPI_COMM_WORLD, &
     f2, cf2, f_hat_sum)

  ! calculate error of PNFFT
  call compare_f(f1, f2, local_M, f_hat_sum, MPI_COMM_WORLD)

  ! Free mem and finalize
  call pnfft_free(cf1)
  call pnfft_free(cf2)
  call pnfft_cleanup()
  call MPI_Finalize(ierror)
end program main



subroutine compare_f(f1, f2, local_M, f_hat_sum, comm)
  use, intrinsic :: iso_c_binding
  use mpi
  implicit none
  include "fftw3-mpi.f03"
  include "pfft.f03"
  include "pnfft.f03"

  integer(C_INTPTR_T),       intent(in) :: local_M
  complex(C_DOUBLE_COMPLEX), intent(in) :: f1(local_M), f2(local_M)
  real(C_DOUBLE),            intent(in) :: f_hat_sum
  integer(C_INTPTR_T) j
  integer, intent(in) :: comm
  double precision local_error, global_error
  integer myrank, ierror

  call MPI_Comm_rank(comm, myrank, ierror)

  local_error = 0
  do j=1,local_M
    if ( cdabs(f1(j) - f2(j)) > local_error ) then
      local_error = cdabs(f1(j) - f2(j))
    endif
  enddo
  call MPI_Allreduce(local_error, global_error, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierror)

  ! Print output Fourier coefficents
  if(myrank .eq. 0) then
    write (*,*) ""
    write (*,*) "absolute error = ", global_error
    write (*,*) "relative error = ", global_error/f_hat_sum
  endif
end subroutine compare_f

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

double precision function func(l,N)
  use, intrinsic :: iso_c_binding
  implicit none

  integer(C_INTPTR_T), intent(in) :: l, N

  func = exp(-real(l*l)/real(N*N))
end function func

subroutine init_f_hat(N, local_N, local_N_start, f_hat)
  use, intrinsic :: iso_c_binding
  implicit none
  interface
    double precision function func(l,N)
      use, intrinsic :: iso_c_binding
      integer(C_INTPTR_T), intent(in) :: l, N
    end function func
  end interface

  integer(C_INTPTR_T),       intent(in) :: N(3), local_N(3), local_N_start(3)
  ! transposed FFT output AND C-like memory order
  complex(C_DOUBLE_COMPLEX), intent(out) :: f_hat(local_N(1),local_N(3),local_N(2))
  integer(C_INTPTR_T) :: l1, l2, l3

  ! use C-like row-major order here  
  do l2 = 1,local_N(2)
    do l3 = 1,local_N(3)
      do l1 = 1,local_N(1)
        f_hat(l1,l3,l2) = &
            func(l1 + local_N_start(1), N(1)) & 
          * func(l2 + local_N_start(2), N(2)) &
          * func(l3 + local_N_start(3), N(3)) 
      enddo
    enddo
  enddo
end subroutine init_f_hat

subroutine pnfft_perform_guru( &
    N, Nos, local_M, m, x_max, window_flag, &
    np, comm, &
    f, cf, f_hat_sum &
  )
  use, intrinsic :: iso_c_binding
  use mpi
  implicit none
  include "fftw3-mpi.f03"
  include "pfft.f03"
  include "pnfft.f03"

  integer(C_INTPTR_T), intent(in)  :: N(3), Nos(3), local_M
  integer,             intent(in)  :: m, window_flag, np(3), comm
  real(C_DOUBLE),      intent(in)  :: x_max(3)
  complex(C_DOUBLE_COMPLEX), pointer, intent(inout) :: f(:)
  type(C_PTR),                        intent(inout) :: cf
  real(C_DOUBLE),      intent(out) :: f_hat_sum

  integer(C_INTPTR_T) :: local_N(3), local_N_start(3), k1, k2, k3, d=3
  complex(C_DOUBLE_COMPLEX), pointer :: f_hat(:,:,:)
  real(C_DOUBLE), pointer :: x(:,:)
  real(C_DOUBLE) :: lower_border(3), upper_border(3)
  real(C_DOUBLE) :: local_f_hat_sum
  type(C_PTR) :: pnfft, cf_hat, cx
  integer ierror, myrank, comm_cart_3d
  double precision time, max_time

  call MPI_Comm_rank(comm, myrank, ierror)

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
      PNFFT_TRANSPOSED_F_HAT, &
      local_N, local_N_start, lower_border, upper_border)

  ! Plan parallel NFFT
  pnfft = pnfft_init_guru(3, N, Nos, x_max, local_M, m, &
      PNFFT_MALLOC_X + PNFFT_MALLOC_F_HAT + PNFFT_MALLOC_F + &
      PNFFT_TRANSPOSED_F_HAT + window_flag, &
      PFFT_ESTIMATE, comm_cart_3d)

  ! Get data pointers in C format
  cf_hat = pnfft_get_f_hat(pnfft)
  cf     = pnfft_get_f(pnfft)
  cx     = pnfft_get_x(pnfft)

  ! Convert data pointers to Fortran format with transposed f_hat
  call c_f_pointer(cf_hat, f_hat, [local_N(1),local_N(3),local_N(2)])
  call c_f_pointer(cf,     f,     [local_M])
  call c_f_pointer(cx,     x,     [d,local_M])

  ! Initialize Fourier coefficients
  call init_f_hat(N, local_N, local_N_start, &
      f_hat)

  ! Initialize nonequispaced nodes
  call init_random_x(lower_border, upper_border, local_M, &
     x)

  ! Execute parallel NFFT
  time = -MPI_Wtime()
  call pnfft_trafo(pnfft)
  time = time + MPI_Wtime()

  call MPI_Reduce(time, max_time, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, comm_cart_3d, ierror)
  if(myrank .eq. 0) then
    write(*,*) "pnfft_trafo needs ", max_time, " s"
  endif

  local_f_hat_sum = 0d0;
  do k2=1,local_N(2)
    do k3=1,local_N(3)
      do k1=1,local_N(1)
        local_f_hat_sum = local_f_hat_sum + cdabs(f_hat(k1,k3,k2))
      enddo
    enddo
  enddo
  call MPI_Allreduce(local_f_hat_sum, f_hat_sum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm_cart_3d, ierror)

  call pnfft_finalize(pnfft, PNFFT_FREE_X + PNFFT_FREE_F_HAT)
  call MPI_Comm_free(comm_cart_3d, ierror)
end subroutine pnfft_perform_guru
