#ifdef HAVE_FCONFIG_H
#include <fconfig.h>
#endif

program test_pepc
  use, intrinsic :: iso_c_binding

#ifdef NEED_USE_IFPORT
  use ifport
#endif
  use FCSError_fh
  use FCSInput_fh
  use FCSOutput_fh
  use FCSParameter_fh
  use fcs_pepc_fh

  implicit none

 include 'mpif.h'
!  include 'mpi.h'

  ! Parameters
  double precision :: pi = 3.141592653589793D0

   ! MPI variables
  integer :: mpi_err
  integer :: mpi_size

  integer :: mpi_rank
 
  double precision :: starttime, endtime

  ! Argument counter
  integer :: argc

  ! Timer variables
  double precision :: t_start, t_stop
  character(len=80) :: timer_fmt

  ! Particle counter and numbers
  integer :: p, n_particles, n_local_particles

  ! Particle data
  real(kind=c_double), dimension(:), allocatable :: fx, fy, fz
  real(kind=c_double), dimension(:), allocatable, target :: q

  ! Total energy
  double precision :: e_sum, e_sum_local

  ! Total charge
  double precision :: q_sum, q_sum_local

  ! Sum of forces
  double precision, dimension(3) :: f_sum_local, f_sum

  ! Maximum absolute values of components
  double precision :: f_max_abs_x_local, f_max_abs_y_local, f_max_abs_z_local
  double precision :: f_max_abs_x, f_max_abs_y, f_max_abs_z

  ! Sum of squares of components of force
  double precision :: f_sum_squared_x_components, f_sum_squared_x_components_local
  double precision :: f_sum_squared_y_components, f_sum_squared_y_components_local
  double precision :: f_sum_squared_z_components, f_sum_squared_z_components_local

  ! Variables for reading data
  double precision :: my_x, my_y, my_z, my_q

  ! Variables for fcs-interface written in C
  type( c_ptr ) :: param_handle
  type( c_ptr ) :: props_handle
  type( c_ptr ) :: error_handle
  type( c_ptr ) :: input_handle
  type( c_ptr ) :: output_handle

  !dummies for correkt initalization of the FCSParameter-object
  real(kind = c_double), dimension(3), target :: bl
  character(kind = c_char), dimension(3), target :: pf 
  type(c_ptr) :: box_length
  type(c_ptr) :: periodic_flags

  type( c_ptr ) :: f
  type( c_ptr ) :: pot
  type( c_ptr ) :: xyz
  type( c_ptr ) :: charges
  real(kind=c_double), dimension(:), target, allocatable :: positions
  real(kind=c_double), dimension(:), pointer :: forces => null()
  real(kind=c_double), dimension(:), pointer :: e => null()
  character( kind=c_char, len = 300) :: parameterstring
  character(kind=c_char), dimension(:), pointer :: ptr => null()
  type( c_ptr ) :: test_c_ptr
  integer(kind=c_int) :: system_dim
  integer(kind=c_int) :: n_local_particles_c
  integer(kind=c_int) :: ret_val
  integer :: DIM, k, str_length

  ! optional parameter
  !real(kind = c_double) :: theta, eps
  double precision :: theta, eps

  ! Format for timer output
  timer_fmt = "(A, T50, F15.5, A)"

  DIM = 3
  system_dim = DIM
  ! Initializing MPI
  call MPI_Init(mpi_err)
  call MPI_Comm_size(MPI_COMM_WORLD, mpi_size, mpi_err)
  call MPI_Comm_rank(MPI_COMM_WORLD, mpi_rank, mpi_err)


  bl(1:3) = 1.0
  pf(1:3) = "111"
  box_length = c_loc(bl)
  periodic_flags = c_loc(pf)
  
  n_local_particles = 20
  n_particles = mpi_size * n_local_particles

  ! create FCSParameter-object for fcsinit-function
  error_handle = fcsParameter_create(param_handle, n_particles, n_local_particles, n_local_particles, &
       system_dim, box_length, periodic_flags, MPI_COMM_WORLD, "1");
  ! check status
  if(fcsError_getErrorCode(error_handle) /= FCS_SUCCESS) then
     call fcsError_printError(error_handle)
     write(*, *)'exit.....'
     call MPI_Finalize(mpi_err)
     call EXIT(-1)
  end if

  error_handle = fcsinit_pepc(param_handle, props_handle)
  ! check status
  if(fcsError_getErrorCode(error_handle) /= FCS_SUCCESS) then
     call fcsError_printError(error_handle)
     write(*, *)'exit.....'
     call MPI_Finalize(mpi_err)
     call EXIT(-1)
  end if

  
  ! allocate data structures 
  allocate(q(1:n_local_particles))

  allocate(fx(1:n_local_particles), fy(1:n_local_particles), fz(1:n_local_particles))

  ! allocate and initialize positions-array for fcs interface
  allocate(positions(1:(n_local_particles*DIM)))
  call srand(mpi_size+10)
  p=1
  do k=1, n_local_particles
     positions(p) = rand()
     positions(p+1) = rand()
     positions(p+2) = rand()
     p = p + 3      
  end do

  ! Initialize forces and charges
  do p=1,n_local_particles
     fx(p) = 0.0D0
     fy(p) = 0.0D0
     fz(p) = 0.0D0
     q(p) = 1e-2
  end do

  ! get addresses of the positions- and charges-array
  xyz = c_loc(positions)
  charges = c_loc(q)

  ! create input structure for fcs interface
  error_handle = fcsInput_create(input_handle, n_particles, n_local_particles, &
       & xyz, system_dim, charges)

  ! check status
  if(fcsError_getErrorCode(error_handle) /= FCS_SUCCESS) then
     call fcsError_printError(error_handle)
     write(*, *)'exit.....'
     call MPI_Finalize(mpi_err)
     call EXIT(-1)
  end if

  ! create output structure for fcs interface
  error_handle = fcsOutput_create(output_handle, n_local_particles, system_dim)

  ! check status
  if(fcsError_getErrorCode(error_handle) /= FCS_SUCCESS) then
     call fcsError_printError(error_handle)
     call MPI_Finalize(mpi_err)
     call EXIT(-1)
  end if

  ! create a parameter string
  parameterstring = pepc_THETA // ',' // pepc_EPS // c_null_char

  theta = 0.65D0
  eps = 1.31e-6

  call MPI_Barrier(MPI_COMM_WORLD,mpi_err)
  if (mpi_rank == 0) then
     starttime = MPI_Wtime()
  end if

  error_handle = call_pepc_with_opt_param(props_handle, &
       input_handle,&
       output_handle, &
       parameterstring, &
       theta, eps)

  ! check status
  if(fcsError_getErrorCode(error_handle) /= pepc_OK) then
     call fcsError_printError(error_handle)
     call MPI_Finalize(mpi_err)
     call EXIT(-1)
  end if

  ! get calculated forces
  f = fcsOutput_getForces(output_handle)
  call c_f_pointer( f, forces, [n_local_particles*DIM] )

  ! get calculated potentials
  pot = fcsOutput_getPotentials(output_handle)
  call c_f_pointer( pot, e, [n_local_particles] )

  call MPI_Barrier(MPI_COMM_WORLD,mpi_err)
  if (mpi_rank == 0) then
     endtime = MPI_Wtime()
     write(*, *) "processors: ", mpi_size
     write(*,timer_fmt) "pepc runtime: ",endtime - starttime," s"
  end if

  ! save forces to the local data structures for a more convenient usage
  p = 1
  do k=1, n_local_particles
     fx(k) = forces(p)
     fy(k) = forces(p+1)
     fz(k) = forces(p+2)
     p = p+3
  end do

  !----------------------------------------------------------------------------
  !
  ! Energy calculation
  !

  e_sum_local = 0.0D0
  do p=1,n_local_particles
     e_sum_local = e_sum_local + e(p)
  end do

  e_sum = 0.0D0
  call MPI_Reduce(e_sum_local,e_sum,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
       & MPI_COMM_WORLD,mpi_err)

  if (mpi_rank == 0) then
!!$     write(*,"(/, A, T60, ES15.8)") &
!!$          & "Self energy:                                  ", &
!!$          & 1.0D0/(4.0D0 * pi) * 14.0D0/(5.0D0*1/(2*ghosts*min(m,n,o)))
     write(*,"(A, T60, ES15.8)") &
          & "Total energy:                                 ",e_sum
     write(*,"(A, T60, ES15.8)") &
          & "Approx. Madelung's constant:                  ", &
          & e_sum / n_particles * 2.0D0 * pi
     write(*,"(A, T60, ES15.8)") &
          & "Relative error of Approx. Madelung's constant:", &
          & (-1.74756459463318219063621D0 - &
          & e_sum / n_particles * 2.0D0 * pi) / &
          & (-1.74756459463318219063621D0)
  end if


  ! Sums and maxima of forces
  f_sum_local = 0.0D0
  f_sum_squared_x_components_local = 0.0D0
  f_sum_squared_y_components_local = 0.0D0
  f_sum_squared_z_components_local = 0.0D0
  f_max_abs_x_local = abs(f_sum_local(1))
  f_max_abs_y_local = abs(f_sum_local(2))
  f_max_abs_z_local = abs(f_sum_local(3))
  do p=1,n_local_particles
     f_sum_local(1) = f_sum_local(1) + fx(p)
     f_sum_local(2) = f_sum_local(2) + fy(p)
     f_sum_local(3) = f_sum_local(3) + fz(p)
     f_sum_squared_x_components_local = f_sum_squared_x_components_local + &
          & fx(p)**2
     f_sum_squared_y_components_local = f_sum_squared_y_components_local + &
          & fy(p)**2
     f_sum_squared_z_components_local = f_sum_squared_z_components_local + &
          & fz(p)**2
     if (abs(fx(p)) > f_max_abs_x_local) then
        f_max_abs_x_local = abs(fx(p))
     end if
     if (abs(fy(p)) > f_max_abs_x_local) then
        f_max_abs_y_local = abs(fy(p))
     end if
     if (abs(fz(p)) > f_max_abs_x_local) then
        f_max_abs_z_local = abs(fz(p))
     end if
  end do
  f_sum = 0.0D0
  call MPI_Reduce(f_sum_local,f_sum,3,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,mpi_err)
  f_sum_squared_x_components = 0.0D0
  call MPI_Reduce(f_sum_squared_x_components_local,f_sum_squared_x_components, &
       & 1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,mpi_err)
  f_sum_squared_y_components = 0.0D0
  call MPI_Reduce(f_sum_squared_y_components_local,f_sum_squared_y_components, &
       & 1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,mpi_err)
  f_sum_squared_z_components = 0.0D0
  call MPI_Reduce(f_sum_squared_z_components_local,f_sum_squared_z_components, &
       & 1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,mpi_err)
  f_max_abs_x = f_max_abs_x_local
  call MPI_Reduce(f_max_abs_x_local,f_max_abs_x,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,mpi_err)
  f_max_abs_y = f_max_abs_y_local
  call MPI_Reduce(f_max_abs_y_local,f_max_abs_y,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,mpi_err)
  f_max_abs_z_local = f_max_abs_z_local
  call MPI_Reduce(f_max_abs_z_local,f_max_abs_z,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,mpi_err)

  if (mpi_rank == 0) then
     write(*,"(A, T60, ES15.8)") 'Norm of sum of forces:', &
          & sqrt(f_sum(1)**2 + f_sum(2)**2 + f_sum(3)**2)
     write(*,"(A, T60, ES15.8)") 'Sqrt. of sum of squares of x-components of forces:', &
          & sqrt(f_sum_squared_x_components)
     write(*,"(A, T60, ES15.8)") 'Sqrt. of sum of squares of y-components of forces:', &
          & sqrt(f_sum_squared_y_components)
     write(*,"(A, T60, ES15.8)") 'Sqrt. of sum of squares of z-components of forces:', &
          & sqrt(f_sum_squared_z_components)
     write(*,"(A, T60, ES15.8)") &
          & 'Sqrt. of sum of squares of all components of forces:', &
          & sqrt(f_sum_squared_x_components + &
          & f_sum_squared_x_components + f_sum_squared_x_components)
     write(*,"(A, T60, ES15.8)") 'Max. x-component of forces:', &
          & f_max_abs_x
     write(*,"(A, T60, ES15.8)") 'Max. y-component of forces:', &
          & f_max_abs_y
     write(*,"(A, T60, ES15.8)") 'Max. z-component of forces:', &
          & f_max_abs_z
  end if


  error_handle = fcsfree_pepc()
  ! check status
  if(fcsError_getErrorCode(error_handle) /= pepc_OK) then
     call fcsError_printError(error_handle)
     call MPI_Finalize(mpi_err)
     call EXIT(-1)
  end if
  

  if (n_local_particles > 0) then
     deallocate(q)
     deallocate(fx, fy, fz)
     deallocate(positions)

     error_handle = fcsInput_destroy(input_handle)
     error_handle = fcsOutput_destroy(output_handle)
     ret_val = fcsError_destroy(error_handle)
  end if

  ! Finalizing MPI
  call MPI_Finalize(mpi_err)

  stop

end program test_pepc


