#ifdef HAVE_FCONFIG_H
#include <fconfig.h>
#endif

program test_pp3mg
  use iso_c_binding
 
  use FCSError_fh
  use FCSInput_fh
  use FCSOutput_fh
  use FCSParameter_fh
  use FCSMethodProperties_fh
  use fcs_pp3mg_fh

  implicit none

  include 'mpif.h'

  
  ! Parameters
  double precision :: pi = 3.141592653589793D0

  integer :: m, n, o
  integer :: maxiter
  integer :: nu1, nu2
  double precision :: tol
  integer :: ghosts
  integer :: max_particles
  character(len=150) :: fname = "particles.txt"

  ! MPI variables
  integer :: mpi_err
  integer :: mpi_size

  integer :: mpi_rank
  integer :: mpi_comm_cart
  integer, dimension(1:3), target :: mpi_dims
  logical, dimension(1:3) :: mpi_periods
  integer, dimension(1:3) :: mpi_coords

  double precision :: starttime, endtime

  ! Argument counter
  integer :: argc

  ! Timer variables
  double precision :: t_start, t_stop
  character(len=80) :: timer_fmt

  ! Particle counter and numbers
  integer :: p, p_local, n_particles, n_local_particles

  ! Particle data
  real(kind=c_double), dimension(:), allocatable :: x, y, z, fx, fy, fz
  real(kind=c_double), dimension(:), allocatable, target :: q
 
  ! Size of local domain
  double precision :: x_start, y_start, z_start
  double precision :: x_end, y_end, z_end

  ! Degree of interpolation polynomial
  integer :: degree
  
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
  real(kind = c_double), dimension(3), target :: bl
  character(kind = c_char), dimension(3), target :: pf 

  ! Format for timer output
  timer_fmt = "(A, T50, F15.5, A)"


  DIM = 3
  system_dim = DIM
  ! Initializing MPI
  call MPI_Init(mpi_err)
  call MPI_Comm_size(MPI_COMM_WORLD, mpi_size, mpi_err)
  call MPI_Comm_rank(MPI_COMM_WORLD, mpi_rank, mpi_err)
  
  ! Get and distribute input arguments
  if (mpi_rank == 0) then
     argc = command_argument_count()
     if (argc >= 1) then
        call get_command_argument(1,fname)
        open(unit=99,file=fname)
        read(unit=99,fmt='(/)')
        read(unit=99,fmt=*) m, n, o
        read(unit=99,fmt='(/)')
        read(unit=99,fmt=*) maxiter
        read(unit=99,fmt='(/)')
        read(unit=99,fmt=*) nu1, nu2
        read(unit=99,fmt='(/)')
        read(unit=99,fmt=*) tol
        read(unit=99,fmt='(/)')
        read(unit=99,fmt=*) ghosts
        read(unit=99,fmt='(/)')
        read(unit=99,fmt=*) max_particles
        read(unit=99,fmt='(/)')
       ! read(unit=99,fmt='(2X)')
        read(unit=99,fmt='(2X,A)') fname
        close(unit=99)
     end if
     write(*,"(/, A)") "Active parameters:"
     write(*,"(A)") "------------------"
     write(*,"(A,I15)") "m              = ",m
     write(*,"(A,I15)") "n              = ",n
     write(*,"(A,I15)") "o              = ",o
     write(*,"(A,I15)") "maxiter        = ",maxiter
     write(*,"(A,I15)") "nu1            = ",nu1
     write(*,"(A,I15)") "nu2            = ",nu2
     write(*,"(A,ES15.5)") "tol            = ",tol
     write(*,"(A,I15)") "ghosts         = ",ghosts
     write(*,"(A,I15)") "max_particles  = ",max_particles
  end if
  call MPI_Bcast(m,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpi_err)
  call MPI_Bcast(n,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpi_err)
  call MPI_Bcast(o,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpi_err)
  call MPI_Bcast(maxiter,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpi_err)
!!$  call MPI_Bcast(nu1,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpi_err)
!!$  call MPI_Bcast(nu2,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpi_err)
!!$  call MPI_Bcast(tol,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpi_err)
  call MPI_Bcast(ghosts,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpi_err)
  call MPI_Bcast(max_particles,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpi_err)
  call MPI_Bcast(fname,80,MPI_CHARACTER,0,MPI_COMM_WORLD,mpi_err)

  degree = 4
  ghosts = max(ghosts,ceiling(0.5 * dble(degree)))

  ! Initialize cartesian process grid
    mpi_dims = 0
    mpi_periods = .TRUE.
    call MPI_Comm_size(MPI_COMM_WORLD, mpi_size, mpi_err)
    call MPI_Dims_create(mpi_size, 3, mpi_dims, mpi_err)
    call MPI_Cart_create(MPI_COMM_WORLD, 3, mpi_dims, mpi_periods, .true., &
         & mpi_comm_cart, mpi_err)
    call MPI_Comm_rank(mpi_comm_cart, mpi_rank, mpi_err)
    call MPI_Cart_coords(mpi_comm_cart, mpi_rank, 3, mpi_coords, mpi_err)

    ! Initialize parameters for particle grid
    x_start = dble(1.0)/dble(mpi_dims(1)) &
         * dble(mpi_coords(1))
    x_end = dble(1.0)/dble(mpi_dims(1)) &
         * dble(mpi_coords(1)+1)
    y_start = dble(1.0)/dble(mpi_dims(2)) &
         * dble(mpi_coords(2))
    y_end = dble(1.0)/dble(mpi_dims(2)) &
         * dble(mpi_coords(2)+1)
    z_start = dble(1.0)/dble(mpi_dims(3)) &
         * dble(mpi_coords(3))
    z_end = dble(1.0)/dble(mpi_dims(3)) &
         * dble(mpi_coords(3)+1)
    
    ! count local particles 
    open (unit = 99, file = fname)
    read(99,*) n_particles
    n_local_particles = 0
    do p=1, n_particles
       read(99,*) my_x, my_y, my_z, my_q
       if (x_start<=my_x .and. my_x<x_end .and. &
            & y_start<=my_y .and. my_y<y_end .and. &
            & z_start<=my_z .and. my_z<z_end) then
          n_local_particles = n_local_particles + 1
       end if
    end do
    close(99)

    ! allocate data structures 
    n_local_particles_c = n_local_particles
    allocate(x(1:n_local_particles), y(1:n_local_particles), &
         & z(1:n_local_particles), q(1:n_local_particles))
    
    ! read particles coordinates from the input file
    open (unit = 99, file = fname)
    read(99,*) n_particles
    p_local = 0
    do p=1, n_particles
       read(99,*) my_x, my_y, my_z, my_q
       if (x_start<=my_x .and. my_x<x_end .and. &
            & y_start<=my_y .and. my_y<y_end .and. &
            & z_start<=my_z .and. my_z<z_end) then
          p_local = p_local + 1
          x(p_local) = my_x
          y(p_local) = my_y
          z(p_local) = my_z
          q(p_local) = my_q
       end if
    end do
    close(99)
    
    allocate(fx(1:n_local_particles), fy(1:n_local_particles), fz(1:n_local_particles))
  
    ! Initialize forces
    do p=1,n_local_particles
       fx(p) = 0.0D0
       fy(p) = 0.0D0
       fz(p) = 0.0D0
    end do

    ! allocate and initialize positions-array for fcs interface
    allocate(positions(1:(n_local_particles*DIM)))
    p=1
    do k=1, n_local_particles
       positions(p) = x(k)
       positions(p+1) = y(k)
       positions(p+2) = z(k)
       p = p + 3      
    end do

    bl(1:3) = 1.0
    pf(1:3) = "111"
    box_length = c_loc(bl)
    periodic_flags = c_loc(pf)


    ! create FCSParameter-object for fcsinit-function
    error_handle = fcsParameter_create(param_handle, n_particles, n_local_particles_c, n_local_particles_c, &
         system_dim, box_length, periodic_flags, mpi_comm_cart, "1");
    ! check status
    if(fcsError_getErrorCode(error_handle) /= FCS_SUCCESS) then
       write(*, *)'exit.....'
       call MPI_Finalize(mpi_err)
       call EXIT(-1)
    end if

    ! create a parameter string
    parameterstring = pp3mg_MPI_DIMS // ',' // pp3mg_MAX_PARTICLES// ',' // pp3mg_CELLS_X // ','// pp3mg_CELLS_Y //',' // &
                      & pp3mg_CELLS_Z //','// pp3mg_GHOST_CELLS // ','// pp3mg_ERR_BOUND // c_null_char

    call MPI_Barrier(mpi_comm_cart,mpi_err)
    if (mpi_rank == 0) then
       starttime = MPI_Wtime()
    end if

    ! 1. step: initialize the pp3mg-method. The init-method of the pp3mg
    ! requires at least a cartesian MPI communicator as an input parameter!
    error_handle = call_fcsinit_pp3mg_with_opt_param(param_handle, props_handle,&
         parameterstring, mpi_dims=c_loc(mpi_dims), max_particles=max_particles,&
         cells_x=m, cells_y=n, cells_z=o, ghost_cells=ghosts, err_bound=0.00001D0)
  
    ! check status
    if(fcsError_getErrorCode(error_handle) /= pp3mg_OK) then
       call fcsError_printError(error_handle)
       call MPI_Finalize(mpi_err)
       call EXIT(-1)
    end if


    ! get addresses of the positions- and charges-array
    xyz = c_loc(positions)
    charges = c_loc(q)
    
    ! create input structure for fcs interface
    error_handle = fcsInput_create(input_handle, n_particles, n_local_particles_c, &
         & xyz, system_dim, charges)

    ! check status
    if(fcsError_getErrorCode(error_handle) /= FCS_SUCCESS) then
       call fcsError_printError(error_handle)
       write(*, *)'exit.....'
       call MPI_Finalize(mpi_err)
       call EXIT(-1)
    end if

     !*** example of the output of c-strings ***
!!$     test_c_ptr = fcsError_getMessage(error_handle)
!!$     if(c_associated(test_c_ptr)) then
!!$        str_length = strlen(test_c_ptr)
!!$        call c_f_pointer(test_c_ptr, ptr, [m])
!!$        write(*, *) "fcsError_getErrorMessage: ", ptr
!!$     end if
     !*******************************************
    
    ! create output structure for fcs interface
    error_handle = fcsOutput_create(output_handle, n_local_particles_c, system_dim)

    ! check status
    if(fcsError_getErrorCode(error_handle) /= FCS_SUCCESS) then
       call fcsError_printError(error_handle)
       call MPI_Finalize(mpi_err)
       call EXIT(-1)
    end if

  
    ! 2. step: run pp3mg. There are no optional parameters
    ! available to run the pp3mg. All possible optional parameters have to be specified
    ! during the call of the init_pp3mg_with_opt_param-function 
    error_handle = fcsrun_pp3mg(props_handle, input_handle, output_handle)

    ! check status
    if(fcsError_getErrorCode(error_handle) /= pp3mg_OK) then
       call fcsError_printError(error_handle)
       call MPI_Finalize(mpi_err)
       call EXIT(-1)
    end if

    ! 3. step: free pp3mg resources
    error_handle = fcsfree_pp3mg()

    ! check status
    if(fcsError_getErrorCode(error_handle) /= pp3mg_OK) then
       call fcsError_printError(error_handle)
    end if


    ! get calculated forces
    f = fcsOutput_getForces(output_handle)
    call c_f_pointer( f, forces, [n_local_particles*DIM] )
   
    ! get calculated potentials
    pot = fcsOutput_getPotentials(output_handle)
    call c_f_pointer( pot, e, [n_local_particles] )

    call MPI_Barrier(mpi_comm_cart,mpi_err)
    if (mpi_rank == 0) then
       endtime = MPI_Wtime()
       write(*, *) "processors: ", mpi_size
       write(*,timer_fmt) "pp3mg runtime: ",endtime - starttime," s"
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
         & mpi_comm_cart,mpi_err)

    if (mpi_rank == 0) then
       write(*,"(/, A, T60, ES15.8)") &
            & "Self energy:                                  ", &
            & 1.0D0/(4.0D0 * pi) * 14.0D0/(5.0D0*1/(2*ghosts*min(m,n,o)))
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
  call MPI_Reduce(f_sum_local,f_sum,3,MPI_DOUBLE_PRECISION,MPI_SUM,0,mpi_comm_cart,mpi_err)
  f_sum_squared_x_components = 0.0D0
  call MPI_Reduce(f_sum_squared_x_components_local,f_sum_squared_x_components, &
                  & 1,MPI_DOUBLE_PRECISION,MPI_SUM,0,mpi_comm_cart,mpi_err)
  f_sum_squared_y_components = 0.0D0
  call MPI_Reduce(f_sum_squared_y_components_local,f_sum_squared_y_components, &
                  & 1,MPI_DOUBLE_PRECISION,MPI_SUM,0,mpi_comm_cart,mpi_err)
  f_sum_squared_z_components = 0.0D0
  call MPI_Reduce(f_sum_squared_z_components_local,f_sum_squared_z_components, &
                  & 1,MPI_DOUBLE_PRECISION,MPI_SUM,0,mpi_comm_cart,mpi_err)
  f_max_abs_x = f_max_abs_x_local
  call MPI_Reduce(f_max_abs_x_local,f_max_abs_x,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,mpi_comm_cart,mpi_err)
  f_max_abs_y = f_max_abs_y_local
  call MPI_Reduce(f_max_abs_y_local,f_max_abs_y,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,mpi_comm_cart,mpi_err)
  f_max_abs_z_local = f_max_abs_z_local
  call MPI_Reduce(f_max_abs_z_local,f_max_abs_z,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,mpi_comm_cart,mpi_err)

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
  

  if (n_local_particles > 0) then
     deallocate(x, y, z, q)
     deallocate(fx, fy, fz)
     deallocate(positions)
  end if
  error_handle = fcsInput_destroy(input_handle)
  error_handle = fcsOutput_destroy(output_handle)
  ret_val = fcsError_destroy(error_handle)

  ! Finalizing MPI
  call MPI_Finalize(mpi_err)
  
  stop
  
end program test_pp3mg


