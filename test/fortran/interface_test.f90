#ifdef HAVE_FCONFIG_H
#include <fconfig.h>
#endif
program test

    use fcs_module
    use iso_fortran_env
    use iso_c_binding
    
    implicit none
    include 'mpif.h'

    type(c_ptr)                                                   ::  handle
    type(c_ptr)                                                   ::  ret
    integer(kind = fcs_integer_kind_isoc)                         ::  return_value

    character(len = 160)                                          ::  source = "test_program"
    character(len = 160)                                          ::  message = "testing message"
    real(kind = fcs_real_kind), dimension(:), allocatable         ::  potentials
    real(kind = fcs_real_kind), dimension(:), allocatable         ::  fields
    real(kind = fcs_real_kind), dimension(9)                      ::  virial
    real(kind = fcs_real_kind), parameter                         ::  BOX_SIZE = 1.00d0
    integer(kind = fcs_integer_kind), parameter                   ::  TEST_XYZ_INTERVALL = 1
    integer(kind = fcs_integer_kind), parameter                   ::  RUN_STEP_INTERVAL = 1

    integer, dimension(3)                                         ::  dims

    integer                                                       ::  communicator = MPI_COMM_WORLD
    integer                                                       ::  cart_comm
    logical(kind = fcs_integer_kind)                              ::  short_range_flag = .true.
    real(kind = fcs_real_kind_isoc), dimension(3)                 ::  box_a
    real(kind = fcs_real_kind_isoc), dimension(3)                 ::  box_b
    real(kind = fcs_real_kind_isoc), dimension(3)                 ::  box_c
    real(kind = fcs_real_kind_isoc), dimension(3)                 ::  offset
    integer(kind = fcs_integer_kind_isoc)                         ::  total_particles = -1
    integer(kind = fcs_integer_kind_isoc)                         ::  total_charge_count = -1
    integer(kind = fcs_integer_kind_isoc)                         ::  local_particle_count = -1 
    integer(kind = fcs_integer_kind_isoc)                         ::  local_max_particles = -1
    logical, dimension(3)                                         ::  periodicity
    logical, dimension(3)                                         ::  periodicity_i
    real(kind = fcs_real_kind_isoc), dimension(:), allocatable    ::  local_particles, local_charges
    real(kind = fcs_real_kind_isoc), dimension(:), allocatable    ::  velocities, masses
    integer(kind = fcs_integer_kind_isoc)                         ::  result_destroy

#if FCS_ENABLE_FMM
    integer(kind = fcs_integer_kind_isoc)                         ::  fmm_absrel = FCS_FMM_CUSTOM_RELATIVE
    integer(kind = fcs_integer_kind_isoc)                         ::  fmm_dipole_correction = FCS_FMM_NO_DIPOLE_CORRECTION
    real(kind = fcs_real_kind_isoc)                               ::  fmm_deltaE = 1.d-3
#endif
#if FCS_ENABLE_PEPC
    real(kind = fcs_real_kind_isoc)                               ::  pepc_theta = 0.5d0
    real(kind = fcs_real_kind_isoc)                               ::  pepc_epsilon = 0.5d0
    integer(kind = fcs_integer_kind_isoc)                         ::  pepc_debuglevel = -1
#endif
#if FCS_ENABLE_DIRECT
    real(kind = fcs_real_kind_isoc)                               ::  direct_cutoff = 0.0
#endif
#if FCS_ENABLE_P2NFFT
    real(kind = fcs_real_kind_isoc)                               ::  p2nfft_accuracy = 0.001
#endif
#if FCS_ENABLE_P3M
    real(kind = fcs_real_kind_isoc)                               ::  p3m_accuracy = 0.001
#endif
#if FCS_ENABLE_VMG
    integer(kind = fcs_integer_kind_isoc)                         ::  vmg_max_level = 6
    integer(kind = fcs_integer_kind_isoc)                         ::  vmg_max_iteration = 6
    integer(kind = fcs_integer_kind_isoc)                         ::  vmg_smooth_steps = 3
    integer(kind = fcs_integer_kind_isoc)                         ::  vmg_cycle_type = 2
    real(kind = fcs_real_kind_isoc)                               ::  vmg_precision = 0.001
    integer(kind = fcs_integer_kind_isoc)                         ::  vmg_near_cells = 8
#endif
    integer                                                       ::  my_rank, comm_size, ierr

    integer(kind = fcs_integer_kind)                              ::  command_count
    character(len=120)                                            ::  filename
    character(len = 8)                                            ::  method
    integer(kind = fcs_integer_kind)                              ::  run_count
    character(len=400)                                            ::  read_dummy
    real(kind = fcs_real_kind)                                    ::  read_real_dummy
    character                                                     ::  read_char_dummy

    integer(kind = fcs_integer_kind)                              ::  i,j,k
    integer(kind = fcs_integer_kind)                              ::  chosen_method
    integer                                                       ::  mpifh
    integer, dimension(MPI_STATUS_SIZE)                           ::  mpistatus
    real(kind = fcs_real_kind),dimension(7)                       ::  mpibuf


    ! mass scale = 1.6e-26 kg
    ! charge scale = 1.602e-19 C
    ! time scale = 1.0e-15 s
    ! length scale = 1.0e-9 m
    real(kind = fcs_real_kind),parameter                          ::  time_step = 1e-9
    real(kind = fcs_real_kind),parameter                          ::  charge_scale = 1.602e-19
    real(kind = fcs_real_kind),parameter                          ::  mass_scale = 1.6e-26
    real(kind = fcs_real_kind),parameter                          ::  scaling_factor = &
                                                                      time_step*charge_scale/mass_scale
    real(kind = fcs_real_kind),dimension(3)                       ::  system_momentum
    real(kind = fcs_real_kind),dimension(3)                       ::  system_momentum_local
    real(kind = fcs_real_kind),dimension(3)                       ::  system_momentum_start


#if FCS_ENABLE_DIRECT
    character(len=256, kind = c_char)                             ::  direct_parameters = &
                                            "direct_cutoff,0.0,direct_periodic_images,1,1,1"
#endif
#if FCS_ENABLE_EWALD
    character(len=256, kind = c_char)                             ::  ewald_parameters = &
                                "ewald_required_accuracy,1e-6"
#endif
#if FCS_ENABLE_FMM
    character(len=256, kind = c_char)                             ::  fmm_parameters = &
 "fmm_absrel,2,fmm_dipole_correction,0,fmm_deltaE,1e-6,fmm_maxdepth,19ll," // &
 "fmm_unroll_limit,50ll,fmm_balanceload,1ll,fmm_internal_tuning,1ll"
    integer(kind = c_long_long)                                   ::  fmm_internal_tuning = 1 
#endif
#if FCS_ENABLE_MMM1D
    character(len=256, kind = c_char)                             ::  mmm1d_parameters = &            
     "mmm1d_bessel_cutoff,3,mmm1d_far_switch_radius,6.0,mmm1d_maxPWerror,1e-6"
#endif
#if FCS_ENABLE_MMM2D
    character(len=256, kind = c_char)                             ::  mmm2d_parameters = &
"mmm2d_max_PWerror,1e-6,mmm2d_far_cutoff,1.73,mmm2d_dielectric_contrasts,3.17,2.13,mmm2d_layer_per_node,10,mmm2d_skin,0.5"
#endif
#if FCS_ENABLE_PEPC
    character(len=256, kind = c_char)                             ::  pepc_parameters = & 
                   "pepc_debuglevel,-1,pepc_epsilon,0.5,pepc_theta,0.5"
    integer :: provided
#endif
#if FCS_ENABLE_P2NFFT
    character(len=256, kind = c_char)                             ::  p2nfft_parameters = & 
                   "p2nfft_required_accuracy,1e-6"
#endif
#if FCS_ENABLE_P3M
    character(len=256, kind = c_char)                             ::  p3m_parameters = &
                   "p3m_required_accuracy,1e-6"
#endif
#if FCS_ENABLE_VMG
    character(len=256, kind = c_char)                             ::  vmg_parameters = & 
"vmg_cycle_type,2,vmg_max_iterations,20,vmg_max_level,6,vmg_near_field_cells,6,vmg_precision,1e-6,vmg_smoothing_steps,3"
#endif
    character(len=256, kind = c_char)                             ::  common_parameters = &
"box_a,1.01,0.0,0.0,box_b,0.0,1.01,0.0,box_c,0.0,0.0,1.01,periodicity,1,1,1,offset,0.0,0.0,0.0,near_field_flag,0"

#if FCS_ENABLE_PEPC
    call MPI_INIT_THREAD(MPI_THREAD_MULTIPLE, provided, ierr)
#else
    call MPI_INIT(ierr)
#endif

    periodicity(3) = .true.
    periodicity(2) = .true.
    periodicity(1) = .true.
    dims = 0

    call MPI_COMM_SIZE(communicator, comm_size, ierr)
    write(*,*) "comm_size: ", comm_size
    call MPI_DIMS_CREATE(comm_size,3,dims,ierr)
    write(*,*) "dims: ", dims
    call MPI_CART_CREATE(communicator,3,dims,periodicity,.false.,cart_comm,ierr)
    communicator = cart_comm
    call MPI_COMM_RANK(communicator, my_rank, ierr)
    call MPI_COMM_SIZE(communicator, comm_size, ierr)
    call MPI_CART_COORDS(communicator,my_rank,3,dims,ierr)
    !write(*,'(a,i7,a,i7,a,3i7)') "rank/size -> coords", my_rank, "/", comm_size, " -> ", dims
    
    command_count = COMMAND_ARGUMENT_COUNT()
    ! check if command is called correctly
    if (command_count /= 3) then
        if (my_rank == 0) then
            write(*,'(a,i3)') "command arguments: ", command_count
            write(*,'(a)') "The test program was called in a wrong way, please use:"
            write(*,'(a)') "test <solver> <input_file_name> <run_count>"
            write(*,'(a)') "(currently only with similar format as files from test_env/data are supported)"
        end if
        call MPI_FINALIZE(ierr)
    endif

    call GET_COMMAND_ARGUMENT(1,method)
    call GET_COMMAND_ARGUMENT(2,filename)
    call GET_COMMAND_ARGUMENT(3,read_dummy)

    read(read_dummy,'(i10)') run_count

    open(1409, FILE=filename, ACTION="read", FORM="formatted")

    read(1409,'(a)') read_dummy
    read(1409,'(a)') read_dummy
    read(1409,'(a)') read_dummy
    read(1409,'(a)') read_dummy
    read(1409,'(a,i20)') read_char_dummy, total_particles
    total_charge_count = total_particles

    if (my_rank /= comm_size-1) then
        local_particle_count = total_particles / comm_size
        allocate(local_particles(3*local_particle_count),local_charges(local_particle_count))
        allocate(fields(3*local_particle_count),potentials(local_particle_count))
        allocate(velocities(3*local_particle_count),masses(local_particle_count))
        j = 0
        do i = 0,total_particles-1
            if (i .ge. local_particle_count * my_rank .and. i .lt. local_particle_count * (my_rank+1)) then
                j = j+1
                read(1409,'(8es14.6)') local_particles(3*j-2:3*j), masses(j), local_charges(j),velocities(3*j-2:3*j)
            else
                read(1409,'(a)') read_dummy
            endif
        end do
    else
        local_particle_count = total_particles - ((total_particles / comm_size ) * (comm_size-1))
        allocate(local_particles(3*local_particle_count),local_charges(local_particle_count))
        allocate(fields(3*local_particle_count),potentials(local_particle_count))
        allocate(velocities(3*local_particle_count),masses(local_particle_count))
        j = 0
        do i = 0,total_particles-1
            if (i .ge. total_particles - local_particle_count) then
                j = j+1
                read(1409,'(8es14.6)') local_particles(3*j-2:3*j), masses(j), local_charges(j),velocities(3*j-2:3*j)
            else
                read(1409,'(a)') read_dummy
            endif
        end do
    endif
    where (masses < 1e-6)
        masses = 1.0d0
    end where
!    write (*,'(a,i7,i7)') 'local particles: ', my_rank, local_particle_count

    local_max_particles = 2*local_particle_count

!    do i = 0, comm_size-1
!        if (i == my_rank) then
!            do j = 1, local_particle_count
!                write(*,'(i7,4es14.6)') my_rank, local_particles(3*j-2:3*j), local_charges(j)
!            end do
!        end if
!        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!    end do


!    allocate(potentials(local_particle_count),field(local_particle_count),virial(local_particle_count))

    box_a(1) = BOX_SIZE
    box_a(2) = 0.0
    box_a(3) = 0.0
    box_b(1) = 0.0
    box_b(2) = BOX_SIZE
    box_b(3) = 0.0
    box_c(1) = 0.0
    box_c(2) = 0.0
    box_c(3) = BOX_SIZE
    offset = 0.0

    

    fields = 0.0
    potentials = 0.0

    if (my_rank == 0) write(*,*) "chosen method: ", method

    if (my_rank == 0) write(*,*) "----------------------------call init---------------------------------"
    ret = fcs_init(handle, trim(adjustl(method)) // c_null_char, communicator)
    if (my_rank == 0) write(*,*) "fcs_init returns: ", fcsResult_getReturnCode(ret)
    if (my_rank == 0) write(*,*) "fcs_init returns: ", trim(adjustl(fcsResult_getErrorMessage(ret)))
    if (my_rank == 0) write(*,*) "fcs_init returns: ", trim(adjustl(fcsResult_getErrorSource(ret)))
    if (my_rank == 0) result_destroy = fcsResult_destroy(ret)
    if (my_rank == 0 .and. result_destroy == -1) write(*,*) "fcsResult_destroy failed!"
    
    call MPI_BARRIER(communicator,ierr)

    if (my_rank == 0) write(*,*) "----------------------------call common setter---------------------------------"
    ret = fcs_parser(handle, trim(adjustl(common_parameters)) // c_null_char, FCS_FALSE)
    !ret =  fcs_set_common(handle, short_range_flag, box_a, box_b, box_c, offset, periodicity, total_particles)
    if (my_rank == 0) write(*,*) "fcs_set_common (parser) returns: ", fcsResult_getReturnCode(ret)
    if (my_rank == 0) write(*,*) "fcs_set_common (parser) returns: ", trim(adjustl(fcsResult_getErrorMessage(ret)))
    if (my_rank == 0) write(*,*) "fcs_set_common (parser) returns: ", trim(adjustl(fcsResult_getErrorSource(ret)))
    if (my_rank == 0) result_destroy = fcsResult_destroy(ret)
    if (my_rank == 0 .and. result_destroy == -1) write(*,*) "fcsResult_destroy failed!"
    ret = fcs_set_total_particles(handle, total_particles)
    if (my_rank == 0) write(*,*) "fcs_set_common (total particles) returns: ", fcsResult_getReturnCode(ret)
    if (my_rank == 0) write(*,*) "fcs_set_common (total particles) returns: ", trim(adjustl(fcsResult_getErrorMessage(ret)))
    if (my_rank == 0) write(*,*) "fcs_set_common (total particles) returns: ", trim(adjustl(fcsResult_getErrorSource(ret)))
    if (my_rank == 0) result_destroy = fcsResult_destroy(ret)
    if (my_rank == 0 .and. result_destroy == -1) write(*,*) "fcsResult_destroy failed!"
    if (fcs_get_method(handle) == FCS_MMM1D) then
      periodicity_i(1) = .false.
      periodicity_i(2) = .false.
      periodicity_i(3) = .true.
      ret = fcs_set_periodicity(handle, periodicity_i)
      if (my_rank == 0) write(*,*) "fcs_set_common (mmm1d) returns: ", fcsResult_getReturnCode(ret)
      if (my_rank == 0) write(*,*) "fcs_set_common (mmm1d) returns: ", trim(adjustl(fcsResult_getErrorMessage(ret)))
      if (my_rank == 0) write(*,*) "fcs_set_common (mmm1d) returns: ", trim(adjustl(fcsResult_getErrorSource(ret)))
      if (my_rank == 0) result_destroy = fcsResult_destroy(ret)
    if (my_rank == 0 .and. result_destroy == -1) write(*,*) "fcsResult_destroy failed!"
    else if (fcs_get_method(handle) == FCS_MMM2D) then
      periodicity_i(1) = .true.
      periodicity_i(2) = .true.
      periodicity_i(3) = .false.
      ret = fcs_set_periodicity(handle, periodicity_i)
      if (my_rank == 0) write(*,*) "fcs_set_common (mmm2d) returns: ", fcsResult_getReturnCode(ret)
      if (my_rank == 0) write(*,*) "fcs_set_common (mmm2d) returns: ", trim(adjustl(fcsResult_getErrorMessage(ret)))
      if (my_rank == 0) write(*,*) "fcs_set_common (mmm2d) returns: ", trim(adjustl(fcsResult_getErrorSource(ret)))    
      if (my_rank == 0) result_destroy = fcsResult_destroy(ret)
    if (my_rank == 0 .and. result_destroy == -1) write(*,*) "fcsResult_destroy failed!"
  end if
    

    call MPI_BARRIER(communicator,ierr)

    if (my_rank == 0) write(*,*) "-----------------------call print content (standard) ----------------------"
    if (my_rank == 0) call fcs_printHandle(handle)
    if (my_rank == 0) write(*,*) "-----------------------call method-specific setter-----------------------------"
    select case (fcs_get_method(handle))
    
#ifdef FCS_ENABLE_DIRECT
      case (FCS_DIRECT)
        !ret = fcs_direct_setup(handle, direct_cutoff)
        ret = fcs_parser(handle, trim(adjustl(direct_parameters)) // c_null_char, FCS_FALSE)
#endif
#ifdef FCS_ENABLE_EWALD
      case (FCS_EWALD)
        ret = fcs_parser(handle, trim(adjustl(ewald_parameters)) // c_null_char, FCS_FALSE)
#endif
#ifdef FCS_ENABLE_FMM
      case (FCS_FMM)
        !ret = fcs_fmm_setup(handle, fmm_absrel, fmm_deltaE, fmm_dipole_correction)
        ret = fcs_parser(handle, trim(adjustl(fmm_parameters)) // c_null_char, FCS_FALSE)
        ret = fcs_fmm_set_internal_tuning(handle,fmm_internal_tuning)
#endif
#ifdef FCS_ENABLE_MEMD
      case (FCS_MEMD)
#endif
#ifdef FCS_ENABLE_MMM1D
      case (FCS_MMM1D)
        ret = fcs_parser(handle, trim(adjustl(mmm1d_parameters)) // c_null_char, FCS_FALSE)
#endif
#ifdef FCS_ENABLE_MMM2D
      case (FCS_MMM2D)
        write(*,*) 'FINBLA!!'
        ret = fcs_parser(handle, trim(adjustl(mmm2d_parameters)) // c_null_char, FCS_FALSE)
#endif
#ifdef FCS_ENABLE_PEPC
      case (FCS_PEPC)
        !ret = fcs_pepc_setup(handle, pepc_epsilon, pepc_theta, pepc_debuglevel)
        ret = fcs_parser(handle, trim(adjustl(pepc_parameters)) // c_null_char, FCS_FALSE)
#endif
#ifdef FCS_ENABLE_P2NFFT
      case (FCS_P2NFFT)
        !ret = fcs_p2nfft_set_required_accuracy(handle, p2nfft_accuracy)
        ret = fcs_parser(handle, trim(adjustl(p2nfft_parameters)) // c_null_char, FCS_FALSE)
#endif
#ifdef FCS_ENABLE_PP3MG
      case (FCS_PP3MG)
#endif
#ifdef FCS_ENABLE_P3M
      case (FCS_P3M)
        !ret = fcs_p3m_set_required_accuracy(handle, p3m_accuracy)
        ret = fcs_parser(handle, trim(adjustl(p3m_parameters)) // c_null_char, FCS_FALSE)
#endif
#ifdef FCS_ENABLE_VMG
      case (FCS_VMG)
        !ret = fcs_vmg_setup(handle, vmg_max_level, vmg_max_iteration, vmg_smooth_steps, vmg_cycle_type, vmg_precision, vmg_near_cells)
        ret = fcs_parser(handle, trim(adjustl(vmg_parameters)) // c_null_char, FCS_FALSE)
#endif
      case (FCS_NO_METHOD_CHOSEN)
        if (my_rank == 0) write(*,*) "unknown / no method chosen, aborting program run ..."
        call MPI_FINALIZE(ierr)
        stop
    end select
    if (my_rank == 0) write(*,*) "method specific setter returns: ", fcsResult_getReturnCode(ret)
    if (my_rank == 0) write(*,*) "method specific setter returns: ", trim(adjustl(fcsResult_getErrorMessage(ret)))
    if (my_rank == 0) write(*,*) "method specific setter returns: ", trim(adjustl(fcsResult_getErrorSource(ret)))
    if (my_rank == 0) result_destroy = fcsResult_destroy(ret)
    if (my_rank == 0 .and. result_destroy == -1) write(*,*) "fcsResult_destroy failed!"

    if (my_rank == 0) write(*,*) "----------------------------call print content---------------------------------"
    if (my_rank == 0) call fcs_printHandle(handle)

    call MPI_Barrier(communicator,ierr)

    if (my_rank == 0) write(*,*) "----------------------------call tune---------------------------------"
!     call fcs_tune(handle, local_particle_count, local_max_particles, local_particles, local_charges, return_value)
!     if (my_rank == 0) write(*,*) "fcs_tune returns: ", return_value
    ret = fcs_tune(handle, local_particle_count, local_max_particles, local_particles, local_charges)
    if (my_rank == 0) write(*,*) "fcs_tune returns: ", fcsResult_getReturnCode(ret)
    if (my_rank == 0) write(*,*) "fcs_tune returns: ", trim(adjustl(fcsResult_getErrorMessage(ret)))
    if (my_rank == 0) write(*,*) "fcs_tune returns: ", trim(adjustl(fcsResult_getErrorSource(ret)))
    if (my_rank == 0) result_destroy = fcsResult_destroy(ret)
    if (my_rank == 0 .and. result_destroy == -1) write(*,*) "fcsResult_destroy failed!"

    call MPI_BARRIER(communicator,ierr)

    if (my_rank == 0) write(*,*) "----------------------------call require virial---------------------------------"
    ret = fcs_require_virial(handle,.true.)
    if (my_rank == 0) write(*,*) "fcs_require_virial returns: ", fcsResult_getReturnCode(ret)
    if (my_rank == 0) write(*,*) "fcs_require_virial returns: ", trim(adjustl(fcsResult_getErrorMessage(ret)))
    if (my_rank == 0) write(*,*) "fcs_require_virial returns: ", trim(adjustl(fcsResult_getErrorSource(ret)))
    if (my_rank == 0) result_destroy = fcsResult_destroy(ret)
    if (my_rank == 0 .and. result_destroy == -1) write(*,*) "fcsResult_destroy failed!"

    !output of starting configuration
    do j = 0,comm_size-1
        if(my_rank == j) then
            if (my_rank == 0 .and. i == TEST_XYZ_INTERVALL) then
                 open(984, FILE="F_test.xyz", ACTION="write")
            else
                 open(984, FILE="F_test.xyz", ACTION="write", POSITION="append")
            end if
            if (my_rank == 0) write(984, '(i10)') total_particles
            if (my_rank == 0) write(984, '(2a,2i12,3L2)') "F ", method, comm_size, run_count, periodicity
            do k = 1,local_particle_count
                if (local_charges(k) > 0.0) then 
                  write (984,'(a,7es17.6)') 'A ',local_particles(3*k-2:3*k), local_charges(k), velocities(3*k-2:3*k)
                else
                  write (984,'(a,7es17.6)') 'B ',local_particles(3*k-2:3*k), local_charges(k), velocities(3*k-2:3*k)
                endif
            end do
            close(984)
        end if
        call MPI_BARRIER(communicator,ierr) 
    end do
    
    system_momentum = 0.0d0
    system_momentum_local = 0.0d0
    do i = 1, local_particle_count
        system_momentum_local = system_momentum_local + masses(i) * velocities(3*i-2:3*i)
    end do
    call MPI_REDUCE(system_momentum_local,system_momentum,3,MPI_DOUBLE_PRECISION,MPI_SUM,0,communicator,ierr)
    if (my_rank == 0) system_momentum_start = system_momentum
    if (my_rank == 0) write(*,'(a,6es16.9)') 'system momentum (absolute value, absolute to start value)', system_momentum, &
                                                                                      system_momentum-system_momentum_start

    do i = 1, run_count
        if (my_rank == 0 .and. modulo(i,RUN_STEP_INTERVAL) == 0) write(*,'(a,i7,a)') "----------------------------call run ", i,&
                                               "---------------------------------"
        ret = fcs_run(handle, local_particle_count, local_max_particles, local_particles, local_charges, fields, &
                      potentials)
        if (my_rank == 0 .and. modulo(i,RUN_STEP_INTERVAL) == 0) write(*,*) "fcs_run returns: ",&
        fcsResult_getReturnCode(ret)
        if (my_rank == 0 .and. modulo(i,RUN_STEP_INTERVAL) == 0) write(*,*) "fcs_run returns: ",&
        trim(adjustl(fcsResult_getErrorMessage(ret)))
        if (my_rank == 0 .and. modulo(i,RUN_STEP_INTERVAL) == 0) write(*,*) "fcs_run returns: ",&
        trim(adjustl(fcsResult_getErrorSource(ret)))
        if (my_rank == 0 .and. modulo(i,RUN_STEP_INTERVAL) == 0) result_destroy = fcsResult_destroy(ret)
        if (my_rank == 0 .and. modulo(i,RUN_STEP_INTERVAL) == 0 .and. result_destroy == -1) &
            write (*,*) "fcsResult_destroy failed!"
    
        !LEAP-FROG step
        do j = 1,local_particle_count
            !v(t+0.5dt) = v(t-0.5dt)+E(t)*q*dt/m
            velocities(3*j-2:3*j) = velocities(3*j-2:3*j) + fields(3*j-2:3*j) *&
                                    local_charges(j) / masses(j) * scaling_factor
            !r(t+dt) = r(t)+v(t+0.5dt)*dt
            local_particles(3*j-2:3*j) = local_particles(3*j-2:3*j) + velocities(3*j-2:3*j) *&
                                         time_step

            if (periodicity(1)) then
                local_particles(3*j-2) = modulo(local_particles(3*j-2),BOX_SIZE)
            end if
            if (periodicity(2)) then
                local_particles(3*j-1) = modulo(local_particles(3*j-1),BOX_SIZE)
            end if
            if (periodicity(3)) then
                local_particles(3*j) = modulo(local_particles(3*j),BOX_SIZE)
            end if
        end do

        if (modulo (i,TEST_XYZ_INTERVALL) == 0) then
        do j = 0,comm_size-1
               if(my_rank == j) then
                    if (my_rank == 0 .and. i == TEST_XYZ_INTERVALL) then
                    open(984, FILE="F_test.xyz", ACTION="write")
                else
                    open(984, FILE="F_test.xyz", ACTION="write", POSITION="append")
                end if
                if (my_rank == 0) write(984, '(i10)') total_particles
                if (my_rank == 0) write(984, '(2a,2i12,3L2)') "F ", method, comm_size, run_count, periodicity
                do k = 1,local_particle_count
                    if (local_charges(k) > 0.0) then 
              write (984,'(a,7es17.6)') 'A ',local_particles(3*k-2:3*k), local_charges(k), velocities(3*k-2:3*k)
            else
              write (984,'(a,7es17.6)') 'B ',local_particles(3*k-2:3*k), local_charges(k), velocities(3*k-2:3*k)
            endif
                end do
                close(984)
            end if
            call MPI_BARRIER(communicator,ierr)
        end do
        end if

    if (modulo (i,TEST_XYZ_INTERVALL) == 0) then
        do j = 0,comm_size-1
               if(my_rank == j) then
                    if (my_rank == 0 .and. i == TEST_XYZ_INTERVALL) then
                    open(984, FILE="F_test.pot", ACTION="write")
                else
                    open(984, FILE="F_test.pot", ACTION="write", POSITION="append")
                end if
                if (my_rank == 0) write(984, '(i10)') total_particles
                if (my_rank == 0) write(984, '(2a,2i12,3L2)') "F ", method, comm_size, run_count, periodicity
                do k = 1,local_particle_count
                    mpibuf(1:3) = local_particles(3*k-2:3*k)
                    mpibuf(4) = local_charges(k)
                    mpibuf(5:7) = velocities(3*k-2:3*k)
                    if (local_charges(k) > 0.0) then 
              write (984,'(a,i7,5es17.6)') 'A ', k, local_charges(k), potentials(k), fields(3*k-2:3*k)
            else
              write (984,'(a,i7,5es17.6)') 'B ', k, local_charges(k), potentials(k), fields(3*k-2:3*k)
            endif
                end do
                close(984)
            end if
            call MPI_BARRIER(communicator,ierr)
        end do
        end if

        if (modulo(i,RUN_STEP_INTERVAL) == 0) then
          system_momentum = 0.0d0
          system_momentum_local = 0.0d0
          do j = 1, local_particle_count
              system_momentum_local = system_momentum_local + masses(j) * velocities(3*j-2:3*j)
          end do
          call MPI_REDUCE(system_momentum_local,system_momentum,3,MPI_DOUBLE_PRECISION,MPI_SUM,0,communicator,ierr)
        end if

        if (my_rank == 0 .and. modulo(i,RUN_STEP_INTERVAL) == 0) then
          write(*,*) "----------------------------call get virial---------------------------------"
          ret = fcs_get_virial(handle,virial)
          write(*,*) "fcs_get_virial returns: ", fcsResult_getReturnCode(ret)
          write(*,*) "fcs_get_virial returns: ", trim(adjustl(fcsResult_getErrorMessage(ret)))
          write(*,*) "fcs_get_virial returns: ", trim(adjustl(fcsResult_getErrorSource(ret)))
          result_destroy = fcsResult_destroy(ret)
          if (result_destroy == -1) write(*,*) "fcsResult_destory failed!"
        
          write(*,'(a,i7)') "system virial in step ", i
          do j = 0,2
            write(*,'(3f20.6)') virial(3*j+1), virial(3*j+2), virial(3*j+3)
          enddo

          write(*,'(a,6es16.9)') 'system momentum (absolute value, absolute to start value)', system_momentum, &
                                                                                            system_momentum-system_momentum_start

        end if


    end do

    close(1409)

    if (my_rank == 0) write(*,*) "----------------------------call destroy---------------------------------"
    ret = fcs_destroy(handle)
    if (my_rank == 0) write(*,*) "fcs_destroy returns: ", fcsResult_getReturnCode(ret)
    if (my_rank == 0) write(*,*) "fcs_destroy returns: ", trim(adjustl(fcsResult_getErrorMessage(ret)))
    if (my_rank == 0) write(*,*) "fcs_destroy returns: ", trim(adjustl(fcsResult_getErrorSource(ret)))
    if (my_rank == 0) result_destroy = fcsResult_destroy(ret)
    if (my_rank == 0 .and. result_destroy == -1) write(*,*) "fcsResult_destroy failed!"
    

    deallocate(local_particles)
    deallocate(velocities)
    deallocate(fields)
    deallocate(potentials)
    deallocate(masses)
    deallocate(local_charges)
    
    call MPI_FINALIZE(ierr)

!    deallocate(potentials,field,virial)

end program test
