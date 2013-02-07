program test
#include <fcs_fconfig.h>

    use fcs_module
    use iso_fortran_env
    use iso_c_binding
    
    implicit none
    include 'mpif.h'

    integer(kind = fcs_integer_kind), parameter     ::  run_count = 1
    type(c_ptr)                                     ::  handle
    type(c_ptr)                                     ::  ret
    real(kind = fcs_real_kind), dimension(8)        ::    potentials
    real(kind = fcs_real_kind), dimension(24)       ::    fields
    real(kind = fcs_real_kind), dimension(9)        ::  virial
    real(kind = fcs_real_kind), parameter           ::  BOX_SIZE = 1.00d0
    integer                                         ::  communicator = MPI_COMM_WORLD
    logical(kind = fcs_integer_kind)                ::  short_range_flag = .true.
    real(kind = fcs_real_kind_isoc), dimension(3)   ::  box_a = (/BOX_SIZE,0.0d0,0.0d0/)
    real(kind = fcs_real_kind_isoc), dimension(3)   ::  box_b = (/0.0d0,BOX_SIZE,0.0d0/)
    real(kind = fcs_real_kind_isoc), dimension(3)   ::  box_c = (/0.0d0,0.0d0,BOX_SIZE/)
    real(kind = fcs_real_kind_isoc), dimension(3)   ::  offset = (/0.0d0,0.0d0,0.0d0/)
    logical, dimension(3)                           ::  periodicity = (/.true.,.true.,.true./)
    logical                                         ::  l_virial = .true.
    integer(kind = fcs_integer_kind_isoc)           ::  total_particles = -1
    integer(kind = fcs_integer_kind_isoc)           ::  local_particle_count = -1 
    integer(kind = fcs_integer_kind_isoc)           ::  local_max_particles = -1
    real(kind = fcs_real_kind_isoc), dimension(8)   ::  local_charges
    real(kind = fcs_real_kind_isoc), dimension(24)  ::  local_particles

#if FCS_ENABLE_FMM
    ! FMM parameters
    integer(kind = fcs_integer_kind_isoc)           ::  fmm_absrel = FCS_FMM_CUSTOM_RELATIVE
    integer(kind = fcs_integer_kind_isoc)           ::  fmm_dipole_correction =  FCS_FMM_NO_DIPOLE_CORRECTION
    real(kind = fcs_real_kind_isoc)                 ::  fmm_deltaE = 1.d-3
    integer(kind = c_long_long)                     ::  fmm_maxdepth = 19
    integer(kind = c_long_long)                     ::  fmm_unroll_limit = 50
    integer(kind = c_long_long)                     ::  fmm_balance = 0
#endif
    real(kind = fcs_real_kind_isoc)                 ::  e_local, e_total
    integer                                         ::  my_rank, comm_size, ierr
    character(len = 8)                              ::  method
    
    integer(kind = fcs_integer_kind)                ::  i,j,k,l

   
    call MPI_INIT(ierr)

    method = "fmm"

    ! set periodicity 
    periodicity(3) = .true.
    periodicity(2) = .true.
    periodicity(1) = .true.

    call MPI_COMM_SIZE(communicator, comm_size, ierr)
    call MPI_COMM_RANK(communicator, my_rank, ierr)
    
    total_particles = 8

    local_particle_count = 8 / comm_size
    l = 0
    ! different cases to be able to run with 1,2,4,8 processes (since no internal sorting is within the fmm to
    ! ensure that at least one particle gets on each process)
    if (comm_size == 8) then
        do i = iand(my_rank,1),iand(my_rank,1)
            do j = iand(my_rank,2)/2,iand(my_rank,2)/2
                do k = iand(my_rank,4)/4,iand(my_rank,4)/4
                    l = l + 1
                    local_particles(3*l-2) = 0.25d0 + 0.5d0*k
                    local_particles(3*l-1) = 0.25d0 + 0.5d0*j
                    local_particles(3*l)   = 0.25d0 + 0.5d0*i
                    local_charges(l)       = (-1.0d0)**(i+j+k)
                end do
            end do
        end do
    else if (comm_size == 4) then
        do i = iand(my_rank,1),iand(my_rank,1)
            do j = iand(my_rank,2)/2,iand(my_rank,2)/2
                do k = 0,1
                    l = l + 1
                    local_particles(3*l-2) = 0.25d0 + 0.5d0*k
                    local_particles(3*l-1) = 0.25d0 + 0.5d0*j
                    local_particles(3*l)   = 0.25d0 + 0.5d0*i
                    local_charges(l)       = (-1.0d0)**(i+j+k)
                end do
            end do
        end do
    else if (comm_size == 2) then
        do i = iand(my_rank,1),iand(my_rank,1)
            do j = 0,1
                do k = 0,1
                    l = l + 1
                    local_particles(3*l-2) = 0.25d0 + 0.5d0*k
                    local_particles(3*l-1) = 0.25d0 + 0.5d0*j
                    local_particles(3*l)   = 0.25d0 + 0.5d0*i
                    local_charges(l)       = (-1.0d0)**(i+j+k)
                end do
            end do
        end do
    else if (comm_size == 1) then
        do i = 0,1
            do j = 0,1
                do k = 0,1
                    l = l + 1
                    local_particles(3*l-2) = 0.25d0 + 0.5d0*k
                    local_particles(3*l-1) = 0.25d0 + 0.5d0*j
                    local_particles(3*l)   = 0.25d0 + 0.5d0*i
                    local_charges(l)       = (-1.0d0)**(i+j+k)
                end do
            end do
        end do
    end if

    write (*,'(a,i7,i7)') 'local particles: ', my_rank, local_particle_count
    do i = 1,local_particle_count
        write(*,'(i7,a,i7,4es14.7)') my_rank, '---->', i, local_particles(3*i-2:3*i), local_charges(i) 
    end do

    local_max_particles = 2*local_particle_count

    fields = 0.0d0
    potentials = 0.0d0

    if (my_rank == 0) write(*,*) "----------------------------call init---------------------------------"
    ret = fcs_init(handle, trim(adjustl(method)) // c_null_char, communicator)
    if (my_rank == 0) write(*,*) "fcs_init returns: ", fcsResult_getReturnCode(ret)
    if (my_rank == 0) write(*,*) "fcs_init returns: ", trim(adjustl(fcsResult_getErrorMessage(ret)))
    if (my_rank == 0) write(*,*) "fcs_init returns: ", trim(adjustl(fcsResult_getErrorSource(ret)))
    
    call MPI_BARRIER(communicator,ierr)

    if (my_rank == 0) write(*,*) "----------------------------call common setter---------------------------------"
    ret =  fcs_common_set(handle, short_range_flag, box_a, box_b, box_c, offset, periodicity, total_particles)
    if (my_rank == 0) write(*,*) "fcs_common_set (parser) returns: ", fcsResult_getReturnCode(ret)
    if (my_rank == 0) write(*,*) "fcs_common_set (parser) returns: ", trim(adjustl(fcsResult_getErrorMessage(ret)))
    if (my_rank == 0) write(*,*) "fcs_common_set (parser) returns: ", trim(adjustl(fcsResult_getErrorSource(ret)))
    
    call MPI_BARRIER(communicator,ierr)

    if (my_rank == 0) write(*,*) "-----------------------call print content (standard) ----------------------"
    if (my_rank == 0) call fcs_printContent(handle)
    if (my_rank == 0) write(*,*) "-----------------------call method-specific setter-----------------------------"

    ret = fcs_fmm_set_absrel(handle, fmm_absrel)
    if (my_rank == 0) write(*,*) "setter (absrel) returns: ", fcsResult_getReturnCode(ret)
    if (my_rank == 0) write(*,*) "setter (absrel) returns: ", trim(adjustl(fcsResult_getErrorMessage(ret)))
    if (my_rank == 0) write(*,*) "setter (absrel) returns: ", trim(adjustl(fcsResult_getErrorSource(ret)))
    ret = fcs_fmm_set_deltaE(handle, fmm_deltaE)
    if (my_rank == 0) write(*,*) "setter (deltaE) returns: ", fcsResult_getReturnCode(ret)
    if (my_rank == 0) write(*,*) "setter (deltaE) returns: ", trim(adjustl(fcsResult_getErrorMessage(ret)))
    if (my_rank == 0) write(*,*) "setter (deltaE) returns: ", trim(adjustl(fcsResult_getErrorSource(ret)))
    ret = fcs_fmm_set_dipole_correction(handle, fmm_dipole_correction)
    if (my_rank == 0) write(*,*) "setter (dipole_correction) returns: ", fcsResult_getReturnCode(ret)
    if (my_rank == 0) write(*,*) "setter (dipole_correction) returns: ", trim(adjustl(fcsResult_getErrorMessage(ret)))
    if (my_rank == 0) write(*,*) "setter (dipole_correction) returns: ", trim(adjustl(fcsResult_getErrorSource(ret)))
    ret = fcs_fmm_set_maxdepth(handle, fmm_maxdepth)
    if (my_rank == 0) write(*,*) "setter (maxdepth) returns: ", fcsResult_getReturnCode(ret)
    if (my_rank == 0) write(*,*) "setter (maxdepth) returns: ", trim(adjustl(fcsResult_getErrorMessage(ret)))
    if (my_rank == 0) write(*,*) "setter (maxdepth) returns: ", trim(adjustl(fcsResult_getErrorSource(ret)))
    ret = fcs_fmm_set_unroll_limit(handle, fmm_unroll_limit)
    if (my_rank == 0) write(*,*) "setter (unroll_limit) returns: ", fcsResult_getReturnCode(ret)
    if (my_rank == 0) write(*,*) "setter (unroll_limit) returns: ", trim(adjustl(fcsResult_getErrorMessage(ret)))
    if (my_rank == 0) write(*,*) "setter (unroll_limit) returns: ", trim(adjustl(fcsResult_getErrorSource(ret)))
    ret = fcs_fmm_set_balanceload(handle, fmm_balance)
    if (my_rank == 0) write(*,*) "setter (balanceload) returns: ", fcsResult_getReturnCode(ret)
    if (my_rank == 0) write(*,*) "setter (balanceload) returns: ", trim(adjustl(fcsResult_getErrorMessage(ret)))
    if (my_rank == 0) write(*,*) "setter (balanceload) returns: ", trim(adjustl(fcsResult_getErrorSource(ret)))

    if (my_rank == 0) write(*,*) "----------------------------call print content---------------------------------"
    if (my_rank == 0) call fcs_printContent(handle)

    call MPI_Barrier(communicator,ierr)

    if (my_rank == 0) write(*,*) "----------------------------call tune---------------------------------"
    ret = fcs_tune(handle, local_particle_count, local_max_particles, local_particles, local_charges)
    if (my_rank == 0) write(*,*) "fcs_tune returns: ", fcsResult_getReturnCode(ret)
    if (my_rank == 0) write(*,*) "fcs_tune returns: ", trim(adjustl(fcsResult_getErrorMessage(ret)))
    if (my_rank == 0) write(*,*) "fcs_tune returns: ", trim(adjustl(fcsResult_getErrorSource(ret)))

    call MPI_BARRIER(communicator,ierr)

    if (my_rank == 0) write(*,*) "----------------------------call require virial---------------------------------"
    ret = fcs_require_virial(handle,l_virial)
    if (my_rank == 0) write(*,*) "fcs_require_virial returns: ", fcsResult_getReturnCode(ret)
    if (my_rank == 0) write(*,*) "fcs_require_virial returns: ", trim(adjustl(fcsResult_getErrorMessage(ret)))
    if (my_rank == 0) write(*,*) "fcs_require_virial returns: ", trim(adjustl(fcsResult_getErrorSource(ret)))

    
    do i = 1, run_count
      if (my_rank == 0) write(*,*) "----------------------------call run---------------------------------"
      ret = fcs_run(handle, local_particle_count, local_max_particles, local_particles, local_charges, fields, &
                    potentials)
      if (my_rank == 0) write(*,*) "fcs_run returns: ", fcsResult_getReturnCode(ret)
      if (my_rank == 0) write(*,*) "fcs_run returns: ", trim(adjustl(fcsResult_getErrorMessage(ret)))
      if (my_rank == 0) write(*,*) "fcs_run returns: ", trim(adjustl(fcsResult_getErrorSource(ret)))

      e_local = 0.0d0
      do j = 1, local_particle_count
        e_local = e_local + local_charges(j) * potentials(j)
      end do
      call MPI_REDUCE(e_local,e_total,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,communicator,ierr)
      if (my_rank == 0) write (*,'(a,es14.7)') "total energy: ", e_total

      if (my_rank == 0) then
        write(*,*) "----------------------------call get virial---------------------------------"
        ret = fcs_get_virial(handle,virial)
        write(*,*) "fcs_get_virial returns: ", fcsResult_getReturnCode(ret)
        write(*,*) "fcs_get_virial returns: ", trim(adjustl(fcsResult_getErrorMessage(ret)))
        write(*,*) "fcs_get_virial returns: ", trim(adjustl(fcsResult_getErrorSource(ret)))
        
        write(*,'(a,i7)') "system virial in step ", i
        do j = 0,2
          write(*,'(3f20.6)') virial(3*j+1), virial(3*j+2), virial(3*j+3)
        enddo
      end if
    end do

    if (my_rank == 0) write(*,*) "----------------------------call destroy---------------------------------"
    ret = fcs_destroy(handle)
    if (my_rank == 0) write(*,*) "fcs_destroy returns: ", fcsResult_getReturnCode(ret)
    if (my_rank == 0) write(*,*) "fcs_destroy returns: ", trim(adjustl(fcsResult_getErrorMessage(ret)))
    if (my_rank == 0) write(*,*) "fcs_destroy returns: ", trim(adjustl(fcsResult_getErrorSource(ret)))
    

    call MPI_FINALIZE(ierr)

!    deallocate(potentials,field,virial)

end program test
