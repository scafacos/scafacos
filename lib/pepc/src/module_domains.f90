! This file is part of PEPC - The Pretty Efficient Parallel Coulomb Solver.
! 
! Copyright (C) 2002-2014 Juelich Supercomputing Centre, 
!                         Forschungszentrum Juelich GmbH,
!                         Germany
! 
! PEPC is free software: you can redistribute it and/or modify
! it under the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
! 
! PEPC is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
! 
! You should have received a copy of the GNU Lesser General Public License
! along with PEPC.  If not, see <http://www.gnu.org/licenses/>.
!

!>
!> Encapsulates domain decomposition and restoration of original particle order
!>
module module_domains
  use module_comm_env, only: t_comm_env
  use module_pepc_types
  implicit none
  save
  private

  type, public :: t_decomposition
    integer(kind_particle) :: npold  !< original particle number
    integer(kind_particle) :: npnew  !< particle number after domain decomposition
    integer(kind_particle) :: nppmax !< maximum number of local particles after domain decomposition
    !> these have to be kind_default as MPI_ALLTOALLV expects default integer kind arguments
    integer(kind_default), allocatable :: indxl(:)  !< permutes `particles(:)` into the MPI buffer to be use for shipping
    integer(kind_default), allocatable :: irnkl(:)  !< permutes receiving MPI buffer into `particles(:)`
    integer(kind_default), allocatable :: fposts(:) !< send displacements
    integer(kind_default), allocatable :: gposts(:) !< receive displacements
    integer(kind_default), allocatable :: islen(:)  !< send counts
    integer(kind_default), allocatable :: irlen(:)  !< receive counts

    type(t_comm_env) :: comm_env !< the communication topology over which the domain is decomposed
  end type t_decomposition

  integer(kind_default), public :: weighted = 1 !< set to 0 to disable load balancing, 1 to enable load balancing

  public decomposition_destroy
  public decomposition_allocated
  public domain_decompose
  public domain_restore

  contains

  subroutine decomposition_create(d, nl, n, c)
    use module_comm_env, only: comm_env_dup
    use module_debug
    implicit none

    type(t_decomposition), intent(inout) :: d
    integer(kind_particle), intent(in) :: nl
    integer(kind_particle), intent(in) :: n
    type(t_comm_env), intent(in) :: c

    call comm_env_dup(c, d%comm_env)
    
    !TODO: make adjustable by user or find a good estimation. Additional Question: Does this value have to be globally constant?
    ! allow 25% fluctuation around average particle number per PE in sorting library for load balancing
    d%nppmax = int(1.25 * max(int(n / d%comm_env%size), 1000))

    if (d%nppmax .lt. nl) then
      DEBUG_WARNING_ALL(*, 'nppmax = ', d%nppmax, 'is smaller than np_local = ', nl, ' - fixing and continuing. This could lead to load balancing issues. See ticket no. 10')
      d%nppmax = max(d%nppmax, nl)
    end if

    ! fields for sorting library results
    allocate(d%indxl(d%nppmax), d%irnkl(d%nppmax), d%fposts(d%comm_env%size + 1), &
      d%gposts(d%comm_env%size + 1), d%islen(d%comm_env%size), d%irlen(d%comm_env%size))
  end subroutine decomposition_create


  subroutine decomposition_destroy(d)
    use module_comm_env, only: comm_env_destroy
    implicit none

    type(t_decomposition), intent(inout) :: d

    d%npnew = 0
    d%npold = 0
    d%nppmax = 0

    deallocate(d%indxl, d%irnkl, d%fposts, d%gposts, d%islen, d%irlen)

    call comm_env_destroy(d%comm_env)
  end subroutine decomposition_destroy


  logical function decomposition_allocated(d)
    implicit none

    type(t_decomposition), intent(in) :: d

    decomposition_allocated = allocated(d%indxl) .or. allocated(d%irnkl) .or. allocated(d%fposts) .or. &
      allocated(d%gposts) .or. allocated(d%islen) .or. allocated(d%irlen)

  end function decomposition_allocated


  !>
  !>  Domain decomposition:
  !>  Share particle keys amoung PEs
  !>  - weighting according to load incurred on previous timestep
  !>
  subroutine domain_decompose(d, b, n, particles, bp, c)
    use module_pepc_types, only: t_particle, mpi_type_particle
    use module_box, only: t_box
    use module_timings
    use module_spacefilling, only: key_to_coord
    use module_debug
    implicit none
    include 'mpif.h'

    type(t_decomposition), intent(inout) :: d
    type(t_box), intent(in) :: b
    integer(kind_particle), intent(in) :: n
    type(t_particle), allocatable, intent(inout) :: particles(:)
    type(t_particle), intent(out) :: bp(2)
    type(t_comm_env), intent(in) :: c

    integer(kind_default) :: ierr
    integer(kind_particle) :: i, j
    real*8 :: imba
    type(t_particle), allocatable :: ship_parts(:), get_parts(:) !< arrays for parallel sort
    integer(kind_key), allocatable :: temp(:), local_keys(:), key_diffs(:)
    real*8, allocatable :: workload(:)
    integer(kind_default), allocatable :: irnkl2(:)

    interface
      subroutine slsort_keys(nin, nmax, keys, workload, balance_weight, max_imbalance, nout, indxl, &
        irnkl, scounts, rcounts, sdispls, rdispls, keys2, irnkl2, size, rank, comm)
        use module_pepc_types
        integer(kind_particle), intent(in) :: nin
        integer(kind_particle), intent(in) :: nmax
        integer(kind_key), intent(inout) :: keys(*)
        real*8,intent(inout) :: workload(*)
        integer(kind_default), intent(in) :: balance_weight
        real*8,intent(in) :: max_imbalance
        integer(kind_particle), intent(out) :: nout
        integer(kind_default), intent(out) :: indxl(*), irnkl(*), scounts(*), rcounts(*), sdispls(*), rdispls(*)
        integer(kind_key), intent(out) :: keys2(*)
        integer(kind_default), intent(out) :: irnkl2(*)
        integer(kind_pe), intent(in) :: size, rank
        integer(kind_default), intent(in) :: comm
      end subroutine slsort_keys
    end interface

    call pepc_status('DOMAIN DECOMPOSITION')
    call timer_start(t_domains)

    d%npold = size(particles, 1)

    imba = 0.01

    ! workload per particle must be nonzero
    do i = 1, d%npold
      particles(i)%work = max(particles(i)%work, 1._8)
    end do

    if (dbg(DBG_DOMAIN)) call print_particle_list(particles, d%npold, &
                                  'Particle list before key sort (see t_particle in module_pepc_types.f90 for meaning of the columns)')

    call decomposition_create(d, d%npold, n, c)
    call timer_start(t_domains_sort)
    call timer_start(t_domains_add_sort)
    call timer_start(t_domains_sort_pure)

    allocate(temp(d%nppmax), workload(d%nppmax), irnkl2(d%nppmax), local_keys(d%nppmax))

    ! start permutation of local key list
    workload(1:d%npold) = particles(1:d%npold)%work
    local_keys(1:d%npold) = particles(1:d%npold)%key

    ! perform index sort on keys
    call slsort_keys(d%npold, d%nppmax, local_keys, workload, weighted, imba, d%npnew, d%indxl, d%irnkl, &
      d%islen, d%irlen, d%fposts, d%gposts, temp, irnkl2, d%comm_env%size, d%comm_env%rank , &
      d%comm_env%comm)

    ! FIXME: every processor has to have at least one particle
    if (d%npnew < 2) then
        DEBUG_ERROR('("rank has less than two particles after sorting (had ", I8, " before). Did you initialise particle field %work correctly? --> aborting")', d%npold)
    end if

    call timer_stop(t_domains_sort_pure)

    call timer_stop(t_domains_sort)
    call timer_start(t_domains_ship)

    ! Now permute particle properties
    ! Set up particle structure
    call timer_start(t_domains_add_pack)
    
    allocate(ship_parts(d%nppmax))
    do i = 1, d%npold
      ship_parts(i) = particles( d%indxl(i) )
    end do

    call timer_stop(t_domains_add_pack)

    deallocate(particles) ! has size npold until here, i.e. npp == npold
    allocate(get_parts(d%nppmax))

    call timer_start(t_domains_add_alltoallv)

    ! perform permute
    call MPI_ALLTOALLV(ship_parts, d%islen, d%fposts, mpi_type_particle, &
      get_parts, d%irlen, d%gposts, MPI_TYPE_particle, d%comm_env%comm, ierr)

    call timer_stop(t_domains_add_alltoallv)

    deallocate(ship_parts)
    allocate(particles(d%npnew))

    call timer_start(t_domains_add_unpack)

    do i = 1, d%npnew
      particles(d%irnkl(i)) = get_parts(i)
    end do
    deallocate(get_parts)

    call timer_stop(t_domains_add_unpack)

    if (d%npnew > d%nppmax) then
      DEBUG_ERROR('("More than nppm particles after sorting: nppm = ", I0, " < npp = ",I0,". All local particle fields are too short. Aborting.")', d%nppmax, d%npnew)
    end if

    call timer_stop(t_domains_ship)
    call timer_stop(t_domains_add_sort)
    call timer_start(t_domains_bound)

    if (dbg(DBG_DOMAIN)) call print_particle_list(particles, d%npnew, &
      'Particle list after key sort (see t_particle in module_pepc_types.f90 for meaning of the columns)')

    ! Each PE now has sorted segment of particles of length npp
    ! Note that now npp /= npart/num_pe, only approx depending on key distribution, or target shape.

    ! check for duplicate keys

    ! first, we exchange boundary particles to avoid duplicate keys across processor boundaries in the following
    call exchange_boundary_particles()

    ! we will first work on a copy of the original keys to be able to only adjust coordinates
    ! of particles where the key really has changed (instead of being shifted up and down only)
    local_keys(1:d%npnew) = particles(:)%key

    ! check whether there is enough space in the local key domain
    if ( ( local_keys(d%npnew) - local_keys(1) + 1) < d%npnew) then
      DEBUG_WARNING_ALL('("There are more particles than available keys in the local domain: npp=",I0,", but upper and lower key (octal) =", 2(x,O0),"; upper-lower (dec) = ",I0,". Usually, particle coordinates are invalid in this case. Printing particle list to diag file.")', d%npnew, local_keys(d%npnew), local_keys(1), local_keys(d%npnew) - local_keys(1))
      DEBUG_INFO(*, '1st column: KEY, folowing columns: t_particle structure')
      do i=1,d%npnew
        DEBUG_INFO(*, local_keys(i), particles(i))
      end do
      call debug_mpi_abort()
    end if

    ! key differences for identifiying duplicates and/or overlap, be aware of the problem, that for me==num_pe-1, key_diffs(npp) is invalid since there is no right neighbour
    allocate(key_diffs(d%npnew))
    key_diffs(1:d%npnew - 1) = local_keys(2:d%npnew) - local_keys(1:d%npnew - 1) ! key_diffs(i) = local_keys(i+1) - local_keys(i)
    if (d%comm_env%last) then
      key_diffs(d%npnew) = 1
    else
      key_diffs(d%npnew) = bp(2)%key - local_keys(d%npnew)
    end if

    do i = 1, d%npnew
      if (key_diffs(i) < 1) then
        DEBUG_INFO('("Identical keys found for i = ", I0, " of ", I0, " and its successor, key = ", O0, ", label(i) = ", I0)', i, d%npnew, local_keys(i), particles(i)%label)

        ! looking upwards and downwards synchronously, we try to find the nearest gap
        do j = 1, max(d%npnew - i, i)
          if (key_diffs(max(1_kind_particle, i - j)) > 1) then
            ! there is a near gap below the current keys --> shift keys downwards
            DEBUG_INFO('("Fixing by shifting keys of section ", I0,":",I0," downwards")', i-j+1, i)
            local_keys(i-j+1:i) = local_keys(i-j+1:i) - 1
            key_diffs(i)        = key_diffs(i)        + 1
            key_diffs(i-j)      = key_diffs(i-j)      - 1
            exit ! from inner loop
          else if (key_diffs(min(d%npnew - 1, i + j)) > 1) then
            ! there is a near gap above the current keys --> shift keys upwards
            DEBUG_INFO('("Fixing by shifting keys of section ", I0,":",I0," upwards")', i + 1, i + j)
            local_keys(i+1:i+j) = local_keys(i+1:i+j) + 1
            key_diffs(i)        = key_diffs(i)        + 1
            key_diffs(i+j)      = key_diffs(i+j)      - 1
            exit ! from inner loop
          end if
        end do
      end if
    end do

    ! adjust particle coordinates to new keys if necessary
    do i = 1, d%npnew
      if (local_keys(i) .ne. particles(i)%key) then
        particles(i)%key = local_keys(i)
        call key_to_coord(b, particles(i)%key, particles(i)%x)
      end if
    end do

    ! since we possibly modified the key of particles(npp), we repeat the boundary exchange
    call exchange_boundary_particles()
    if (dbg(DBG_DOMAIN)) call print_particle_list(particles, d%npnew, &
                                  'Particle list after boundary swap (see t_particle in module_pepc_types.f90 for meaning of the columns)')

    call timer_stop(t_domains_bound)
    deallocate(temp, workload, irnkl2, local_keys, key_diffs)

    ! reset work load, do not need a (overflowing) history of all my sins...
    particles(:)%work = 1.
    call timer_stop(t_domains)

    contains
    !>
    !> Copy boundary particles to adjacent PEs to ensure proper tree construction
    !>  - if we do not do this, can get two particles on separate PEs 'sharing' a leaf
    !> after calling this routine,
    !> if (me /= 0)
    !>   bp(1) = particles[me - 1](npp)
    !> if (me /= num_pe)
    !>   bp(2) = particles[me + 1](1)
    !>
    subroutine exchange_boundary_particles()
      implicit none
      include 'mpif.h'

      integer :: state(MPI_STATUS_SIZE)
      integer(kind_particle) :: npp
      integer(kind_pe) :: prev, next

      npp = d%npnew

      ! Define neighbours for non-circular shift
      if (d%comm_env%first) then
        prev = MPI_PROC_NULL
      else
        prev = d%comm_env%rank - 1_kind_pe
      end if
      if (d%comm_env%last) then
        next = MPI_PROC_NULL
      else
        next = d%comm_env%rank + 1_kind_pe
      end if

      ! Copy boundary particles to adjacent PEs to ensure proper tree construction
      !  - if we do not do this, can get two particles on separate PEs 'sharing' a leaf

      ! Ship 1st particle data to LH neighbour PE
      call MPI_SENDRECV(particles(1), 1, mpi_type_particle, prev, 1990, &
        bp(2), 1, mpi_type_particle, next, 1990, &
        d%comm_env%comm, state, ierr)

      ! Ship end particle data to RH neighbour PE
      call MPI_SENDRECV(particles(npp), 1, mpi_type_particle, next, 2990, &
        bp(1), 1, mpi_type_particle, prev, 2990, &
        d%comm_env%comm, state, ierr)

    end subroutine exchange_boundary_particles
  end subroutine domain_decompose


  !>
  !>  Restore initial particle order
  !>
  subroutine domain_restore(d, p)
      use module_pepc_types, only: t_particle, mpi_type_particle
      use module_debug, only : pepc_status
      implicit none
      include 'mpif.h'

      type(t_decomposition), intent(inout) :: d
      type(t_particle), intent(inout), allocatable :: p(:)

      integer(kind_particle) :: i
      integer(kind_default) :: ierr
      type (t_particle), allocatable :: get_parts(:), ship_parts(:)

      call pepc_status('RESTORE DOMAINS')

      allocate(ship_parts(d%npnew))
      do i = 1, d%npnew
        ship_parts(i) = p(d%irnkl(i))
      end do

      deallocate(p) ! had size npnew
      allocate(get_parts(d%npold))

      ! perform permute
      call MPI_ALLTOALLV(ship_parts, d%irlen, d%gposts, MPI_TYPE_particle, &
            get_parts, d%islen, d%fposts, MPI_TYPE_particle, &
            d%comm_env%comm, ierr)

      deallocate(ship_parts)
      allocate(p(d%npold))
      do i = 1, d%npold
          p(d%indxl(i)) = get_parts(i)
      end do

      deallocate(get_parts)
      call decomposition_destroy(d)
  end subroutine domain_restore


  subroutine print_particle_list(particles, npart, callinfo)
    use module_pepc_types, only: t_particle
    use module_debug
    implicit none
    type(t_particle), intent(in) :: particles(:)
    integer(kind_particle), intent(in) :: npart
    character(*), intent(in) :: callinfo

    integer(kind_particle) :: j

    call debug_ipefile_open()
    write (debug_ipefile,'(/a/)') callinfo
    do j=1,npart
      write(debug_ipefile,'(i10)',advance='no') j
      write(debug_ipefile,*)                     particles(j)
    end do
    write(debug_ipefile,'(/)')
    call debug_ipefile_close()
  end subroutine
end module module_domains
