! This file is part of PEPC - The Pretty Efficient Parallel Coulomb Solver.
! 
! Copyright (C) 2002-2013 Juelich Supercomputing Centre, 
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!> Encapsulates domain decomposition and restoration of original particle order
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module module_domains
    implicit none
    save
    private

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!  public variable declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    integer, public :: weighted = 1 !< set to 0 to disable load balancing, 1 to enable load balancing
    logical, public :: force_cubic_domain = .false. !< if set to .true., pepc uses an overall cubic enclosure of the particle cloud instead of the cuboid (closer) one



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!  public subroutine declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    public tree_domains
    public restore


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!  private variable declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    contains

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !>  Domain decomposition:
    !>  Share particle keys amoung PEs
    !>  - weighting according to load incurred on previous timestep
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine tree_domains(particles, nppm, indxl,irnkl,islen,irlen,fposts,gposts,npnew,npold,neighbour_pe_particles)

        use treevars
        use module_interaction_specific
        use module_utils
        use module_timings
        use module_spacefilling
        use module_branching
        use module_debug
        implicit none
        include 'mpif.h'

        type(t_particle), allocatable, intent(inout) :: particles(:)

        integer, intent(in) :: nppm !< maximum allowed number of particles on this PE
        integer, intent(out) :: indxl(nppm),irnkl(nppm)
        integer, intent(out) :: islen(num_pe),irlen(num_pe)
        integer, intent(out) :: fposts(num_pe+1),gposts(num_pe+1)
        integer :: npnew,npold
        integer, intent(out) :: neighbour_pe_particles !< number of particles that have been fetched from neighbouring PEs - they are stored in particles(npp+1:npp+neighbour_pe_particles)

        integer :: i, j
        real*8 :: min_local(3), max_local(3)
        integer :: ierr
        type (t_particle) :: ship_parts(nppm), get_parts(nppm) !< arrays for parallel sort
        integer*8 :: w1(nppm)
        real*8 imba

        interface
            subroutine slsort_keys(nin,nmax,keys,workload,balance_weight,max_imbalance,nout,indxl,irnkl,scounts,rcounts,sdispls,rdispls,keys2,irnkl2,size,rank,comm)
                integer,intent(in) :: nin,nmax,balance_weight,size,rank,comm
                real*8,intent(in) :: max_imbalance
                integer,intent(out) :: nout,indxl(*),irnkl(*),scounts(*),rcounts(*),sdispls(*),rdispls(*),irnkl2(*)
                integer*8,intent(out) :: keys2(*)
                integer*8,intent(inout) :: keys(*)
                real*8,intent(inout) :: workload(*)
            end subroutine slsort_keys
        end interface

        real*8 :: work2(nppm)
        integer :: irnkl2(nppm)
        integer*8 :: local_keys(nppm)
        integer*8, allocatable :: key_diffs(:)

        call timer_start(t_domains)
        call timer_start(t_domains_keys)

        call pepc_status('DOMAIN DECOMPOSITION')

        ! Find limits of local simulation region
        min_local(1) = minval(particles(1:npp)%x(1))
        max_local(1) = maxval(particles(1:npp)%x(1))
        min_local(2) = minval(particles(1:npp)%x(2))
        max_local(2) = maxval(particles(1:npp)%x(2))
        min_local(3) = minval(particles(1:npp)%x(3))
        max_local(3) = maxval(particles(1:npp)%x(3))

        if (force_cubic_domain) then
           min_local(1:3) = minval(min_local(1:3))
           max_local(1:3) = maxval(max_local(1:3))
        endif

        ! Find global limits
        call MPI_ALLREDUCE(min_local, boxmin, 3, MPI_REAL8, MPI_MIN,  MPI_COMM_lpepc, ierr )
        call MPI_ALLREDUCE(max_local, boxmax, 3, MPI_REAL8, MPI_MAX,  MPI_COMM_lpepc, ierr )

        ! Safety margin - put buffer region around particles
        boxsize = boxmax - boxmin
        boxmin  = boxmin - boxsize/10000.0
        boxmax  = boxmax + boxsize/10000.0
        boxsize = boxmax - boxmin


        if (dbg(DBG_DOMAIN)) then
          DEBUG_WARNING('(4(a15,f12.4/))',
            'xmin = ',boxmin(1),'xmax = ',boxmax(1),
            'ymin = ',boxmin(2),'ymax = ',boxmax(2),
            'zmin = ',boxmin(3),'zmax = ',boxmax(3),
            'boxsize = ',boxsize )
        endif

        call compute_particle_keys(particles)

        if (dbg(DBG_DOMAIN)) call print_particle_list(particles, npp, &
                                     'Particle list before key sort (see t_particle in module_pepc_types.f90 for meaning of the columns)')

        call timer_stop(t_domains_keys)
        call timer_start(t_domains_sort)

        imba = 0.01

        npold = npp
        npnew = npp

        call timer_start(t_domains_add_sort)

        ! start permutation of local key list
        work2(1:npp) = particles(1:npp)%work

        call timer_start(t_domains_sort_pure)

        local_keys(1:npold) = particles(1:npold)%key
        ! perform index sort on keys !TODO: remove the "-2", compare other cases with "+2" and "npp+1" etc.
        call slsort_keys(npold,nppm-2,local_keys,work2,weighted,imba,npnew,indxl,irnkl,islen,irlen,fposts,gposts,w1,irnkl2,num_pe,me,MPI_COMM_lpepc)

        ! FIXME: every processor has to have at least one particle
        if (npnew < 2) then
            DEBUG_ERROR('("rank less than two particles after sorting (had ", I8, " before) - currently this can lead to errors --> aborting")', npold)
        endif

        call timer_stop(t_domains_sort_pure)

        call timer_stop(t_domains_sort)
        call timer_start(t_domains_ship)

        ! Now permute particle properties
        ! Set up particle structure
        call timer_start(t_domains_add_pack)

        do i=1,npold
            ship_parts(i) = particles( indxl(i) )
        enddo

        call timer_stop(t_domains_add_pack)

        deallocate(particles) ! has size npold until here, i.e. npp == npold

        call timer_start(t_domains_add_alltoallv)

        ! perform permute
        call MPI_alltoallv(  ship_parts, islen, fposts, mpi_type_particle, &
        get_parts, irlen, gposts, mpi_type_particle, &
        MPI_COMM_lpepc,ierr)

        call timer_stop(t_domains_add_alltoallv)

        allocate(particles(npnew+2)) ! TODO: the limit particles from neighbouring PEs are put into the final two places - this is for branching and correct insertion into the tree and should be done there with local variables instead
        npp = npnew

        call timer_start(t_domains_add_unpack)

        do i=1,npp
            particles( irnkl(i) ) = get_parts(i)
        enddo

        call timer_stop(t_domains_add_unpack)

        if (npp > nppm) then
            DEBUG_ERROR('("More than nppm particles after sorting: nppm = ", I0, " < npp = ",I0,". All local particle fields are too shirt. Aborting.")', nppm, npp)
        endif

        call timer_stop(t_domains_ship)
        call timer_stop(t_domains_add_sort)
        call timer_start(t_domains_bound)

        if (dbg(DBG_DOMAIN)) call print_particle_list(particles, npp, &
                                     'Particle list after key sort (see t_particle in module_pepc_types.f90 for meaning of the columns)')

        particles(1:npp)%pid = me  ! new owner

        ! Each PE now has sorted segment of particles of length npp
        ! Note that now npp /= npart/num_pe, only approx depending on key distribution, or target shape.

        ! check for duplicate keys

        ! first, we exchange boundary particles to avoid duplicate keys across processor boundaries in the following
        call exchange_boundary_particles(particles, npp, neighbour_pe_particles, me, num_pe)

        ! we will first work on a copy of the original keys to be able to only adjust coordinates
        ! of particles where the key really has changed (instead of being shifted up and down only)
        local_keys(1:npp+1) = particles(1:npp+1)%key

        ! check whether there is enough space in the local key domain
        if ( ( local_keys(npp) - local_keys(1) + 1) < npp) then
          DEBUG_ERROR('("There are more particles than available keys in the local domain: npp=",I0,", but upper and lower key (octal) =", 2(x,O0),"; upper-lower (dec) = ",I0)',
                           npp, local_keys(npp), local_keys(1), local_keys(npp) - local_keys(1))
        endif

        ! key differences for identifiying duplicates and/or overlap, be aware of the problem, that for me==num_pe-1, key_diffs(npp) is invalid since there is no right neighbour
        allocate(key_diffs(npp))
        key_diffs(1:npp) = local_keys(2:npp+1) - local_keys(1:npp) ! key_diffs(i) = local_keys(i+1) - local_keys(i)
        if (me==num_pe-1) key_diffs(npp) = 1

        do i=1,npp
            if (key_diffs(i) < 1) then
                DEBUG_INFO('("Identical keys found for i = ", I0, " and its successor, key = ", O0, ", labels = ", I0,x,I0)', i, local_keys(i), particles(i)%label, particles(i+1)%label)

                ! looking upwards and downwards synchronously, we try to find the nearest gap
                do j=1,max(npp-i,i)
                  if (key_diffs(max(1,i-j)) > 1) then
                    ! there is a near gap below the current keys --> shift keys downwards
                    DEBUG_INFO('("Fixing by shifting keys of section ", I0,":",I0," downwards")', i-j+1, i)
                    local_keys(i-j+1:i) = local_keys(i-j+1:i) - 1
                    key_diffs(i)        = key_diffs(i)        + 1
                    key_diffs(i-j)      = key_diffs(i-j)      - 1
                    exit ! from inner loop
                  else if (key_diffs(min(npp-1,i+j)) > 1) then
                    ! there is a near gap above the current keys --> shift keys upwards
                    DEBUG_INFO('("Fixing by shifting keys of section ", I0,":",I0," upwards")', i+1, i+j)
                    local_keys(i+1:i+j) = local_keys(i+1:i+j) + 1
                    key_diffs(i)        = key_diffs(i)        + 1
                    key_diffs(i+j)      = key_diffs(i+j)      - 1
                    exit ! from inner loop
                  end if
                end do
            endif
        end do

        deallocate(key_diffs)

        ! adjust particle coordinates to new keys if necessary
        do i=1,npp
          if (local_keys(i) .ne. particles(i)%key) then
            particles(i)%key = local_keys(i)
            call key_to_coord(particles(i)%key, particles(i)%x)
          endif
        end do

        ! since we possibly modified the key of particles(npp), we repeat the boundary exchange
        call exchange_boundary_particles(particles, npp, neighbour_pe_particles, me, num_pe)

        ! Initialize VLD-stuff
        call branches_initialize_VLD(particles)

        if (dbg(DBG_DOMAIN)) call print_particle_list(particles, npp + neighbour_pe_particles, &
                                     'Particle list after boundary swap (see t_particle in module_pepc_types.f90 for meaning of the columns)')

        call timer_stop(t_domains_bound)
        call timer_stop(t_domains)

    end subroutine tree_domains



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !>  Copy boundary particles to adjacent PEs to ensure proper tree construction
    !>  - if we do not do this, can get two particles on separate PEs 'sharing' a leaf
    !> after calling this routine,
    !> if (me == 0)
    !>    particles[me](npp+1)   == particles[me+1](1)
    !>    neighbour_pe_particles == 1
    !>  elseif (me == num_pe)
    !>    particles[me](npp+1)   == particles[me-1](npp)
    !>    neighbour_pe_particles == 1
    !>  else
    !>    particles[me](npp+1)   == particles[me+1](1)
    !>    particles[me](npp+2)   == particles[me-1](npp)
    !>    neighbour_pe_particles == 2
    !>  endif
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine exchange_boundary_particles(particles, npp, neighbour_pe_particles, me, num_pe)
      use treevars, only : MPI_COMM_lpepc
      use module_pepc_types
      implicit none
      include 'mpif.h'
      integer, intent(out) :: neighbour_pe_particles
      integer, intent(in) :: me, num_pe, npp
      type(t_particle), intent(inout) :: particles(1:npp+2)

      integer :: prev, next, ierr, state(MPI_STATUS_SIZE)

        ! Define neighbours for non-circular shift
        prev = me - 1
        next = me + 1
        if (me == 0) then
            prev = MPI_PROC_NULL
        end if
        if (me == num_pe - 1) then
            next = MPI_PROC_NULL
        end if

        ! Copy boundary particles to adjacent PEs to ensure proper tree construction
        !  - if we do not do this, can get two particles on separate PEs 'sharing' a leaf
        neighbour_pe_particles = 0

        ! Ship 1st particle data to end of list of LH neighbour PE
        if ( me /= num_pe-1) then
            neighbour_pe_particles = neighbour_pe_particles + 1
        end if
        call MPI_SENDRECV(particles(1), 1, mpi_type_particle, prev, 1990, &
          particles(npp + neighbour_pe_particles), 1, mpi_type_particle, next, 1990, &
          MPI_COMM_lpepc, state, ierr)

        ! Ship end particle data to end of list of RH neighbour PE
        if (me /= 0) then
            neighbour_pe_particles = neighbour_pe_particles + 1
        end if
        call MPI_SENDRECV(particles(npp), 1, mpi_type_particle, next, 2990, &
          particles(npp + neighbour_pe_particles), 1, mpi_type_particle, prev, 2990, &
          MPI_COMM_lpepc, state, ierr)

    end subroutine


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !>  Restore initial particle order
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine restore(npnew,npold,nppm_ori,indxl,irnkl,islen,irlen,fposts,gposts,&
                             particles)
        use module_interaction_specific
        use treevars, only : num_pe, npp, MPI_COMM_lpepc
        use module_pepc_types
        use module_debug, only : pepc_status
        implicit none
        include 'mpif.h'

        integer, intent(in) :: npnew,npold,nppm_ori
        integer, intent(in) :: indxl(nppm_ori),irnkl(nppm_ori)
        integer, intent(in) :: islen(num_pe),irlen(num_pe)
        integer, intent(in) :: fposts(num_pe+1),gposts(num_pe+1)
        type(t_particle),         intent(inout), allocatable :: particles(:)

        integer :: i, ierr

        type (t_particle)         :: get_parts(npold), ship_parts(npnew)

        call pepc_status('RESTORE DOMAINS')

        do i=1,npnew
          ship_parts(i) = particles(indxl(i))
        enddo

        deallocate(particles) ! had size npnew

        ! perform permute
        call MPI_alltoallv(  ship_parts, islen, fposts, MPI_TYPE_particle, &
              get_parts, irlen, gposts, MPI_TYPE_particle, &
              MPI_COMM_lpepc, ierr )

        allocate(particles(npold))
        npp = npold

        do i=1,npold
            particles(irnkl(i)) = get_parts(i)
        enddo

    end subroutine restore


    subroutine print_particle_list(particles, npart, callinfo)
      use module_pepc_types
      use module_debug
      implicit none
      type(t_particle), intent(in) :: particles(:)
      integer, intent(in) :: npart
      character(*), intent(in) :: callinfo

      integer :: j

      call debug_ipefile_open()
      write (debug_ipefile,'(/a/)') callinfo
      do j=1,npart
        write(debug_ipefile,'(i10)',advance='no') j
        write(debug_ipefile,*)                     particles(j)
        write(debug_ipefile,'(/)')
      end do
      write(debug_ipefile,'(/)')
      call debug_ipefile_close()

   end subroutine

end module module_domains
