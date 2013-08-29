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

!>
!> Contains all routines that should be callable from the frontend
!> in most cases
!>
module module_pepc
    use module_tree, only: t_tree
    use module_pepc_types
    implicit none
    private

    public pepc_initialize                !< mandatory, once per simulation

    ! for intializing internal parameters, you might want to consider calling one of the following functions at least once per simulation
    public pepc_read_parameters
    public pepc_read_parameters_from_file_name
    public pepc_read_parameters_from_first_argument

    public pepc_prepare                   !< mandatory, once per simulation or after changing internal parameters

    public pepc_particleresults_clear     ! usually once per timestep

    public pepc_grow_tree                 !< mandatory, once per timestep
    public pepc_traverse_tree             !< mandatory, several times per timestep with different particles possible
    public pepc_statistics                !< once or never per timestep
    public pepc_check_sanity              !< as often as necessary
    public pepc_restore_particles         !< once or never per timestep
    public pepc_timber_tree               !< mandatory, once per timestep

    public pepc_grow_and_traverse         !< once per timestep, calls pepc_grow_tree, pepc_traverse_tree, pepc_statistics, pepc_restore_particles, pepc_timber_tree

    public pepc_finalize                  !< mandatory, once per simulation

    public pepc_get_para_file
    public pepc_write_parameters

    logical :: pepc_initializes_mpi !< is set to .true., if pepc has to care for MPI_INIT and MPI_FINALIZE; otherwise, the frontend must care for that
    type(t_tree), public, target, save :: global_tree

    contains

    !>
    !> Initializes MPI library and data structures for treecode kernel,
    !> reads several parameters from file, that is given as first parameter
    !> to actual executable, initializes submodules
    !> 
    !> Call this function at program startup before any MPI calls
    !>
    subroutine pepc_initialize(frontendname, my_rank, n_cpu, init_mpi, db_level_in, comm)
      use treevars, only : np_mult, me, num_pe, MPI_COMM_lpepc, main_thread_processor_id
      use module_pepc_types, only : register_lpepc_mpi_types
      use module_utils, only : create_directory, MPI_IN_PLACE_test
      use module_walk
      use module_domains
      use module_debug
      use pthreads_stuff, only: get_my_core
      implicit none
      include 'mpif.h'
      character(*), intent(in) :: frontendname !< name of the program that uses the treecode (only for output purposes)
      integer(kind_pe), intent(out) :: my_rank !< MPI rank of this instance as returned from MPI
      integer(kind_pe), intent(out) :: n_cpu !< number of MPI ranks as returned from MPI
      logical, intent(in) :: init_mpi !< if set to .true., if pepc has to care for MPI_INIT and MPI_FINALIZE; otherwise, the frontend must care for that
      integer, intent(in), optional :: db_level_in !< sets debug level for treecode kernel (overrides settings, that may be read from libpepc-section in input file)
      integer, intent(inout), optional :: comm !< communicator. if pepc initializes MPI, it returns an MPI_COMM_DUP-copy of its own communicator (the frontend is responsible for calling MPI_COMM_FREE(comm) prior to calling pepc_finalize() or supply it as argument to pepc_finalize() in this case); otherwise, it uses an MPI_COMM_DUP copy of the given comm
      integer(kind_default) :: ierr, provided

      ! "Multiple threads may call MPI, with no restrictions." - MPI-2.2, p. 385
      integer(kind_default), parameter :: MPI_THREAD_LEVEL = MPI_THREAD_MULTIPLE

      call pepc_status('SETUP')

      pepc_initializes_mpi = init_mpi

      if (pepc_initializes_mpi) then
        ! Initialize the MPI system (thread safe version, will fallback automatically if thread safety cannot be guaranteed)
        call MPI_INIT_THREAD(MPI_THREAD_LEVEL, provided, ierr)
        call MPI_COMM_DUP(MPI_COMM_WORLD, MPI_COMM_lpepc, ierr)
        if (present(comm)) call MPI_COMM_DUP(MPI_COMM_lpepc, comm, ierr)
      else
        ! check if MPI was initialized with sufficient thread support
        call MPI_QUERY_THREAD(provided, ierr)
      
        if (present(comm)) then
           call MPI_COMM_DUP(comm, MPI_COMM_lpepc, ierr)
        else
          call MPI_COMM_DUP(MPI_COMM_WORLD, MPI_COMM_lpepc, ierr)
        endif
      endif

      call MPI_IN_PLACE_test(MPI_COMM_lpepc)

      ! Get the id number of the current task
      call MPI_COMM_RANK(MPI_COMM_lpepc, my_rank, ierr)
      ! Get the number of MPI tasks
      call MPI_COMM_SIZE(MPI_COMM_lpepc, n_cpu, ierr)
      ! Get the processor ID of the main thread
      main_thread_processor_id = get_my_core()

      call MPI_FILE_SET_ERRHANDLER(MPI_FILE_NULL, MPI_ERRORS_ARE_FATAL, ierr)

      if (my_rank == 0 .and. pepc_initializes_mpi) then
        ! verbose startup-output
        write(*,'(a)') "   ____    ____    ____    ____        "
        write(*,'(a)') "  /\  _`\ /\  _`\ /\  _`\ /\  _`\      "
        write(*,'(a)') "  \ \ \L\ \ \ \L\_\ \ \L\ \ \ \/\_\      The Pretty Efficient"
        write(*,'(a)') "   \ \ ,__/\ \  _\L\ \ ,__/\ \ \/_/_           Parallel Coulomb Solver"
        write(*,'(a)') "    \ \ \/  \ \ \L\ \ \ \/  \ \ \L\ \  "
        write(*,'(a)') "     \ \_\   \ \____/\ \_\   \ \____/           pepc@fz-juelich.de"
        write(*,'(a)') "      \/_/    \/___/  \/_/    \/___/   "
        write(*,'(/"Starting PEPC, svn revision [",a,"] with frontend {", a, "} on ", I0, " MPI ranks."//)') &
                       SVNREVISION, frontendname, n_cpu
      endif

      if (my_rank == 0 .and. provided < MPI_THREAD_LEVEL) then 
        !inform the user about possible issues concerning MPI thread safety
        if (pepc_initializes_mpi) then
          write(*,'("Call to MPI_INIT_THREAD failed. Requested/provided level of multithreading:", I2, "/" ,I2)') &
                         MPI_THREAD_LEVEL, provided
          write(*,'(a/)') 'Initialized with provided level of multithreading. This can lead to incorrect results or crashes.'
        else
          write(*,'("Frontend application did not call to MPI_INIT_THREAD correctly. Needed/provided level of multithreading:", I2, "/" ,I2)') &
                         MPI_THREAD_LEVEL, provided
          write(*,'(a/)') 'Trying to run with provided level of multithreading. This can lead to incorrect results or crashes.'
        endif
      end if

      ! copy call parameters to treevars module
      me     = my_rank
      num_pe = n_cpu

      if (present(db_level_in)) then
          debug_level = db_level_in
      else
          debug_level = 0
      endif
      np_mult         = -45.
      weighted        =   1

      ! create and register mpi types
      call register_lpepc_mpi_types()

      call pepc_prepare()
    end subroutine


    !>
    !> Initializes internal variables by reading them from a file
    !> that is given as first argument to the actual executable
    !> the optional parameters file_available and filename return
    !> whether and which file was used
    !>
    subroutine pepc_read_parameters_from_first_argument(file_available, filename)
      use module_debug, only : pepc_status
      use treevars, only : me
      implicit none
      character(*), optional :: filename
      logical, optional :: file_available

      integer, parameter :: para_file_id = 10
      character(len=255) :: para_file_name
      logical :: para_file_available

      ! read in parameter file
      call pepc_get_para_file(para_file_available, para_file_name, me)

      if (present(file_available)) file_available = para_file_available

      if (para_file_available) then
        call pepc_status("INIT FROM FILE "//para_file_name)
        call pepc_read_parameters_from_file_name(para_file_name)
        if (present(filename)) filename = para_file_name
      else
        call pepc_status("INIT WITH DEFAULT PARAMETERS")
        if (present(filename)) filename = ''
      end if
    end subroutine


    !>
    !> Initializes internal variables by reading them from the given
    !> filename using several namelists and initializes derived variables
    !> if you just want to pass a handle of an open file,
    !> consider using pepc_read_parameters()
    !>
    subroutine pepc_read_parameters_from_file_name(filename)
      implicit none
      character(*), intent(in) :: filename
      integer, parameter :: para_file_id = 10

      open(para_file_id,file=trim(filename),action='read')
      call pepc_read_parameters(para_file_id)
      close(para_file_id)
    end subroutine


    !>
    !> Initializes internal variables by reading them from the given
    !> filehandle using several namelists and initializes derived variables
    !> if you just want to pass a filename, consider using pepc_read_parameters_from_file_name()
    !>
    subroutine pepc_read_parameters(filehandle)
      use module_debug, only : pepc_status
      use module_interaction_specific, only : calc_force_read_parameters
      use module_walk, only: tree_walk_read_parameters
      use module_libpepc_main, only: libpepc_read_parameters
      implicit none
      integer, intent(in) :: filehandle

      call pepc_status("READ PARAMETERS")
      rewind(filehandle)
      call libpepc_read_parameters(filehandle)
      rewind(filehandle)
      call calc_force_read_parameters(filehandle)
      rewind(filehandle)
      call tree_walk_read_parameters(filehandle)
    end subroutine


    !>
    !> Writes internal variables to the given filehandle
    !>
    subroutine pepc_write_parameters(filehandle)
      use module_debug, only : pepc_status
      use module_interaction_specific, only : calc_force_write_parameters
      use module_walk, only: tree_walk_write_parameters
      use module_libpepc_main, only: libpepc_write_parameters
      implicit none
      integer, intent(in) :: filehandle

      call pepc_status("WRITE PARAMETERS")
      call calc_force_write_parameters(filehandle)
      call tree_walk_write_parameters(filehandle)
      call libpepc_write_parameters(filehandle)
    end subroutine


    !>
    !> Initializes derived variables in pepc, walk, and calc_force
    !> should be called after changing those module`s variables and
    !> before performing tree buildup
    !>
    subroutine pepc_prepare(idim)
      use treevars, only : treevars_prepare
      use module_walk, only: tree_walk_prepare
      use module_interaction_specific, only: calc_force_prepare
      use module_mirror_boxes, only: calc_neighbour_boxes
      use pthreads_stuff, only: pthreads_init, set_prefetching
      use module_tree_communicator, only: tree_communicator_prepare
      use module_debug
      implicit none

      integer(kind_dim), optional, intent(in) :: idim

      ERROR_ON_FAIL(pthreads_init())
      ERROR_ON_FAIL(set_prefetching())
      call treevars_prepare(idim)
      call calc_neighbour_boxes() ! initialize mirror boxes
      call calc_force_prepare() ! prepare interaction-specific routines
      call tree_walk_prepare()
      call tree_communicator_prepare()
    end subroutine


    !>
    !> Finalizes MPI library and reverses all initialization from pepc_initialize
    !> Call this function at program termination after all MPI calls
    !>
    subroutine pepc_finalize(comm)
      use module_debug
      use module_pepc_types, only : free_lpepc_mpi_types
      use module_walk, only : tree_walk_finalize 
      use module_interaction_specific, only : calc_force_finalize
      use treevars, only : treevars_finalize, MPI_COMM_lpepc
      use pthreads_stuff, only: pthreads_uninit
      use module_tree_communicator, only: tree_communicator_finalize
      implicit none
      include 'mpif.h'
      integer(kind_default) :: ierr

      integer, intent(inout), optional :: comm !< communicator. if pepc_initialize() initializes MPI, it returns an MPI_COMM_DUP-copy of its own communicator in comm, that can be given here to be freed automatically

      call pepc_status('FINALIZE')
      ! finalize internal data structures
      call calc_force_finalize()
      call tree_communicator_finalize()
      call tree_walk_finalize()
      ! deregister mpi types
      call free_lpepc_mpi_types()

      call treevars_finalize()
      if (0 /= pthreads_uninit()) then
        DEBUG_INFO(*, "pthreads_uninit() failed!")
      end if

      call MPI_COMM_FREE(MPI_COMM_lpepc, ierr)
      if (pepc_initializes_mpi) then
        if (present(comm)) then; call MPI_COMM_FREE(comm, ierr); end if 
        call MPI_FINALIZE(ierr)
      end if
    end subroutine


    !>
    !> checks if the first application argument was set
    !> broadcasts the filename to all mpi processes
    !>
    subroutine pepc_get_para_file(available, file_name, my_rank, comm)
        use treevars, only : MPI_COMM_lpepc
        implicit none
        include 'mpif.h'

        logical, intent(out)          :: available
        character(len=255), intent(out) :: file_name
        integer(kind_pe), intent(in)           :: my_rank
        integer, intent(in), optional :: comm

        integer :: ierr, MPI_COMM_local

        if (present(comm)) then
            MPI_COMM_local = comm
        else
            MPI_COMM_local = MPI_COMM_lpepc
        end if

        ! rank 0 reads in first command line argument
        available = .false.
        if (my_rank .eq. 0) then
            if( COMMAND_ARGUMENT_COUNT() .ne. 0 ) then
                call GET_COMMAND_ARGUMENT(1, file_name)
                available = .true.
                if(my_rank .eq. 0) write(*,*) "found parameter file: ", file_name
            end if
        end if

        call MPI_BCAST( available, 1, MPI_LOGICAL, 0, MPI_COMM_local, ierr )

        ! broadcast file name, read actual inputs from namelist file
        if (available) then
            call MPI_BCAST( file_name, 255, MPI_CHARACTER, 0, MPI_COMM_local, ierr )
        end if
    end subroutine pepc_get_para_file


    !>
    !> Builds the tree from the given particles, redistributes particles
    !> to other MPI ranks if necessary (i.e. reallocates particles and changes size(particles))
    !>
    subroutine pepc_grow_and_traverse(particles, itime, no_dealloc, no_restore)
      use module_pepc_types, only: t_particle
      use module_debug
      use module_tree_communicator, only : tree_communicator_stop
      implicit none
      type(t_particle), allocatable, intent(inout) :: particles(:) !< input particle data, initializes %x, %data, %work appropriately (and optionally set %label) before calling this function
      integer, intent(in) :: itime !> current timestep (used as filename suffix for statistics output)
      logical, optional, intent(in) :: no_dealloc ! if .true., the internal data structures are not deallocated (e.g. for a-posteriori diagnostics)
      logical, optional, intent(in) :: no_restore ! if .true., the particles are not backsorted to their pre-domain-decomposition order

      logical :: restore, dealloc

      restore = .true.
      dealloc = .true.

      if (present(no_dealloc)) dealloc = .not. no_dealloc
      if (present(no_restore)) restore = .not. no_restore

      call pepc_grow_tree(particles)
      call pepc_traverse_tree(particles)

      if (dbg(DBG_STATS)) call pepc_statistics(itime)
      if (restore) then
        ! for better thread-safety we have to kill the communicator thread before trying to perform any other mpi stuff
        call tree_communicator_stop(global_tree)
        call pepc_restore_particles(particles)
      endif
      
      if (dealloc) call pepc_timber_tree()
      
    end subroutine


    !>
    !> Builds the tree from the given particles, redistributes particles
    !> to other MPI ranks if necessary (i.e. reallocates particles and changes np_local)
    !>
    subroutine pepc_grow_tree(particles)
      use module_pepc_types, only: t_particle
      use module_libpepc_main, only: libpepc_grow_tree
      implicit none
      type(t_particle), allocatable, intent(inout) :: particles(:) !< input particle data, initializes %x, %data, %work appropriately (and optionally set %label) before calling this function

      call libpepc_grow_tree(global_tree, particles)
      
    end subroutine


    !>
    !> Traverses the complete tree for the given particles, i.e. computes
    !> the field values at their positions. Although missing information
    !> is automatically requested from remote MPI ranks, it is important
    !> that the particle coordinates fit to the local MPI ranks domain
    !> to avoid excessive communication
    !> If field values on some regular grid are needed, they can be
    !> generated using pepc_prepare_local_grid() [TODO: provide this function]
    !> Otherwise, it makes sense to provide the same particles as given/returned
    !> from to pepc_grow_tree()
    !>
    subroutine pepc_traverse_tree(particles)
      use module_pepc_types
      use module_libpepc_main
      implicit none

      type(t_particle), target, intent(inout) :: particles(:) !< input particle data, initialize %x, %data, %work appropriately (and optionally set %label) before calling this function
      
      call libpepc_traverse_tree(global_tree, particles)

    end subroutine


    !>
    !> Writes detailed statistics on the treecode into stats/stats.ITIME.
    !>
    subroutine pepc_statistics(itime)
        use module_tree, only: tree_stats
        use module_walk, only: tree_walk_statistics, interactions_local, mac_evaluations_local
        use module_utils, only: create_directory
        use treevars, only: me, stats_u
        use module_debug
        use module_timings
        implicit none
        integer, intent(in) :: itime !< current timestep (used as file suffix)

        logical, save :: firstcall = .true.
        character(30) :: cfile

        call timer_start(t_fields_stats)
        if (firstcall) then
          call create_directory("stats")
          if (dbg(DBG_LOADFILE)) then
            call create_directory("load")
          end if
          firstcall = .false.
        end if

        write (cfile, '("stats/stats.",i6.6)') itime
        if (0 == me) then; open (stats_u, file = trim(cfile)); end if
        call tree_stats(global_tree, stats_u)
        call tree_walk_statistics(stats_u)
        if (0 == me) then; close (stats_u); end if

        if( dbg(DBG_LOADFILE) ) then
            write(cfile,'("load/load_",i6.6,".dat")') me
            open(60, file=trim(cfile),STATUS='UNKNOWN', POSITION = 'APPEND')
            write(60,'(i5,2f20.10, i12)') itime, interactions_local, mac_evaluations_local, global_tree%npart_me
            close(60)
        end if

        call timer_stop(t_fields_stats)
    end subroutine


    !>
    !> Checks the internal state of PEPC and dumps the hash table if requested.
    !>
    subroutine pepc_check_sanity(caller, dump, particles)
      use module_pepc_types, only: t_particle
      use module_tree, only: tree_check, tree_dump
      use module_debug
      implicit none

      character(*), intent(in) :: caller !< describes the caller
      logical, optional, intent(in) :: dump !< whether to dump the hash table
      type(t_particle), optional, intent(in) :: particles(:) !< list of particles to dump along with the hash table

      if (.not. tree_check(global_tree, caller)) then
        call tree_dump(global_tree, particles)
        DEBUG_ERROR(*, "Sanity check failed, aborting!")
      else if (dump) then
        call tree_dump(global_tree, particles)
      end if
    end subroutine pepc_check_sanity


    !>
    !> Restores the initial particle distribution (before calling pepc_grow_tree() ).
    !>
    subroutine pepc_restore_particles(particles)
      use module_pepc_types
      use module_libpepc_main
      implicit none
      type(t_particle), allocatable, intent(inout) :: particles(:) !< input particle data on local MPI rank - is replaced by original particle data that was given before calling pepc_grow_tree()

      call libpepc_restore_particles(global_tree, particles)

    end subroutine


    !>
    !> Frees all tree specific data fields that were allocated in pepc_grow_tree().
    !>
    subroutine pepc_timber_tree()
      use module_libpepc_main, only: libpepc_timber_tree
      implicit none

      call libpepc_timber_tree(global_tree)
    end subroutine


    !>
    !> clears result in t_particle datatype
    !>
    subroutine pepc_particleresults_clear(particles)
      use module_pepc_types
      use module_interaction_specific
      implicit none
      type(t_particle), intent(inout) :: particles(:)
      
      call particleresults_clear(particles)
    end subroutine

end module module_pepc
