! This file is part of PEPC - The Pretty Efficient Parallel Coulomb Solver.
! 
! Copyright (C) 2002-2012 Juelich Supercomputing Centre, 
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
!> Contains all routines that should be callable from the frontend
!> in most cases
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module module_pepc
    use module_debug, only : debug_level
    use treevars, only : np_mult, interaction_list_length_factor
    use module_spacefilling, only : curve_type
    use module_domains, only : weighted
    implicit none
    private

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!  public variable declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!  public subroutine declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
    public pepc_restore_particles         !< once or never per timestep
    public pepc_timber_tree               !< once or never per timestep

    public pepc_grow_and_traverse         !< once per timestep, calls pepc_grow_tree, pepc_traverse_tree, pepc_statistics, pepc_restore_particles, pepc_timber_tree

    public pepc_finalize                  !< mandatory, once per simulation

    public pepc_get_para_file
    public pepc_write_parameters

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!  private variable declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    logical :: pepc_initializes_mpi !< is set to .true., if pepc has to care for MPI_INIT and MPI_FINALIZE; otherwise, the frontend must care for that
    namelist /libpepc/ debug_level, np_mult, curve_type, weighted, interaction_list_length_factor

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!  subroutine-implementation  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    contains


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> Initializes MPI library and data structures for treecode kernel,
    !> reads several parameters from file, that is given as first parameter
    !> to actual executable, initializes submodules
    !> 
    !> Call this function at program startup before any MPI calls
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine pepc_initialize(frontendname, my_rank,n_cpu,init_mpi, db_level_in, comm)
      use treevars, only : me, num_pe, treevars_idim => idim, MPI_COMM_lpepc
      use module_pepc_types, only : register_lpepc_mpi_types
      use module_walk
      use module_domains
      use module_debug, only : pepc_status, debug_level, dbg, DBG_LOADFILE, DBG_TIMINGFILE
      use module_utils, only : create_directory, MPI_IN_PLACE_test
      implicit none
      include 'mpif.h'
      character(*), intent(in) :: frontendname !< name of the program that uses the treecode (only for output purposes)
      integer, intent(out) :: my_rank !< MPI rank of this instance as returned from MPI
      integer, intent(out) :: n_cpu !< number of MPI ranks as returned from MPI
      logical, intent(in) :: init_mpi !< if set to .true., if pepc has to care for MPI_INIT and MPI_FINALIZE; otherwise, the frontend must care for that
      integer, intent(in), optional :: db_level_in !< sets debug level for treecode kernel (overrides settings, that may be read from libpepc-section in input file)
      integer, intent(inout), optional :: comm !< communicator. if pepc initializes MPI, it returns an MPI_COMM_DUP-copy of its own communicator; otherwise, it uses an MPI_COMM_DUP copy of the given comm
      integer :: ierr, provided

      integer, parameter :: MPI_THREAD_LEVEL = MPI_THREAD_FUNNELED ! "The process may be multi-threaded, but the application
                                                                   !  must ensure that only the main thread makes MPI calls."

      call pepc_status('SETUP')

      pepc_initializes_mpi = init_mpi

      if (pepc_initializes_mpi) then
        ! Initialize the MPI system (thread safe version, will fallback automatically if thread safety cannot be guaranteed)
        call MPI_INIT_THREAD(MPI_THREAD_LEVEL, provided, ierr)
        call MPI_COMM_DUP(MPI_COMM_WORLD, MPI_COMM_lpepc, ierr)
        if (present(comm)) call MPI_COMM_DUP(MPI_COMM_lpepc, comm, ierr)
      else
        if (present(comm)) then
           call MPI_COMM_DUP(comm, MPI_COMM_lpepc, ierr)
        else
          call MPI_COMM_DUP(MPI_COMM_WORLD, MPI_COMM_lpepc, ierr)
        endif
      endif

      call MPI_IN_PLACE_test()


      ! Get the id number of the current task
      call MPI_COMM_RANK(MPI_COMM_lpepc, my_rank, ierr)
      ! Get the number of MPI tasks
      call MPI_COMM_size(MPI_COMM_lpepc, n_cpu, ierr)

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

        if ((pepc_initializes_mpi) .and. (provided < MPI_THREAD_LEVEL)) then
          !inform the user about possible issues concerning MPI thread safety
          write(*,'("Call to MPI_INIT_THREAD failed. Requested/provided level of multithreading:", I2, "/" ,I2)') &
                         MPI_THREAD_LEVEL, provided
          write(*,'(a/)') "Initializing with provided level of multithreading. Usually, this is no problem."
        end if
      endif

      ! create all necessary directories
      if (my_rank == 0) then
        call create_directory("diag")
        if( dbg(DBG_LOADFILE) )    call create_directory("load")
        if ( dbg(DBG_TIMINGFILE) ) call create_directory("timing")
      endif

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

      call pepc_prepare(treevars_idim)

    end subroutine


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> Initializes internal variables by reading them from the given a file
    !> that is given as first argument to the actual executable
    !> the optional parameters file_available and filename return
    !> whether and which file was used
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> Initializes internal variables by reading them from the given
    !> filename using several namelists and initializes derived variables
    !> if you just want to pass a handle of an open file,
    !> consider using pepc_read_parameters()
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine pepc_read_parameters_from_file_name(filename)
      implicit none
      character(*), intent(in) :: filename
      integer, parameter :: para_file_id = 10

      open(para_file_id,file=trim(filename),action='read')
      call pepc_read_parameters(para_file_id)
      close(para_file_id)
    end subroutine


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> Initializes internal variables by reading them from the given
    !> filehandle using several namelists and initializes derived variables
    !> if you just want to pass a filename, consider using pepc_read_parameters_from_file_name()
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine pepc_read_parameters(filehandle)
      use module_debug, only : pepc_status
      use module_interaction_specific, only : calc_force_read_parameters
      use module_walk, only: tree_walk_read_parameters
      implicit none
      integer, intent(in) :: filehandle

      call pepc_status("READ PARAMETERS")
      rewind(filehandle)
      call calc_force_read_parameters(filehandle)
      rewind(filehandle)
      call tree_walk_read_parameters(filehandle)
      rewind(filehandle)
      call pepc_status("READ PARAMETERS, section libpepc")
      read(filehandle,NML=libpepc)

    end subroutine

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> Writes internal variables to the given filehandle
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine pepc_write_parameters(filehandle)
      use module_debug, only : pepc_status
      use module_interaction_specific, only : calc_force_write_parameters
      use module_walk, only: tree_walk_write_parameters
      implicit none
      integer, intent(in) :: filehandle

      call pepc_status("WRITE PARAMETERS")
      call calc_force_write_parameters(filehandle)
      call tree_walk_write_parameters(filehandle)
      write(filehandle,NML=libpepc)

    end subroutine


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> Initializes derived variables in pepc, walk, and calc_force
    !> should be called after changing those module`s variables and
    !> before performing tree buildup
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine pepc_prepare(idim)
      use treevars, only : treevars_idim => idim
      use module_walk
      use module_branching, only : branches_initialize
      use module_interaction_specific
      use module_mirror_boxes
      implicit none
      integer, intent(in) :: idim

      treevars_idim = idim

      ! initialize mirror boxes
      call calc_neighbour_boxes()
      ! prepare interaction-specific routines
      call calc_force_prepare()

      call tree_walk_prepare()
      ! initialize data structures in module_branches
      call branches_initialize()

    end subroutine


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> Finalizes MPI library and reverses all initialization from pepc_initialize
    !> Call this function at program termination after all MPI calls
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine pepc_finalize()
      use module_branching
      use module_debug, only : pepc_status
      use module_pepc_types, only : free_lpepc_mpi_types
      use module_walk, only : tree_walk_finalize 
      use module_interaction_specific, only : calc_force_finalize
      implicit none
      include 'mpif.h'
      integer :: ierr

      call pepc_status('FINALIZE')
      ! finalize internal data structures
      call branches_finalize()
      call calc_force_finalize()
      call tree_walk_finalize()
      ! deregister mpi types
      call free_lpepc_mpi_types()

      if (pepc_initializes_mpi) call MPI_FINALIZE(ierr)
    end subroutine


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> checks if the first application argument was set
    !> broadcasts the filename to all mpi processes
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine pepc_get_para_file(available, file_name, my_rank, comm)
        use treevars, only : MPI_COMM_lpepc
        implicit none
        include 'mpif.h'

        logical,   intent(out)          :: available
        character(len=255), intent(out) :: file_name
        integer,   intent(in)           :: my_rank
        integer,   intent(in), optional :: comm

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


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> Builds the tree from the given particles, redistributes particles
    !> to other MPI ranks if necessary (i.e. reallocates particles and changes np_local)
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine pepc_grow_and_traverse(np_local, npart_total, particles, itime, no_dealloc, no_restore)
      use module_pepc_types
      use module_libpepc_main
      use module_debug
      implicit none
      integer, intent(inout) :: np_local    !< number of particles on this CPU, i.e. number of particles in particles-array
      integer, intent(in) :: npart_total !< total number of simulation particles (sum over np_local over all MPI ranks)
      type(t_particle), allocatable, intent(inout) :: particles(:) !< input particle data, initializes %x, %data, %work appropriately (and optionally set %label) before calling this function
      integer, intent(in) :: itime !> current timestep (used as filename suffix for statistics output)
      logical, optional, intent(in) :: no_dealloc ! if .true., the internal data structures are not deallocated (e.g. for a-posteriori diagnostics)
      logical, optional, intent(in) :: no_restore ! if .true., the particles are not backsorted to their pre-domain-decomposition order

      logical :: restore, dealloc

      restore = .true.
      dealloc = .true.

      if (present(no_dealloc)) dealloc = .not. no_dealloc
      if (present(no_restore)) restore = .not. no_restore

      call pepc_grow_tree(np_local, npart_total, particles)
      call pepc_traverse_tree(np_local, particles)
      if (dbg(DBG_STATS)) call pepc_statistics(itime)
      if (restore)        call pepc_restore_particles(np_local, particles)
      if (dealloc)        call pepc_timber_tree()

    end subroutine



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> Builds the tree from the given particles, redistributes particles
    !> to other MPI ranks if necessary (i.e. reallocates particles and changes np_local)
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine pepc_grow_tree(np_local, npart_total, particles)
      use module_pepc_types
      use module_libpepc_main
      implicit none
      integer, intent(inout) :: np_local    !< number of particles on this CPU, i.e. number of particles in particles-array
      integer, intent(in) :: npart_total !< total number of simulation particles (sum over np_local over all MPI ranks)
      type(t_particle), allocatable, intent(inout) :: particles(:) !< input particle data, initializes %x, %data, %work appropriately (and optionally set %label) before calling this function

      call libpepc_grow_tree(np_local, npart_total, particles)

    end subroutine


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine pepc_traverse_tree(nparticles, particles)
      use module_pepc_types
      use module_libpepc_main
      implicit none
      integer, intent(in) :: nparticles    !< number of particles on this CPU, i.e. number of particles in particles-array
      type(t_particle), allocatable, intent(inout) :: particles(:) !< input particle data, initializes %x, %data, %work appropriately (and optionally set %label) before calling this function

      call libpepc_traverse_tree(nparticles, particles)

    end subroutine



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> Writes detailed statistics on the treecode into stats/stats.ITIME.
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine pepc_statistics(itime)
        use module_debug
        use treevars
        use module_timings
        implicit none
        integer, intent(in) :: itime !< current timestep (used as file suffix)
        character(30) :: cfile

        call timer_start(t_fields_stats)
        call tree_stats(itime)

        if( dbg(DBG_LOADFILE) ) then
            write(cfile,'("load/load_",i6.6,".dat")') me
            open(60, file=trim(cfile),STATUS='UNKNOWN', POSITION = 'APPEND')
            write(60,'(i5,2f20.10, i12)') itime,interactions_local, mac_evaluations_local,npp
            close(60)
        end if

        call timer_stop(t_fields_stats)

    end subroutine


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> Restores the initial particle distribution (before calling pepc_grow_tree() ).
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine pepc_restore_particles(np_local, particles)
      use module_pepc_types
      use module_libpepc_main
      implicit none
      integer, intent(inout) :: np_local    !< number of particles on this CPU, i.e. number of particles in particles-array
      type(t_particle), allocatable, intent(inout) :: particles(:) !< input particle data on local MPI rank - is replaced by original particle data that was given before calling pepc_grow_tree()

      call libpepc_restore_particles(np_local, particles)

    end subroutine



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> Frees all tree_specific data fields that were allocated in pepc_groe_tree().
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine pepc_timber_tree()
      use module_timings
      use module_debug, only : pepc_status
      use module_allocation, only : deallocate_tree
      implicit none

      call pepc_status('TIMBER TREE')

     ! deallocate particle and result arrays
      call timer_start(t_deallocate)
      call deallocate_tree()
      call timer_stop(t_deallocate)
      call timer_stop(t_all)

      call pepc_status('TREE HAS FALLEN')

    end subroutine

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> clears result in t_particle datatype
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine pepc_particleresults_clear(particles, nparticles)
      use module_pepc_types
      use module_interaction_specific
      implicit none
      type(t_particle), intent(inout) :: particles(nparticles)
      integer, intent(in) :: nparticles

      call particleresults_clear(particles, nparticles)

    end subroutine


end module module_pepc
