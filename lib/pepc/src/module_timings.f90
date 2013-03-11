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
!> All stuff concerning timing: timings are contained in a single
!> array. certain entries therein are addressed via integer
!> parameters, eg.
!>
!>      tim(t_allocate) = 0.1234
!>
!> you can use own (frontend-defined) timer constants in the range
!> `t_userdefined_first` .. `t_userdefined_last`, e.g.
!>
!>      call timer_start(t_userdefined + 0)
!>      call timer_start(t_userdefined + 7)
!>
!> etc.
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! The following is a diagram of the hierarchy of a subset of the timers 
! provided by this module.
!
!
!  t_all
!  +
!  |
!  +-> t_fields_tree
!  |   +
!  |   |
!  |   +-> t_domains
!  |   |   +
!  |   |   |
!  |   |   +-> t_domains_keys
!  |   |   |
!  |   |   +-> t_domains_add_sort
!  |   |   |   +
!  |   |   |   |
!  |   |   |   +-> t_domains_sort
!  |   |   |   |   +
!  |   |   |   |   |
!  |   |   |   |   +-> t_domains_sort_pure
!  |   |   |   |
!  |   |   |   +-> t_domains_ship
!  |   |   |       +
!  |   |   |       |
!  |   |   |       +-> t_domains_add_pack
!  |   |   |       |
!  |   |   |       +-> t_domains_add_alltoallv
!  |   |   |       |
!  |   |   |       +-> t_domains_add_unpack
!  |   |   |
!  |   |   +-> t_domains_bound
!  |   |
!  |   +-> t_allocate
!  |   |
!  |   +-> t_local
!  |   |   +
!  |   |   |
!  |   |   +-> t_build_pure
!  |   |   |
!  |   |   +-> t_props_leafs
!  |   |
!  |   +-> t_branches_find
!  |   |
!  |   +-> t_exchange_branches
!  |   |
!  |   +-> t_global
!  |
!  +-> t_fields_passes
!      +
!      |
!      +-> t_walk
!      |   +
!      |   |
!      |   +-> t_comm_total
!      |       +
!      |       |
!      |       +-> t_comm_sendreqs
!      |       |
!      |       +-> t_comm_recv
!      |
!      +-> t_lattice
!

module module_timings
  implicit none

    ! global timings
    integer, parameter :: t_domains            =  1
    integer, parameter :: t_allocate           =  2
    integer, parameter :: t_unused_3            =  3
    integer, parameter :: t_exchange_branches_pack        =  4
    integer, parameter :: t_exchange_branches_allgatherv  =  5
    integer, parameter :: t_exchange_branches_integrate   =  6
    integer, parameter :: t_restore            =  7
    integer, parameter :: t_walk               =  8
    integer, parameter :: t_walk_local         =  9
    integer, parameter :: t_unused10           = 10
    integer, parameter :: t_deallocate         = 11
    integer, parameter :: t_all                = 12
    integer, parameter :: t_local              = 13
    integer, parameter :: t_exchange_branches  = 14
    integer, parameter :: t_global             = 15
    integer, parameter :: t_lattice            = 16
    ! fields internal
    integer, parameter :: t_unused_17          = 17
    integer, parameter :: t_fields_tree        = 18
    integer, parameter :: t_unused19           = 19
    integer, parameter :: t_fields_passes      = 20
    integer, parameter :: t_fields_stats       = 21
    integer, parameter :: t_unused22           = 22
    ! tree_domains
    integer, parameter :: t_domains_keys       = 23
    integer, parameter :: t_domains_sort       = 24
    integer, parameter :: t_domains_sort_pure  = 25
    integer, parameter :: t_domains_ship       = 26
    integer, parameter :: t_domains_bound      = 27
    ! tree_allocate
    integer, parameter :: t_unused28           = 28
    ! tree_build
    integer, parameter :: t_unused_29          = 29
    integer, parameter :: t_unused_30          = 30
    integer, parameter :: t_unused_31          = 31
    ! tree_branches
    integer, parameter :: t_branches_find      = 32
    integer, parameter :: t_exchange_branches_admininstrative = 33
    integer, parameter :: t_build_pure         = 34
    ! pepc_grid_fields
    integer, parameter :: t_walk_grid          = 35
    integer, parameter :: t_lattice_grid       = 36
    ! tree_props
    integer, parameter :: t_props_leafs        = 37
    integer, parameter :: t_unused_38          = 38
    integer, parameter :: t_unused_39          = 39
    integer, parameter :: t_unused_40          = 40
    ! timings for outside fields()
    integer, parameter :: t_tot                = 41
    ! timings for tree_walk_communicator
    integer, parameter :: t_comm_total         = 42
    integer, parameter :: t_comm_recv          = 43
    integer, parameter :: t_comm_sendreqs      = 44
    ! tree_domains, additional timings
    integer, parameter :: t_domains_add_sort   = 45
    integer, parameter :: t_domains_add_pack   = 46
    integer, parameter :: t_domains_add_unpack = 47
    integer, parameter :: t_domains_add_alltoallv = 48

    integer, parameter :: t_direct_force       = 49
    integer, parameter :: t_direct_comm        = 50
    integer, parameter :: t_unused_51          = 51
    integer, parameter :: t_unused_52          = 52
    integer, parameter :: t_unused_53          = 53

    integer, parameter :: t_userdefined_first  = 60
    integer, parameter :: t_userdefined_last   = 90

    !> number of timing entries - dont forget to adjust if you add some timing variables
    integer, private, parameter :: numtimings = t_userdefined_last

    !> array for local timings
    real*8, private, dimension(1:numtimings) :: tim = 0.

  contains

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> Give me a specified timer value
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function timer_read(id)
      implicit none
      integer, intent(in) :: id !< the affected timer address
      real*8 :: timer_read

      timer_read = tim(id)

    end function timer_read

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> Resets a certain timer to zero
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine timer_reset(id)
      implicit none
      integer, intent(in) :: id !< the affected timer address

      tim(id) = 0.

    end subroutine timer_reset

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> Resets all timers to zero
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine timer_reset_all()
      implicit none

      tim(1:numtimings) = 0.

    end subroutine timer_reset_all

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> Starts a timer, i.e. sets
    !>
    !>      tim(id) = - MPI_WTIME()
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine timer_start(id)
      implicit none
      include 'mpif.h'
      integer, intent(in) :: id !< the affected timer address

      tim(id) = - MPI_WTIME()

    end subroutine timer_start


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> Stops a timer, i.e. sets
    !>
    !>      tim(id) = tim(id) + MPI_WTIME()
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine timer_stop(id)
      implicit none
      include 'mpif.h'
      integer, intent(in) :: id !< the affected timer address

      tim(id) = tim(id) + MPI_WTIME()

    end subroutine timer_stop


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> Resumes a timer, i.e. sets
    !>
    !>      tim(id) = tim(id) - MPI_WTIME()
    !>
    !> This can be used to accumulate the durations of a task that is
    !> performed repeatedly between accesses to the timer, e.g.
    !> multiple particles are treated per timestep. Before use, the
    !> timer has to be reset.
    !>
    !>      call timer_reset(t_example)
    !>      do i=1,n
    !>        ...
    !>        call timer_resume(t_example)
    !>        call subroutine_to_be_timed()
    !>        call timer_stop(t_example)
    !>        ...
    !>      end do
    !>      accumulated_time = timer_read(t_example)
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine timer_resume(id)
      implicit none
      include 'mpif.h'
      integer, intent(in) :: id !< the affected timer address

      tim(id) = tim(id) - MPI_WTIME()

    end subroutine timer_resume


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> Adds the given value to a timer for cumulative measurements
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine timer_add(id, val)
      implicit none
      integer, intent(in) :: id !< the affected timer address
      real*8, intent(in) :: val !< value to be added

      tim(id) = tim(id) + val

    end subroutine timer_add


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> Outputs the given timing array to a file with the given filename
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine timings_ToFile(itime, iuserflag, tdata, filename, printheader)
      implicit none
      integer, intent(in) :: itime !< current timestep
      integer, intent(in) :: iuserflag !< frontend-defined flag that is passed through and output to the second column
      character(*), intent(in) :: filename !< output filename
      real*8, dimension(1:numtimings), intent(in) :: tdata !< array with timing data
      logical, optional, intent(in) :: printheader !< if set to true, a header with column numbers is output

      integer, parameter :: ifile = 60
      integer :: i
      character(30) :: formatstring

      open  (ifile, file=filename,STATUS='UNKNOWN', POSITION = 'APPEND')


      if (printheader) then
        write(formatstring,'(a,i5,a)' ) '(a1,2(1x,a20),', numtimings, '(1x,i20))'
        write (ifile,formatstring) "#", "timestep", "userflag", [ (i,i=1,numtimings) ]
      endif

      write(formatstring,'(a,i5,a)' ) '(x,2(1x,i20),', numtimings, '(1x,e20.5))'
      write (ifile,formatstring) itime, iuserflag, tdata
      close (ifile)

    end subroutine timings_ToFile




    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> Outputs all local timing data to timing_XXXX.dat
    !> if `itime <=1`, an additional header is printed
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine timings_LocalOutput(itime, iuserflag)
      use module_debug
      use treevars

      implicit none
      integer, intent(in) :: itime !< current timestep
      integer, optional, intent(in) :: iuserflag !< frontend-defined flag that is passed through and output to the second column
      character(30) :: cfile
      integer :: flag

      if (present(iuserflag)) then
        flag = iuserflag
      else
        flag = 0
      endif

      if ( dbg(DBG_TIMINGFILE) ) then
         write(cfile,'("timing/timing_",i6.6,".dat")') me
         call timings_ToFile(itime, flag, tim, cfile, itime<=1)
      end if

    end subroutine timings_LocalOutput


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> Gathers global timing data and outputs to
    !> timing_avg.dat, timing_min.dat, timing_max.dat.
    !> If `printheader` is present, its value controls the output of
    !> descriptive column headers, otherwise, if `itime <=1`, headers
    !> are printed.
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine timings_GatherAndOutput(itime, iuserflag, printheader)
      use treevars
      implicit none
      include 'mpif.h'
      integer, intent(in) :: itime !< current timestep
      integer, optional, intent(in) :: iuserflag !< frontend-defined flag that is passed through and output to the second column
      logical, optional, intent(in) :: printheader
      integer :: ierr

      integer :: flag
      logical :: do_printheader

      real*8, dimension(1:numtimings) :: tim_max
      real*8, dimension(1:numtimings) :: tim_avg
      real*8, dimension(1:numtimings) :: tim_min
      real*8, dimension(1:numtimings) :: tim_dev

      if (present(iuserflag)) then
        flag = iuserflag
      else
        flag = 0
      endif

      if (present(printheader)) then
        do_printheader = printheader
      else
        do_printheader = itime <= 1
      end if

      call MPI_REDUCE(tim, tim_max, numtimings, MPI_REAL8, MPI_MAX, 0, MPI_COMM_lpepc,ierr);
      call MPI_REDUCE(tim, tim_min, numtimings, MPI_REAL8, MPI_MIN, 0, MPI_COMM_lpepc,ierr);
      call MPI_REDUCE(tim, tim_avg, numtimings, MPI_REAL8, MPI_SUM, 0, MPI_COMM_lpepc,ierr);

     if (me==0) then
        tim_avg = tim_avg / num_pe
        tim_dev = tim_max - tim_min
        call timings_ToFile(itime, flag, tim_max, 'timing_max.dat', do_printheader)
        call timings_ToFile(itime, flag, tim_avg, 'timing_avg.dat', do_printheader)
        call timings_ToFile(itime, flag, tim_min, 'timing_min.dat', do_printheader)
        call timings_ToFile(itime, flag, tim_dev, 'timing_dev_abs.dat', do_printheader)
        tim_dev = tim_dev / tim_min
        call timings_ToFile(itime, flag, tim_dev, 'timing_dev_rel.dat', do_printheader)

        write(*,'(a20,f16.10," s")') "t_all = ",       tim(t_all)
        write(*,'(a20,f16.10," s")') "t_tot = ",       tim(t_tot)
     endif

    end subroutine timings_GatherAndOutput

end module module_timings
