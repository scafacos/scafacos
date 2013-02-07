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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!>  Encapsulates some debugging and i/o specific routines
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module module_debug
     implicit none
     save
     private

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  public variable declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      integer, public :: debug_ipefile = 21
      integer, public :: debug_my_rank = -1
      integer, public :: debug_stdout  =  6

      integer, public :: debug_level = 0 !< or-combination of the bitmasks below
      ! set to the following values to get the old behaviour:
      !             db_level = 0      --> debug_level = 0
      !             db_level = 1      --> debug_level = 1 + 2 + 32 = 35
      !             db_level = 2      --> debug_level = 1 + 2 + 32 + 64 = 99
      !             db_level = 3      --> debug_level = 1 + 2 + 32 + 64 + 8 + 2048 = 2155
      !             db_level = 4      --> debug_level = 1 + 2 + 32 + 64 + 8 + 2048 + 4 + 16 = 2175
      !             db_level = 5      --> debug_level = 128
      !             db_level = 6      --> debug_level = 1 + 2 + 32 + 64 + 8 + 2048 + 4 + 16 + 128 + 256 + 512 = 3071
      !
      integer, parameter, public :: DBG_STATUS      = B'0000000000000001'    ! 1
      integer, parameter, public :: DBG_TREE        = B'0000000000000010'    ! 2
      integer, parameter, public :: DBG_BUILD       = B'0000000000000100'    ! 4
      integer, parameter, public :: DBG_DOMAIN      = B'0000000000001000'    ! 8
      integer, parameter, public :: DBG_BRANCH      = B'0000000000010000'    ! 16
      integer, parameter, public :: DBG_STATS       = B'0000000000100000'    ! 32
      integer, parameter, public :: DBG_WALKSUMMARY = B'0000000001000000'    ! 64
      integer, parameter, public :: DBG_DUMPTREE    = B'0000000010000000'    ! 128
      integer, parameter, public :: DBG_TIMINGFILE  = B'0000000100000000'    ! 256
      integer, parameter, public :: DBG_LOADFILE    = B'0000001000000000'    ! 512
      integer, parameter, public :: DBG_WALK        = B'0000010000000000'    ! 1024
      integer, parameter, public :: DBG_PERIODIC    = B'0000100000000000'    ! 2048

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  private variable declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      logical, private       :: debug_initialized = .false.
      character(30), private :: debug_ipefile_name

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  public subroutine declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      public dbg
      public pepc_status

      public debug_ipefile_open
      public debug_ipefile_close
      public debug_mpi_abort
      public debug_barrier


   contains

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !>
     !>  lpepc status output
     !>
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     subroutine pepc_status(stat)
       use treevars, only : me
       implicit none
       character(*), intent(in) :: stat

       if (dbg(DBG_STATUS)) then
          if (debug_initialized) then ! output to ipefile only if it already has been created
            call debug_ipefile_open()
            write(debug_ipefile,'("LPEPC | ", a)') stat
            call debug_ipefile_close()
          endif

          if (me==0) write(*,'("LPEPC | ", a)') stat
       endif

     end subroutine

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !>
     !>  module initialization (is called automatically on first call to debug_ipefile_open)
     !>
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     subroutine debug_initialize()
       use treevars
       implicit none
       include 'mpif.h'

       character(MPI_MAX_PROCESSOR_NAME) :: procname
       integer :: resultlen, ierr

       debug_my_rank = me
       call MPI_GET_PROCESSOR_NAME( procname, resultlen, ierr )

       write(debug_ipefile_name,'("diag/diag_",i6.6,".dat")') me

       open(debug_ipefile, file=trim(debug_ipefile_name),STATUS='UNKNOWN', POSITION = 'REWIND')
       call timstamp(debug_ipefile, "PEPC on ["//procname(1:resultlen)//"]")
       close(debug_ipefile)

       debug_initialized = .true.

     end subroutine debug_initialize

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !>
     !>  prints timestamp and reason to istream
     !>
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     subroutine timstamp(istream,reason)
       implicit none

       character :: cdate*8, ctime*10, czone*5
       integer, intent(in) :: istream
       character(*), intent(in) :: reason

       call DATE_AND_TIME(cdate,ctime,czone)

       write(istream, '(2(a2,"/"),a4,", ",2(a2,":"),a2, " (GMT", a5") - ", a, //)') &
             cdate(7:8), cdate(5:6), cdate(1:4), &
             ctime(1:2), ctime(3:4), ctime(5:6), &
             czone, &
             reason

     end subroutine


     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !>
     !>  opens the processor diagnostic file, afterwards, the application may
     !>  write to file stream debug_ipefile
     !>
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     subroutine debug_ipefile_open()
       implicit none

       if (.not. debug_initialized) call debug_initialize()

       open(debug_ipefile, file=trim(debug_ipefile_name),STATUS='UNKNOWN', POSITION = 'APPEND')
     end subroutine


     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !>
     !>  calls MPI_ABORT(MPI_COMM_lpepc, 1, ierr)
     !>
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     subroutine debug_mpi_abort()
       use treevars, only : MPI_COMM_lpepc
       implicit none
       include 'mpif.h'
       integer :: ierr

       call MPI_ABORT(MPI_COMM_lpepc, 1, ierr)

     end subroutine

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !>
     !>  calls MPI_BARRIER(MPI_COMM_lpepc, ierr)
     !>
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     subroutine debug_barrier()
       use treevars, only : MPI_COMM_lpepc
       implicit none
       include 'mpif.h'
       integer :: ierr

       call MPI_BARRIER(MPI_COMM_lpepc, ierr)

     end subroutine

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !>
     !>  closes the processor diagnostic file, afterwards, the application may not
     !>  write to file stream debug_ipefile
     !>
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     subroutine debug_ipefile_close()
       implicit none
       close(debug_ipefile)
     end subroutine


     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !>
     !>  Debug flag query function
     !>
     !>  Usage:
     !>      if (dbg(DBG_DOMAIN)) call do_some_debug()
     !>
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     function dbg(flag)
       implicit none
       logical :: dbg
       integer, intent(in) :: flag

       dbg = (iand(debug_level, flag) .ne. 0)
     end function

end module module_debug
