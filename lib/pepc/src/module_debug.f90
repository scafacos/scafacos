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
!>  Encapsulates some debugging and i/o specific routines
!>
module module_debug
     use module_pepc_types
     implicit none
     save
     private

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
      ! deprecated: integer, parameter, public :: DBG_LOADFILE    = B'0000001000000000'    ! 512
      integer, parameter, public :: DBG_WALK        = B'0000010000000000'    ! 1024
      integer, parameter, public :: DBG_PERIODIC    = B'0000100000000000'    ! 2048

      logical, private       :: debug_initialized = .false.
      character(30), private :: debug_ipefile_name

      public dbg
      public pepc_status
      public debug_ipefile_open
      public debug_ipefile_close
      public debug_mpi_abort
      public debug_barrier
      public debug_print_timestamp

   contains

     !>
     !>  lpepc status output
     !>
     subroutine pepc_status(stat)
       use, intrinsic :: iso_fortran_env, only: output_unit
       use treevars, only : me
       implicit none
       character(*), intent(in) :: stat

       if (dbg(DBG_STATUS)) then
          if (debug_initialized) then ! output to ipefile only if it already has been created
            call debug_ipefile_open()
            call debug_print_timestamp(debug_ipefile)
            write(debug_ipefile,'(" LPEPC | ", a)') stat
            call debug_ipefile_close()
          endif

          if (me==0) then
            call debug_print_timestamp(output_unit)
            write(*,'(" LPEPC | ", a)') stat
          end if
       endif
     end subroutine


     !>
     !>  module initialization (is called automatically on first call to debug_ipefile_open)
     !>
     subroutine debug_initialize()
       use treevars
       use module_utils, only: create_directory
       implicit none
       include 'mpif.h'

       character(MPI_MAX_PROCESSOR_NAME) :: procname
       integer(kind_default) :: resultlen, ierr

       debug_my_rank = me
       call MPI_GET_PROCESSOR_NAME( procname, resultlen, ierr )

       write(debug_ipefile_name,'("diag/diag_",i6.6,".dat")') me
       call create_directory("diag")

       open(debug_ipefile, file=trim(debug_ipefile_name),STATUS='UNKNOWN', POSITION = 'REWIND')
       call debug_print_timestamp(debug_ipefile)
       write (debug_ipefile, '(a)') " PEPC on ["//procname(1:resultlen)//"]"
       close(debug_ipefile)

       debug_initialized = .true.
     end subroutine debug_initialize


     !>
     !> Writes the current date and time into the argument in the format YYYY-MM-DD hh:mm:ss
     !>
     subroutine debug_timestamp(str)
       implicit none
       character(*), intent(out) :: str

       ! [
       !  year,
       !  month of year,
       !  day of month,
       !  difference to UTC in minutes,
       !  hour of day,
       !  minutes of hour,
       !  seconds of minute,
       !  milliseconds of second
       ! ]
       integer :: v(8)

       if (len(str) >= 19) then
         call date_and_time(values = v)
         write (str, '(i4, "-", i2.2, "-", i2.2, " ", i2.2, ":", i2.2, ":", i2.2)') v(1), v(2), v(3), v(5), v(6), v(7)
       end if
     end subroutine debug_timestamp


     !>
     !> Writes the output of `debug_timestamp` surrounded by angled brackets to the I/O unit specified in the argument
     !> without advancing to the next line.
     !>
     subroutine debug_print_timestamp(iounit)
       implicit none

       integer(kind_default), intent(in) :: iounit

       character(19) :: str

       call debug_timestamp(str)

       write (iounit, '("<", a, ">")', advance = 'no') str
    end subroutine debug_print_timestamp


     !>
     !>  opens the processor diagnostic file, afterwards, the application may
     !>  write to file stream debug_ipefile
     !>
     subroutine debug_ipefile_open()
       implicit none

       if (.not. debug_initialized) call debug_initialize()

       open(debug_ipefile, file=trim(debug_ipefile_name),STATUS='UNKNOWN', POSITION = 'APPEND')
     end subroutine


     !>
     !>  calls MPI_ABORT(MPI_COMM_lpepc, 1, ierr)
     !>
     subroutine debug_mpi_abort()
       use treevars, only : MPI_COMM_lpepc
       #if defined(__ICC) || defined(__INTEL_COMPILER)
         use ifcore
       #endif
       implicit none
       include 'mpif.h'
       integer(kind_default) :: ierr

       ! starting from GCC version 4.8, a backtrace() subroutine is provided by gfortran
       #if defined(__GNUC__)
         #define GCC_VERSION (__GNUC__ * 10000 \
                            + __GNUC_MINOR__ * 100 \
                            + __GNUC_PATCHLEVEL__)
         #if GCC_VERSION >= 40800
           ! http://gcc.gnu.org/onlinedocs/gfortran/BACKTRACE.html
           call backtrace()
         #endif
       #elif defined(__ICC) || defined(__INTEL_COMPILER)
         ! http://software.intel.com/sites/products/documentation/hpc/composerxe/en-us/2011Update/fortran/lin/lref_for/source_files/rftrace.htm
         call tracebackqq("stacktrace", -1)
       #elif defined(__IBMC__) || defined(__IBMCPP__)
         ! http://publib.boulder.ibm.com/infocenter/comphelp/v8v101/index.jsp?topic=%2Fcom.ibm.xlf101a.doc%2Fxlflr%2Fsup-xltrbk.htm
         call xl__trbk()
       #endif

       call MPI_ABORT(MPI_COMM_lpepc, 1, ierr)
     end subroutine


     !>
     !>  calls MPI_BARRIER(MPI_COMM_lpepc, ierr)
     !>
     subroutine debug_barrier()
       use treevars, only : MPI_COMM_lpepc
       implicit none
       include 'mpif.h'
       integer(kind_default) :: ierr

       call MPI_BARRIER(MPI_COMM_lpepc, ierr)
     end subroutine


     !>
     !>  closes the processor diagnostic file, afterwards, the application may not
     !>  write to file stream debug_ipefile
     !>
     subroutine debug_ipefile_close()
       implicit none
       close(debug_ipefile)
     end subroutine


     !>
     !>  Debug flag query function
     !>
     !>  Usage:
     !>      if (dbg(DBG_DOMAIN)) call do_some_debug()
     !>
     function dbg(flag)
       implicit none
       logical :: dbg
       integer, intent(in) :: flag

       dbg = (iand(debug_level, flag) .ne. 0)
     end function
end module module_debug
