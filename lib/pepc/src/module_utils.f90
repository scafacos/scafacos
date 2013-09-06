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
!> Utility functions.
!>
module module_utils
  implicit none
  private

  interface
    subroutine create_directory_c(dirname) bind(c)
      use, intrinsic :: iso_c_binding, only : c_char
      implicit none
      character(kind=c_char), dimension(*) :: dirname
    end subroutine
  end interface

  public create_directory
  public file_exists
  public MPI_IN_PLACE_test

  contains

  !> creates the directory with te relative pathname dirname
  subroutine create_directory(dirname)
    use, intrinsic :: iso_c_binding, only : c_char, c_null_char
    implicit none
    character(*), intent(in) :: dirname

    call create_directory_c(dirname//c_null_char)
  end subroutine
  
  !> check if file exists
  logical function file_exists(filename)
    implicit none
    character(*), intent(in) :: filename

    inquire(file=trim(filename), exist=file_exists)
  end function
    


  !> checks if MPI_IN_PLACE might be damaged and aborts the application if necessary
  subroutine MPI_IN_PLACE_test(comm)
    implicit none
    include 'mpif.h'

    integer, intent(inout) :: comm

    integer, parameter :: initval = 47
    integer :: ierr, n_cpu, data = initval

    ! Get the number of MPI tasks
    call MPI_COMM_SIZE(comm, n_cpu, ierr)

    call MPI_ALLREDUCE(MPI_IN_PLACE, data, 1, MPI_INTEGER, MPI_SUM, comm, ierr)

    if (initval*n_cpu .ne. data) then
      print '(3(a,/))', 'Serious Issue: MPI_IN_PLACE is not working in your configuration of MPI distribution, compiler(flags) and compiler-optimization.', &
                      'If you are using GCC, you might want to deactivate link time optimization (flags -flto, -fwhole-program, etc.).', &
                      'If running on OSX, please update to at least OpenMPI 1.5.5 or MPICH2 (see also https://svn.open-mpi.org/trac/ompi/ticket/1982).'
      call MPI_ABORT(comm, 1, ierr)
    end if
  end subroutine
end module module_utils
