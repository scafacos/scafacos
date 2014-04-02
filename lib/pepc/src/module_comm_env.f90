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
!> Defines a derived type that describes a communication environment and
!> associated procedures.
!>
module module_comm_env
  use module_pepc_types
  implicit none
  private

  !>
  !> Derived type representing a communication environment.
  !>
  type, public :: t_comm_env
    integer(kind_default) :: comm !< an MPI communicator
    integer(kind_pe) :: size !< size of the communicator
    integer(kind_pe) :: rank !< rank within the communicator
    logical :: first !< whether this process is the first process
    logical :: last !< whether this process is the last process
  end type t_comm_env

  interface comm_env_mirror
    module procedure comm_env_mirror
    module procedure comm_env_mirror_from_mpi
  end interface comm_env_mirror

  interface comm_env_dup
    module procedure comm_env_dup
    module procedure comm_env_dup_from_mpi
  end interface comm_env_dup

  public comm_env_mirror
  public comm_env_dup
  public comm_env_destroy

  contains

  !>
  !> Initialize a communication environment `c` from an MPI communicator `comm`.
  !>
  !> The MPI communicator `comm` is duplicated via `mpi_comm_dup` in the
  !> process.
  !>
  subroutine comm_env_mirror_from_mpi(comm, c)
    implicit none
    include 'mpif.h'

    integer, intent(in) :: comm !< the MPI communicator used for communication
    type(t_comm_env), intent(out) :: c !< the communication environment to initialize
    integer(kind_default) :: ierr

    c%comm = comm
    call MPI_COMM_SIZE(c%comm, c%size, ierr)
    call MPI_COMM_RANK(c%comm, c%rank, ierr)
    c%first = c%rank == 0
    c%last = c%rank == (c%size - 1)
  end subroutine comm_env_mirror_from_mpi


  !>
  !> Mirrors a communication environment `c1` to another, `c2`.
  !>
  !> Both environments share the same communicator.
  !>
  subroutine comm_env_mirror(c1, c2)
    implicit none

    type(t_comm_env), intent(in) :: c1 !< environment to mirror
    type(t_comm_env), intent(out) :: c2 !< result of mirroring

    c2 = c1
  end subroutine comm_env_mirror


  !>
  !> Duplicate MPI communicator `comm` and encapsulate.
  !>
  subroutine comm_env_dup_from_mpi(comm, c)
    implicit none

    integer, intent(in) :: comm
    type(t_comm_env), intent(out) :: c

    integer(kind_default) :: cdup, ierr

    call mpi_comm_dup(comm, cdup, ierr)
    call comm_env_mirror_from_mpi(cdup, c)
  end subroutine comm_env_dup_from_mpi


  !>
  !> Duplicate an existing communication environment.
  !>
  !> The underlying MPI communicator is duplicated via `mpi_comm_dup`.
  !>
  subroutine comm_env_dup(c1, c2)
    implicit none

    type(t_comm_env), intent(in) :: c1 !< environment to duplicate
    type(t_comm_env), intent(out) :: c2 !< duplicate

    call comm_env_dup_from_mpi(c1%comm, c2)
  end subroutine comm_env_dup


  !>
  !> Destroy a communication environment.
  !>
  !> Calls `mpi_comm_free` on the MPI communicator contained in `c`.
  !>
  subroutine comm_env_destroy(c)
    implicit none
    include 'mpif.h'

    type(t_comm_env), intent(inout) :: c !< environment to destroy
    integer(kind_default) :: ierr
    
    call mpi_comm_free(c%comm, ierr)
    c%size = 0
    c%rank = 0
    c%first = .false.
    c%last = .false.
  end subroutine comm_env_destroy

end module module_comm_env
