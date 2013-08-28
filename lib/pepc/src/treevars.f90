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
!>  Encapsulates all global variables for lpepc
!>
module treevars
  use module_pepc_types
  implicit none

  ! I/O units
  integer, parameter :: stats_u = 60 !< statistics

  !  Associated parallelization stuff
  integer(kind_pe) :: me       !< Rank of current task
  integer(kind_pe) :: num_pe   !< # cpus used by program
  integer(kind_default) :: MPI_COMM_lpepc !< communicator that has been supplied to or created by pepc_initialize
  integer :: num_threads = 3 !< number of threads to be used for hybrid parallelization (Pthreads, OpenMP, etc.), for compatibility, we set it to num_walk_threads in tree_walk_read_parameters() for now
  integer :: main_thread_processor_id !< id of processor that runs the applications main thread

  integer(kind_level) :: nlev !< max refinement level
  integer(kind_dim)   :: idim = 3_kind_dim !< dimension of the system

! Memory control
  real    :: np_mult = 1.5
  integer :: interaction_list_length_factor = 1 !< factor for increasing todo_list_length and defer_list_length in case of respective warning (e.g. for very inhomogeneous or 2D cases set to 2..8)

  contains

  subroutine treevars_prepare(dim)
    implicit none

    integer(kind_dim), optional, intent(in) :: dim

    if (present(dim)) then; idim = dim; end if

                                ! we do not use the first (sign) bit and need space for one additional placeholder-bit
    nlev = int(bit_size(1_kind_key) - 2, kind_level) / idim
  end subroutine treevars_prepare


  subroutine treevars_finalize()
    implicit none
  end subroutine treevars_finalize
end module treevars



