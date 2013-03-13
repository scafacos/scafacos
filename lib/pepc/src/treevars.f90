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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!>  Encapsulates all global variables for lpepc
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module treevars
  
  use module_interaction_specific_types
  use module_pepc_types

  implicit none

  !  Associated parallelization stuff
  integer :: me       !< Rank of current task
  integer :: num_pe   !< # cpus used by program
  integer :: MPI_COMM_lpepc !< communicator that has been supplied to or created by pepc_initialize
  integer :: num_threads = 3 !< number of threads to be used for hybrid parallelization (Pthreads, OpenMP, etc.), for compatibility, we set it to num_walk_threads in tree_walk_read_parameters() for now

  !  tree variables
  integer*8, allocatable :: &
                                branch_key(:), &    !< keys of branch nodes covering all domains
                                pebranch(:)         !< keys of branch nodes covering local domain

  integer, allocatable :: &
                                nbranches(:), &       !< # branches in local domain
                                branch_owner(:)       !< owners of branch nodes covering all domains

  type(t_tree_node_interaction_data), target, allocatable  :: tree_nodes(:)                 !< Tree node properties TODO: move to module_tree

  integer   :: nlev !< max refinement level
  integer*8 :: iplace !< value of place holder bit = 2^(idim*nlev)

  integer :: &
             nleaf, &          ! total # leaf nodes in local #table 
             ntwig, &          ! total # twig nodes in local #table
             nleaf_me, &       ! total # leaves in local domain
             ntwig_me, &       ! total # twigs in local domain
             nlist, &          ! # particles/PE + boundary bodies (1 or 2)
             nbranch, &        ! min # branch nodes covering local domain
             nbranch_sum, &    ! total # branch nodes covering all domains
             nintmax, &        ! max # terms allowed in interaction list
             maxleaf, &        ! max leaf allowed in #table
             maxtwig, &        ! max twig allowed in #table
             maxships, &       ! max # multipole ships per traversal 
             sum_ships, &      ! total # multipole ships per iteration  
             sum_fetches, &    ! total # key fetches  per iteration  
             npart, &          ! actual # particles (total)
             npp, &            ! actual  # particles/PE
             idim              ! dimension of the system

  real*8 :: boxmin(3)  ! box min limits
  real*8 :: boxmax(3)  ! box max limits
  real*8 :: boxsize(3) ! box extension
  real*8 :: interactions_local = 0. !< number of interactions that have been processed locally
  real*8 :: mac_evaluations_local = 0.!< number of mac evaluations that have been processed locally

! Memory control
  real    :: np_mult = 1.5
  integer :: interaction_list_length_factor = 1 !< factor for increasing todo_list_length and defer_list_length in case of respective warning (e.g. for very inhomogeneous or 2D cases set to 2..8)

  contains

  subroutine treevars_prepare(dim)
    implicit none
    integer, intent(in) :: dim

    idim = dim

    nlev = 60 / idim

    iplace = 2_8**(idim * nlev)

  end subroutine treevars_prepare

  subroutine treevars_finalize()
    implicit none

  end subroutine treevars_finalize

end module treevars



