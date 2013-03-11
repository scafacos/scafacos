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

!
!                          Tree sort utility module
!                       *   Parallel Sorting
!

module module_utils

  interface
    subroutine create_directory_c(dirname) bind(c)
      use, intrinsic :: iso_c_binding, only : c_char
      implicit none
      character(kind=c_char), dimension(*) :: dirname
    end subroutine
  end interface


  interface sort
     module procedure sort_i
  end interface

  interface swap
     module procedure swap_ab
  end interface

contains

  !> creates the directory with te relative pathname dirname
  subroutine create_directory(dirname)
    use, intrinsic :: iso_c_binding, only : c_char, c_null_char
    implicit none
    character(*), intent(in) :: dirname

    call create_directory_c(dirname//c_null_char)
  end subroutine


  !> checks if MPI_IN_PLACE might be damaged and aborts the application if necessary
  subroutine MPI_IN_PLACE_test()
    use module_debug
    use treevars, only : MPI_COMM_lpepc
    implicit none
    include 'mpif.h'

    integer :: ierr, n_cpu
    integer, parameter :: initval = 47
    integer :: data = initval


    ! Get the number of MPI tasks
    call MPI_COMM_size(MPI_COMM_lpepc, n_cpu, ierr)

    call MPI_ALLREDUCE(MPI_IN_PLACE, data, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_lpepc, ierr)

    if (initval*n_cpu .ne. data) then
      DEBUG_ERROR_NO_DIAGFILE('(a,/)','Serious Issue: MPI_IN_PLACE is not working in your configuration of MPI distribution, compiler(flags) and compiler-optimization.',
                          'If you are using GCC, you might want to deactivate link time optimization (flags -flto, -fwhole-program, etc.).',
                          'If running on OSX, please update to at least OpenMPI 1.5.5 or MPICH2 (see also https://svn.open-mpi.org/trac/ompi/ticket/1982).')
    endif

  end subroutine


  !  ================================
  !
  !         SORT_HEAP
  !
  !     Sort 64-bit integer array into ascending order
  !     using heap-sort algorithm
  !    (Numerical Recipes f90, p1171)
  !    Modified 2011.12.02 by Andreas Breslau
  !    to optionally return an integer array containing the original positions of the sorted array (iarr).
  !    This can be used by the calling code to sort several other arrays according to iarr.
  !
  !  ================================

  subroutine sort_i(iarr, map)
    use module_debug
    implicit none
    integer*8, intent(inout) :: iarr(:)
    integer, optional, intent(out) :: map(:)           !< the optional integer array filled with positions for sorting other arrays according to iarr
    integer :: i,n

    if (present(map)) then
      if (size(map) .ne. size(iarr)) then
         DEBUG_ERROR(*, "Error in tree_utils sort_i: optional second parameter has not the same size as first argument. Ignoring second argument.")
      endif
    end if

    n = size(iarr)

    ! only sort map, if map present
    if(present(map)) map = [ (i, i=1,n) ]

    do i=n/2,1,-1                      ! Left range - hiring phase (heap creation)
       call sift_down(i,n)
    end do

    do i=n,2,-1                        ! Right range of sift-down is decr. from n-1 ->1
       ! during retirement/promotion (heap selection) phase.
       call sort_swap_ab( 1,i )      ! Clear space at end of array and retire top of heap into it
       call sift_down( 1,i-1)
    end do

  contains
    subroutine sift_down(l,r)        ! Carry out the sift-down on element arr(l) to maintain the heap structure
      ! Modified 2011.12.02 by Andreas Breslau to sort the map array according to iarr (see above description)
      integer, intent(in) :: l,r
      integer :: j,jold    ! index
      integer*8 :: a
      integer :: b
      
      a = iarr(l)

      ! only sort map, if map present
      if(present(map)) b= map(l)

      jold = l
      j = l + l
      do                   ! do while j <= r
         if (j > r) exit
         if (j < r) then
            if (iarr(j) < iarr(j+1)) j = j+1
         endif
         if (a >= iarr(j)) exit       ! Found a`s level, so terminate sift-down
         iarr(jold) = iarr(j)

         ! only sort map, if map present
         if(present(map)) map(jold) = map(j)

         jold = j                    ! Demote a and continue
         j = j+j
      end do
      iarr(jold) = a                  ! Put a into its slot

      ! only sort map, if map present
      if(present(map)) map(jold) = b

    end subroutine sift_down

    subroutine sort_swap_ab(p,q)
      integer :: p,q, dum
      integer*8 dum8

      dum8 = iarr(p)
      iarr(p)=iarr(q)
      iarr(q) = dum8

      ! only sort map, if map present
      if(present(map)) then
         dum = map(p)
         map(p)=map(q)
         map(q) = dum
      end if

    end subroutine sort_swap_ab

  end subroutine sort_i


  subroutine swap_ab(p,q)
    integer*8 :: p,q, dum
    dum = p
    p=q
    q = dum
  end subroutine swap_ab

end module module_utils
