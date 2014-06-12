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
!> Provides `sort()`, a local heap sort for integer arrays.
!>
module module_sort
  implicit none
  private

  interface sort
     module procedure sort_i8
     module procedure sort_i8i4
     module procedure sort_i4
  end interface

  interface swap
    module procedure swap4
    module procedure swap8
  end interface

  public sort

  contains

  !>
  !> Sort 64-bit integer array into ascending order using heap-sort algorithm. (Numerical Recipes f90, p1171)
  !>
  !> @note to optionally return an integer array containing the original positions of the sorted array (iarr).
  !> This can be used by the calling code to sort several other arrays according to iarr.
  !>
  subroutine sort_i8(iarr, map)
    use module_debug
    implicit none
    integer*8, intent(inout) :: iarr(:)
    integer*8, optional, intent(out) :: map(:) !< the optional integer array filled with positions for sorting other arrays according to iarr
    integer*8 :: i,n

    if (present(map)) then
      DEBUG_ASSERT_MSG(size(map)==size(iarr), *, "Error in tree_utils sort_i: optional second parameter has not the same size as first argument. Ignoring second argument.")
    end if

    n = size(iarr, kind=kind(n))

    ! only sort map, if map present
    if(present(map)) map = [ (i, i=1,n) ]

    do i=n/2,1,-1                      ! Left range - hiring phase (heap creation)
       call sift_down(i,n)
    end do

    do i=n,2,-1                        ! Right range of sift-down is decr. from n-1 ->1
       ! during retirement/promotion (heap selection) phase.
       call swap( iarr(1), iarr(i) )      ! Clear space at end of array and retire top of heap into it
       if (present(map)) call swap(map(1), map(i))
       call sift_down(    1_8, i-1_8)
    end do

    contains

    subroutine sift_down(l,r)        ! Carry out the sift-down on element arr(l) to maintain the heap structure
      integer*8, intent(in) :: l,r
      integer*8 :: j,jold    ! index
      integer*8 :: a, b
      
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

  end subroutine sort_i8

  !>
  !> Sort 64/32-bit integer array into ascending order using heap-sort algorithm. (Numerical Recipes f90, p1171)
  !>
  !> @note to optionally return an integer array containing the original positions of the sorted array (iarr).
  !> This can be used by the calling code to sort several other arrays according to iarr.
  !>
  subroutine sort_i8i4(iarr, map)
    use module_debug
    implicit none
    integer*8, intent(inout) :: iarr(:)
    integer*4, intent(out) :: map(:) !< the optional integer array filled with positions for sorting other arrays according to iarr
    integer*4 :: i,n

    DEBUG_ASSERT_MSG(size(map)==size(iarr), *, "Error in tree_utils sort_i: optional second parameter has not the same size as first argument. Ignoring second argument.")

    n = size(iarr, kind=kind(n))

    map = [ (i, i=1,n) ]

    do i=n/2,1,-1                      ! Left range - hiring phase (heap creation)
       call sift_down(i,n)
    end do

    do i=n,2,-1                        ! Right range of sift-down is decr. from n-1 ->1
       ! during retirement/promotion (heap selection) phase.
       call swap( iarr(1), iarr(i) )      ! Clear space at end of array and retire top of heap into it
       call swap(map(1), map(i))
       call sift_down(    1_4, i-1_4)
    end do

    contains

    subroutine sift_down(l,r)        ! Carry out the sift-down on element arr(l) to maintain the heap structure
      integer*4, intent(in) :: l,r
      integer*4 :: j,jold    ! index
      integer*8 :: a
      integer*4 :: b
      
      a = iarr(l)
      b= map(l)

      jold = l
      j = l + l
      do                   ! do while j <= r
         if (j > r) exit
         if (j < r) then
            if (iarr(j) < iarr(j+1)) j = j+1
         endif
         if (a >= iarr(j)) exit       ! Found a`s level, so terminate sift-down
         iarr(jold) = iarr(j)

         map(jold) = map(j)

         jold = j                    ! Demote a and continue
         j = j+j
      end do
      iarr(jold) = a                  ! Put a into its slot

      map(jold) = b
    end subroutine sift_down

  end subroutine sort_i8i4

  !>
  !> Sort 32-bit integer array into ascending order using heap-sort algorithm. (Numerical Recipes f90, p1171)
  !>
  !> @note to optionally return an integer array containing the original positions of the sorted array (iarr).
  !> This can be used by the calling code to sort several other arrays according to iarr.
  !>
  subroutine sort_i4(iarr, map)
    use module_debug
    implicit none
    integer*4, intent(inout) :: iarr(:)
    integer*4, optional, intent(out) :: map(:) !< the optional integer array filled with positions for sorting other arrays according to iarr
    integer*4 :: i,n

    if (present(map)) then
      DEBUG_ASSERT_MSG(size(map)==size(iarr), *, "Error in tree_utils sort_i: optional second parameter has not the same size as first argument. Ignoring second argument.")
    end if

    n = size(iarr, kind=kind(n))

    ! only sort map, if map present
    if(present(map)) map = [ (i, i=1,n) ]

    do i=n/2,1,-1                      ! Left range - hiring phase (heap creation)
       call sift_down(i,n)
    end do

    do i=n,2,-1                        ! Right range of sift-down is decr. from n-1 ->1
       ! during retirement/promotion (heap selection) phase.
       call swap( iarr(1), iarr(i) )      ! Clear space at end of array and retire top of heap into it
       if (present(map)) call swap(map(1), map(i))
       call sift_down(    1_4, i-1_4)
    end do

    contains

    subroutine sift_down(l,r)        ! Carry out the sift-down on element arr(l) to maintain the heap structure
      integer*4, intent(in) :: l,r
      integer*4 :: j,jold    ! index
      integer*4 :: a, b
      
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

  end subroutine sort_i4
  
  subroutine swap4(p,q)
    integer*4 :: p,q,dum

    dum = p
    p   = q
    q   = dum
  end subroutine

  subroutine swap8(p,q)
    integer*8 :: p,q,dum

    dum = p
    p   = q
    q   = dum
  end subroutine

end module module_sort
