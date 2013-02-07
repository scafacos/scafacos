module allocation
use iso_c_binding
use fmmkinds
implicit none
 public
 interface
  type(c_ptr) function dummy_malloc_aligned(bsize,ierr) bind(c,Name='dummy_malloc_aligned')
  use iso_c_binding, only : c_ptr,c_int
  implicit none
  type(c_ptr) :: ptr
  integer(c_int), value :: bsize
  integer(c_int) :: ierr
  end function dummy_malloc_aligned

  type(c_ptr) function dummy_malloc(bsize,ierr) bind(c,Name='dummy_malloc')
  use iso_c_binding, only : c_ptr,c_int
  implicit none
  type(c_ptr) :: ptr
  integer(c_int), value :: bsize
  integer(c_int) :: ierr
  end function dummy_malloc

  subroutine dummy_free(ptr) bind(c,Name='dummy_free')
  use iso_c_binding, only : c_ptr
  implicit none
  type(c_ptr), value :: ptr
  end subroutine dummy_free
 end interface

 interface allocate_aligned
  module procedure allocate_aligned_1d
  module procedure allocate_aligned_2d
 end interface allocate_aligned

 contains

 subroutine allocate_aligned_1d(cptr,lo,hi,ierr)
 implicit none
 type(c_ptr) :: cptr
 integer(kind=fmm_integer) :: lo,hi
 integer(kind=fmm_integer) :: ierr
 integer(c_int) :: bsize, err

  bsize = (hi-lo+1)*fmm_real
  cptr = dummy_malloc_aligned(bsize,err)

  ierr = err

 end subroutine allocate_aligned_1d

 subroutine allocate_aligned_2d(cptr,lo1,hi1,lo2,hi2,ierr)
 implicit none
 type(c_ptr) :: cptr
 integer(kind=fmm_integer) :: lo1,hi1,lo2,hi2
 integer(kind=fmm_integer) :: ierr
 integer(c_int) :: bsize, err

  bsize = (hi1-lo1+1)*(hi2-lo2+1)*fmm_real
  cptr = dummy_malloc_aligned(bsize,err)

  ierr = err

 end subroutine allocate_aligned_2d

 subroutine deallocate_aligned(cptr,ierr)
 implicit none
 type (c_ptr):: cptr
 integer(kind=fmm_integer) :: ierr

  call dummy_free(cptr)
  ierr = 0

 end subroutine deallocate_aligned

 end module allocation
