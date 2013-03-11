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

module module_atomic_ops
  use, intrinsic :: iso_c_binding
  implicit none

  private

  public t_atomic_int
  public atomic_allocate_int
  public atomic_deallocate_int
  public atomic_store_int
  public atomic_load_int
  public atomic_fetch_and_increment_int
  public atomic_mod_increment_and_fetch_int
  public atomic_write_barrier
  public atomic_read_barrier
  public atomic_read_write_barrier

  type t_atomic_int
    type(c_ptr) :: p
  end type t_atomic_int

  interface

    type(c_ptr) function c_atomic_alloc_int() bind(C, name='_atomic_alloc_int')
      use, intrinsic :: iso_c_binding
      implicit none
    end function c_atomic_alloc_int

    subroutine c_atomic_free_int(storage) bind(C, name='_atomic_free_int')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), intent(in), value :: storage
    end subroutine c_atomic_free_int

    integer(kind=c_int) function c_atomic_load_int(storage) bind(C, name='_atomic_load_int')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), intent(in), value :: storage
    end function c_atomic_load_int

    subroutine c_atomic_store_int(storage, val) bind(C, name='_atomic_store_int')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), intent(in), value :: storage
      integer(kind=c_int), intent(in), value :: val
    end subroutine c_atomic_store_int

    integer(kind=c_int) function c_atomic_fetch_and_increment_int(storage) bind(C, name='_atomic_fetch_and_increment_int')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), intent(in), value :: storage
    end function c_atomic_fetch_and_increment_int

    integer(kind=c_int) function c_atomic_mod_increment_and_fetch_int(storage, mod) bind(C, name='_atomic_mod_increment_and_fetch_int')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), intent(in), value :: storage
      integer(kind=c_int), intent(in), value :: mod
    end function c_atomic_mod_increment_and_fetch_int    

    subroutine atomic_write_barrier() bind(C, name='_atomic_write_barrier')
      use, intrinsic :: iso_c_binding
      implicit none
    end subroutine atomic_write_barrier

    subroutine atomic_read_barrier() bind(C, name='_atomic_read_barrier')
      use, intrinsic :: iso_c_binding
      implicit none
    end subroutine atomic_read_barrier

    subroutine atomic_read_write_barrier() bind(C, name='_atomic_read_write_barrier')
      use, intrinsic :: iso_c_binding
      implicit none
    end subroutine atomic_read_write_barrier

  end interface

  contains


  subroutine atomic_allocate_int(storage)
    implicit none

    type(t_atomic_int), pointer, intent(out) :: storage

    type(t_atomic_int), pointer :: tmp_f
    type(c_ptr) :: tmp_c

    allocate(tmp_f)
    tmp_c = c_atomic_alloc_int()

    if (.not. (associated(tmp_f) .and. c_associated(tmp_c))) then
      call c_atomic_free_int(tmp_c)
      deallocate(tmp_f)
      storage => null()
    else
      tmp_f%p = tmp_c
      storage => tmp_f
    end if
  end subroutine atomic_allocate_int


  subroutine atomic_deallocate_int(storage)
    implicit none

    type(t_atomic_int), pointer, intent(inout) :: storage

    call c_atomic_free_int(storage%p)
    deallocate(storage)
    storage => null()
  end subroutine atomic_deallocate_int


  integer function atomic_load_int(storage)
    implicit none

    type(t_atomic_int), intent(in) :: storage

    atomic_load_int = c_atomic_load_int(storage%p)
  end function atomic_load_int


  subroutine atomic_store_int(storage, val)
    implicit none

    type(t_atomic_int), intent(in) :: storage
    integer, intent(in) :: val

    call c_atomic_store_int(storage%p, val)
  end subroutine atomic_store_int


  integer function atomic_fetch_and_increment_int(storage)
    implicit none

    type(t_atomic_int), intent(in) :: storage

    atomic_fetch_and_increment_int = c_atomic_fetch_and_increment_int(storage%p)
  end function atomic_fetch_and_increment_int


  integer function atomic_mod_increment_and_fetch_int(storage, mod)
    implicit none

    type(t_atomic_int), intent(in) :: storage
    integer, intent(in) :: mod

    atomic_mod_increment_and_fetch_int = c_atomic_mod_increment_and_fetch_int(storage%p, mod)
  end function atomic_mod_increment_and_fetch_int
  
end module module_atomic_ops
