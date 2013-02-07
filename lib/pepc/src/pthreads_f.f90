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

module pthreads_stuff
  implicit none

  interface
    integer(c_int) function get_my_core() bind(C, name='get_my_core')
      use, intrinsic :: iso_c_binding
      implicit none
    end function
  end interface


  interface

    integer(c_int) function pthreads_init(numthreads) bind(C, name='pthreads_init')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), intent(in), value :: numthreads
    end function

     integer(c_int) function pthreads_uninit() bind(C, name='pthreads_uninit')
      use, intrinsic :: iso_c_binding
      implicit none
    end function
  end interface

  interface
    integer(c_int) function pthreads_createthread(id, start_routine, arg) bind(C, name='pthreads_createthread')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int),   intent(in), value :: id
      type( c_funptr ), intent(in), value :: start_routine
      type( c_ptr ),    intent(in), value :: arg
    end function

    integer(c_int) function pthreads_jointhread(id) bind(C, name='pthreads_jointhread')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), intent(in), value :: id
    end function

    integer(c_int) function pthreads_exitthread() bind(C, name='pthreads_exitthread')
      use, intrinsic :: iso_c_binding
      implicit none
    end function

    integer(c_int) function pthreads_sched_yield() bind(C, name='pthreads_sched_yield')
      use, intrinsic :: iso_c_binding
      implicit none
    end function

  end interface




  interface

    integer(c_int) function rwlocks_init(numlocks) bind(C, name='rwlocks_init')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), intent(in), value :: numlocks
    end function

    integer(c_int) function rwlocks_uninit() bind(C, name='rwlocks_uninit')
      use, intrinsic :: iso_c_binding
      implicit none
    end function

    integer(c_int) function rwlocks_wrlock(id) bind(C, name='rwlocks_wrlock')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), intent(in), value :: id
    end function

    integer(c_int) function rwlocks_rdlock(id) bind(C, name='rwlocks_rdlock')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), intent(in), value :: id
    end function

    integer(c_int) function rwlocks_unlock(id) bind(C, name='rwlocks_unlock')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), intent(in), value :: id
    end function

  end interface


  interface

    integer(c_int) function get_my_tid() bind(C, name='get_my_tid')
      use, intrinsic :: iso_c_binding
      implicit none
    end function

    integer(c_int) function get_my_pid() bind(C, name='get_my_pid')
      use, intrinsic :: iso_c_binding
      implicit none
    end function

  end interface


  contains

    function getfullid()
      implicit none
      character(20) :: getfullid

      write(getfullid,'("{", I8, ".", I8, "}")') get_my_pid(), get_my_tid()
   end function

end module pthreads_stuff
