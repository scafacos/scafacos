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

module pthreads_stuff
  use, intrinsic :: iso_c_binding
  implicit none

  integer, parameter :: THREAD_TYPE_DEFAULT = 0
  integer, parameter :: THREAD_TYPE_COMMUNICATOR = 1
  integer, parameter :: THREAD_TYPE_WORKER = 2

  type, bind(C) :: t_pthread_with_type
    private
    type(c_ptr) :: thread
    type(c_funptr) :: start_routine
    type(c_ptr) :: arg
    integer(c_int) :: thread_type
    integer(c_int) :: counter
  end type

  interface
    integer(c_int) function get_my_core() bind(C, name='get_my_core')
      use, intrinsic :: iso_c_binding
      implicit none
    end function
  end interface

  interface  
    integer(c_int) function set_prefetching() bind(C, name='set_prefetching')
      use, intrinsic :: iso_c_binding
      implicit none
    end function
  end interface


  interface
    integer(c_int) function pthreads_init() bind(C, name='pthreads_init')
      use, intrinsic :: iso_c_binding
      implicit none
    end function

    integer(c_int) function pthreads_uninit() bind(C, name='pthreads_uninit')
      use, intrinsic :: iso_c_binding
      implicit none
    end function
  end interface


  interface
    integer(c_int) function pthreads_createthread_c(thread) bind(C, name='pthreads_createthread_c')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), intent(in), value :: thread
    end function

    integer(c_int) function pthreads_jointhread(thread) bind(C, name='pthreads_jointhread')
      use, intrinsic :: iso_c_binding
      import :: t_pthread_with_type
      implicit none
      type(t_pthread_with_type), intent(in), value :: thread
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

  integer(c_int) function pthreads_createthread(thread, start_routine, arg, thread_type, counter)
    use, intrinsic :: iso_c_binding
    implicit none

    type(t_pthread_with_type), target, intent(inout) :: thread
    type(c_funptr), intent(in) :: start_routine
    type(c_ptr), intent(in) :: arg
    integer, optional, intent(in) :: thread_type
    integer, optional, intent(in) :: counter

    thread%start_routine = start_routine
    thread%arg = arg

    thread%thread_type = THREAD_TYPE_DEFAULT
    if (present(thread_type)) thread%thread_type = thread_type

    thread%counter = 0
    if (present(counter)) thread%counter = counter

    pthreads_createthread = pthreads_createthread_c(c_loc(thread))
  end function


  function getfullid()
    implicit none
    character(20) :: getfullid

    write(getfullid,'("{", I8, ".", I8, "}")') get_my_pid(), get_my_tid()
  end function
end module pthreads_stuff
