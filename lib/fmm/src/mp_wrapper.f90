#include "fmm.h"
#include "fmmkinds.h"
#include "mp_wrapper.h"
!==================================================================
!ARMCI Bindings
!==================================================================
module armci_types
use iso_c_binding, only : c_long,c_ptr,c_double,c_int
implicit none
  type, bind(c) :: armci_hdl_t
    integer (c_long), dimension(4) :: mydata
    real (c_double), dimension(72) :: dummy
  end type armci_hdl_t

  type, bind(c) :: armci_giov_t
    type(c_ptr) :: src_ptr_array
    type(c_ptr) :: dst_ptr_array
    integer (c_long) :: ptr_array_length
    integer (c_long) :: bytes
  end type armci_giov_t

  ! armci and amrci-mpi define it as c_long, A1 as c_int
# if (FMM_MP == FMM_MP_ARMCI) || (FMM_MP == FMM_MP_ARMCIMPI)
  integer, parameter :: armci_size_t = c_long
# elif (FMM_MP == FMM_MP_A1)
  integer, parameter :: armci_size_t = c_int
# else
# error "communication library not supported"
# endif

end module armci_types

module armci_wrapper
use iso_c_binding, only : c_int,c_long,c_long_long,c_float,c_double,c_ptr,c_char,c_null_char,c_loc
use myarmci_constants
implicit none
  interface
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

    integer (c_int) function armci_init() bind(c,Name='ARMCI_Init')
    use iso_c_binding, only : c_int
    implicit none
    end function armci_init

    subroutine armci_finalize() bind(c,Name='ARMCI_Finalize')
    implicit none
    end subroutine armci_finalize

    integer (c_int) function armci_malloc(ptr,bsize) bind(c,Name='ARMCI_Malloc')
    use iso_c_binding, only : c_int,c_ptr
    use armci_types, only : armci_size_t
    implicit none
    type (c_ptr),dimension(*),target :: ptr
    integer (armci_size_t), value :: bsize
    end function armci_malloc

    integer (c_int) function armci_free(ptr) bind(c,Name='ARMCI_Free')
    use iso_c_binding, only : c_int,c_ptr
    implicit none
    type (c_ptr), value :: ptr
    end function armci_free

    integer (c_int) function armci_put(src,dst,bsize,proc) bind(c,name='ARMCI_Put')
    use iso_c_binding, only : c_int,c_ptr
    implicit none
    type (c_ptr), value :: src,dst
    integer (c_int), value :: bsize,proc
    end function armci_put

    subroutine armci_fence(proc) bind(c,name='ARMCI_Fence')
    use iso_c_binding, only : c_int
    implicit none
    integer (c_int), value :: proc
    end subroutine armci_fence

    subroutine armci_allfence() bind(c,name='ARMCI_AllFence')
    implicit none
    end subroutine armci_allfence

    subroutine armci_barrier() bind(c,name='ARMCI_Barrier')
    implicit none
    end subroutine armci_barrier

#   if FMM_MP == FMM_MP_ARMCI
    integer (c_int) function armci_notify(proc) bind(c,name='armci_notify')
    use iso_c_binding, only : c_int
    implicit none
    integer (c_int), value :: proc
    end function armci_notify

    integer (c_int) function armci_notify_wait(proc,pval) bind(c,name='armci_notify_wait')
    use iso_c_binding, only : c_int
    implicit none
    integer (c_int), value :: proc
    integer (c_int) :: pval
    end function armci_notify_wait
#   endif

    subroutine armci_cleanup() bind(c,name='ARMCI_Cleanup')
    implicit none
    end subroutine armci_cleanup

    subroutine armci_error(msg,ierr) bind(c,name='ARMCI_Error')
    use iso_c_binding, only : c_char, c_int
    implicit none
    character (c_char), dimension(*) :: msg
    integer (c_int), value :: ierr
    end subroutine armci_error
  end interface

  ! put overloading 
  interface mp_put
    module procedure mp_put_int
    module procedure mp_put_int_1d
    module procedure mp_put_int_2d
    module procedure mp_put_real4
    module procedure mp_put_real4_1d
    module procedure mp_put_real4_2d
    module procedure mp_put_real8
    module procedure mp_put_real8_1d
    module procedure mp_put_real8_2d
  end interface mp_put

contains
  subroutine mp_fence(proc)
  implicit none
  integer (FMM_INTEGER) :: proc
  integer (MyARMCI_Proc) :: proc_tmp
    proc_tmp = proc
    call armci_fence(proc_tmp)
  end subroutine mp_fence

  subroutine mp_allfence()
  implicit none
    call armci_allfence()
  end subroutine mp_allfence

# if FMM_MP == FMM_MP_ARMCI
  subroutine mp_notify(proc)
  implicit none
  integer (FMM_INTEGER) :: proc
  integer (MyARMCI_Proc) :: proc_tmp
  integer (MyARMCI_Waitcount) :: waitcount
    proc_tmp = proc
    waitcount = armci_notify(proc_tmp)
  end subroutine mp_notify

  subroutine mp_notifywait(proc,waitcount)
  implicit none
  integer (FMM_INTEGER) :: proc,waitcount
  integer (MyARMCI_Proc) :: proc_tmp
  integer (MyARMCI_Waitcount) :: waitcount_tmp, stat
    proc_tmp = proc
    waitcount_tmp = waitcount
    stat = armci_notify_wait(proc_tmp,waitcount_tmp)
    waitcount = stat
  end subroutine mp_notifywait
# endif

!==================================================================
!ARMCI Put Integer
!==================================================================

  subroutine mp_put_int(src,rptr,proc)
  implicit none
  integer (FMM_INTEGER), target :: src
  type (c_ptr), target :: rptr
  integer (FMM_INTEGER) :: proc
  integer (MyARMCI_Sendsize) :: sendsize
  integer (MyARMCI_Proc) :: proc_tmp
  integer (MyARMCI_Errorcode) :: ierr
  character (kind=c_char,len=40) :: msg
	
    proc_tmp = proc
    sendsize = kind(src)

#   ifdef FMM_CHECKSIZE
    if(sendsize.le.0) then
      msg = "FMM: mp_wrapper: put (FMM_INTEGER) sendsize <= 0" // c_null_char
      call armci_error(msg,ierr)
    endif
#   endif

    ierr = armci_put(c_loc(src),rptr,sendsize,proc_tmp)

    if (ierr.ne.MyARMCI_SUCCESS) then
      msg = "FMM: ARMCI: put (FMM_INTEGER) error" // c_null_char
      call armci_error(msg,ierr)
    end if
  end subroutine mp_put_int

  subroutine mp_put_int_1d(src,elem,rptr,proc)
  implicit none
  integer (FMM_INTEGER), target, dimension(*) :: src
  type (c_ptr), target :: rptr
  integer (FMM_INTEGER) :: elem,proc
  integer (MyARMCI_Sendsize) :: sendsize
  integer (MyARMCI_Proc) :: proc_tmp
  integer (MyARMCI_Errorcode) :: ierr
  character (kind=c_char,len=40) :: msg

    proc_tmp = proc
    sendsize = elem*kind(src)

#   ifdef FMM_CHECKSIZE
    if(sendsize.le.0) then
      msg = "FMM: mp_wrapper: put (FMM_INTEGER) sendsize <= 0" // c_null_char
      call armci_error(msg,ierr)
    endif
#endif

    ierr = armci_put(c_loc(src),rptr,sendsize,proc_tmp)

    if (ierr.ne.MyARMCI_SUCCESS) then
      msg = "FMM: ARMCI: put (FMM_INTEGER) 1d error" // c_null_char
      call armci_error(msg,ierr)
    end if
  end subroutine mp_put_int_1d
	
  subroutine mp_put_int_2d(src,elem1,elem2,rptr,proc)
  implicit none
  integer (FMM_INTEGER) :: elem1,elem2,proc
  integer (FMM_INTEGER), target, dimension(elem1,elem2) :: src
  type (c_ptr), target :: rptr
  integer (MyARMCI_Sendsize) :: sendsize
  integer (MyARMCI_Proc) :: proc_tmp	
  integer (MyARMCI_Errorcode) :: ierr
  character (kind=c_char,len=40) :: msg

    proc_tmp = proc
    sendsize = elem1*elem2*kind(src)

#   ifdef FMM_CHECKSIZE
    if(sendsize.le.0) then
      msg = "FMM: mp_wrapper: put (FMM_INTEGER) sendsize <= 0" // c_null_char
      call armci_error(msg,ierr)
    endif
#   endif

    ierr = armci_put(c_loc(src),rptr,sendsize,proc_tmp)

    if (ierr.ne.MyARMCI_SUCCESS) then
      msg = "FMM: ARMCI: put (FMM_INTEGER) 2d error" // c_null_char 
      call armci_error(msg,ierr)
    end if
  end subroutine mp_put_int_2d

!==================================================================
!ARMCI Put Real4
!==================================================================

	subroutine mp_put_real4(src,rptr,proc)
	implicit none
	real (c_float), target :: src
	type (c_ptr), target :: rptr
	integer (FMM_INTEGER) :: proc
	integer (MyARMCI_Sendsize) :: sendsize
	integer (MyARMCI_Proc) :: proc_tmp	
	integer (MyARMCI_Errorcode) :: ierr
	character (kind=c_char,len=40) :: msg

		proc_tmp = proc
		sendsize = kind(src)

#ifdef FMM_CHECKSIZE
			if(sendsize.le.0) then
			msg = "FMM: mp_wrapper: put (real4) sendsize <= 0" // c_null_char 
				call armci_error(msg,ierr)
			endif		
#endif		
		
		ierr = armci_put(c_loc(src),rptr,sendsize,proc_tmp)				

		if (ierr.ne.MyARMCI_SUCCESS) then
			msg = "FMM: ARMCI: put (real4) error" // c_null_char 
			call armci_error(msg,ierr)
		end if
	end subroutine mp_put_real4
	
	subroutine mp_put_real4_1d(src,elem,rptr,proc)
	implicit none
	real (c_float), target, dimension(*) :: src
	type (c_ptr), target :: rptr
	integer (FMM_INTEGER) :: elem,proc
	integer (MyARMCI_Sendsize) :: sendsize
	integer (MyARMCI_Proc) :: proc_tmp	
	integer (MyARMCI_Errorcode) :: ierr
	character (kind=c_char,len=40) :: msg
	
		proc_tmp = proc
		sendsize = elem*kind(src)
		
#ifdef FMM_CHECKSIZE
			if(sendsize.le.0) then
			msg = "FMM: mp_wrapper: put (real4) sendsize <= 0" // c_null_char 
				call armci_error(msg,ierr)
			endif		
#endif
		
		ierr = armci_put(c_loc(src),rptr,sendsize,proc_tmp)
		
		if (ierr.ne.MyARMCI_SUCCESS) then
			msg = "FMM: ARMCI: put (real4) 1d error" // c_null_char 
			call armci_error(msg,ierr)
		end if
	end subroutine mp_put_real4_1d

	subroutine mp_put_real4_2d(src,lb1,ub1,lb2,ub2,rptr,proc)
	implicit none
	integer (FMM_INTEGER) :: lb1,ub1,lb2,ub2,proc
	real (c_float), target, dimension(lb1:ub1,lb2:ub2) :: src
	real (c_float), pointer :: src_tmp
	type (c_ptr), target :: rptr
	integer (MyARMCI_Sendsize) :: sendsize
	integer (MyARMCI_Proc) :: proc_tmp	
	integer (MyARMCI_Errorcode) :: ierr
	character (kind=c_char,len=40) :: msg	
	
		proc_tmp = proc	
		sendsize = (ub1-lb1+1)*(ub2-lb2+1)*kind(src)

#ifdef FMM_CHECKSIZE
			if(sendsize.le.0) then
			msg = "FMM: mp_wrapper: put (real4) sendsize <= 0" // c_null_char 
				call armci_error(msg,ierr)
			endif		
#endif
				
                src_tmp => src(lb1,lb2)
		
		ierr = armci_put(c_loc(src_tmp),rptr,sendsize,proc_tmp)

		if (ierr.ne.MyARMCI_SUCCESS) then
			msg = "FMM: ARMCI: put (real4) 2d error" // c_null_char 
			call armci_error(msg,ierr)
		end if
	end subroutine mp_put_real4_2d		

!==================================================================
!ARMCI Put Real8
!==================================================================

	subroutine mp_put_real8(src,rptr,proc)
	implicit none
	real (c_double), target :: src
	type (c_ptr), target :: rptr
	integer (FMM_INTEGER) :: proc
	integer (MyARMCI_Sendsize) :: sendsize
	integer (MyARMCI_Proc) :: proc_tmp	
	integer (MyARMCI_Errorcode) :: ierr
	character (kind=c_char,len=40) :: msg
	
		proc_tmp = proc
		sendsize = kind(src)
		
#ifdef FMM_CHECKSIZE
			if(sendsize.le.0) then
			msg = "FMM: mp_wrapper: put (real8) sendsize <= 0" // c_null_char 
				call armci_error(msg,ierr)
			endif		
#endif				
		
		ierr = armci_put(c_loc(src),rptr,sendsize,proc_tmp)
		
		if (ierr.ne.MyARMCI_SUCCESS) then
			msg = "FMM: ARMCI: put (real8) error" // c_null_char 
			call armci_error(msg,ierr)
		end if
	end subroutine mp_put_real8
	
	subroutine mp_put_real8_1d(src,elem,rptr,proc)
	implicit none
	integer (FMM_INTEGER) :: elem,proc	
	real (c_double), target, dimension(1:elem) :: src	
	type (c_ptr), target :: rptr
	integer (MyARMCI_Sendsize) :: sendsize
	integer (MyARMCI_Proc) :: proc_tmp	
	integer (MyARMCI_Errorcode) :: ierr
	character (kind=c_char,len=40) :: msg
	
		proc_tmp = proc
		sendsize = elem*kind(src)

#ifdef FMM_CHECKSIZE
			if(sendsize.le.0) then
			msg = "FMM: mp_wrapper: put (real8) sendsize <= 0" // c_null_char 
				call armci_error(msg,ierr)
			endif		
#endif		
		
		ierr = armci_put(c_loc(src),rptr,sendsize,proc_tmp)
		
		if (ierr.ne.MyARMCI_SUCCESS) then
			msg = "FMM: ARMCI: put (real8) 1d error" // c_null_char 
			call armci_error(msg,ierr)
		end if
	end subroutine mp_put_real8_1d
	
	subroutine mp_put_real8_2d(src,lb1,ub1,lb2,ub2,rptr,proc)
	implicit none
	integer (FMM_INTEGER) :: lb1,ub1,lb2,ub2,proc
	real (c_double), target, dimension(lb1:ub1,lb2:ub2) :: src
	real (c_double), pointer :: src_tmp
	type (c_ptr), target :: rptr
	integer (MyARMCI_Sendsize) :: sendsize
	integer (MyARMCI_Proc) :: proc_tmp	
	integer (MyARMCI_Errorcode) :: ierr
	character (kind=c_char,len=40) :: msg
		
		proc_tmp = proc	
		sendsize = (ub1-lb1+1)*(ub2-lb2+1)*kind(src)
		
#ifdef FMM_CHECKSIZE
			if(sendsize.le.0) then
			msg = "FMM: mp_wrapper: put (real8) sendsize <= 0" // c_null_char 
				call armci_error(msg,ierr)
			endif		
#endif		
		
                src_tmp => src(lb1,lb2)
		
		ierr = armci_put(c_loc(src_tmp),rptr,sendsize,proc_tmp)
		
		if (ierr.ne.MyARMCI_SUCCESS) then
			msg = "FMM: ARMCI: put (real8) 2d error" // c_null_char 
			call armci_error(msg,ierr)
		end if
	end subroutine mp_put_real8_2d
end module armci_wrapper

!==================================================================
!MPI Bindings
!==================================================================
module mpi_wrapper
use mympi_constants
use armci_wrapper, only : armci_cleanup
implicit none
	!ToDo initialization interface, do not use save variables
	integer(8), dimension(10042) :: notifyerbuffer
	integer(FMM_INTEGER) :: attached = 0
	
	interface mp_allreduce
		module procedure mp_allreduce_int8	
		module procedure mp_allreduce_int8_1d
		module procedure mp_allreduce_real4		
		module procedure mp_allreduce_real4_1d
                module procedure mp_allreduce_real4_2d
		module procedure mp_allreduce_real8		
		module procedure mp_allreduce_real8_1d		
		module procedure mp_allreduce_real8_2d		
	end interface

	interface mp_allgather
		module procedure mp_allgather_byte_2d
		module procedure mp_allgather_byte_3d
		module procedure mp_allgather_real4
		module procedure mp_allgather_int8		
		module procedure mp_allgather_int8_2d		
	end interface
	
	interface mp_walltime
	 module procedure mp_walltime_real4
	 module procedure mp_walltime_real8
#ifdef FMM_TIMER16
	 module procedure mp_walltime_real16
#endif
	end interface mp_walltime
	
	contains
	
!==================================================================
!MPI Allreduce
!==================================================================	

	subroutine mp_allreduce_int8(dst,op,comm)
	implicit none
	integer(kind=8)			:: dst
	integer(MyMPI_Op)		:: op
	integer(MyMPI_Comm) 		:: comm
	integer(MyMPI_Errorcode)	:: ierr,ierr2
	integer(MyMPI_Entries)		:: elem_tmp
	
		elem_tmp = 1
		
		call mpi_allreduce(MPI_IN_PLACE,dst,elem_tmp,MPI_INTEGER8,op,comm,ierr)
		
		if (ierr.ne.MPI_SUCCESS) then
			call armci_cleanup()
			call mpi_abort(comm,ierr,ierr2)
		endif
	end subroutine mp_allreduce_int8
	
	subroutine mp_allreduce_int8_1d(dst,elem,op,comm)
	implicit none
	integer(kind=8), dimension(*)	:: dst
	integer(FMM_INTEGER)		:: elem
	integer(MyMPI_Op)		:: op
	integer(MyMPI_Comm) 		:: comm
	integer(MyMPI_Errorcode)	:: ierr,ierr2
	integer(MyMPI_Entries)		:: elem_tmp
	
		elem_tmp = elem
		
		call mpi_allreduce(MPI_IN_PLACE,dst,elem_tmp,MPI_INTEGER8,op,comm,ierr)
		
		if (ierr.ne.MPI_SUCCESS) then
			call armci_cleanup()
			call mpi_abort(comm,ierr,ierr2)
		endif
	end subroutine mp_allreduce_int8_1d
	
	subroutine mp_allreduce_real4(dst,op,comm)
	implicit none
	real(kind=4)		:: dst
	integer(MyMPI_Op)	:: op
	integer(MyMPI_Comm) 	:: comm
	integer(MyMPI_Errorcode):: ierr,ierr2
	integer(MyMPI_Entries)	:: elem_tmp
	
		elem_tmp = 1
		
		call mpi_allreduce(MPI_IN_PLACE,dst,elem_tmp,MPI_REAL4,op,comm,ierr)
		
		if (ierr.ne.MPI_SUCCESS) then
			call armci_cleanup()
			call mpi_abort(comm,ierr,ierr2)
		endif
	end subroutine mp_allreduce_real4
	
	subroutine mp_allreduce_real4_1d(dst,elem,op,comm)
	implicit none
	real(kind=4), dimension(*)	:: dst
	integer(FMM_INTEGER)		:: elem
	integer(MyMPI_Op)		:: op
	integer(MyMPI_Comm) 		:: comm
	integer(MyMPI_Errorcode)	:: ierr,ierr2
	integer(MyMPI_Entries)		:: elem_tmp
	
		elem_tmp = elem
		
		call mpi_allreduce(MPI_IN_PLACE,dst,elem_tmp,MPI_REAL4,op,comm,ierr)
		
		if (ierr.ne.MPI_SUCCESS) then
			call armci_cleanup()
			call mpi_abort(comm,ierr,ierr2)
		endif
	end subroutine mp_allreduce_real4_1d

        subroutine mp_allreduce_real4_2d(dst,elem,op,comm)
        implicit none
        real(kind=4), dimension(:,:)    :: dst
        integer(FMM_INTEGER)            :: elem,lo,hi
        integer(MyMPI_Op)               :: op
        integer(MyMPI_Comm)             :: comm
        integer(MyMPI_Errorcode)        :: ierr,ierr2
        integer(MyMPI_Entries)          :: elem_tmp
        
                lo = lbound(dst,1)
                hi = ubound(dst,1)
                elem_tmp = elem*(hi-lo+1)
                
                call mpi_allreduce(MPI_IN_PLACE,dst,elem_tmp,MPI_REAL4,op,comm,ierr)
                
                if (ierr.ne.MPI_SUCCESS) then
                        call armci_cleanup()
                        call mpi_abort(comm,ierr,ierr2)
                endif
        end subroutine mp_allreduce_real4_2d    

	subroutine mp_allreduce_real8(dst,op,comm)
	implicit none
	real(kind=8)		:: dst
	integer(MyMPI_Op)	:: op
	integer(MyMPI_Comm) 	:: comm
	integer(MyMPI_Errorcode):: ierr,ierr2
	integer(MyMPI_Entries)	:: elem_tmp
	
		elem_tmp = 1
		
		call mpi_allreduce(MPI_IN_PLACE,dst,elem_tmp,MPI_REAL8,op,comm,ierr)
	
		if (ierr.ne.MPI_SUCCESS) then
			call armci_cleanup()
			call mpi_abort(comm,ierr,ierr2)
		endif
	end subroutine mp_allreduce_real8
		
	subroutine mp_allreduce_real8_1d(dst,elem,op,comm)
	implicit none
	real(kind=8), dimension(*)	:: dst
	integer(FMM_INTEGER)		:: elem
	integer(MyMPI_Op)		:: op
	integer(MyMPI_Comm) 		:: comm
	integer(MyMPI_Errorcode)	:: ierr,ierr2
	integer(MyMPI_Entries)		:: elem_tmp
	
		elem_tmp = elem
		
		call mpi_allreduce(MPI_IN_PLACE,dst,elem_tmp,MPI_REAL8,op,comm,ierr)
		
		if (ierr.ne.MPI_SUCCESS) then
			call armci_cleanup()
			call mpi_abort(comm,ierr,ierr2)
		endif
	end subroutine mp_allreduce_real8_1d
	
	subroutine mp_allreduce_real8_2d(dst,elem,op,comm)
	implicit none
	real(kind=8), dimension(:,:)	:: dst
	integer(FMM_INTEGER)		:: elem,lo,hi
	integer(MyMPI_Op)		:: op
	integer(MyMPI_Comm) 		:: comm
	integer(MyMPI_Errorcode)	:: ierr,ierr2
	integer(MyMPI_Entries)		:: elem_tmp
	
		lo = lbound(dst,1)
		hi = ubound(dst,1)
		elem_tmp = elem*(hi-lo+1)
		
		call mpi_allreduce(MPI_IN_PLACE,dst,elem_tmp,MPI_REAL8,op,comm,ierr)
		
		if (ierr.ne.MPI_SUCCESS) then
			call armci_cleanup()
			call mpi_abort(comm,ierr,ierr2)
		endif
	end subroutine mp_allreduce_real8_2d	
	
!==================================================================
!MPI Allgather
!==================================================================		

	subroutine mp_allgather_byte_2d(dst,elem,comm)
	implicit none
	byte, dimension(:,:)	:: dst
	integer(FMM_INTEGER)	:: elem,lo,hi
	integer(MyMPI_Comm) 	:: comm
	integer(MyMPI_Errorcode):: ierr,ierr2
	integer(MyMPI_Entries)	:: elem_tmp
	
		lo = lbound(dst,1)
		hi = ubound(dst,1)	
		elem_tmp = elem*(hi-lo+1)
		
		call mpi_allgather(MPI_IN_PLACE,elem_tmp,MPI_BYTE,dst,elem_tmp,MPI_BYTE,comm,ierr)
		
		if (ierr.ne.MPI_SUCCESS) then
			call armci_cleanup()
			call mpi_abort(comm,ierr,ierr2)
		endif
	end subroutine mp_allgather_byte_2d

	subroutine mp_allgather_byte_3d(dst,elem,comm)
	implicit none
	byte, dimension(:,:,:)	:: dst
	integer(FMM_INTEGER)	:: elem,lo1,hi1,lo2,hi2
	integer(MyMPI_Comm) 	:: comm
	integer(MyMPI_Errorcode):: ierr,ierr2
	integer(MyMPI_Entries)	:: elem_tmp
	
		lo1 = lbound(dst,1)
		hi1 = ubound(dst,1)
		lo2 = lbound(dst,2)
		hi2 = ubound(dst,2)			
		elem_tmp = elem*(hi1-lo1+1)*(hi2-lo2+1)
		
		call mpi_allgather(MPI_IN_PLACE,elem_tmp,MPI_BYTE,dst,elem_tmp,MPI_BYTE,comm,ierr)
		
		if (ierr.ne.MPI_SUCCESS) then
			call armci_cleanup()
			call mpi_abort(comm,ierr,ierr2)
		endif
	end subroutine mp_allgather_byte_3d	
	
	subroutine mp_allgather_real4(dst,elem,comm)
	implicit none
	real(kind=4), dimension(*)	:: dst
	integer(FMM_INTEGER)		:: elem
	integer(MyMPI_Comm) 		:: comm
	integer(MyMPI_Errorcode)	:: ierr,ierr2
	integer(MyMPI_Entries)		:: elem_tmp
	
		elem_tmp = elem
		
		call mpi_allgather(MPI_IN_PLACE,elem_tmp,MPI_REAL4,dst,elem_tmp,MPI_REAL4,comm,ierr)
		
		if (ierr.ne.MPI_SUCCESS) then
			call armci_cleanup()
			call mpi_abort(comm,ierr,ierr2)
		endif
	end subroutine mp_allgather_real4
	
	subroutine mp_allgather_int8(dst,elem,comm)
	implicit none
	integer(kind=8), dimension(*)	:: dst
	integer(FMM_INTEGER)		:: elem
	integer(MyMPI_Comm) 		:: comm
	integer(MyMPI_Errorcode)	:: ierr,ierr2
	integer(MyMPI_Entries)		:: elem_tmp
	
		elem_tmp = elem
		
		call mpi_allgather(MPI_IN_PLACE,elem_tmp,MPI_INTEGER8,dst,elem_tmp,MPI_INTEGER8,comm,ierr)
		
		if (ierr.ne.MPI_SUCCESS) then
			call armci_cleanup()
			call mpi_abort(comm,ierr,ierr2)
		endif
	end subroutine mp_allgather_int8
	
	subroutine mp_allgather_int8_2d(dst,elem,comm)
	implicit none
	integer(kind=8), dimension(:,:)	:: dst
	integer(FMM_INTEGER)		:: elem,lo,hi
	integer(MyMPI_Comm) 		:: comm
	integer(MyMPI_Errorcode)	:: ierr,ierr2
	integer(MyMPI_Entries)		:: elem_tmp
	
		lo = lbound(dst,1)
		hi = ubound(dst,1)	
		elem_tmp = elem*(hi-lo+1)
		
		call mpi_allgather(MPI_IN_PLACE,elem_tmp,MPI_INTEGER8,dst,elem_tmp,MPI_INTEGER8,comm,ierr)
		
		if (ierr.ne.MPI_SUCCESS) then
			call armci_cleanup()
			call mpi_abort(comm,ierr,ierr2)
		endif
	end subroutine mp_allgather_int8_2d						

!==================================================================
!MPI WTime
!==================================================================
	
	subroutine mp_walltime_real4(mytime)
	implicit none
	real (kind=4) :: mytime
	real (MyMPI_Timer) :: mytime_tmp
		mytime_tmp = mpi_wtime()
		mytime = mytime_tmp
	end subroutine mp_walltime_real4
		subroutine mp_walltime_real8(mytime)
	implicit none
	real (kind=8) :: mytime
	real (MyMPI_Timer) :: mytime_tmp			
		mytime_tmp = mpi_wtime()
		mytime = mytime_tmp			
	end subroutine mp_walltime_real8
#ifdef FMM_TIMER16
	subroutine mp_walltime_real16(mytime)
	implicit none
	real (kind=16) :: mytime
	real (MyMPI_Timer) :: mytime_tmp
		mytime_tmp = mpi_wtime()
		mytime = mytime_tmp
	end subroutine mp_walltime_real16
#endif

!==================================================================
!MPI Notify/NotifyWait on BGP
!==================================================================

#if FMM_MP == FMM_MP_A1 || FMM_MP == FMM_MP_ARMCIMPI
	subroutine mp_notify(proc)
	implicit none
	integer (FMM_INTEGER) :: proc
	integer (MyMPI_Rank) :: proc_tmp
	integer (MyMPI_Entries) :: elem
	integer (MyMPI_Errorcode) :: ierr
	integer (kind=4),save :: notifyersnd = 0
	integer (kind=4) :: mytag
		mytag = 42
		elem = 1
		proc_tmp = proc
		notifyersnd = notifyersnd + 1
		call mpi_bsend(notifyersnd,elem,MPI_INTEGER4,proc_tmp,mytag,MPI_COMM_WORLD,ierr)
	end subroutine mp_notify
	
	subroutine mp_notifywait(proc,waitcount)
	implicit none
	integer (FMM_INTEGER) :: proc,waitcount
	integer (MyMPI_Rank) :: proc_tmp
	integer (MyMPI_Entries) :: elem
	integer(4), dimension(MPI_STATUS_SIZE) :: rcvstatus
	integer (MyMPI_Errorcode) :: ierr
	integer (kind=4) :: notifyerrcv
	integer (kind=4) :: mytag
		mytag = 42
		elem = 1
		proc_tmp = proc
		call mpi_recv(notifyerrcv,elem,MPI_INTEGER4,proc_tmp,mytag,MPI_COMM_WORLD,rcvstatus,ierr)
		waitcount = notifyerrcv
	end subroutine mp_notifywait
#endif				
end module mpi_wrapper

!==================================================================
!FMM Bindings
!==================================================================
module mp_wrapper
use iso_c_binding
use armci_wrapper
use armci_types
use mpi_wrapper
implicit none

	private

	! module mpi_wrapper -> module mpi
    integer(MyMPI_Comm), public, parameter :: MP_COMM = 8!kind(MyMPI_Comm)
	integer(MyMPI_Comm), public 	:: MP_ALLNODES = MPI_COMM_WORLD
	integer(MyMPI_Op), public, parameter 	:: MP_MIN = MPI_MIN
	integer(MyMPI_Op), public, parameter 	:: MP_MAX = MPI_MAX	
	integer(MyMPI_Op), public, parameter 	:: MP_SUM = MPI_SUM

	integer(MyMPI_Datatype), parameter 	:: MP_REAL4 = MPI_REAL
	integer(MyMPI_Datatype), parameter 	:: MP_REAL8 = MPI_DOUBLE_PRECISION
	integer(MyMPI_Datatype), parameter 	:: MP_INT4 = MPI_INTEGER
	integer(MyMPI_Datatype), parameter	:: MP_INT8 = MPI_INTEGER8

	!Communicator
	public :: MPI_COMM_WORLD
	public :: MPI_COMM_SELF
		
	!module armci_types
!	public :: armci_hdl_t
!	public :: armci_size_t	

	!mp_wrapper public functions
	public :: mp_init
	public :: mp_error
	public :: mp_rank
	public :: mp_nnodes
	public :: mp_allocate
	public :: mp_deallocate
	public :: mp_finalize
	public :: mp_notify,mp_notifywait
	public :: mp_barrier
    public :: mp_fence
    public :: mp_allfence
	public :: mp_put
	public :: mp_allgather
	public :: mp_allreduce
	public :: mp_walltime

	! c pointer size function
	public :: ctob
	! difference of two c pointers
	public :: diffcpointers

	! casting of external pointer arithmetic functions
	interface
		! offset 4 byte integer
		pure type (c_ptr) function ptroffset_pi(ptr,offset) bind(c,name='ptroffset_pi')
		use iso_c_binding
		implicit none
		type (c_ptr), value, intent(in) :: ptr
		integer (c_int), value, intent(in) :: offset 
		end function ptroffset_pi

		! offset 8 byte integer
		pure type (c_ptr) function ptroffset_pll(ptr,offset) bind(c,name='ptroffset_pll')
		use iso_c_binding
		implicit none
		type (c_ptr), value, intent(in) :: ptr
		integer (c_long_long), value, intent(in) :: offset 
		end function ptroffset_pll
	end interface

	! pointer arithmetic overloading into .add. operator
	interface operator (.add.)
		! offset -> c_int
		module procedure fort_ptroffset_pi
		! offset -> c_long_long
		module procedure fort_ptroffset_pll
	end interface operator (.add.)

	! mp_wrapper public operators
	public :: operator (.add.)
		
	interface
		pure integer (c_int) function cptrsize_i() bind(c,name='cptrsize_i')
		use iso_c_binding
		implicit none
		end function cptrsize_i
		
		pure integer (c_long_long) function cptrsize_ll() bind(c,name='cptrsize_ll')
		use iso_c_binding	
		implicit none
		end function cptrsize_ll		
	end interface
	
	interface ctob
		module procedure cptrsize_int
		module procedure cptrsize_longlong
	end interface ctob
	
	interface
		pure integer (c_int) function diffcpointers_i(ptr1,ptr2) bind(c,name='diffcpointers_i')
		use iso_c_binding
		implicit none
		type (c_ptr), value, intent(in) :: ptr1,ptr2
		end function diffcpointers_i
		
		pure integer (c_long_long) function diffcpointers_ll(ptr1,ptr2) bind(c,name='diffcpointers_ll')
		use iso_c_binding	
		implicit none
		type (c_ptr), value, intent(in) :: ptr1,ptr2		
		end function diffcpointers_ll
	end interface
	
	interface diffcpointers
		module procedure diffcpointers_int
		module procedure diffcpointers_longlong
	end interface diffcpointers
	
	interface mp_error
		module procedure mp_error_int4
		module procedure mp_error_int8		
	end interface mp_error		
	
contains

	type (c_ptr) function fort_ptroffset_pi(ptr,offset)
	use iso_c_binding
	implicit none
	type (c_ptr),intent(in) :: ptr
	type (c_ptr) :: tmpptr
	integer (c_int), value, intent(in) :: offset
		tmpptr = ptr
		fort_ptroffset_pi = ptroffset_pi(tmpptr,offset)
	end function fort_ptroffset_pi	

	type (c_ptr) function fort_ptroffset_pll(ptr,offset)
	use iso_c_binding
	implicit none
	type (c_ptr),intent(in) :: ptr
	type (c_ptr) :: tmpptr
	integer (c_long_long), value, intent(in) :: offset
		tmpptr = ptr
		fort_ptroffset_pll = ptroffset_pll(tmpptr,offset)
	end function fort_ptroffset_pll		

	subroutine diffcpointers_int(ptr1,ptr2,diff)
	use iso_c_binding
	implicit none
	type (c_ptr) :: ptr1,ptr2
	integer (c_int) :: diff
		diff = diffcpointers_i(ptr1,ptr2)
	end subroutine diffcpointers_int

	subroutine diffcpointers_longlong(ptr1,ptr2,diff)
	use iso_c_binding
	implicit none
	type (c_ptr) :: ptr1,ptr2
	integer (c_long_long) :: diff
		diff = diffcpointers_ll(ptr1,ptr2)
	end subroutine diffcpointers_longlong
	
	subroutine cptrsize_int(ptrsize)
	use iso_c_binding
	implicit none
	integer (c_int) :: ptrsize
		ptrsize = cptrsize_i()
	end subroutine cptrsize_int
	
	subroutine cptrsize_longlong(ptrsize)
	use iso_c_binding
	implicit none
	integer (c_long_long) :: ptrsize
		ptrsize = cptrsize_ll()
	end subroutine cptrsize_longlong	
	
	subroutine mp_init()
	implicit none
	integer (MyMPI_Errorcode) :: ierr1
	integer (MyARMCI_Errorcode) ::ierr2
	logical (MyMPI_Initialized) :: initialized
	integer (4) :: elems
!		call mpi_initialized(initialized,ierr1)
!		if (ierr1.ne.MPI_SUCCESS) call mp_error(ierr1)
!		if (initialized) then		
!			ierr2 = armci_init()
!			if (ierr2.ne.MyARMCI_SUCCESS) call mp_error(ierr2)		
!		else
!			call mpi_init(ierr1)
!			if (ierr1.ne.MPI_SUCCESS) call mp_error(ierr1)
				ierr2 = armci_init()
!			if (ierr2.ne.0) call mp_error(ierr2)
!		endif
		
#if FMM_MP == FMM_MP_A1 || FMM_MP == FMM_MP_ARMCIMPI
		elems = 10042
		if(attached.eq.0) then
			call mpi_buffer_attach(notifyerbuffer,elems,ierr1)
			attached = 1
		endif
#endif
	end subroutine mp_init


	subroutine mp_error_int4(ierr)
	implicit none
	integer (kind=4) :: ierr
	integer (MyMPI_Errorcode) :: ierr_tmp,ierr_tmp2
	integer (MyMPI_Rank) :: me
	character(len=100) :: errmsg
		ierr_tmp = ierr
		
		errmsg = "FMM: Error in parallel environment"
		call mpi_comm_rank(MPI_COMM_WORLD,me,ierr_tmp2)
		write(6,*) "Node:",me,trim(errmsg),ierr_tmp
		
		call armci_cleanup()
		call mpi_abort(MPI_COMM_WORLD,ierr_tmp,ierr_tmp2)
	end subroutine mp_error_int4

	subroutine mp_error_int8(ierr)
	implicit none
	integer (kind=8) :: ierr
	integer (MyMPI_Errorcode) :: ierr_tmp,ierr_tmp2
	integer (MyMPI_Rank) :: me
	character(len=100) :: errmsg
		ierr_tmp = ierr
		
		errmsg = "FMM: Error in parallel environment"
		call mpi_comm_rank(MPI_COMM_WORLD,me,ierr_tmp2)
		write(6,*) "Node:",me,trim(errmsg),ierr_tmp
		
		call armci_cleanup()
		call mpi_abort(MPI_COMM_WORLD,ierr_tmp,ierr_tmp2)
	end subroutine mp_error_int8		

	subroutine mp_rank(comm,me)
	implicit none
	integer (MyMPI_Comm) :: comm
	integer (FMM_INTEGER) :: me
	integer (MyMPI_Errorcode) :: ierr
	integer (MyMPI_Rank) :: me_tmp
	
		call mpi_comm_rank(comm,me_tmp,ierr)
		if (ierr.ne.MPI_SUCCESS) call mp_error(ierr)
		me = me_tmp
	end subroutine mp_rank

	subroutine mp_nnodes(comm,nnodes)
	implicit none
	integer (MyMPI_Comm) :: comm	
	integer (FMM_INTEGER) :: nnodes
	integer (MyMPI_Errorcode) :: ierr
        integer (MyMPI_Rank) :: nnodes_tmp
			
		call mpi_comm_size(comm,nnodes_tmp,ierr)
		if (ierr.ne.MPI_SUCCESS) call mp_error(ierr)
		nnodes = nnodes_tmp
	end subroutine mp_nnodes
	
	subroutine mp_allocate(ptr,bsize,ierr)
	implicit none
	type (c_ptr), dimension(*), target :: ptr
	integer(FMM_INTEGER) :: bsize,ierr
	integer (MyARMCI_Errorcode) :: ierr_tmp
	integer (armci_size_t) :: bsize_tmp
	integer (c_int) :: bsize_tmp2,ierr_tmp2

		if(MPI_COMM_WORLD.eq.MP_ALLNODES) then
  	  	    bsize_tmp = abs(bsize)
		    ierr_tmp = armci_malloc(ptr,bsize_tmp)
		    ierr = ierr_tmp
                elseif(MPI_COMM_SELF.eq.MP_ALLNODES) then
		 if(bsize.ne.0) then
		  bsize_tmp2 = abs(bsize)
  		  ptr(1) = dummy_malloc(bsize_tmp2,ierr_tmp2)
		  ierr = ierr_tmp2
		 else
		  ierr = 0 
		 endif 
                endif
	end subroutine mp_allocate

	subroutine mp_deallocate(ptr,ierr)
	implicit none
	type (c_ptr), value :: ptr
	integer (FMM_INTEGER) :: ierr,i
	integer (MyARMCI_Errorcode) :: ierr_tmp

	  if(MPI_COMM_WORLD.eq.MP_ALLNODES) then
#if FMM_MP == FMM_MP_A1 || FMM_MP == FMM_MP_ARMCIMPI
      ierr_tmp = armci_free(ptr)
#elif FMM_MP == FMM_MP_ARMCI
        if(ierr.eq.0) then
          ierr_tmp = 0
        else
	      ierr_tmp = armci_free(ptr)
        	  endif
#endif
		elseif(MPI_COMM_SELF.eq.MP_ALLNODES) then
		 if(ierr.eq.0) then
              ierr_tmp = 0
         else
		   call dummy_free(ptr)
		   ierr_tmp = 0
		 endif  
	  endif

		ierr = ierr_tmp
	end subroutine mp_deallocate

	subroutine mp_finalize()
	implicit none
	integer (MyMPI_Errorcode) :: ierr
	integer (4) :: elems	
#if FMM_MP == FMM_MP_A1 || FMM_MP == FMM_MP_ARMCIMPI
		elems = 10042
		if (attached.eq.1) then
			call mpi_buffer_detach(notifyerbuffer,elems,ierr)
			attached = 0
		endif
#endif
		
		call armci_finalize()
!		call mpi_finalize(ierr)
	end subroutine mp_finalize

	subroutine mp_barrier(nnodes)
	implicit none
	integer(FMM_INTEGER) :: nnodes
		! armci_barrier calls mpi_barrier intrinsicly
                if(nnodes.gt.1) then
                   call armci_barrier()
                endif
	end subroutine mp_barrier
	
end module mp_wrapper
