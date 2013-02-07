#include "fmm.h"

#ifdef FMM_PARALLEL

module pvlist
use fmmkinds
use fmmalloc
use mp_info
implicit none
private

type trootelem
 type(telem), pointer :: first,last,current
 integer(kind=fmm_integer) :: lastidx,currentidx,defaultval
 integer(kind=fmm_integer), allocatable, dimension(:) :: mylist
end type trootelem

type(telem), target :: dummyelem

public :: trootelem
public :: initlist
public :: destroylist
public :: setlist
public :: setelem
public :: getelem
public :: addelem

contains

 subroutine initlist(rootelem,myval)
 implicit none
 type(trootelem), intent(inout) :: rootelem
 integer(kind=fmm_integer), intent(in) :: myval
 integer(kind=fmm_integer) :: i

 call fmmallocate(rootelem%mylist,0,nnodes-1,i)
 if(i.ne.0) call bummer('fmm: initlist error, allocation stat =',i)

 rootelem%mylist = myval

 end subroutine initlist

 subroutine destroylist(rootelem)
 implicit none
 type(trootelem), intent(inout) :: rootelem
 integer(kind=fmm_integer) :: loopidx,i

 if(allocated(rootelem%mylist)) then
   call fmmdeallocate(rootelem%mylist,i)
   if(i.ne.0) call bummer('fmm: destroylist error, deallocation stat =',i)
 endif
 end subroutine destroylist 

 subroutine setlist(rootelem,myval)
 implicit none
 type(trootelem), intent(inout) :: rootelem
 integer(kind=fmm_integer), intent(in) :: myval
 integer(kind=fmm_integer) :: loopidx

 if(allocated(rootelem%mylist)) then
   rootelem%mylist = myval
 endif
 end subroutine setlist 
 
 subroutine setelem(rootelem,mypos,myval)
 implicit none
 type(trootelem), intent(inout) :: rootelem
 integer(kind=fmm_integer), intent(in) :: mypos,myval
 type(telem), pointer :: foundelem,newelem
 integer(kind=fmm_integer) :: found,i

 rootelem%mylist(mypos) = myval

 end subroutine setelem

 subroutine addelem(rootelem,mypos,myval)
 implicit none
 type(trootelem), intent(inout) :: rootelem
 integer(kind=fmm_integer), intent(in) :: mypos,myval
 type(telem), pointer :: foundelem,newelem
 integer(kind=fmm_integer) :: found,i

 rootelem%mylist(mypos) = rootelem%mylist(mypos)+myval

 end subroutine addelem

 function getelem(rootelem,mypos) result(myval)
 implicit none
 type(trootelem), intent(inout) :: rootelem
 integer(kind=fmm_integer), intent(in) :: mypos
 type(telem), pointer :: foundelem,newelem
 integer(kind=fmm_integer) :: myval,found

 myval = rootelem%mylist(mypos)

 end function getelem
end module pvlist

#endif
