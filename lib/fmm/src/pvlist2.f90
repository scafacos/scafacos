#include "fmm.h"

#ifdef FMM_PARALLEL

module pvlist
use fmmkinds
use fmmalloc
implicit none
private

type trootelem
 type(telem), pointer :: first,last,current
 integer(kind=fmm_integer) :: lastidx,currentidx,defaultval
#ifdef DEBUG
 integer(kind=fmm_integer), allocatable, dimension(:) :: mylist
#endif
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

#ifdef DEBUG
 call fmmallocate(rootelem%mylist,0,100,i)
 if(i.ne.0) call bummer('fmm: initlist error, allocation stat =',i)
 rootelem%mylist = myval
#endif

 rootelem%lastidx = 0
 rootelem%currentidx = 0
 rootelem%first => null()
 rootelem%last => null()
 rootelem%current => null()

 rootelem%defaultval = myval
 end subroutine initlist

 subroutine destroylist(rootelem)
 implicit none
 type(trootelem), intent(inout) :: rootelem
 integer(kind=fmm_integer) :: loopidx,i

 select case (rootelem%lastidx)
   case (1:) ! one or more elems in list
     if(associated(rootelem%last%next)) then
       call bummer('fmm: destroylist error, lastidx =',rootelem%lastidx)
     endif
     do loopidx = rootelem%lastidx,1,-1
       if(associated(rootelem%last)) then
         rootelem%current => rootelem%last%prev       
         call fmmdeallocate(rootelem%last,i)
         if(i.ne.0) then
           call bummer('fmm: destroylist error, deallocation stat =',i)
         endif
         rootelem%last => rootelem%current
       else
          call bummer('fmm: destroylist error, idx =',loopidx)
       endif  
     end do

     ! check if list is empty
     if(associated(rootelem%last)) call bummer('fmm: destroylist error, not empty =',1)

     ! reset counter and pointer
     rootelem%currentidx = 0
     rootelem%lastidx = 0
     rootelem%first => null()
   case (:-1) ! invalid number of elems
     call bummer('fmm: destroylist error, lastidx =',rootelem%lastidx)
 end select

#ifdef DEBUG
 if(allocated(rootelem%mylist)) then
   fmmdeallocate(rootelem%mylist,i)
   if(i.ne.0) call bummer('fmm: destroylist error, deallocation stat =',i)
 endif
#endif
 end subroutine destroylist 

 subroutine setlist(rootelem,myval)
 implicit none
 type(trootelem), intent(inout) :: rootelem
 integer(kind=fmm_integer), intent(in) :: myval
 integer(kind=fmm_integer) :: loopidx

 select case (rootelem%lastidx)
   case (1:) ! one or more elems in list
     if(associated(rootelem%last%next)) then
       call bummer('fmm: setlist error, lastidx =',rootelem%lastidx)
     endif
     rootelem%current => rootelem%last
     do loopidx = rootelem%lastidx,1,-1
       if(associated(rootelem%current)) then
         rootelem%current%val = myval
         rootelem%current => rootelem%current%prev       
       else
          call bummer('fmm: setlist error, idx =',loopidx)
       endif  
     end do

     ! check if list is empty
     if(associated(rootelem%current)) call bummer('fmm: setlist error, not empty =',1)

     ! reset counter
     rootelem%currentidx = 1
     rootelem%current => rootelem%first
   case (:-1) ! invalid number of elems
     call bummer('fmm: setlist error, lastidx =',rootelem%lastidx)
 end select

#ifdef DEBUG
 if(allocated(rootelem%mylist)) then
   mylist = myval
 endif
#endif
 end subroutine setlist 
 
 subroutine setelem(rootelem,mypos,myval)
 implicit none
 type(trootelem), intent(inout) :: rootelem
 integer(kind=fmm_integer), intent(in) :: mypos,myval
 type(telem), pointer :: foundelem,newelem
 integer(kind=fmm_integer) :: found,i

#ifdef DEBUG
 rootelem%mylist(mypos) = myval
#endif
 
 select case (rootelem%lastidx)
   case (1:)
   
!    rootelem%currentidx = 1
!    rootelem%current => rootelem%first
      
    if(rootelem%current%pos.eq.mypos) then
      rootelem%current%val = myval
      return
    else
      foundelem => dummyelem
      found = findelem(rootelem,mypos,foundelem)
    endif
   case (0)
     foundelem => dummyelem
     foundelem%prev => null()
     foundelem%next => null()
     found = 0
   case (:-1)
     call bummer('fmm: addelem error, lastidx',rootelem%lastidx)
 end select
 select case (found)
   case (:0) ! did not find mypos in list
     call fmmallocate(newelem,i)
     if(i.ne.0) call bummer('fmm: setelem allocation error, stat =',i)
     newelem%pos = mypos
     newelem%val = myval
     newelem%prev => foundelem%prev
     newelem%next => foundelem%next

     if(associated(foundelem%prev)) then 
       foundelem%prev%next => newelem
     else
       rootelem%first => newelem
     endif
     if(associated(foundelem%next)) then
       foundelem%next%prev => newelem
     else
       rootelem%last => newelem
     endif

     rootelem%lastidx = rootelem%lastidx+1
     rootelem%current => newelem
     rootelem%currentidx = abs(found)+1
   case (1:) ! found mypos in list
     foundelem%val = myval

     rootelem%current => foundelem
     rootelem%currentidx = found
 end select
 end subroutine setelem

 subroutine addelem(rootelem,mypos,myval)
 implicit none
 type(trootelem), intent(inout) :: rootelem
 integer(kind=fmm_integer), intent(in) :: mypos,myval
 type(telem), pointer :: foundelem,newelem
 integer(kind=fmm_integer) :: found,i

#ifdef DEBUG
 rootelem%mylist(mypos) = rootelem%mylist(mypos)+myval
#endif

 select case (rootelem%lastidx)
   case (1:)
   
!    rootelem%currentidx = 1
!    rootelem%current => rootelem%first
    
    if(rootelem%current%pos.eq.mypos) then
      rootelem%current%val = rootelem%current%val + myval
#ifdef DEBUG
      if(rootelem%current%val.ne.rootelem%mylist(mypos)) call bummer('fmm: addelem error',1)      
#endif      
      return
    else
      foundelem => dummyelem
      found = findelem(rootelem,mypos,foundelem)
    endif
   case (0)
     foundelem => dummyelem
     foundelem%prev => null()
     foundelem%next => null()
     found = 0
   case (:-1)
     call bummer('fmm: addelem error, lastidx',rootelem%lastidx)
 end select
 select case (found)
   case (:0) ! did not found mypos in list
     call fmmallocate(newelem,i)
     if(i.ne.0) call bummer('fmm: addelem allocation error, stat =',i)
     newelem%pos = mypos
     newelem%val = myval
     newelem%prev => foundelem%prev
     newelem%next => foundelem%next

     if(associated(foundelem%prev)) then 
       foundelem%prev%next => newelem
     else
       rootelem%first => newelem
     endif
     if(associated(foundelem%next)) then
       foundelem%next%prev => newelem
     else
       rootelem%last => newelem
     endif

     rootelem%lastidx = rootelem%lastidx+1
     rootelem%current => newelem
     rootelem%currentidx = abs(found)+1
   case (1:) ! found mypos in list
     foundelem%val = foundelem%val + myval
     rootelem%current => foundelem
     rootelem%currentidx = found
 end select
 end subroutine addelem

 function getelem(rootelem,mypos) result(myval)
 implicit none
 type(trootelem), intent(inout) :: rootelem
 integer(kind=fmm_integer), intent(in) :: mypos
 type(telem), pointer :: foundelem,newelem
 integer(kind=fmm_integer) :: myval,found

 select case (rootelem%lastidx)
   case (1:)
     foundelem => dummyelem
     found = findelem(rootelem,mypos,foundelem)
     myval = foundelem%val 
   case (0)
     myval = rootelem%defaultval
   case (:-1)  
     call bummer('fmm: getelem error, lastidx',rootelem%lastidx)
 end select 
    
#ifdef DEBUG
 if(myval.ne.rootelem%mylist(mypos)) call bummer('fmm: getelem, myval-rootelem%mylist(mypos) =',myval-rootelem%mylist(mypos))
#endif
 end function getelem

 function findelem(rootelem,mypos,foundelem) result(found)
 implicit none
 type(trootelem), intent(inout) :: rootelem
 integer(kind=fmm_integer), intent(in) :: mypos
 type(telem), pointer, intent(out) :: foundelem ! >0 found (currentidx) <0 not found currentidx
 type(telem), pointer :: thiselem
 integer(kind=fmm_integer) :: loopidx,found

 select case (rootelem%lastidx)
   case (1:) ! one or more elements in list
! if enabled, search starts from first elem
!     rootelem%current => rootelem%first
!     rootelem%currentidx = 1

     thiselem => rootelem%current
     select case (mypos-thiselem%pos)
       case (0) ! elem found => return ptr to thiselem
         foundelem => thiselem
         found = rootelem%currentidx
         return
       case (1:) ! go to the right in list
         select case (rootelem%lastidx-rootelem%currentidx)
          case (0) ! current elem is last elem => mypos would create new last elem
            foundelem%val = rootelem%defaultval
	     foundelem%prev => rootelem%last
	     foundelem%next => null()
            found = -rootelem%lastidx
	     return
	   case (1:) ! current elem is not last elem
	     thiselem => thiselem%next
             do loopidx = rootelem%currentidx+1,rootelem%lastidx
	       select case (mypos-thiselem%pos)
		  case(:-1) ! elem not found => mypos would create elem left to thiselem
                  foundelem%val = rootelem%defaultval
                  foundelem%prev => thiselem%prev
                  foundelem%next => thiselem
                  found = -loopidx+1
		    return
		  case(0) ! elem found => return ptr to thiselem
                 foundelem => thiselem
                 found = loopidx
		   return
	         case(1:)
		   thiselem => thiselem%next
	       end select
	     end do
            ! current elem is beyond last elem => mypos would create new last elem
            foundelem%val = rootelem%defaultval
            foundelem%prev => rootelem%last
            foundelem%next => null()
            found = -rootelem%lastidx
             return   
	   case default ! current elem is beyond last elem
             call bummer('fmm: findelem error, lastidx =',rootelem%lastidx)  
	 end select
       case(:-1) ! go to the left in list
         select case (rootelem%currentidx-1)
	   case (0) ! current elem is first elem => mypos would create new first elem
            foundelem%val = rootelem%defaultval
            foundelem%prev => null()
            foundelem%next => rootelem%first
            found = 0
	     return
	   case (1:) ! current elem is not first elem
	     thiselem => thiselem%prev
             do loopidx = rootelem%currentidx-1,1,-1
	       select case (mypos-thiselem%pos)
	         case(:-1)
		   thiselem => thiselem%prev
		 case(0) ! elem found => return ptr to thiselem
                 foundelem => thiselem
                 found = loopidx
		   return
		 case(1:) ! elem not found => mypos would create elem right to thiselem
                 foundelem%val = rootelem%defaultval
                 foundelem%prev => thiselem
                 foundelem%next => thiselem%next
                 found = -loopidx
		   return
	       end select
	     end do
            ! current elem is ahead of first elem => mypos would create new first elem
            foundelem%val = rootelem%defaultval
            foundelem%prev => null()
            foundelem%next => rootelem%first
            found = 0
            return   
	   case default
             call bummer('fmm: findelem error, lastidx =',rootelem%lastidx)  
	 end select       
     end select
   case (0) ! empty list => mypos would create new first and last elem
     foundelem%val = rootelem%defaultval
     foundelem%prev => null()
     foundelem%next => null()
     found = 0
     return
   case default ! error condition <0 elems in list
     call bummer('fmm: findelem error, lastidx =',rootelem%lastidx)
 end select
 end function findelem
end module pvlist

#endif
