      subroutine getneighbors(int3x,int3y,int3z, &
				 & vedge,casejump,boxaddr, &
				 & px,py,pz,ix,iy,iz, &
				 & child,fwdstart,vcboxaddr,dist,fwdonly,&
				 & fmmpass,int3exist)
      !version 2: including int3x, int3y, int3z, vedge independent box generation with additinonal parameter int3exist      
      ! getnbrs -> getneighbors
      use getneighbors_vars, only:fvedge,fmm_integer
      implicit none
      integer(kind=fmm_integer) int3x(0:*),int3y(0:*),int3z(0:*)
      integer(kind=fmm_integer) :: child !position of child 1-8
      integer(kind=fmm_integer) :: boxaddr !box adresse centerbox
      integer(kind=fmm_integer) :: ix,iy,iz !boxaddr 1D index
      integer(kind=fmm_integer) :: px,py,pz !parent box 1D index
      integer(kind=fmm_integer) :: fwdstart !first boxnumber larger than boxaddr
      integer(kind=fmm_integer) :: sboxaddr !short addr max 8,64,512,4096 
      integer(kind=fmm_integer) :: eix,eiy,eiz !edge ix,iy,iz temps
      integer(kind=fmm_integer) :: fix,fiy,fiz !child box temps
      integer(kind=fmm_integer),dimension(27) :: vcboxaddr ! vector of neighborbox adresses
      integer(kind=fmm_integer) :: dist !1D distance from neighborbox to centerbox: entry getdist, input:
      integer(kind=fmm_integer) :: fwdonly ! full (27) neighbor calculation or righthandside only
      integer(kind=fmm_integer) :: fmmpass ! pass2 fmmpass =1, pass5 fmmpass=2
      integer(kind=fmm_integer) :: int3exist ! 1= int3x,y,z available and allocated, 2= int3x,y,z not available and not allocated
      integer(kind=fmm_integer),dimension(0:*) :: casejump !(0:8) (0:64) (0:512) (0:4096) entry array to find correct case select
      integer(kind=fmm_integer),dimension(0:*) :: vedge !corresponding edge
      select case(fmmpass) !holds bitmask 7,63,511,4095
        case(1)
        sboxaddr=iand(ishft(boxaddr,-3),casejump(0))+1 ! boxaddr-1 done in main fmm call used in casejump
        case(2)
        sboxaddr=iand(boxaddr,casejump(0))+1 ! pass 5 neighbors! 
      end select
      select case(child)
        case(1)
          fix = ix
          fiy = iy
          fiz = iz
        case(2)
          fix = ix-1
          fiy = iy
          fiz = iz
        case(3)
          fix = ix
          fiy = iy-1
          fiz = iz
        case(4)
          fix = ix-1
          fiy = iy-1
          fiz = iz
        case(5)
          fix = ix
          fiy = iy
          fiz = iz-1
        case(6)
          fix = ix-1
          fiy = iy
          fiz = iz-1
        case(7)
          fix = ix
          fiy = iy-1
          fiz = iz-1
        case(8)
          fix = ix-1
          fiy = iy-1
          fiz = iz-1
      end select
      select case(casejump(sboxaddr))
      case(  1)
          call caladdr4(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
          return
      case(  2)
          call caladdr6(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
          return
      case(  3)
          call caladdr7(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
          return
      case(  4)
          call caladdr8(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
          return
      case(  5)
          call caladdr9(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
          return
      case(  6)
          call caladdr10(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
          return
      case(  7)
          call caladdr11(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
          return
      case(  8)
          call caladdr12(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
          return
      case(  9)
          call caladdr13(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
          return
      case( 10)
          call caladdr14(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
          return
      case( 11)
          call caladdr2(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
          return
      case( 12)
          call caladdr16(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
          return
      case( 13)
          call caladdr17(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
          return
      case( 14)
          call caladdr18(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
          return
      case( 15)
          call caladdr19(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
          return
      case( 16)
          call caladdr1(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
          return
      case( 17)
          call caladdr20(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
          return
      case( 18)
          call caladdr21(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
          return
      case( 19)
          call caladdr22(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
          return
      case( 20)
          call caladdr23(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
          return
      case( 21)
          call caladdr24(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
          return
      case( 22)
          call caladdr26(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
          return
      case( 23)
          call caladdr27(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
          return
      case( 24)
          call caladdr28(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
          return
      case( 25)
          call caladdr29(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
          return
      case( 26)
          call caladdr30(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
          return
      case( 27)
          call caladdr31(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
          return
      case( 28)
          call caladdr32(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
          return
      case( 29)
          call caladdr33(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
          return
      case( 30)
          call caladdr34(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
          return
      case( 31)
          call caladdr35(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
          return
      case( 32)
          call caladdr36(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
          return
      case( 33)
          call caladdr37(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
          return
      case( 34)
          call caladdr38(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
          return
      case( 35)
          call caladdr3(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
          return
      case( 36)
          call caladdr5(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
          return
      case( 37)
          call caladdr15(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
          return
      case( 38)
          call caladdr40(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
          return
      case( 39)
          call caladdr41(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
          return
      case( 40)
          call caladdr25(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
          return
      case( 41)
          call caladdr42(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
          return
      case( 42)
          call caladdr44(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
          return
      case( 43)
          call caladdr45(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
          return
      case( 44)
          call caladdr46(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
          return
      case( 45)
          call caladdr47(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
          return
      case( 46)
          call caladdr48(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
          return
      case( 47)
          call caladdr39(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
          return
      case( 48)
          call caladdr43(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
          return
      case( 49)
! possible calls  2  2 26  2  2  2 26  2
! coordinates  1  0  0
        select case(int3exist)
	case(1)
          eiy = vedge(py)
          eiz = vedge(pz)
	case(2)
          eiy = fvedge(py)
          eiz = fvedge(pz)
        end select
        if(eiy.le.eiz) then
!        if((vedge(py+ 0).le.vedge(pz+ 0))) then
          call caladdr2(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
          return
        else
          call caladdr26(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
          return
        endif
      case( 50)
! possible calls  3  3  3  3 27  3 27  3
! coordinates  0  1  0
        select case(int3exist)
	case(1)
          eix = vedge(px)
          eiz = vedge(pz)
	case(2)
          eix = fvedge(px)
          eiz = fvedge(pz)
        end select
        if(eix.le.eiz) then	
!        if((vedge(px+ 0).le.vedge(pz+ 0))) then
          call caladdr3(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
          return
        else
          call caladdr27(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
          return
        endif
      case( 51)
! possible calls  5  5  5  5 47 47  5  5
! coordinates  0  0  1
        select case(int3exist)
	case(1)
          eix = vedge(px)
          eiy = vedge(py)
	case(2)
          eix = fvedge(px)
          eiy = fvedge(py)
        end select	
        if(eix.le.eiy) then
!        if((vedge(px+ 0).le.vedge(py+ 0))) then
          call caladdr5(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
          return
        else
          call caladdr47(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
          return
        endif
      case( 52)
! possible calls  1  1 28  1  1  1 28  1
! coordinates  2  0  0
        select case(int3exist)
	case(1)
          eiy = vedge(py)
          eiz = vedge(pz)
	case(2)
          eiy = fvedge(py)
          eiz = fvedge(pz)
	end select	
        if(eiy.le.eiz) then
!        if((vedge(py+ 0).le.vedge(pz+ 0))) then
          call caladdr1(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
          return
        else
          call caladdr28(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
          return
        endif
      case( 53)
! possible calls 15 15 15 15 32 15 32 15
! coordinates  0  2  0
        select case(int3exist)
	case(1)
          eix = vedge(px)
          eiz = vedge(pz)
	case(2)
          eix = fvedge(px)
          eiz = fvedge(pz)
        end select
        if(eix.le.eiz) then
!        if((vedge(px+ 0).le.vedge(pz+ 0))) then
          call caladdr15(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
          return
        else
          call caladdr32(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
          return
        endif
      case( 54)
! possible calls 25 25 25 25 48 48 25 25
! coordinates  0  0  2
        select case(int3exist)
	case(1)
          eix = vedge(px)
	  eiy = vedge(py)
	case(2)
          eix = fvedge(px)
	  eiy = fvedge(py)
        end select
        if(eix.le.eiy) then
!        if((vedge(px+ 0).le.vedge(py+ 0))) then
          call caladdr25(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
          return
        else
          call caladdr48(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
          return
        endif
      case( 55)
! possible calls 30 10 30 10 10 10 10 10
! coordinates  7  1  0
        select case(int3exist)
	case(1)
          eix = vedge(px+1)
          eiz = vedge(pz)
	case(2)
          eix = fvedge(px+1)
          eiz = fvedge(pz)
        end select
        if(eix.le.eiz) then
!        if((vedge(px+ 1).le.vedge(pz+ 0))) then
          call caladdr10(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
          return
        else
          call caladdr30(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
          return
        endif
      case( 56)
! possible calls 39 39 12 12 12 12 12 12
! coordinates  7  0  1
        select case(int3exist)
	case(1)
          eix = vedge(px+1)
	  eiy = vedge(py)
	case(2)
          eix = fvedge(px+1)
	  eiy = fvedge(py)	  
        end select
        if(eix.le.eiy) then
!        if((vedge(px+ 1).le.vedge(py+ 0))) then
          call caladdr12(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
          return
        else
          call caladdr39(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
          return
        endif
      case( 57)
! possible calls 35 20 35 20 20 20 20 20
! coordinates  7  2  0
        select case(int3exist)
	case(1)
          eix = vedge(px+1)
          eiz = vedge(pz)
	case(2)
          eix = fvedge(px+1)
          eiz = fvedge(pz)
        end select
        if(eix.le.eiz) then
!        if((vedge(px+ 1).le.vedge(pz+ 0))) then
          call caladdr20(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
          return
        else
          call caladdr35(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
          return
        endif
      case( 58)
! possible calls 43 43 29 29 29 29 29 29
! coordinates  7  0  2
        select case(int3exist)
	case(1)
          eix = vedge(px+1)
	  eiy = vedge(py)
	case(2)
          eix = fvedge(px+1)
	  eiy = fvedge(py)
        end select
        if(eix.le.eiy) then
!        if((vedge(px+ 1).le.vedge(py+ 0))) then
          call caladdr29(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
          return
        else
          call caladdr43(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
          return
        endif
      case( 59)
! possible calls 34  4  4  4 34  4  4  4
! coordinates  1  7  0
        select case(int3exist)
	case(1)
          eiy = vedge(py+1)
          eiz = vedge(pz)
	case(2)
          eiy = fvedge(py+1)
          eiz = fvedge(pz)
        end select
        if(eiy.le.eiz) then
!        if((vedge(py+ 1).le.vedge(pz+ 0))) then
          call caladdr4(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
          return
        else
          call caladdr34(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
          return
        endif
      case( 60)
! possible calls 40 40 40 40 40 40 18 18
! coordinates  0  7  1
        select case(int3exist)
	case(1)
          eix = vedge(px)
	  eiy = vedge(py+1)
	case(2)
          eix = fvedge(px)
	  eiy = fvedge(py+1)
        end select
        if(eix.le.eiy) then
!        if((vedge(px+ 0).le.vedge(py+ 1))) then
          call caladdr40(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
          return
        else
          call caladdr18(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
          return
        endif
      case( 61)
! possible calls 36  9  9  9 36  9  9  9
! coordinates  2  7  0
        select case(int3exist)
	case(1)
          eiy = vedge(py+1)
          eiz = vedge(pz)
	case(2)
          eiy = fvedge(py+1)
          eiz = fvedge(pz)
        end select
        if(eiy.le.eiz) then
!        if((vedge(py+ 1).le.vedge(pz+ 0))) then
          call caladdr9(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
          return
        else
          call caladdr36(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
          return
        endif
      case( 62)
! possible calls 44 44 44 44 44 44 33 33
! coordinates  0  7  2
        select case(int3exist)
	case(1)
          eix = vedge(px)
          eiy = vedge(py+1)
  	case(2)
          eix = fvedge(px)
          eiy = fvedge(py+1)	  
        end select
        if(eix.le.eiy) then
!        if((vedge(px+ 0).le.vedge(py+ 1))) then
          call caladdr44(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
          return
        else
          call caladdr33(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
          return
        endif
      case( 63)
! possible calls 24 24 41 41 24 24 24 24
! coordinates  7  7  1
        select case(int3exist)
	case(1)
          eix = vedge(px+1)
          eiy = vedge(py+1)
	case(2)
          eix = fvedge(px+1)
          eiy = fvedge(py+1)	  
        end select
        if(eix.le.eiy) then
!        if((vedge(px+ 1).le.vedge(py+ 1))) then
          call caladdr24(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
          return
        else
          call caladdr41(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
          return
        endif
      case( 64)
! possible calls 37 37 46 46 37 37 37 37
! coordinates  7  7  2
        select case(int3exist)
	case(1)
          eix = vedge(px+1)
	  eiy = vedge(py+1)
	case(2)
          eix = fvedge(px+1)
	  eiy = fvedge(py+1)
        end select
        if(eix.le.eiy) then
!        if((vedge(px+ 1).le.vedge(py+ 1))) then
          call caladdr37(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
          return
        else
          call caladdr46(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
          return
        endif
      case( 65)
! possible calls 17 17 17  6 17 17 17  6
! coordinates  1  0  7
        select case(int3exist)
	case(1)
          eiy = vedge(py)
          eiz = vedge(pz+1)
	case(2)
          eiy = fvedge(py)
          eiz = fvedge(pz+1)
        end select
        if(eiy.le.eiz) then
!        if((vedge(py+ 0).le.vedge(pz+ 1))) then
          call caladdr17(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
          return
        else
          call caladdr6(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
          return
        endif
      case( 66)
! possible calls 42 42 42 42 42  7 42  7
! coordinates  0  1  7
        select case(int3exist)
	case(1)
          eix = vedge(px)
          eiz = vedge(pz+1)
	case(2)
          eix = fvedge(px)
          eiz = fvedge(pz+1)	  
        end select
        if(eix.le.eiz) then
!        if((vedge(px+ 0).le.vedge(pz+ 1))) then
          call caladdr42(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
          return
        else
          call caladdr7(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
          return
        endif
      case( 67)
! possible calls 21 21 21 11 21 21 21 11
! coordinates  2  0  7
        select case(int3exist)
	case(1)
          eiy = vedge(py)
          eiz = vedge(pz+1)
	case(2)
          eiy = fvedge(py)
          eiz = fvedge(pz+1)
        end select
        if(eiy.le.eiz) then
!        if((vedge(py+ 0).le.vedge(pz+ 1))) then
          call caladdr21(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
          return
        else
          call caladdr11(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
          return
        endif
      case( 68)
! possible calls 45 45 45 45 45 16 45 16
! coordinates  0  2  7
        select case(int3exist)
	case(1)
          eix = vedge(px)
          eiz = vedge(pz+1)
	case(2)
          eix = fvedge(px)
          eiz = fvedge(pz+1)
        end select
        if(eix.le.eiz) then
!        if((vedge(px+ 0).le.vedge(pz+ 1))) then
          call caladdr45(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
          return
        else
          call caladdr16(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
          return
        endif
      case( 69)
! possible calls 31 14 31 14 31 31 31 31
! coordinates  7  1  7
        select case(int3exist)
	case(1)
          eix = vedge(px+1)
          eiz = vedge(pz+1)
	case(2)
          eix = fvedge(px+1)
          eiz = fvedge(pz+1)
        end select
        if(eix.le.eiz) then
!        if((vedge(px+ 1).le.vedge(pz+ 1))) then
          call caladdr31(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
          return
        else
          call caladdr14(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
          return
        endif
      case( 70)
! possible calls 38 22 38 22 38 38 38 38
! coordinates  7  2  7
        select case(int3exist)
	case(1)
          eix = vedge(px+1)
          eiz = vedge(pz+1)
	case(2)
          eix = fvedge(px+1)
          eiz = fvedge(pz+1)
        end select
        if(eix.le.eiz) then
!        if((vedge(px+ 1).le.vedge(pz+ 1))) then
          call caladdr38(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
          return
        else
          call caladdr22(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
          return
        endif
      case( 71)
! possible calls  8 19  8  8  8 19  8  8
! coordinates  1  7  7
        select case(int3exist)
	case(1)
          eiy = vedge(py+1)
          eiz = vedge(pz+1)
	case(2)
          eiy = fvedge(py+1)
          eiz = fvedge(pz+1)
        end select
        if(eiy.le.eiz) then
!        if((vedge(py+ 1).le.vedge(pz+ 1))) then
          call caladdr8(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
          return
        else
          call caladdr19(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
          return
        endif
      case( 72)
! possible calls 13 23 13 13 13 23 13 13
! coordinates  2  7  7
        select case(int3exist)
	case(1)
          eiy = vedge(py+1)
          eiz = vedge(pz+1)
	case(2)
          eiy = fvedge(py+1)
          eiz = fvedge(pz+1)
        end select
        if(eiy.le.eiz) then
!        if((vedge(py+ 1).le.vedge(pz+ 1))) then
          call caladdr13(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
          return
        else
          call caladdr23(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
          return
        endif
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------	
      case( 73)
        select case(int3exist)
	case(1)
          eix=vedge(px)
          eiy=vedge(py)
          eiz=vedge(pz)
	case(2)
          eix=fvedge(px)
          eiy=fvedge(py)
          eiz=fvedge(pz)
        end select      
        if(eix.le.eiy) then
          if(eiy.le.eiz) then
            call caladdr1(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
          else
            if(eix.le.eiz) then
              call caladdr28(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
            else
              call caladdr25(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
            endif
          endif
        else
          if(eiy.le.eiz) then
            if(eix.le.eiz) then
              call caladdr15(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
            else
              call caladdr32(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
            endif
          else
            call caladdr48(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
          endif
        endif
        return
      case( 74)
        select case(int3exist)
	case(1)      
          eix=vedge(px+1)
          eiy=vedge(py)
          eiz=vedge(pz)
	case(2)      
          eix=fvedge(px+1)
          eiy=fvedge(py)
          eiz=fvedge(pz)
        end select      
        if(eix.le.eiy) then
          if(eiy.le.eiz) then
            call caladdr2(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
          else
            if(eix.le.eiz) then
              call caladdr26(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
            else
              call caladdr29(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
            endif
          endif
        else
          if(eiy.le.eiz) then
            if(eix.le.eiz) then
              call caladdr20(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
            else
              call caladdr35(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
            endif
          else
            call caladdr43(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
          endif
        endif
        return
      case( 75)
        select case(int3exist)
	case(1)
          eix=vedge(px)
          eiy=vedge(py+1)
          eiz=vedge(pz)
	case(2)
          eix=fvedge(px)
          eiy=fvedge(py+1)
          eiz=fvedge(pz)
        end select      
        if(eix.le.eiy) then
          if(eiy.le.eiz) then
            call caladdr9(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
          else
            if(eix.le.eiz) then
              call caladdr36(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
            else
              call caladdr44(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
            endif
          endif
        else
          if(eiy.le.eiz) then
            if(eix.le.eiz) then
              call caladdr3(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
            else
              call caladdr27(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
            endif
          else
            call caladdr33(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
          endif
        endif
        return
      case( 76)
        select case(int3exist)
	case(1)      
          eix=vedge(px+1)
          eiy=vedge(py+1)
          eiz=vedge(pz)
	case(2)      
          eix=fvedge(px+1)
          eiy=fvedge(py+1)
          eiz=fvedge(pz)
        end select      
        if(eix.le.eiy) then
          if(eiy.le.eiz) then
            call caladdr4(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
          else
            if(eix.le.eiz) then
              call caladdr34(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
            else
              call caladdr37(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
            endif
          endif
        else
          if(eiy.le.eiz) then
            if(eix.le.eiz) then
              call caladdr10(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
            else
              call caladdr30(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
            endif
          else
            call caladdr46(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
          endif
        endif
        return
      case( 77)
        select case(int3exist)
	case(1)      
          eix=vedge(px)
          eiy=vedge(py)
          eiz=vedge(pz+1)
	case(2)      
          eix=fvedge(px)
          eiy=fvedge(py)
          eiz=fvedge(pz+1)
        end select      
        if(eix.le.eiy) then
          if(eiy.le.eiz) then
            call caladdr21(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
          else
            if(eix.le.eiz) then
              call caladdr11(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
            else
              call caladdr5(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
            endif
          endif
        else
          if(eiy.le.eiz) then
            if(eix.le.eiz) then
              call caladdr45(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
            else
              call caladdr16(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
            endif
          else
            call caladdr47(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
          endif
        endif
        return
      case( 78)
        select case(int3exist)
	case(1)      
          eix=vedge(px+1)
          eiy=vedge(py)
          eiz=vedge(pz+1)
	case(2)      
          eix=fvedge(px+1)
          eiy=fvedge(py)
          eiz=fvedge(pz+1)
        end select      
        if(eix.le.eiy) then
          if(eiy.le.eiz) then
            call caladdr17(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
          else
            if(eix.le.eiz) then
              call caladdr6(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
            else
              call caladdr12(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
            endif
          endif
        else
          if(eiy.le.eiz) then
            if(eix.le.eiz) then
              call caladdr38(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
            else
              call caladdr22(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
            endif
          else
            call caladdr39(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
          endif
        endif
        return
      case( 79)
        select case(int3exist)
	case(1)      
          eix=vedge(px)
          eiy=vedge(py+1)
          eiz=vedge(pz+1)
	case(2)      
          eix=fvedge(px)
          eiy=fvedge(py+1)
          eiz=fvedge(pz+1)
        end select	  
        if(eix.le.eiy) then
          if(eiy.le.eiz) then
            call caladdr13(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
          else
            if(eix.le.eiz) then
              call caladdr23(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
            else
              call caladdr40(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
            endif
          endif
        else
          if(eiy.le.eiz) then
            if(eix.le.eiz) then
              call caladdr42(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
            else
              call caladdr7(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
            endif
          else
            call caladdr18(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
          endif
        endif
        return
      case( 80)
        select case(int3exist)
	case(1)      
          eix=vedge(px+1)
          eiy=vedge(py+1)
          eiz=vedge(pz+1)
	case(2)      
          eix=fvedge(px+1)
          eiy=fvedge(py+1)
          eiz=fvedge(pz+1)
        end select      
        if(eix.le.eiy) then
          if(eiy.le.eiz) then
            call caladdr8(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
          else
            if(eix.le.eiz) then
              call caladdr19(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
            else
              call caladdr24(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
            endif
          endif
        else
          if(eiy.le.eiz) then
            if(eix.le.eiz) then
              call caladdr31(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
            else
              call caladdr14(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
            endif
          else
            call caladdr41(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
          endif
        endif
        return
      end select
      stop 'error getneighbors'
      end subroutine getneighbors
