!
! Copyright (C) 2012 Ivo Kabadshow, Holger Dachsel
!
!  This file is part of ScaFaCoS.
!
!  ScaFaCoS is free software: you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.
!
!  ScaFaCoS is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with this program.  If not, see <http://www.gnu.org/licenses/>.
!

#include <fconfig.h>
!#include "../../../../fmm.h"

#ifndef FMM_ALLOCALIGNED 
#  error FMM_ALLOCALIGNED is not defined, simd intrinsics cannot be used.
#endif

      subroutine pass2trfrqdcach(nmultipoles,&
                                 nsqmultipoles,&
				 jaddress,&
                                 jacc,&
                                 romega1,&
				 iomega1,&
				 romega2,&
				 iomega2,&
				 romega1a,&
				 iomega1a,&
				 romega2a,&
				 iomega2a,&
				 romega1b,&
				 iomega1b,&
				 romega2b,&
				 iomega2b,&
				 romega1c,&
				 iomega1c,&
				 romega2c,&
				 iomega2c,&
				 mu1,&
				 mu2,&
				 mu3,&
				 mu4,&
				 mu5,&
				 mu6,&
				 mu7,&
				 mu8,&
				 mu9,&
				 mu10,&
				 mu11,&
                                 mu12,&
				 mu13,&
				 mu14,&
				 mu15,&
				 mu16,&
				 cmphi,&
				 smphi,&
				 cmphipi,&
				 smphipi,&
				 csmphi,&
                                 csmphipi,&
				 sg,&
				 fr,&
				 d2,&
				 d3f,&
				 scr1,&
				 scr2)

      use fmmkinds

      implicit none

      real(kind=fmm_real) :: romega1(*),iomega1(*)
      real(kind=fmm_real) :: romega2(*),iomega2(*)
      real(kind=fmm_real) :: romega1a(*),iomega1a(*)
      real(kind=fmm_real) :: romega2a(*),iomega2a(*)
      real(kind=fmm_real) :: romega1b(*),iomega1b(*)
      real(kind=fmm_real) :: romega2b(*),iomega2b(*)
      real(kind=fmm_real) :: romega1c(*),iomega1c(*)
      real(kind=fmm_real) :: romega2c(*),iomega2c(*)
      
      real(kind=fmm_real) :: mu1(*),mu2(*),mu3(*),mu4(*)
      real(kind=fmm_real) :: mu5(*),mu6(*),mu7(*),mu8(*)
      real(kind=fmm_real) :: mu9(*),mu10(*),mu11(*),mu12(*)
      real(kind=fmm_real) :: mu13(*),mu14(*),mu15(*),mu16(*)
      
      real(kind=fmm_real) :: cmphi(0:*),smphi(0:*)
      real(kind=fmm_real) :: cmphipi(0:*),smphipi(0:*)
      real(kind=fmm_real) :: csmphi(*),csmphipi(*)
      
      real(kind=fmm_real) :: sg(0:*),fr(0:*)
      real(kind=fmm_real) :: d2(*),d3f(*)
      real(kind=fmm_real) :: scr1(*),scr2(*)

      integer(kind=fmm_integer) nmultipoles,nsqmultipoles,jaddress(*)
      integer(kind=4) :: i,j,k,l,m,mmmm,mmm,n,nn,mm,mmmmm,mmmmmm,nnn

      logical(kind=fmm_logical) jacc(*)

      real(kind=fmm_real), parameter :: zero = 0.e0_fmm_real
      real(kind=fmm_real), parameter :: one = 1.e0_fmm_real

      complex(kind=8) :: reg00,reg01,reg02,reg03,reg04,reg05,reg06,reg07
      complex(kind=8) :: reg08,reg09,reg10,reg11,reg12,reg13,reg14,reg15

      complex(kind=8) :: reg16,reg17

      interface
       subroutine rotate_around_z(nmultipoles,&
                                  romega11,romega12,romega13,romega14,&
                                  iomega11,iomega12,iomega13,iomega14,&
                                  romega21,romega22,romega23,romega24,&
                                  iomega21,iomega22,iomega23,iomega24,&
                                  csmphi,csmphipi,scr1,scr2,d2) &
				  bind(c,Name='rotate_around_z')
         use iso_c_binding, only : c_double,c_long_long
         implicit none
         real(c_double), dimension(*) :: romega11,romega12,romega13,romega14
         real(c_double), dimension(*) :: iomega11,iomega12,iomega13,iomega14
         real(c_double), dimension(*) :: romega21,romega22,romega23,romega24
         real(c_double), dimension(*) :: iomega21,iomega22,iomega23,iomega24
         real(c_double), dimension(*) :: csmphi,csmphipi
         real(c_double), dimension(*) :: scr1,scr2,d2
         integer(c_long_long), value :: nmultipoles
       end subroutine rotate_around_z

       subroutine rotate_around_y(nmultipoles,scr1,scr2,d2) bind(c,Name='rotate_around_y')
         use iso_c_binding, only : c_double,c_long_long
         implicit none
         real(c_double), dimension(*) :: scr1,scr2
	 real(c_double), dimension(*) :: d2
         integer(c_long_long), value :: nmultipoles
       end subroutine rotate_around_y

       subroutine m2l_along_z(nmultipoles,scr1,scr2,d2,fr,sg) bind(c,Name='m2l_along_z')
         use iso_c_binding, only : c_double,c_long_long
         implicit none
         real(c_double), dimension(*) :: scr1,scr2
	 real(c_double), dimension(*) :: d2,fr,sg
         integer(c_long_long), value :: nmultipoles
       end subroutine m2l_along_z

       subroutine rotate_back_around_yz(nmultipoles,csmphi,csmphipi,scr1,scr2,d3f,jaddress) bind(c,Name='rotate_back_around_yz')
         use iso_c_binding, only : c_double,c_long_long
         implicit none
         real(c_double), dimension(*) :: scr1,scr2
	 real(c_double), dimension(*) :: d3f,csmphi,csmphipi
         integer(c_long_long), dimension(*) :: jaddress
         integer(c_long_long), value :: nmultipoles
       end subroutine rotate_back_around_yz
      end interface

!     rotate around z
      call rotate_around_z(nmultipoles,&
                           romega1,romega1a,romega1b,romega1c,&
                           iomega1,iomega1a,iomega1b,iomega1c,&
                           romega2,romega2a,romega2b,romega2c,&
                           iomega2,iomega2a,iomega2b,iomega2c,&
                           csmphi,csmphipi,scr1,scr2,d2)

!     rotate around y
      call rotate_around_y(nmultipoles,scr1,scr2,d2)

!     perform shift
      call m2l_along_z(nmultipoles,scr1,scr2,d2,fr,sg)

!     rotate back around yz
      call rotate_back_around_yz(nmultipoles,csmphi,csmphipi,scr1,scr2,d3f,jaddress)

!     store translations
      i = 0

      if(jacc(1)) then
         do j = 1,nsqmultipoles
            i = i+1
            mu1(j) = mu1(j)+scr2(i)
         enddo
      else
         j = i+nsqmultipoles
         k = j+nsqmultipoles
         j = j+1
         do l = j,k
            i = i+1
            scr2(l) = scr2(l)+scr2(i)
         enddo
      endif

      if(jacc(2)) then
         do j = 1,nsqmultipoles
            i = i+1
            mu2(j) = mu2(j)+scr2(i)
         enddo
      else
         j = i+nsqmultipoles
         k = j+nsqmultipoles
         j = j+1
         do l = j,k
            i = i+1
            scr2(l) = scr2(l)+scr2(i)
         enddo
      endif

      if(jacc(3)) then
         do j = 1,nsqmultipoles
            i = i+1
            mu3(j) = mu3(j)+scr2(i)
         enddo
      else
         j = i+nsqmultipoles
         k = j+nsqmultipoles
         j = j+1
         do l = j,k
            i = i+1
            scr2(l) = scr2(l)+scr2(i)
         enddo
      endif

      if(jacc(4)) then
         do j = 1,nsqmultipoles
            i = i+1
            mu4(j) = mu4(j)+scr2(i)
         enddo
      else
         j = i+nsqmultipoles
         k = j+nsqmultipoles
         j = j+1
         do l = j,k
            i = i+1
            scr2(l) = scr2(l)+scr2(i)
         enddo
      endif

      if(jacc(5)) then
         do j = 1,nsqmultipoles
            i = i+1
            mu5(j) = mu5(j)+scr2(i)
         enddo
      else
         j = i+nsqmultipoles
         k = j+nsqmultipoles
         j = j+1
         do l = j,k
            i = i+1
            scr2(l) = scr2(l)+scr2(i)
         enddo
      endif

      if(jacc(6)) then
         do j = 1,nsqmultipoles
            i = i+1
            mu6(j) = mu6(j)+scr2(i)
         enddo
      else
         j = i+nsqmultipoles
         k = j+nsqmultipoles
         j = j+1
         do l = j,k
            i = i+1
            scr2(l) = scr2(l)+scr2(i)
         enddo
      endif

      if(jacc(7)) then
         do j = 1,nsqmultipoles
            i = i+1
            mu7(j) = mu7(j)+scr2(i)
         enddo
      else
         j = i+nsqmultipoles
         k = j+nsqmultipoles
         j = j+1
         do l = j,k
            i = i+1
            scr2(l) = scr2(l)+scr2(i)
         enddo
      endif

      if(jacc(8)) then
         do j = 1,nsqmultipoles
            i = i+1
            mu8(j) = mu8(j)+scr2(i)
         enddo
      else
         j = i+nsqmultipoles
         k = j+nsqmultipoles
         j = j+1
         do l = j,k
            i = i+1
            scr2(l) = scr2(l)+scr2(i)
         enddo
      endif
      if(jacc(9)) then
         do j = 1,nsqmultipoles
            i = i+1
            mu9(j) = mu9(j)+scr2(i)
         enddo
      else
         j = i+nsqmultipoles
         k = j+nsqmultipoles
         j = j+1
         do l = j,k
            i = i+1
            scr2(l) = scr2(l)+scr2(i)
         enddo
      endif

      if(jacc(10)) then
         do j = 1,nsqmultipoles
            i = i+1
            mu10(j) = mu10(j)+scr2(i)
      enddo
      else
         j = i+nsqmultipoles
         k = j+nsqmultipoles
         j = j+1
         do l = j,k
            i = i+1
            scr2(l) = scr2(l)+scr2(i)
         enddo
      endif

      if(jacc(11)) then
         do j = 1,nsqmultipoles
            i = i+1
            mu11(j) = mu11(j)+scr2(i)
         enddo
      else
         j = i+nsqmultipoles
         k = j+nsqmultipoles
         j = j+1
         do l = j,k
            i = i+1
            scr2(l) = scr2(l)+scr2(i)
         enddo
      endif

      if(jacc(12)) then
         do j = 1,nsqmultipoles
            i = i+1
            mu12(j) = mu12(j)+scr2(i)
         enddo
      else
         j = i+nsqmultipoles
         k = j+nsqmultipoles
         j = j+1
         do l = j,k
            i = i+1
            scr2(l) = scr2(l)+scr2(i)
         enddo
      endif

      if(jacc(13)) then
         do j = 1,nsqmultipoles
            i = i+1
            mu13(j) = mu13(j)+scr2(i)
         enddo
      else
         j = i+nsqmultipoles
         k = j+nsqmultipoles
         j = j+1
         do l = j,k
            i = i+1
            scr2(l) = scr2(l)+scr2(i)
         enddo
      endif

      if(jacc(14)) then
         do j = 1,nsqmultipoles
            i = i+1
            mu14(j) = mu14(j)+scr2(i)
         enddo
      else
         j = i+nsqmultipoles
         k = j+nsqmultipoles
         j = j+1
         do l = j,k
            i = i+1
            scr2(l) = scr2(l)+scr2(i)
         enddo
      endif

      if(jacc(15)) then
         do j = 1,nsqmultipoles
            i = i+1
            mu15(j) = mu15(j)+scr2(i)
         enddo
      else
         j = i+nsqmultipoles
         k = j+nsqmultipoles
         j = j+1
         do l = j,k
            i = i+1
            scr2(l) = scr2(l)+scr2(i)
         enddo
      endif

      if(jacc(16)) then
         do j = 1,nsqmultipoles
            i = i+1
            mu16(j) = mu16(j)+scr2(i)
         enddo
      else
         j = i+nsqmultipoles
         k = j+nsqmultipoles
         j = j+1
         do l = j,k
            i = i+1
            scr2(l) = scr2(l)+scr2(i)
         enddo
      endif
      end subroutine pass2trfrqdcach
