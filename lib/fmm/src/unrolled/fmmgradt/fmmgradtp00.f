#include "../../fmm.h"
#if FMM_MAXNMULTIPOLES < 0 || FMM_MAXNMULTIPOLES > 50 || !defined(FMM_UNROLLED)
      subroutine fmmgrd(nmultipoles,rooperator,iooperator,rmu,imu,
     .rmu00,rmu10,rmu11,imu11)
c
      use fmmkinds
c
      implicit none
c
      real(kind=fmm_real) rooperator(*),iooperator(*),rmu(*),imu(*),
     .rmu00,rmu10,
     .rmu11,imu11,zero,a,b
c
      integer(kind=fmm_integer) nmultipoles,ioo,joo,koo,moo,jmu,i,j,k
c
      rmu00 = rmu(1)
c
      if(nmultipoles.gt.2) then
       rmu10 = rmu(2)
       rmu00 = rmu00+rooperator(2)*rmu10
       zero = iooperator(2)
       rmu11 = rmu(3)
       imu11 = imu(3)
       a = rooperator(3)*rmu11-iooperator(3)*imu11
       rmu00 = rmu00+(a+a)
c
       rmu10 = rmu10+rooperator(2)*rmu(4)
       rmu11 = rmu11-rooperator(3)*rmu(4)
       imu11 = imu11+iooperator(3)*rmu(4)
       rmu00 = rmu00+rooperator(4)*rmu(4)
c
       rmu11 = rmu11+rooperator(2)*rmu(5)
       imu11 = imu11+rooperator(2)*imu(5)
       a = rooperator(3)*rmu(5)-iooperator(3)*imu(5)
       b = rooperator(5)*rmu(5)-iooperator(5)*imu(5)
       rmu11 = rmu11+(rooperator(3)*rmu(6)-iooperator(3)*imu(6))
       imu11 = imu11+(rooperator(3)*imu(6)+iooperator(3)*rmu(6))
       b = b+(rooperator(6)*rmu(6)-iooperator(6)*imu(6))
c
       rmu10 = rmu10+(a+a)
       rmu00 = rmu00+(b+b)
c
       ioo = 6
       joo = 3
       koo = 3
       moo = 3
       jmu = 4
       i = 0
c
       do 1 j = 3,nmultipoles
        ioo = ioo+1
        joo = joo+1
        moo = moo+2
        jmu = jmu+3
        i = i+1
        k = koo+2
c
        rmu10 = rmu10+rooperator(joo)*rmu(ioo)
        rmu11 = rmu11-rooperator(k)*rmu(ioo)
        imu11 = imu11+iooperator(k)*rmu(ioo)
        rmu00 = rmu00+rooperator(ioo)*rmu(ioo)
c
        a = zero
        b = zero
c
        do 2 k = 1,i
         ioo = ioo+1
         joo = joo+1
         koo = koo+1
         moo = moo+1
         jmu = jmu+1
c
         rmu11=rmu11+(rooperator(koo)*rmu(ioo)-iooperator(koo)*imu(ioo))
         imu11=imu11+(rooperator(koo)*imu(ioo)+iooperator(koo)*rmu(ioo))
         a = a+(rooperator(joo)*rmu(ioo)-iooperator(joo)*imu(ioo))
         rmu11=rmu11-(rooperator(moo)*rmu(jmu)-iooperator(moo)*imu(jmu))
         imu11=imu11+(rooperator(moo)*imu(jmu)+iooperator(moo)*rmu(jmu))
         b = b+(rooperator(ioo)*rmu(ioo)-iooperator(ioo)*imu(ioo))
 2      continue
c
        ioo = ioo+1
        joo = joo+1
        koo = koo+1
c
        rmu11=rmu11+(rooperator(koo)*rmu(ioo)-iooperator(koo)*imu(ioo))
        imu11=imu11+(rooperator(koo)*imu(ioo)+iooperator(koo)*rmu(ioo))
        a = a+(rooperator(joo)*rmu(ioo)-iooperator(joo)*imu(ioo))
        b = b+(rooperator(ioo)*rmu(ioo)-iooperator(ioo)*imu(ioo))
c
        ioo = ioo+1
        koo = koo+1
c
        rmu11=rmu11+(rooperator(koo)*rmu(ioo)-iooperator(koo)*imu(ioo))
        imu11=imu11+(rooperator(koo)*imu(ioo)+iooperator(koo)*rmu(ioo))
        b = b+(rooperator(ioo)*rmu(ioo)-iooperator(ioo)*imu(ioo))
c
        rmu10 = rmu10+(a+a)
        rmu00 = rmu00+(b+b)
 1     continue
      elseif(nmultipoles.eq.2) then
       rmu10 = rmu(2)
       rmu00 = rmu00+rooperator(2)*rmu10
       rmu11 = rmu(3)
       imu11 = imu(3)
       a = rooperator(3)*rmu11-iooperator(3)*imu11
       rmu00 = rmu00+(a+a)
c
       rmu10 = rmu10+rooperator(2)*rmu(4)
       rmu11 = rmu11-rooperator(3)*rmu(4)
       imu11 = imu11+iooperator(3)*rmu(4)
       rmu00 = rmu00+rooperator(4)*rmu(4)
c
       rmu11 = rmu11+rooperator(2)*rmu(5)
       imu11 = imu11+rooperator(2)*imu(5)
       a = rooperator(3)*rmu(5)-iooperator(3)*imu(5)
       b = rooperator(5)*rmu(5)-iooperator(5)*imu(5)
       rmu11 = rmu11+(rooperator(3)*rmu(6)-iooperator(3)*imu(6))
       imu11 = imu11+(rooperator(3)*imu(6)+iooperator(3)*rmu(6))
       b = b+(rooperator(6)*rmu(6)-iooperator(6)*imu(6))
c
       rmu10 = rmu10+(a+a)
       rmu00 = rmu00+(b+b)
      elseif(nmultipoles.eq.1) then
       rmu10 = rmu(2)
       rmu00 = rmu00+rooperator(2)*rmu10
       rmu11 = rmu(3)
       imu11 = imu(3)
       a = rooperator(3)*rmu11-iooperator(3)*imu11
       rmu00 = rmu00+(a+a)
      elseif(nmultipoles.lt.0) then
       call bummer('fmmgrd: error, nmultipoles = ',nmultipoles)
      endif
      return
      end subroutine fmmgrd
#endif
