#include "../../fmm.h"
c
#if FMM_MAXNMULTIPOLES > 2 && defined(FMM_UNROLLED)
      subroutine fmmgrd3(roop,ioop,rmu,imu,rmu00,rmu10,rmu11,imu11)
c
      use fmmkinds
c
      implicit none
c
      real(kind=fmm_real) roop(*),ioop(*),rmu(*),imu(*),rmu00,rmu10,
     .rmu11,imu11,a,b
c
      rmu00 = rmu(1)
      rmu10 = rmu(2)
      rmu00 = rmu00+roop(2)*rmu10
      rmu11 = rmu(3)
      imu11 = imu(3)
      a = roop(3)*rmu11-ioop(3)*imu11
      rmu00 = rmu00+(a+a)
      rmu10 = rmu10+roop(2)*rmu(4)
      rmu11 = rmu11-roop(3)*rmu(4)
      imu11 = imu11+ioop(3)*rmu(4)
      rmu00 = rmu00+roop(4)*rmu(4)
      rmu11 = rmu11+roop(2)*rmu(5)
      imu11 = imu11+roop(2)*imu(5)
      a = roop(3)*rmu(5)-ioop(3)*imu(5)
      b = roop(5)*rmu(5)-ioop(5)*imu(5)
      rmu11 = rmu11+(roop(3)*rmu(6)-ioop(3)*imu(6))
      imu11 = imu11+(roop(3)*imu(6)+ioop(3)*rmu(6))
      b = b+(roop(6)*rmu(6)-ioop(6)*imu(6))
      rmu10 = rmu10+(a+a)
      rmu00 = rmu00+(b+b)
      rmu10 = rmu10+roop(   4)*rmu(   7)
      rmu11 = rmu11-roop(   5)*rmu(   7)
      imu11 = imu11+ioop(   5)*rmu(   7)
      rmu00 = rmu00+roop(   7)*rmu(   7)
      rmu11 = rmu11+roop(   4)*rmu(   8)
      imu11 = imu11+roop(   4)*imu(   8)
      a = roop(   5)*rmu(   8)-ioop(   5)*imu(   8)
      rmu11 = rmu11-(roop(   6)*rmu(   8)-ioop(   6)*imu(   8))
      imu11 = imu11+(roop(   6)*imu(   8)+ioop(   6)*rmu(   8))
      b = roop(   8)*rmu(   8)-ioop(   8)*imu(   8)
      rmu11 = rmu11+(roop(   5)*rmu(   9)-ioop(   5)*imu(   9))
      imu11 = imu11+(roop(   5)*imu(   9)+ioop(   5)*rmu(   9))
      a = a+(roop(   6)*rmu(   9)-ioop(   6)*imu(   9))
      b = b+(roop(   9)*rmu(   9)-ioop(   9)*imu(   9))
      rmu11 = rmu11+(roop(   6)*rmu(  10)-ioop(   6)*imu(  10))
      imu11 = imu11+(roop(   6)*imu(  10)+ioop(   6)*rmu(  10))
      b = b+(roop(  10)*rmu(  10)-ioop(  10)*imu(  10))
      rmu10 = rmu10+(a+a)
      rmu00 = rmu00+(b+b)
      return
      end subroutine fmmgrd3
#endif
