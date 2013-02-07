#include "../../fmm.h"
c
#if FMM_MAXNMULTIPOLES > 3 && defined(FMM_UNROLLED)
      subroutine fmmgrd4(roop,ioop,rmu,imu,rmu00,rmu10,rmu11,imu11)
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
      rmu10 = rmu10+roop(   7)*rmu(  11)
      rmu11 = rmu11-roop(   8)*rmu(  11)
      imu11 = imu11+ioop(   8)*rmu(  11)
      rmu00 = rmu00+roop(  11)*rmu(  11)
      rmu11 = rmu11+roop(   7)*rmu(  12)
      imu11 = imu11+roop(   7)*imu(  12)
      a = roop(   8)*rmu(  12)-ioop(   8)*imu(  12)
      rmu11 = rmu11-(roop(   9)*rmu(  12)-ioop(   9)*imu(  12))
      imu11 = imu11+(roop(   9)*imu(  12)+ioop(   9)*rmu(  12))
      b = roop(  12)*rmu(  12)-ioop(  12)*imu(  12)
      rmu11 = rmu11+(roop(   8)*rmu(  13)-ioop(   8)*imu(  13))
      imu11 = imu11+(roop(   8)*imu(  13)+ioop(   8)*rmu(  13))
      a = a+(roop(   9)*rmu(  13)-ioop(   9)*imu(  13))
      rmu11 = rmu11-(roop(  10)*rmu(  13)-ioop(  10)*imu(  13))
      imu11 = imu11+(roop(  10)*imu(  13)+ioop(  10)*rmu(  13))
      b = b+(roop(  13)*rmu(  13)-ioop(  13)*imu(  13))
      rmu11 = rmu11+(roop(   9)*rmu(  14)-ioop(   9)*imu(  14))
      imu11 = imu11+(roop(   9)*imu(  14)+ioop(   9)*rmu(  14))
      a = a+(roop(  10)*rmu(  14)-ioop(  10)*imu(  14))
      b = b+(roop(  14)*rmu(  14)-ioop(  14)*imu(  14))
      rmu11 = rmu11+(roop(  10)*rmu(  15)-ioop(  10)*imu(  15))
      imu11 = imu11+(roop(  10)*imu(  15)+ioop(  10)*rmu(  15))
      b = b+(roop(  15)*rmu(  15)-ioop(  15)*imu(  15))
      rmu10 = rmu10+(a+a)
      rmu00 = rmu00+(b+b)
      return
      end subroutine fmmgrd4
#endif
