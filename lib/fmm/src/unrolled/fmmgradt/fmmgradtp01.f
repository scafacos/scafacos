#include "../../fmm.h"
c
#if FMM_MAXNMULTIPOLES > 0 && defined(FMM_UNROLLED)
      subroutine fmmgrd1(roop,ioop,rmu,imu,rmu00,rmu10,rmu11,imu11)
c
      use fmmkinds
c
      implicit none
c
      real(kind=fmm_real) roop(*),ioop(*),rmu(*),imu(*),rmu00,rmu10,
     .rmu11,imu11,a
c
      rmu00 = rmu(1)
      rmu10 = rmu(2)
      rmu00 = rmu00+roop(2)*rmu10
      rmu11 = rmu(3)
      imu11 = imu(3)
      a = roop(3)*rmu11-ioop(3)*imu11
      rmu00 = rmu00+(a+a)
      return
      end subroutine fmmgrd1
#endif
