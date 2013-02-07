#include "../../fmm.h"
c
#if FMM_MAXNMULTIPOLES > 0 && defined(FMM_UNROLLED)
      subroutine fmmoop1(xyzfrom,xyzto,f,g,roop,ioop)
c
      use fmmkinds
c
      implicit none
c
      real(kind=fmm_real) xyzfrom(*),xyzto(*),f(0:*),g(0:*),roop(*),
     .ioop(*),
     .one,zero,x,y,z
c
      one = f(0)
      zero = g(0)
      roop(   1) = one
      ioop(   1) = zero
      roop(   2) = (xyzfrom(3)-xyzto(3))*roop(1)
      ioop(   2) = zero
      roop(   3) = f(1)*((xyzfrom(1)-xyzto(1))*roop(1))
      ioop(   3) = f(1)*((xyzto(2)-xyzfrom(2))*roop(1))
      return
      end subroutine fmmoop1
#endif
