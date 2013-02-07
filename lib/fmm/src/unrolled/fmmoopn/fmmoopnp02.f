#include "../../fmm.h"
c
#if FMM_MAXNMULTIPOLES > 1 && defined(FMM_UNROLLED)
      subroutine fmmoop2(xyzfrom,xyzto,f,g,h,roop,ioop)
c
      use fmmkinds
c
      implicit none
c
      real(kind=fmm_real) xyzfrom(*),xyzto(*),f(0:*),g(0:*),h(0:*),
     .roop(*),ioop(*),
     .one,zero,x,y,z,rsq
c
      one = f(0)
      zero = g(0)
      x = xyzfrom(1)-xyzto(1)
      y = xyzfrom(2)-xyzto(2)
      z = xyzfrom(3)-xyzto(3)
      rsq = x*x+y*y+z*z
      roop(   1) = one
      ioop(   1) = zero
      roop(   2) = z*roop(1)
      ioop(   2) = zero
      roop(   3) = f(1)*(x*roop(1))
      ioop(   3) = -f(1)*(y*roop(1))
      roop(   4) = g(   1)*((h(   0)*z)*roop(   2)-rsq*roop(   1))
      ioop(   4) = zero
      roop(   5) = z*roop(   3)
      ioop(   5) = z*ioop(   3)
      roop(   6) = f(   2)*(x*roop(   3)+y*ioop(   3))
      ioop(   6) = f(   2)*(x*ioop(   3)-y*roop(   3))
      return
      end subroutine fmmoop2
#endif
