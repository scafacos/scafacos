#include "../../fmm.h"
c
#if FMM_MAXNMULTIPOLES > 5 && defined(FMM_UNROLLED)
      subroutine fmmoop6(xyzfrom,xyzto,f,g,h,roop,ioop)
c
      use fmmkinds
c
      implicit none
c
      real(kind=fmm_real) xyzfrom(*),xyzto(*),f(0:*),g(0:*),h(0:*),
     .roop(*),ioop(*),
     .one,zero,x,y,z,rsq,s
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
      s = h(   1)*z
      roop(   7) = g(   2)*(s*roop(   4)-rsq*roop(   2))
      ioop(   7) = zero
      roop(   8) = g(   3)*(s*roop(   5)-rsq*roop(   3))
      ioop(   8) = g(   3)*(s*ioop(   5)-rsq*ioop(   3))
      roop(   9) = z*roop(   6)
      ioop(   9) = z*ioop(   6)
      roop(  10) = f(   3)*(x*roop(   6)+y*ioop(   6))
      ioop(  10) = f(   3)*(x*ioop(   6)-y*roop(   6))
      s = h(   2)*z
      roop(  11) = g(   4)*(s*roop(   7)-rsq*roop(   4))
      ioop(  11) = zero
      roop(  12) = g(   5)*(s*roop(   8)-rsq*roop(   5))
      ioop(  12) = g(   5)*(s*ioop(   8)-rsq*ioop(   5))
      roop(  13) = g(   6)*(s*roop(   9)-rsq*roop(   6))
      ioop(  13) = g(   6)*(s*ioop(   9)-rsq*ioop(   6))
      roop(  14) = z*roop(  10)
      ioop(  14) = z*ioop(  10)
      roop(  15) = f(   4)*(x*roop(  10)+y*ioop(  10))
      ioop(  15) = f(   4)*(x*ioop(  10)-y*roop(  10))
      s = h(   3)*z
      roop(  16) = g(   7)*(s*roop(  11)-rsq*roop(   7))
      ioop(  16) = zero
      roop(  17) = g(   8)*(s*roop(  12)-rsq*roop(   8))
      ioop(  17) = g(   8)*(s*ioop(  12)-rsq*ioop(   8))
      roop(  18) = g(   9)*(s*roop(  13)-rsq*roop(   9))
      ioop(  18) = g(   9)*(s*ioop(  13)-rsq*ioop(   9))
      roop(  19) = g(  10)*(s*roop(  14)-rsq*roop(  10))
      ioop(  19) = g(  10)*(s*ioop(  14)-rsq*ioop(  10))
      roop(  20) = z*roop(  15)
      ioop(  20) = z*ioop(  15)
      roop(  21) = f(   5)*(x*roop(  15)+y*ioop(  15))
      ioop(  21) = f(   5)*(x*ioop(  15)-y*roop(  15))
      s = h(   4)*z
      roop(  22) = g(  11)*(s*roop(  16)-rsq*roop(  11))
      ioop(  22) = zero
      roop(  23) = g(  12)*(s*roop(  17)-rsq*roop(  12))
      ioop(  23) = g(  12)*(s*ioop(  17)-rsq*ioop(  12))
      roop(  24) = g(  13)*(s*roop(  18)-rsq*roop(  13))
      ioop(  24) = g(  13)*(s*ioop(  18)-rsq*ioop(  13))
      roop(  25) = g(  14)*(s*roop(  19)-rsq*roop(  14))
      ioop(  25) = g(  14)*(s*ioop(  19)-rsq*ioop(  14))
      roop(  26) = g(  15)*(s*roop(  20)-rsq*roop(  15))
      ioop(  26) = g(  15)*(s*ioop(  20)-rsq*ioop(  15))
      roop(  27) = z*roop(  21)
      ioop(  27) = z*ioop(  21)
      roop(  28) = f(   6)*(x*roop(  21)+y*ioop(  21))
      ioop(  28) = f(   6)*(x*ioop(  21)-y*roop(  21))
      return
      end subroutine fmmoop6
#endif
