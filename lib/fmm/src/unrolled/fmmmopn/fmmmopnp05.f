#include "../../fmm.h"
c
#if FMM_MAXNMULTIPOLES > 4 && defined(FMM_UNROLLED)
      subroutine fmmmop5(q,x,y,z,rsq,g,h,rmop,imop)
c
      use fmmkinds
c
      implicit none
c
      real(kind=fmm_real) q,x,y,z,rsq,g(0:*),h(0:*),rmop(*),imop(*),
     .zero,s
c
      rsq = rsq*rsq
      x = x*rsq
      y = y*rsq
      z = z*rsq
      zero = g(0)
      rmop(1) = q
      imop(1) = zero
      rmop(2) = z*rmop(1)
      imop(2) = zero
      rmop(3) = x*rmop(1)
      imop(3) = y*rmop(1)
      s = h(   0)*z
      rmop(   4) = s*rmop(   2)-g(   1)*rsq*rmop(   1)
      imop(   4) = zero
      rmop(   5) = s*rmop(   3)
      imop(   5) = s*imop(   3)
      rmop(   6) = h(   0)*(x*rmop(   3)-y*imop(   3))
      imop(   6) = h(   0)*(x*imop(   3)+y*rmop(   3))
      s = h(   1)*z
      rmop(   7) = s*rmop(   4)-g(   2)*rsq*rmop(   2)
      imop(   7) = zero
      rmop(   8) = s*rmop(   5)-g(   3)*rsq*rmop(   3)
      imop(   8) = s*imop(   5)-g(   3)*rsq*imop(   3)
      rmop(   9) = s*rmop(   6)
      imop(   9) = s*imop(   6)
      rmop(  10) = h(   1)*(x*rmop(   6)-y*imop(   6))
      imop(  10) = h(   1)*(x*imop(   6)+y*rmop(   6))
      s = h(   2)*z
      rmop(  11) = s*rmop(   7)-g(   4)*rsq*rmop(   4)
      imop(  11) = zero
      rmop(  12) = s*rmop(   8)-g(   5)*rsq*rmop(   5)
      imop(  12) = s*imop(   8)-g(   5)*rsq*imop(   5)
      rmop(  13) = s*rmop(   9)-g(   6)*rsq*rmop(   6)
      imop(  13) = s*imop(   9)-g(   6)*rsq*imop(   6)
      rmop(  14) = s*rmop(  10)
      imop(  14) = s*imop(  10)
      rmop(  15) = h(   2)*(x*rmop(  10)-y*imop(  10))
      imop(  15) = h(   2)*(x*imop(  10)+y*rmop(  10))
      s = h(   3)*z
      rmop(  16) = s*rmop(  11)-g(   7)*rsq*rmop(   7)
      imop(  16) = zero
      rmop(  17) = s*rmop(  12)-g(   8)*rsq*rmop(   8)
      imop(  17) = s*imop(  12)-g(   8)*rsq*imop(   8)
      rmop(  18) = s*rmop(  13)-g(   9)*rsq*rmop(   9)
      imop(  18) = s*imop(  13)-g(   9)*rsq*imop(   9)
      rmop(  19) = s*rmop(  14)-g(  10)*rsq*rmop(  10)
      imop(  19) = s*imop(  14)-g(  10)*rsq*imop(  10)
      rmop(  20) = s*rmop(  15)
      imop(  20) = s*imop(  15)
      rmop(  21) = h(   3)*(x*rmop(  15)-y*imop(  15))
      imop(  21) = h(   3)*(x*imop(  15)+y*rmop(  15))
      return
      end subroutine fmmmop5
#endif
