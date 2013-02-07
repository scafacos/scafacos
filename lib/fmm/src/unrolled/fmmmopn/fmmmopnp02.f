#include "../../fmm.h"
c
#if FMM_MAXNMULTIPOLES > 1 && defined(FMM_UNROLLED)
      subroutine fmmmop2(q,x,y,z,rsq,g,h,rmop,imop)
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
      return
      end subroutine fmmmop2
#endif
