#include "../../fmm.h"
c
#if FMM_MAXNMULTIPOLES > 0 && defined(FMM_UNROLLED)
      subroutine fmmmop1(q,x,y,z,rsq,g,rmop,imop)
c
      use fmmkinds
c
      implicit none
c
      real(kind=fmm_real) q,x,y,z,rsq,g(0:*),rmop(*),imop(*),zero,s
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
      return
      end subroutine fmmmop1
#endif
