#include "../../fmm.h"
#if FMM_MAXNMULTIPOLES < 0 || FMM_MAXNMULTIPOLES > 50 || !defined(FMM_UNROLLED)
      subroutine fmmmoperator(nmultipoles,q,x,y,z,rsq,g,h,rmoperator,
     .imoperator)
c
      use fmmkinds
c
      implicit none
c
      real(kind=fmm_real) q,x,y,z,rsq,g(0:*),h(0:*),rmoperator(*),
     .imoperator(*),
     .zero,s
c
      integer(kind=fmm_integer) nmultipoles,i,j,k,n,l,m
c
      rsq = rsq*rsq
c
      x = x*rsq
      y = y*rsq
      z = z*rsq
c
      zero = g(0)
c
      rmoperator(1) = q
      imoperator(1) = zero
c
      if(nmultipoles.gt.1) then
         rmoperator(2) = z*rmoperator(1)
         imoperator(2) = zero
c
         rmoperator(3) = x*rmoperator(1)
         imoperator(3) = y*rmoperator(1)
c
         i = -1
         j = 3
         k = 1
         n = 0
c
         do 1 l = 2,nmultipoles
            i = i+1
            s = h(i)*z
c
            do 2 m = 0,i
               j = j+1
               k = k+1
               n = n+1
               rmoperator(j) = s*rmoperator(k)-g(n)*rsq*rmoperator(n)
               imoperator(j) = s*imoperator(k)-g(n)*rsq*imoperator(n)
 2          continue
c
            j = j+1
            k = k+1
c
            rmoperator(j) = s*rmoperator(k)
            imoperator(j) = s*imoperator(k)
c
            j = j+1
c
            rmoperator(j) = h(i)*(x*rmoperator(k)-y*imoperator(k))
            imoperator(j) = h(i)*(x*imoperator(k)+y*rmoperator(k))
 1       continue
      elseif(nmultipoles.eq.1) then
         rmoperator(2) = z*rmoperator(1)
         imoperator(2) = zero
c
         rmoperator(3) = x*rmoperator(1)
         imoperator(3) = y*rmoperator(1)
      elseif(nmultipoles.lt.0) then
         call bummer('fmmmoperator: error, nmultipoles = ',nmultipoles)
      endif
      return
      end subroutine fmmmoperator
#endif
