#include "../../fmm.h"
#if FMM_MAXNMULTIPOLES < 0 || FMM_MAXNMULTIPOLES > 50 || !defined(FMM_UNROLLED)
      subroutine fmmooperator(nmultipoles,xyzfrom,xyzto,f,g,h,
     .rooperator,iooperator)
c
      use fmmkinds
c
      implicit none
c
      real(kind=fmm_real) xyzfrom(*),xyzto(*),f(0:*),g(0:*),h(0:*),
     .rooperator(*),
     .iooperator(*),one,zero,x,y,z,rsq,s
c
      integer(kind=fmm_integer) nmultipoles,i,j,k,n,l,m
c
      one = f(0)
      zero = g(0)
c
      rooperator(1) = one
      iooperator(1) = zero
c
      if(nmultipoles.gt.1) then
         x = xyzfrom(1)-xyzto(1)
         y = xyzfrom(2)-xyzto(2)
         z = xyzfrom(3)-xyzto(3)
c
         rsq = x*x+y*y+z*z
c
         rooperator(2) = z*rooperator(1)
         iooperator(2) = zero
c
         rooperator(3) = f(1)*(x*rooperator(1))
         iooperator(3) = -f(1)*(y*rooperator(1))
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
               rooperator(j) = g(n)*(s*rooperator(k)-rsq*rooperator(n))
               iooperator(j) = g(n)*(s*iooperator(k)-rsq*iooperator(n))
 2          continue
c
            j = j+1
            k = k+1
c
            rooperator(j) = z*rooperator(k)
            iooperator(j) = z*iooperator(k)
c
            j = j+1
c
            rooperator(j) = f(l)*(x*rooperator(k)+y*iooperator(k))
            iooperator(j) = f(l)*(x*iooperator(k)-y*rooperator(k))
 1       continue
      elseif(nmultipoles.eq.1) then
         rooperator(2) = (xyzfrom(3)-xyzto(3))*rooperator(1)
         iooperator(2) = zero
c
         rooperator(3) = f(1)*((xyzfrom(1)-xyzto(1))*rooperator(1))
         iooperator(3) = f(1)*((xyzto(2)-xyzfrom(2))*rooperator(1))
      elseif(nmultipoles.lt.0) then
         call bummer('fmmooperator: error, nmultipoles = ',nmultipoles)
      endif
      return
      end subroutine fmmooperator
#endif
