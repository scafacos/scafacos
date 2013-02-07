#include "../../fmm.h"
c
      subroutine fmmoopn(n,xyzfrom,xyzto,f,g,h,roop,ioop)
c
      use fmmkinds
c
      implicit none
c
      real(kind=fmm_real) xyzfrom(*),xyzto(*),f(0:*),g(0:*),h(0:*),
     .roop(*),ioop(*)
c
      integer(kind=fmm_integer) n
c
#if FMM_MAXNMULTIPOLES == 0 && defined(FMM_UNROLLED)
      go to 1
#elif FMM_MAXNMULTIPOLES == 1 && defined(FMM_UNROLLED)
      go to(1,2) n
#elif FMM_MAXNMULTIPOLES == 2 && defined(FMM_UNROLLED)
      go to(1,2,3) n
#elif FMM_MAXNMULTIPOLES == 3 && defined(FMM_UNROLLED)
      go to(1,2,3,4) n
#elif FMM_MAXNMULTIPOLES == 4 && defined(FMM_UNROLLED)
      go to(1,2,3,4,5) n
#elif FMM_MAXNMULTIPOLES == 5 && defined(FMM_UNROLLED)
      go to(1,2,3,4,5,6) n
#elif FMM_MAXNMULTIPOLES == 6 && defined(FMM_UNROLLED)
      go to(1,2,3,4,5,6,7) n
#elif FMM_MAXNMULTIPOLES == 7 && defined(FMM_UNROLLED)
      go to(1,2,3,4,5,6,7,8) n
#elif FMM_MAXNMULTIPOLES == 8 && defined(FMM_UNROLLED)
      go to(1,2,3,4,5,6,7,8,9) n
#elif FMM_MAXNMULTIPOLES == 9 && defined(FMM_UNROLLED)
      go to(1,2,3,4,5,6,7,8,9,10) n
#elif FMM_MAXNMULTIPOLES == 10 && defined(FMM_UNROLLED)
      go to(1,2,3,4,5,6,7,8,9,10,11) n
#elif FMM_MAXNMULTIPOLES == 11 && defined(FMM_UNROLLED)
      go to(1,2,3,4,5,6,7,8,9,10,11,12) n
#elif FMM_MAXNMULTIPOLES == 12 && defined(FMM_UNROLLED)
      go to(1,2,3,4,5,6,7,8,9,10,11,12,13) n
#elif FMM_MAXNMULTIPOLES == 13 && defined(FMM_UNROLLED)
      go to(1,2,3,4,5,6,7,8,9,10,11,12,13,14) n
#elif FMM_MAXNMULTIPOLES == 14 && defined(FMM_UNROLLED)
      go to(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15) n
#elif FMM_MAXNMULTIPOLES == 15 && defined(FMM_UNROLLED)
      go to(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16) n
#elif FMM_MAXNMULTIPOLES == 16 && defined(FMM_UNROLLED)
      go to(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17) n
#elif FMM_MAXNMULTIPOLES == 17 && defined(FMM_UNROLLED)
      go to(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18) n
#elif FMM_MAXNMULTIPOLES == 18 && defined(FMM_UNROLLED)
      go to(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19) n
#elif FMM_MAXNMULTIPOLES == 19 && defined(FMM_UNROLLED)
      go to(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20) n
#elif FMM_MAXNMULTIPOLES == 20 && defined(FMM_UNROLLED)
      go to(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21) n
#elif FMM_MAXNMULTIPOLES == 21 && defined(FMM_UNROLLED)
      go to(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22) n
#elif FMM_MAXNMULTIPOLES == 22 && defined(FMM_UNROLLED)
      go to(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,
     .23) n
#elif FMM_MAXNMULTIPOLES == 23 && defined(FMM_UNROLLED)
      go to(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,
     .23,24) n
#elif FMM_MAXNMULTIPOLES == 24 && defined(FMM_UNROLLED)
      go to(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,
     .23,24,25) n
#elif FMM_MAXNMULTIPOLES == 25 && defined(FMM_UNROLLED)
      go to(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,
     .23,24,25,26) n
#elif FMM_MAXNMULTIPOLES == 26 && defined(FMM_UNROLLED)
      go to(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,
     .23,24,25,26,27) n
#elif FMM_MAXNMULTIPOLES == 27 && defined(FMM_UNROLLED)
      go to(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,
     .23,24,25,26,27,28) n
#elif FMM_MAXNMULTIPOLES == 28 && defined(FMM_UNROLLED)
      go to(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,
     .23,24,25,26,27,28,29) n
#elif FMM_MAXNMULTIPOLES == 29 && defined(FMM_UNROLLED)
      go to(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,
     .23,24,25,26,27,28,29,30) n
#elif FMM_MAXNMULTIPOLES == 30 && defined(FMM_UNROLLED)
      go to(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,
     .23,24,25,26,27,28,29,30,31) n
#elif FMM_MAXNMULTIPOLES == 31 && defined(FMM_UNROLLED)
      go to(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,
     .23,24,25,26,27,28,29,30,31,32) n
#elif FMM_MAXNMULTIPOLES == 32 && defined(FMM_UNROLLED)
      go to(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,
     .23,24,25,26,27,28,29,30,31,32,33) n
#elif FMM_MAXNMULTIPOLES == 33 && defined(FMM_UNROLLED)
      go to(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,
     .23,24,25,26,27,28,29,30,31,32,33,34) n
#elif FMM_MAXNMULTIPOLES == 34 && defined(FMM_UNROLLED)
      go to(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,
     .23,24,25,26,27,28,29,30,31,32,33,34,35) n
#elif FMM_MAXNMULTIPOLES == 35 && defined(FMM_UNROLLED)
      go to(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,
     .23,24,25,26,27,28,29,30,31,32,33,34,35,36) n
#elif FMM_MAXNMULTIPOLES == 36 && defined(FMM_UNROLLED)
      go to(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,
     .23,24,25,26,27,28,29,30,31,32,33,34,35,36,37) n
#elif FMM_MAXNMULTIPOLES == 37 && defined(FMM_UNROLLED)
      go to(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,
     .23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38) n
#elif FMM_MAXNMULTIPOLES == 38 && defined(FMM_UNROLLED)
      go to(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,
     .23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39) n
#elif FMM_MAXNMULTIPOLES == 39 && defined(FMM_UNROLLED)
      go to(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,
     .23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40) n
#elif FMM_MAXNMULTIPOLES == 40 && defined(FMM_UNROLLED)
      go to(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,
     .23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41) n
#elif FMM_MAXNMULTIPOLES == 41 && defined(FMM_UNROLLED)
      go to(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,
     .23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42) n
#elif FMM_MAXNMULTIPOLES == 42 && defined(FMM_UNROLLED)
      go to(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,
     .23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43) n
#elif FMM_MAXNMULTIPOLES == 43 && defined(FMM_UNROLLED)
      go to(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,
     .23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,
     .44) n
#elif FMM_MAXNMULTIPOLES == 44 && defined(FMM_UNROLLED)
      go to(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,
     .23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,
     .44,45) n
#elif FMM_MAXNMULTIPOLES == 45 && defined(FMM_UNROLLED)
      go to(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,
     .23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,
     .44,45,46) n
#elif FMM_MAXNMULTIPOLES == 46 && defined(FMM_UNROLLED)
      go to(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,
     .23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,
     .44,45,46,47) n
#elif FMM_MAXNMULTIPOLES == 47 && defined(FMM_UNROLLED)
      go to(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,
     .23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,
     .44,45,46,47,48) n
#elif FMM_MAXNMULTIPOLES == 48 && defined(FMM_UNROLLED)
      go to(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,
     .23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,
     .44,45,46,47,48,49) n
#elif FMM_MAXNMULTIPOLES == 49 && defined(FMM_UNROLLED)
      go to(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,
     .23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,
     .44,45,46,47,48,49,50) n
#elif FMM_MAXNMULTIPOLES == 50 && defined(FMM_UNROLLED)
      go to(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,
     .23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,
     .44,45,46,47,48,49,50,51) n
#else
      call fmmooperator((n-1),xyzfrom,xyzto,f,g,h,roop,ioop)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > -1 && defined(FMM_UNROLLED)
 1    roop(1) = f(0)
      ioop(1) = g(0)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 0 && defined(FMM_UNROLLED)
 2    call fmmoop1(xyzfrom,xyzto,f,g,roop,ioop)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 1 && defined(FMM_UNROLLED)
 3    call fmmoop2(xyzfrom,xyzto,f,g,h,roop,ioop)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 2 && defined(FMM_UNROLLED)
 4    call fmmoop3(xyzfrom,xyzto,f,g,h,roop,ioop)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 3 && defined(FMM_UNROLLED)
 5    call fmmoop4(xyzfrom,xyzto,f,g,h,roop,ioop)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 4 && defined(FMM_UNROLLED)
 6    call fmmoop5(xyzfrom,xyzto,f,g,h,roop,ioop)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 5 && defined(FMM_UNROLLED)
 7    call fmmoop6(xyzfrom,xyzto,f,g,h,roop,ioop)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 6 && defined(FMM_UNROLLED)
 8    call fmmoop7(xyzfrom,xyzto,f,g,h,roop,ioop)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 7 && defined(FMM_UNROLLED)
 9    call fmmoop8(xyzfrom,xyzto,f,g,h,roop,ioop)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 8 && defined(FMM_UNROLLED)
 10   call fmmoop9(xyzfrom,xyzto,f,g,h,roop,ioop)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 9 && defined(FMM_UNROLLED)
 11   call fmmoop10(xyzfrom,xyzto,f,g,h,roop,ioop)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 10 && defined(FMM_UNROLLED)
 12   call fmmoop11(xyzfrom,xyzto,f,g,h,roop,ioop)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 11 && defined(FMM_UNROLLED)
 13   call fmmoop12(xyzfrom,xyzto,f,g,h,roop,ioop)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 12 && defined(FMM_UNROLLED)
 14   call fmmoop13(xyzfrom,xyzto,f,g,h,roop,ioop)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 13 && defined(FMM_UNROLLED)
 15   call fmmoop14(xyzfrom,xyzto,f,g,h,roop,ioop)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 14 && defined(FMM_UNROLLED)
 16   call fmmoop15(xyzfrom,xyzto,f,g,h,roop,ioop)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 15 && defined(FMM_UNROLLED)
 17   call fmmoop16(xyzfrom,xyzto,f,g,h,roop,ioop)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 16 && defined(FMM_UNROLLED)
 18   call fmmoop17(xyzfrom,xyzto,f,g,h,roop,ioop)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 17 && defined(FMM_UNROLLED)
 19   call fmmoop18(xyzfrom,xyzto,f,g,h,roop,ioop)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 18 && defined(FMM_UNROLLED)
 20   call fmmoop19(xyzfrom,xyzto,f,g,h,roop,ioop)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 19 && defined(FMM_UNROLLED)
 21   call fmmoop20(xyzfrom,xyzto,f,g,h,roop,ioop)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 20 && defined(FMM_UNROLLED)
 22   call fmmoop21(xyzfrom,xyzto,f,g,h,roop,ioop)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 21 && defined(FMM_UNROLLED)
 23   call fmmoop22(xyzfrom,xyzto,f,g,h,roop,ioop)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 22 && defined(FMM_UNROLLED)
 24   call fmmoop23(xyzfrom,xyzto,f,g,h,roop,ioop)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 23 && defined(FMM_UNROLLED)
 25   call fmmoop24(xyzfrom,xyzto,f,g,h,roop,ioop)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 24 && defined(FMM_UNROLLED)
 26   call fmmoop25(xyzfrom,xyzto,f,g,h,roop,ioop)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 25 && defined(FMM_UNROLLED)
 27   call fmmoop26(xyzfrom,xyzto,f,g,h,roop,ioop)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 26 && defined(FMM_UNROLLED)
 28   call fmmoop27(xyzfrom,xyzto,f,g,h,roop,ioop)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 27 && defined(FMM_UNROLLED)
 29   call fmmoop28(xyzfrom,xyzto,f,g,h,roop,ioop)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 28 && defined(FMM_UNROLLED)
 30   call fmmoop29(xyzfrom,xyzto,f,g,h,roop,ioop)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 29 && defined(FMM_UNROLLED)
 31   call fmmoop30(xyzfrom,xyzto,f,g,h,roop,ioop)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 30 && defined(FMM_UNROLLED)
 32   call fmmoop31(xyzfrom,xyzto,f,g,h,roop,ioop)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 31 && defined(FMM_UNROLLED)
 33   call fmmoop32(xyzfrom,xyzto,f,g,h,roop,ioop)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 32 && defined(FMM_UNROLLED)
 34   call fmmoop33(xyzfrom,xyzto,f,g,h,roop,ioop)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 33 && defined(FMM_UNROLLED)
 35   call fmmoop34(xyzfrom,xyzto,f,g,h,roop,ioop)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 34 && defined(FMM_UNROLLED)
 36   call fmmoop35(xyzfrom,xyzto,f,g,h,roop,ioop)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 35 && defined(FMM_UNROLLED)
 37   call fmmoop36(xyzfrom,xyzto,f,g,h,roop,ioop)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 36 && defined(FMM_UNROLLED)
 38   call fmmoop37(xyzfrom,xyzto,f,g,h,roop,ioop)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 37 && defined(FMM_UNROLLED)
 39   call fmmoop38(xyzfrom,xyzto,f,g,h,roop,ioop)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 38 && defined(FMM_UNROLLED)
 40   call fmmoop39(xyzfrom,xyzto,f,g,h,roop,ioop)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 39 && defined(FMM_UNROLLED)
 41   call fmmoop40(xyzfrom,xyzto,f,g,h,roop,ioop)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 40 && defined(FMM_UNROLLED)
 42   call fmmoop41(xyzfrom,xyzto,f,g,h,roop,ioop)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 41 && defined(FMM_UNROLLED)
 43   call fmmoop42(xyzfrom,xyzto,f,g,h,roop,ioop)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 42 && defined(FMM_UNROLLED)
 44   call fmmoop43(xyzfrom,xyzto,f,g,h,roop,ioop)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 43 && defined(FMM_UNROLLED)
 45   call fmmoop44(xyzfrom,xyzto,f,g,h,roop,ioop)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 44 && defined(FMM_UNROLLED)
 46   call fmmoop45(xyzfrom,xyzto,f,g,h,roop,ioop)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 45 && defined(FMM_UNROLLED)
 47   call fmmoop46(xyzfrom,xyzto,f,g,h,roop,ioop)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 46 && defined(FMM_UNROLLED)
 48   call fmmoop47(xyzfrom,xyzto,f,g,h,roop,ioop)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 47 && defined(FMM_UNROLLED)
 49   call fmmoop48(xyzfrom,xyzto,f,g,h,roop,ioop)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 48 && defined(FMM_UNROLLED)
 50   call fmmoop49(xyzfrom,xyzto,f,g,h,roop,ioop)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 49 && defined(FMM_UNROLLED)
 51   call fmmoop50(xyzfrom,xyzto,f,g,h,roop,ioop)
c
      return
#endif
      end subroutine fmmoopn
