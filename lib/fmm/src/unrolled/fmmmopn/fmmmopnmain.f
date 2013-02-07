#include "../../fmm.h"
c
      subroutine fmmmopn(n,q,x,y,z,rsq,g,h,rmop,imop)
c
      use fmmkinds
c
      implicit none
c
      real(kind=fmm_real) q,x,y,z,rsq,g(0:*),h(0:*),rmop(*),imop(*)
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
      call fmmmoperator((n-1),q,x,y,z,rsq,g,h,rmop,imop)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > -1 && defined(FMM_UNROLLED)
 1    rmop(1) = q
      imop(1) = g(0)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 0 && defined(FMM_UNROLLED)
 2    call fmmmop1(q,x,y,z,rsq,g,rmop,imop)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 1 && defined(FMM_UNROLLED)
 3    call fmmmop2(q,x,y,z,rsq,g,h,rmop,imop)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 2 && defined(FMM_UNROLLED)
 4    call fmmmop3(q,x,y,z,rsq,g,h,rmop,imop)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 3 && defined(FMM_UNROLLED)
 5    call fmmmop4(q,x,y,z,rsq,g,h,rmop,imop)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 4 && defined(FMM_UNROLLED)
 6    call fmmmop5(q,x,y,z,rsq,g,h,rmop,imop)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 5 && defined(FMM_UNROLLED)
 7    call fmmmop6(q,x,y,z,rsq,g,h,rmop,imop)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 6 && defined(FMM_UNROLLED)
 8    call fmmmop7(q,x,y,z,rsq,g,h,rmop,imop)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 7 && defined(FMM_UNROLLED)
 9    call fmmmop8(q,x,y,z,rsq,g,h,rmop,imop)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 8 && defined(FMM_UNROLLED)
 10   call fmmmop9(q,x,y,z,rsq,g,h,rmop,imop)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 9 && defined(FMM_UNROLLED)
 11   call fmmmop10(q,x,y,z,rsq,g,h,rmop,imop)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 10 && defined(FMM_UNROLLED)
 12   call fmmmop11(q,x,y,z,rsq,g,h,rmop,imop)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 11 && defined(FMM_UNROLLED)
 13   call fmmmop12(q,x,y,z,rsq,g,h,rmop,imop)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 12 && defined(FMM_UNROLLED)
 14   call fmmmop13(q,x,y,z,rsq,g,h,rmop,imop)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 13 && defined(FMM_UNROLLED)
 15   call fmmmop14(q,x,y,z,rsq,g,h,rmop,imop)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 14 && defined(FMM_UNROLLED)
 16   call fmmmop15(q,x,y,z,rsq,g,h,rmop,imop)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 15 && defined(FMM_UNROLLED)
 17   call fmmmop16(q,x,y,z,rsq,g,h,rmop,imop)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 16 && defined(FMM_UNROLLED)
 18   call fmmmop17(q,x,y,z,rsq,g,h,rmop,imop)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 17 && defined(FMM_UNROLLED)
 19   call fmmmop18(q,x,y,z,rsq,g,h,rmop,imop)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 18 && defined(FMM_UNROLLED)
 20   call fmmmop19(q,x,y,z,rsq,g,h,rmop,imop)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 19 && defined(FMM_UNROLLED)
 21   call fmmmop20(q,x,y,z,rsq,g,h,rmop,imop)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 20 && defined(FMM_UNROLLED)
 22   call fmmmop21(q,x,y,z,rsq,g,h,rmop,imop)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 21 && defined(FMM_UNROLLED)
 23   call fmmmop22(q,x,y,z,rsq,g,h,rmop,imop)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 22 && defined(FMM_UNROLLED)
 24   call fmmmop23(q,x,y,z,rsq,g,h,rmop,imop)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 23 && defined(FMM_UNROLLED)
 25   call fmmmop24(q,x,y,z,rsq,g,h,rmop,imop)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 24 && defined(FMM_UNROLLED)
 26   call fmmmop25(q,x,y,z,rsq,g,h,rmop,imop)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 25 && defined(FMM_UNROLLED)
 27   call fmmmop26(q,x,y,z,rsq,g,h,rmop,imop)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 26 && defined(FMM_UNROLLED)
 28   call fmmmop27(q,x,y,z,rsq,g,h,rmop,imop)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 27 && defined(FMM_UNROLLED)
 29   call fmmmop28(q,x,y,z,rsq,g,h,rmop,imop)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 28 && defined(FMM_UNROLLED)
 30   call fmmmop29(q,x,y,z,rsq,g,h,rmop,imop)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 29 && defined(FMM_UNROLLED)
 31   call fmmmop30(q,x,y,z,rsq,g,h,rmop,imop)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 30 && defined(FMM_UNROLLED)
 32   call fmmmop31(q,x,y,z,rsq,g,h,rmop,imop)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 31 && defined(FMM_UNROLLED)
 33   call fmmmop32(q,x,y,z,rsq,g,h,rmop,imop)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 32 && defined(FMM_UNROLLED)
 34   call fmmmop33(q,x,y,z,rsq,g,h,rmop,imop)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 33 && defined(FMM_UNROLLED)
 35   call fmmmop34(q,x,y,z,rsq,g,h,rmop,imop)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 34 && defined(FMM_UNROLLED)
 36   call fmmmop35(q,x,y,z,rsq,g,h,rmop,imop)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 35 && defined(FMM_UNROLLED)
 37   call fmmmop36(q,x,y,z,rsq,g,h,rmop,imop)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 36 && defined(FMM_UNROLLED)
 38   call fmmmop37(q,x,y,z,rsq,g,h,rmop,imop)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 37 && defined(FMM_UNROLLED)
 39   call fmmmop38(q,x,y,z,rsq,g,h,rmop,imop)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 38 && defined(FMM_UNROLLED)
 40   call fmmmop39(q,x,y,z,rsq,g,h,rmop,imop)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 39 && defined(FMM_UNROLLED)
 41   call fmmmop40(q,x,y,z,rsq,g,h,rmop,imop)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 40 && defined(FMM_UNROLLED)
 42   call fmmmop41(q,x,y,z,rsq,g,h,rmop,imop)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 41 && defined(FMM_UNROLLED)
 43   call fmmmop42(q,x,y,z,rsq,g,h,rmop,imop)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 42 && defined(FMM_UNROLLED)
 44   call fmmmop43(q,x,y,z,rsq,g,h,rmop,imop)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 43 && defined(FMM_UNROLLED)
 45   call fmmmop44(q,x,y,z,rsq,g,h,rmop,imop)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 44 && defined(FMM_UNROLLED)
 46   call fmmmop45(q,x,y,z,rsq,g,h,rmop,imop)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 45 && defined(FMM_UNROLLED)
 47   call fmmmop46(q,x,y,z,rsq,g,h,rmop,imop)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 46 && defined(FMM_UNROLLED)
 48   call fmmmop47(q,x,y,z,rsq,g,h,rmop,imop)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 47 && defined(FMM_UNROLLED)
 49   call fmmmop48(q,x,y,z,rsq,g,h,rmop,imop)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 48 && defined(FMM_UNROLLED)
 50   call fmmmop49(q,x,y,z,rsq,g,h,rmop,imop)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 49 && defined(FMM_UNROLLED)
 51   call fmmmop50(q,x,y,z,rsq,g,h,rmop,imop)
c
      return
#endif
      end subroutine fmmmopn
