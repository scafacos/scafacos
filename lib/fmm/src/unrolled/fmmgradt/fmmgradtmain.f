#include "../../fmm.h"
c
      subroutine fmmgradt(n,roop,ioop,rmu,imu,rmu00,rmu10,rmu11,imu11)
c
      use fmmkinds
c
      implicit none
c
      real(kind=fmm_real) roop(*),ioop(*),rmu(*),imu(*),rmu00,rmu10,
     .rmu11,imu11
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
      call fmmgrd((n-1),roop,ioop,rmu,imu,rmu00,rmu10,rmu11,imu11)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > -1 && defined(FMM_UNROLLED)
 1    rmu00 = rmu(1)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 0 && defined(FMM_UNROLLED)
 2    call fmmgrd1(roop,ioop,rmu,imu,rmu00,rmu10,rmu11,imu11)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 1 && defined(FMM_UNROLLED)
 3    call fmmgrd2(roop,ioop,rmu,imu,rmu00,rmu10,rmu11,imu11)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 2 && defined(FMM_UNROLLED)
 4    call fmmgrd3(roop,ioop,rmu,imu,rmu00,rmu10,rmu11,imu11)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 3 && defined(FMM_UNROLLED)
 5    call fmmgrd4(roop,ioop,rmu,imu,rmu00,rmu10,rmu11,imu11)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 4 && defined(FMM_UNROLLED)
 6    call fmmgrd5(roop,ioop,rmu,imu,rmu00,rmu10,rmu11,imu11)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 5 && defined(FMM_UNROLLED)
 7    call fmmgrd6(roop,ioop,rmu,imu,rmu00,rmu10,rmu11,imu11)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 6 && defined(FMM_UNROLLED)
 8    call fmmgrd7(roop,ioop,rmu,imu,rmu00,rmu10,rmu11,imu11)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 7 && defined(FMM_UNROLLED)
 9    call fmmgrd8(roop,ioop,rmu,imu,rmu00,rmu10,rmu11,imu11)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 8 && defined(FMM_UNROLLED)
 10   call fmmgrd9(roop,ioop,rmu,imu,rmu00,rmu10,rmu11,imu11)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 9 && defined(FMM_UNROLLED)
 11   call fmmgrd10(roop,ioop,rmu,imu,rmu00,rmu10,rmu11,imu11)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 10 && defined(FMM_UNROLLED)
 12   call fmmgrd11(roop,ioop,rmu,imu,rmu00,rmu10,rmu11,imu11)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 11 && defined(FMM_UNROLLED)
 13   call fmmgrd12(roop,ioop,rmu,imu,rmu00,rmu10,rmu11,imu11)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 12 && defined(FMM_UNROLLED)
 14   call fmmgrd13(roop,ioop,rmu,imu,rmu00,rmu10,rmu11,imu11)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 13 && defined(FMM_UNROLLED)
 15   call fmmgrd14(roop,ioop,rmu,imu,rmu00,rmu10,rmu11,imu11)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 14 && defined(FMM_UNROLLED)
 16   call fmmgrd15(roop,ioop,rmu,imu,rmu00,rmu10,rmu11,imu11)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 15 && defined(FMM_UNROLLED)
 17   call fmmgrd16(roop,ioop,rmu,imu,rmu00,rmu10,rmu11,imu11)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 16 && defined(FMM_UNROLLED)
 18   call fmmgrd17(roop,ioop,rmu,imu,rmu00,rmu10,rmu11,imu11)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 17 && defined(FMM_UNROLLED)
 19   call fmmgrd18(roop,ioop,rmu,imu,rmu00,rmu10,rmu11,imu11)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 18 && defined(FMM_UNROLLED)
 20   call fmmgrd19(roop,ioop,rmu,imu,rmu00,rmu10,rmu11,imu11)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 19 && defined(FMM_UNROLLED)
 21   call fmmgrd20(roop,ioop,rmu,imu,rmu00,rmu10,rmu11,imu11)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 20 && defined(FMM_UNROLLED)
 22   call fmmgrd21(roop,ioop,rmu,imu,rmu00,rmu10,rmu11,imu11)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 21 && defined(FMM_UNROLLED)
 23   call fmmgrd22(roop,ioop,rmu,imu,rmu00,rmu10,rmu11,imu11)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 22 && defined(FMM_UNROLLED)
 24   call fmmgrd23(roop,ioop,rmu,imu,rmu00,rmu10,rmu11,imu11)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 23 && defined(FMM_UNROLLED)
 25   call fmmgrd24(roop,ioop,rmu,imu,rmu00,rmu10,rmu11,imu11)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 24 && defined(FMM_UNROLLED)
 26   call fmmgrd25(roop,ioop,rmu,imu,rmu00,rmu10,rmu11,imu11)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 25 && defined(FMM_UNROLLED)
 27   call fmmgrd26(roop,ioop,rmu,imu,rmu00,rmu10,rmu11,imu11)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 26 && defined(FMM_UNROLLED)
 28   call fmmgrd27(roop,ioop,rmu,imu,rmu00,rmu10,rmu11,imu11)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 27 && defined(FMM_UNROLLED)
 29   call fmmgrd28(roop,ioop,rmu,imu,rmu00,rmu10,rmu11,imu11)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 28 && defined(FMM_UNROLLED)
 30   call fmmgrd29(roop,ioop,rmu,imu,rmu00,rmu10,rmu11,imu11)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 29 && defined(FMM_UNROLLED)
 31   call fmmgrd30(roop,ioop,rmu,imu,rmu00,rmu10,rmu11,imu11)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 30 && defined(FMM_UNROLLED)
 32   call fmmgrd31(roop,ioop,rmu,imu,rmu00,rmu10,rmu11,imu11)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 31 && defined(FMM_UNROLLED)
 33   call fmmgrd32(roop,ioop,rmu,imu,rmu00,rmu10,rmu11,imu11)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 32 && defined(FMM_UNROLLED)
 34   call fmmgrd33(roop,ioop,rmu,imu,rmu00,rmu10,rmu11,imu11)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 33 && defined(FMM_UNROLLED)
 35   call fmmgrd34(roop,ioop,rmu,imu,rmu00,rmu10,rmu11,imu11)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 34 && defined(FMM_UNROLLED)
 36   call fmmgrd35(roop,ioop,rmu,imu,rmu00,rmu10,rmu11,imu11)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 35 && defined(FMM_UNROLLED)
 37   call fmmgrd36(roop,ioop,rmu,imu,rmu00,rmu10,rmu11,imu11)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 36 && defined(FMM_UNROLLED)
 38   call fmmgrd37(roop,ioop,rmu,imu,rmu00,rmu10,rmu11,imu11)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 37 && defined(FMM_UNROLLED)
 39   call fmmgrd38(roop,ioop,rmu,imu,rmu00,rmu10,rmu11,imu11)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 38 && defined(FMM_UNROLLED)
 40   call fmmgrd39(roop,ioop,rmu,imu,rmu00,rmu10,rmu11,imu11)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 39 && defined(FMM_UNROLLED)
 41   call fmmgrd40(roop,ioop,rmu,imu,rmu00,rmu10,rmu11,imu11)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 40 && defined(FMM_UNROLLED)
 42   call fmmgrd41(roop,ioop,rmu,imu,rmu00,rmu10,rmu11,imu11)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 41 && defined(FMM_UNROLLED)
 43   call fmmgrd42(roop,ioop,rmu,imu,rmu00,rmu10,rmu11,imu11)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 42 && defined(FMM_UNROLLED)
 44   call fmmgrd43(roop,ioop,rmu,imu,rmu00,rmu10,rmu11,imu11)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 43 && defined(FMM_UNROLLED)
 45   call fmmgrd44(roop,ioop,rmu,imu,rmu00,rmu10,rmu11,imu11)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 44 && defined(FMM_UNROLLED)
 46   call fmmgrd45(roop,ioop,rmu,imu,rmu00,rmu10,rmu11,imu11)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 45 && defined(FMM_UNROLLED)
 47   call fmmgrd46(roop,ioop,rmu,imu,rmu00,rmu10,rmu11,imu11)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 46 && defined(FMM_UNROLLED)
 48   call fmmgrd47(roop,ioop,rmu,imu,rmu00,rmu10,rmu11,imu11)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 47 && defined(FMM_UNROLLED)
 49   call fmmgrd48(roop,ioop,rmu,imu,rmu00,rmu10,rmu11,imu11)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 48 && defined(FMM_UNROLLED)
 50   call fmmgrd49(roop,ioop,rmu,imu,rmu00,rmu10,rmu11,imu11)
c
      return
#endif
c
#if FMM_MAXNMULTIPOLES > 49 && defined(FMM_UNROLLED)
 51   call fmmgrd50(roop,ioop,rmu,imu,rmu00,rmu10,rmu11,imu11)
c
      return
#endif
      end subroutine fmmgradt
