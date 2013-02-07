#include "../../fmm.h"

subroutine get_fmmec_1d(pole,fmmec_ptr)
  use fmmkinds
  implicit none
  type(c_ptr),target :: fmmec_ptr 
  integer(kind=fmm_integer) :: pole

  interface

#if FMM_MAXNMULTIPOLES > -1
    function get_fmmec_1d_p00() bind(c)
      use, intrinsic :: iso_c_binding
      type(c_ptr) :: get_fmmec_1d_p00
    end function get_fmmec_1d_p00
#endif

#if FMM_MAXNMULTIPOLES > 0
    function get_fmmec_1d_p01() bind(c)
      use, intrinsic :: iso_c_binding
      type(c_ptr) :: get_fmmec_1d_p01
    end function get_fmmec_1d_p01
#endif

#if FMM_MAXNMULTIPOLES > 1
    function get_fmmec_1d_p02() bind(c)
      use, intrinsic :: iso_c_binding
      type(c_ptr) :: get_fmmec_1d_p02
    end function get_fmmec_1d_p02
#endif

#if FMM_MAXNMULTIPOLES > 2
    function get_fmmec_1d_p03() bind(c)
      use, intrinsic :: iso_c_binding
      type(c_ptr) :: get_fmmec_1d_p03
    end function get_fmmec_1d_p03
#endif

#if FMM_MAXNMULTIPOLES > 3
    function get_fmmec_1d_p04() bind(c)
      use, intrinsic :: iso_c_binding
      type(c_ptr) :: get_fmmec_1d_p04
    end function get_fmmec_1d_p04
#endif

#if FMM_MAXNMULTIPOLES > 4
    function get_fmmec_1d_p05() bind(c)
      use, intrinsic :: iso_c_binding
      type(c_ptr) :: get_fmmec_1d_p05
    end function get_fmmec_1d_p05
#endif

#if FMM_MAXNMULTIPOLES > 5
    function get_fmmec_1d_p06() bind(c)
      use, intrinsic :: iso_c_binding
      type(c_ptr) :: get_fmmec_1d_p06
    end function get_fmmec_1d_p06
#endif

#if FMM_MAXNMULTIPOLES > 6
    function get_fmmec_1d_p07() bind(c)
      use, intrinsic :: iso_c_binding
      type(c_ptr) :: get_fmmec_1d_p07
    end function get_fmmec_1d_p07
#endif

#if FMM_MAXNMULTIPOLES > 7
    function get_fmmec_1d_p08() bind(c)
      use, intrinsic :: iso_c_binding
      type(c_ptr) :: get_fmmec_1d_p08
    end function get_fmmec_1d_p08
#endif

#if FMM_MAXNMULTIPOLES > 8
    function get_fmmec_1d_p09() bind(c)
      use, intrinsic :: iso_c_binding
      type(c_ptr) :: get_fmmec_1d_p09
    end function get_fmmec_1d_p09
#endif

#if FMM_MAXNMULTIPOLES > 9
    function get_fmmec_1d_p10() bind(c)
      use, intrinsic :: iso_c_binding
      type(c_ptr) :: get_fmmec_1d_p10
    end function get_fmmec_1d_p10
#endif

#if FMM_MAXNMULTIPOLES > 10
    function get_fmmec_1d_p11() bind(c)
      use, intrinsic :: iso_c_binding
      type(c_ptr) :: get_fmmec_1d_p11
    end function get_fmmec_1d_p11
#endif

#if FMM_MAXNMULTIPOLES > 11
    function get_fmmec_1d_p12() bind(c)
      use, intrinsic :: iso_c_binding
      type(c_ptr) :: get_fmmec_1d_p12
    end function get_fmmec_1d_p12
#endif

#if FMM_MAXNMULTIPOLES > 12
    function get_fmmec_1d_p13() bind(c)
      use, intrinsic :: iso_c_binding
      type(c_ptr) :: get_fmmec_1d_p13
    end function get_fmmec_1d_p13
#endif

#if FMM_MAXNMULTIPOLES > 13
    function get_fmmec_1d_p14() bind(c)
      use, intrinsic :: iso_c_binding
      type(c_ptr) :: get_fmmec_1d_p14
    end function get_fmmec_1d_p14
#endif

#if FMM_MAXNMULTIPOLES > 14
    function get_fmmec_1d_p15() bind(c)
      use, intrinsic :: iso_c_binding
      type(c_ptr) :: get_fmmec_1d_p15
    end function get_fmmec_1d_p15
#endif

#if FMM_MAXNMULTIPOLES > 15
    function get_fmmec_1d_p16() bind(c)
      use, intrinsic :: iso_c_binding
      type(c_ptr) :: get_fmmec_1d_p16
    end function get_fmmec_1d_p16
#endif

#if FMM_MAXNMULTIPOLES > 16
    function get_fmmec_1d_p17() bind(c)
      use, intrinsic :: iso_c_binding
      type(c_ptr) :: get_fmmec_1d_p17
    end function get_fmmec_1d_p17
#endif

#if FMM_MAXNMULTIPOLES > 17
    function get_fmmec_1d_p18() bind(c)
      use, intrinsic :: iso_c_binding
      type(c_ptr) :: get_fmmec_1d_p18
    end function get_fmmec_1d_p18
#endif

#if FMM_MAXNMULTIPOLES > 18
    function get_fmmec_1d_p19() bind(c)
      use, intrinsic :: iso_c_binding
      type(c_ptr) :: get_fmmec_1d_p19
    end function get_fmmec_1d_p19
#endif

#if FMM_MAXNMULTIPOLES > 19
    function get_fmmec_1d_p20() bind(c)
      use, intrinsic :: iso_c_binding
      type(c_ptr) :: get_fmmec_1d_p20
    end function get_fmmec_1d_p20
#endif

#if FMM_MAXNMULTIPOLES > 20
    function get_fmmec_1d_p21() bind(c)
      use, intrinsic :: iso_c_binding
      type(c_ptr) :: get_fmmec_1d_p21
    end function get_fmmec_1d_p21
#endif

#if FMM_MAXNMULTIPOLES > 21
    function get_fmmec_1d_p22() bind(c)
      use, intrinsic :: iso_c_binding
      type(c_ptr) :: get_fmmec_1d_p22
    end function get_fmmec_1d_p22
#endif

#if FMM_MAXNMULTIPOLES > 22
    function get_fmmec_1d_p23() bind(c)
      use, intrinsic :: iso_c_binding
      type(c_ptr) :: get_fmmec_1d_p23
    end function get_fmmec_1d_p23
#endif

#if FMM_MAXNMULTIPOLES > 23
    function get_fmmec_1d_p24() bind(c)
      use, intrinsic :: iso_c_binding
      type(c_ptr) :: get_fmmec_1d_p24
    end function get_fmmec_1d_p24
#endif

#if FMM_MAXNMULTIPOLES > 24
    function get_fmmec_1d_p25() bind(c)
      use, intrinsic :: iso_c_binding
      type(c_ptr) :: get_fmmec_1d_p25
    end function get_fmmec_1d_p25
#endif

#if FMM_MAXNMULTIPOLES > 25
    function get_fmmec_1d_p26() bind(c)
      use, intrinsic :: iso_c_binding
      type(c_ptr) :: get_fmmec_1d_p26
    end function get_fmmec_1d_p26
#endif

#if FMM_MAXNMULTIPOLES > 26
    function get_fmmec_1d_p27() bind(c)
      use, intrinsic :: iso_c_binding
      type(c_ptr) :: get_fmmec_1d_p27
    end function get_fmmec_1d_p27
#endif

#if FMM_MAXNMULTIPOLES > 27
    function get_fmmec_1d_p28() bind(c)
      use, intrinsic :: iso_c_binding
      type(c_ptr) :: get_fmmec_1d_p28
    end function get_fmmec_1d_p28
#endif

#if FMM_MAXNMULTIPOLES > 28
    function get_fmmec_1d_p29() bind(c)
      use, intrinsic :: iso_c_binding
      type(c_ptr) :: get_fmmec_1d_p29
    end function get_fmmec_1d_p29
#endif

#if FMM_MAXNMULTIPOLES > 29
    function get_fmmec_1d_p30() bind(c)
      use, intrinsic :: iso_c_binding
      type(c_ptr) :: get_fmmec_1d_p30
    end function get_fmmec_1d_p30
#endif

#if FMM_MAXNMULTIPOLES > 30
    function get_fmmec_1d_p31() bind(c)
      use, intrinsic :: iso_c_binding
      type(c_ptr) :: get_fmmec_1d_p31
    end function get_fmmec_1d_p31
#endif

#if FMM_MAXNMULTIPOLES > 31
    function get_fmmec_1d_p32() bind(c)
      use, intrinsic :: iso_c_binding
      type(c_ptr) :: get_fmmec_1d_p32
    end function get_fmmec_1d_p32
#endif

#if FMM_MAXNMULTIPOLES > 32
    function get_fmmec_1d_p33() bind(c)
      use, intrinsic :: iso_c_binding
      type(c_ptr) :: get_fmmec_1d_p33
    end function get_fmmec_1d_p33
#endif

#if FMM_MAXNMULTIPOLES > 33
    function get_fmmec_1d_p34() bind(c)
      use, intrinsic :: iso_c_binding
      type(c_ptr) :: get_fmmec_1d_p34
    end function get_fmmec_1d_p34
#endif

#if FMM_MAXNMULTIPOLES > 34
    function get_fmmec_1d_p35() bind(c)
      use, intrinsic :: iso_c_binding
      type(c_ptr) :: get_fmmec_1d_p35
    end function get_fmmec_1d_p35
#endif

#if FMM_MAXNMULTIPOLES > 35
    function get_fmmec_1d_p36() bind(c)
      use, intrinsic :: iso_c_binding
      type(c_ptr) :: get_fmmec_1d_p36
    end function get_fmmec_1d_p36
#endif

#if FMM_MAXNMULTIPOLES > 36
    function get_fmmec_1d_p37() bind(c)
      use, intrinsic :: iso_c_binding
      type(c_ptr) :: get_fmmec_1d_p37
    end function get_fmmec_1d_p37
#endif

#if FMM_MAXNMULTIPOLES > 37
    function get_fmmec_1d_p38() bind(c)
      use, intrinsic :: iso_c_binding
      type(c_ptr) :: get_fmmec_1d_p38
    end function get_fmmec_1d_p38
#endif

#if FMM_MAXNMULTIPOLES > 38
    function get_fmmec_1d_p39() bind(c)
      use, intrinsic :: iso_c_binding
      type(c_ptr) :: get_fmmec_1d_p39
    end function get_fmmec_1d_p39
#endif

#if FMM_MAXNMULTIPOLES > 39
    function get_fmmec_1d_p40() bind(c)
      use, intrinsic :: iso_c_binding
      type(c_ptr) :: get_fmmec_1d_p40
    end function get_fmmec_1d_p40
#endif

#if FMM_MAXNMULTIPOLES > 40
    function get_fmmec_1d_p41() bind(c)
      use, intrinsic :: iso_c_binding
      type(c_ptr) :: get_fmmec_1d_p41
    end function get_fmmec_1d_p41
#endif

#if FMM_MAXNMULTIPOLES > 41
    function get_fmmec_1d_p42() bind(c)
      use, intrinsic :: iso_c_binding
      type(c_ptr) :: get_fmmec_1d_p42
    end function get_fmmec_1d_p42
#endif

#if FMM_MAXNMULTIPOLES > 42
    function get_fmmec_1d_p43() bind(c)
      use, intrinsic :: iso_c_binding
      type(c_ptr) :: get_fmmec_1d_p43
    end function get_fmmec_1d_p43
#endif

#if FMM_MAXNMULTIPOLES > 43
    function get_fmmec_1d_p44() bind(c)
      use, intrinsic :: iso_c_binding
      type(c_ptr) :: get_fmmec_1d_p44
    end function get_fmmec_1d_p44
#endif

#if FMM_MAXNMULTIPOLES > 44
    function get_fmmec_1d_p45() bind(c)
      use, intrinsic :: iso_c_binding
      type(c_ptr) :: get_fmmec_1d_p45
    end function get_fmmec_1d_p45
#endif

#if FMM_MAXNMULTIPOLES > 45
    function get_fmmec_1d_p46() bind(c)
      use, intrinsic :: iso_c_binding
      type(c_ptr) :: get_fmmec_1d_p46
    end function get_fmmec_1d_p46
#endif

#if FMM_MAXNMULTIPOLES > 46
    function get_fmmec_1d_p47() bind(c)
      use, intrinsic :: iso_c_binding
      type(c_ptr) :: get_fmmec_1d_p47
    end function get_fmmec_1d_p47
#endif

#if FMM_MAXNMULTIPOLES > 47
    function get_fmmec_1d_p48() bind(c)
      use, intrinsic :: iso_c_binding
      type(c_ptr) :: get_fmmec_1d_p48
    end function get_fmmec_1d_p48
#endif

#if FMM_MAXNMULTIPOLES > 48
    function get_fmmec_1d_p49() bind(c)
      use, intrinsic :: iso_c_binding
      type(c_ptr) :: get_fmmec_1d_p49
    end function get_fmmec_1d_p49
#endif

#if FMM_MAXNMULTIPOLES > 49
    function get_fmmec_1d_p50() bind(c)
      use, intrinsic :: iso_c_binding
      type(c_ptr) :: get_fmmec_1d_p50
    end function get_fmmec_1d_p50
#endif

  end interface

    select case (pole+1)
#if FMM_MAXNMULTIPOLES > -1 
    case ( 1)
      fmmec_ptr = get_fmmec_1d_p00()
#endif

#if FMM_MAXNMULTIPOLES > 0 
    case ( 2)
      fmmec_ptr = get_fmmec_1d_p01()
#endif

#if FMM_MAXNMULTIPOLES > 1 
    case ( 3)
      fmmec_ptr = get_fmmec_1d_p02()
#endif

#if FMM_MAXNMULTIPOLES > 2 
    case ( 4)
      fmmec_ptr = get_fmmec_1d_p03()
#endif

#if FMM_MAXNMULTIPOLES > 3 
    case ( 5)
      fmmec_ptr = get_fmmec_1d_p04()
#endif

#if FMM_MAXNMULTIPOLES > 4 
    case ( 6)
      fmmec_ptr = get_fmmec_1d_p05()
#endif

#if FMM_MAXNMULTIPOLES > 5 
    case ( 7)
      fmmec_ptr = get_fmmec_1d_p06()
#endif

#if FMM_MAXNMULTIPOLES > 6 
    case ( 8)
      fmmec_ptr = get_fmmec_1d_p07()
#endif

#if FMM_MAXNMULTIPOLES > 7 
    case ( 9)
      fmmec_ptr = get_fmmec_1d_p08()
#endif

#if FMM_MAXNMULTIPOLES > 8 
    case ( 10)
      fmmec_ptr = get_fmmec_1d_p09()
#endif

#if FMM_MAXNMULTIPOLES > 9 
    case ( 11)
      fmmec_ptr = get_fmmec_1d_p10()
#endif

#if FMM_MAXNMULTIPOLES > 10 
    case ( 12)
      fmmec_ptr = get_fmmec_1d_p11()
#endif

#if FMM_MAXNMULTIPOLES > 11 
    case ( 13)
      fmmec_ptr = get_fmmec_1d_p12()
#endif

#if FMM_MAXNMULTIPOLES > 12 
    case ( 14)
      fmmec_ptr = get_fmmec_1d_p13()
#endif

#if FMM_MAXNMULTIPOLES > 13 
    case ( 15)
      fmmec_ptr = get_fmmec_1d_p14()
#endif

#if FMM_MAXNMULTIPOLES > 14 
    case ( 16)
      fmmec_ptr = get_fmmec_1d_p15()
#endif

#if FMM_MAXNMULTIPOLES > 15 
    case ( 17)
      fmmec_ptr = get_fmmec_1d_p16()
#endif

#if FMM_MAXNMULTIPOLES > 16 
    case ( 18)
      fmmec_ptr = get_fmmec_1d_p17()
#endif

#if FMM_MAXNMULTIPOLES > 17 
    case ( 19)
      fmmec_ptr = get_fmmec_1d_p18()
#endif

#if FMM_MAXNMULTIPOLES > 18 
    case ( 20)
      fmmec_ptr = get_fmmec_1d_p19()
#endif

#if FMM_MAXNMULTIPOLES > 19 
    case ( 21)
      fmmec_ptr = get_fmmec_1d_p20()
#endif

#if FMM_MAXNMULTIPOLES > 20 
    case ( 22)
      fmmec_ptr = get_fmmec_1d_p21()
#endif

#if FMM_MAXNMULTIPOLES > 21 
    case ( 23)
      fmmec_ptr = get_fmmec_1d_p22()
#endif

#if FMM_MAXNMULTIPOLES > 22 
    case ( 24)
      fmmec_ptr = get_fmmec_1d_p23()
#endif

#if FMM_MAXNMULTIPOLES > 23 
    case ( 25)
      fmmec_ptr = get_fmmec_1d_p24()
#endif

#if FMM_MAXNMULTIPOLES > 24 
    case ( 26)
      fmmec_ptr = get_fmmec_1d_p25()
#endif

#if FMM_MAXNMULTIPOLES > 25 
    case ( 27)
      fmmec_ptr = get_fmmec_1d_p26()
#endif

#if FMM_MAXNMULTIPOLES > 26 
    case ( 28)
      fmmec_ptr = get_fmmec_1d_p27()
#endif

#if FMM_MAXNMULTIPOLES > 27 
    case ( 29)
      fmmec_ptr = get_fmmec_1d_p28()
#endif

#if FMM_MAXNMULTIPOLES > 28 
    case ( 30)
      fmmec_ptr = get_fmmec_1d_p29()
#endif

#if FMM_MAXNMULTIPOLES > 29 
    case ( 31)
      fmmec_ptr = get_fmmec_1d_p30()
#endif

#if FMM_MAXNMULTIPOLES > 30 
    case ( 32)
      fmmec_ptr = get_fmmec_1d_p31()
#endif

#if FMM_MAXNMULTIPOLES > 31 
    case ( 33)
      fmmec_ptr = get_fmmec_1d_p32()
#endif

#if FMM_MAXNMULTIPOLES > 32 
    case ( 34)
      fmmec_ptr = get_fmmec_1d_p33()
#endif

#if FMM_MAXNMULTIPOLES > 33 
    case ( 35)
      fmmec_ptr = get_fmmec_1d_p34()
#endif

#if FMM_MAXNMULTIPOLES > 34 
    case ( 36)
      fmmec_ptr = get_fmmec_1d_p35()
#endif

#if FMM_MAXNMULTIPOLES > 35 
    case ( 37)
      fmmec_ptr = get_fmmec_1d_p36()
#endif

#if FMM_MAXNMULTIPOLES > 36 
    case ( 38)
      fmmec_ptr = get_fmmec_1d_p37()
#endif

#if FMM_MAXNMULTIPOLES > 37 
    case ( 39)
      fmmec_ptr = get_fmmec_1d_p38()
#endif

#if FMM_MAXNMULTIPOLES > 38 
    case ( 40)
      fmmec_ptr = get_fmmec_1d_p39()
#endif

#if FMM_MAXNMULTIPOLES > 39 
    case ( 41)
      fmmec_ptr = get_fmmec_1d_p40()
#endif

#if FMM_MAXNMULTIPOLES > 40 
    case ( 42)
      fmmec_ptr = get_fmmec_1d_p41()
#endif

#if FMM_MAXNMULTIPOLES > 41 
    case ( 43)
      fmmec_ptr = get_fmmec_1d_p42()
#endif

#if FMM_MAXNMULTIPOLES > 42 
    case ( 44)
      fmmec_ptr = get_fmmec_1d_p43()
#endif

#if FMM_MAXNMULTIPOLES > 43 
    case ( 45)
      fmmec_ptr = get_fmmec_1d_p44()
#endif

#if FMM_MAXNMULTIPOLES > 44 
    case ( 46)
      fmmec_ptr = get_fmmec_1d_p45()
#endif

#if FMM_MAXNMULTIPOLES > 45 
    case ( 47)
      fmmec_ptr = get_fmmec_1d_p46()
#endif

#if FMM_MAXNMULTIPOLES > 46 
    case ( 48)
      fmmec_ptr = get_fmmec_1d_p47()
#endif

#if FMM_MAXNMULTIPOLES > 47 
    case ( 49)
      fmmec_ptr = get_fmmec_1d_p48()
#endif

#if FMM_MAXNMULTIPOLES > 48 
    case ( 50)
      fmmec_ptr = get_fmmec_1d_p49()
#endif

#if FMM_MAXNMULTIPOLES > 49 
    case ( 51)
      fmmec_ptr = get_fmmec_1d_p50()
#endif

    case default
!   TODO implement error call
!      call bummer()
    end select 
end subroutine get_fmmec_1d
