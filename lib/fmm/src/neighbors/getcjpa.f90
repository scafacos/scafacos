module getneighbors_vars
use fmmkinds
integer(kind=fmm_integer),pointer :: gn_nbits
integer(kind=fmm_integer), dimension(:),pointer :: gn_bitpos,gn_mbitpos

contains
subroutine getneighborsinit(nbits,bitpos,mbitpos)
	implicit none
	integer(kind=fmm_integer),target :: nbits
	integer(kind=fmm_integer), dimension(0:),target :: bitpos,mbitpos
	
	gn_nbits => nbits
	gn_bitpos => bitpos
	gn_mbitpos => mbitpos

end subroutine getneighborsinit

function fvedge(myindex)
      implicit none
      integer(kind=fmm_integer) :: fvedge
      integer(kind=fmm_integer) :: myindex
      integer(kind=fmm_integer) :: i
      fvedge = 0
      do i = 0,gn_nbits
        if(btest(myindex,i)) then
          fvedge = i + 1
          return
        end if
      end do
end function fvedge
      
      function fint3x(m)

      implicit none

      integer(kind=fmm_integer) fint3x,m,i,j,k

      i = 0
      j = -2
      k = -1
 1    k = k+1
      if(iand(m,gn_mbitpos(k)).gt.0) then
         j = j+2
         i = ior(i,iand(ishft(iand(m,gn_bitpos(k)),j),gn_mbitpos(j)))
         go to 1
      endif
      fint3x = i
      return
      end function fint3x

      function fint3y(m)
      implicit none

      integer(kind=fmm_integer) fint3y,m,i,j,k

      i = 0
      j = -1
      k = -1
 1    k = k+1
      if(iand(m,gn_mbitpos(k)).gt.0) then
         j = j+2
         i = ior(i,iand(ishft(iand(m,gn_bitpos(k)),j),gn_mbitpos(j)))
         go to 1
      endif
      fint3y = i
      return
      end function fint3y

      function fint3z(m)
      implicit none

      integer(kind=fmm_integer) fint3z,m,i,j,k

      i = 0
      j = 0
      k = -1
 1    k = k+1
      if(iand(m,gn_mbitpos(k)).gt.0) then
         j = j+2
         i = ior(i,iand(ishft(iand(m,gn_bitpos(k)),j),gn_mbitpos(j)))
         go to 1
      endif
      fint3z = i
      return
      end function fint3z

end module getneighbors_vars

subroutine getcjpinit(cjpsize)
       use fmmkinds
	! getcjpinit -> getcasejumpinit cjpsize -> casejumpsize
	implicit none
	integer(kind=fmm_integer) :: cjpsize ! should be 8,64,512,4096

	cjpsize = 512
end subroutine getcjpinit



subroutine getcjp(casejump)
       use fmmkinds
	implicit none
	integer(kind=fmm_integer), dimension(0:*) :: casejump ! must have same size as casejump(0)+1 starting from zero (0..8) (0..64) (0..512) (0..4096) casejump(0) stores casejump size

	select case(casejump(0))
		case(7)
!			call getcjp8(casejump)
			return
		case(63)
!			call getcjp64(casejump)
			return
		case(511)
			call getcjp512(casejump)
			return
		case(4095)
!			call getcjp4096(casejump)
			return
		case default
			stop "error: getcjp"
			return
	end select
end subroutine getcjp

subroutine getcjp512(casejump)
      use fmmkinds
      ! getcjp -> getcasejump
      implicit none
      integer(kind=fmm_integer), dimension(0:512) :: casejump

      casejump(  1) =  73
      casejump(  2) =  49
      casejump(  3) =  50
      casejump(  4) =   1
      casejump(  5) =  51
      casejump(  6) =   2
      casejump(  7) =   3
      casejump(  8) =   4
      casejump(  9) =  52
      casejump( 10) =  49
      casejump( 11) =   5
      casejump( 12) =   6
      casejump( 13) =   7
      casejump( 14) =   8
      casejump( 15) =   9
      casejump( 16) =  10
      casejump( 17) =  53
      casejump( 18) =  11
      casejump( 19) =  50
      casejump( 20) =   1
      casejump( 21) =  12
      casejump( 22) =  13
      casejump( 23) =  14
      casejump( 24) =  15
      casejump( 25) =  16
      casejump( 26) =  17
      casejump( 27) =   5
      casejump( 28) =   1
      casejump( 29) =  18
      casejump( 30) =  19
      casejump( 31) =  20
      casejump( 32) =  21
      casejump( 33) =  54
      casejump( 34) =  22
      casejump( 35) =  23
      casejump( 36) =   1
      casejump( 37) =  51
      casejump( 38) =   2
      casejump( 39) =   3
      casejump( 40) =   4
      casejump( 41) =  24
      casejump( 42) =  25
      casejump( 43) =   5
      casejump( 44) =  26
      casejump( 45) =   7
      casejump( 46) =   2
      casejump( 47) =   9
      casejump( 48) =  27
      casejump( 49) =  28
      casejump( 50) =  11
      casejump( 51) =  29
      casejump( 52) =  30
      casejump( 53) =  12
      casejump( 54) =  13
      casejump( 55) =   3
      casejump( 56) =   4
      casejump( 57) =  16
      casejump( 58) =  31
      casejump( 59) =  32
      casejump( 60) =  33
      casejump( 61) =  18
      casejump( 62) =  34
      casejump( 63) =   9
      casejump( 64) =   4
      casejump( 65) =  52
      casejump( 66) =  49
      casejump( 67) =  35
      casejump( 68) =   1
      casejump( 69) =  36
      casejump( 70) =   2
      casejump( 71) =   3
      casejump( 72) =   4
      casejump( 73) =  52
      casejump( 74) =  74
      casejump( 75) =   5
      casejump( 76) =  55
      casejump( 77) =   7
      casejump( 78) =  56
      casejump( 79) =   9
      casejump( 80) =  10
      casejump( 81) =  37
      casejump( 82) =  11
      casejump( 83) =   5
      casejump( 84) =   1
      casejump( 85) =  12
      casejump( 86) =  13
      casejump( 87) =  38
      casejump( 88) =  15
      casejump( 89) =  16
      casejump( 90) =  57
      casejump( 91) =   5
      casejump( 92) =  55
      casejump( 93) =  18
      casejump( 94) =  19
      casejump( 95) =  20
      casejump( 96) =  39
      casejump( 97) =  40
      casejump( 98) =  22
      casejump( 99) =  23
      casejump(100) =   1
      casejump(101) =   7
      casejump(102) =   2
      casejump(103) =  41
      casejump(104) =   4
      casejump(105) =  24
      casejump(106) =  58
      casejump(107) =   5
      casejump(108) =  26
      casejump(109) =   7
      casejump(110) =  56
      casejump(111) =   9
      casejump(112) =  10
      casejump(113) =  28
      casejump(114) =  11
      casejump(115) =  42
      casejump(116) =  30
      casejump(117) =  43
      casejump(118) =  13
      casejump(119) =   9
      casejump(120) =   4
      casejump(121) =  16
      casejump(122) =  31
      casejump(123) =  32
      casejump(124) =  44
      casejump(125) =  18
      casejump(126) =  19
      casejump(127) =   9
      casejump(128) =  10
      casejump(129) =  53
      casejump(130) =  11
      casejump(131) =  50
      casejump(132) =   1
      casejump(133) =  45
      casejump(134) =   2
      casejump(135) =   3
      casejump(136) =   4
      casejump(137) =  16
      casejump(138) =  11
      casejump(139) =   5
      casejump(140) =   6
      casejump(141) =   7
      casejump(142) =   8
      casejump(143) =   9
      casejump(144) =  10
      casejump(145) =  53
      casejump(146) =  11
      casejump(147) =  75
      casejump(148) =  59
      casejump(149) =  12
      casejump(150) =  13
      casejump(151) =  60
      casejump(152) =  15
      casejump(153) =  16
      casejump(154) =  17
      casejump(155) =  61
      casejump(156) =  59
      casejump(157) =  18
      casejump(158) =  19
      casejump(159) =  20
      casejump(160) =  21
      casejump(161) =  46
      casejump(162) =  22
      casejump(163) =  23
      casejump(164) =   1
      casejump(165) =  12
      casejump(166) =  13
      casejump(167) =   3
      casejump(168) =   4
      casejump(169) =  24
      casejump(170) =  25
      casejump(171) =   5
      casejump(172) =  26
      casejump(173) =  18
      casejump(174) =  13
      casejump(175) =   9
      casejump(176) =  27
      casejump(177) =  28
      casejump(178) =  11
      casejump(179) =  62
      casejump(180) =  30
      casejump(181) =  12
      casejump(182) =  13
      casejump(183) =  60
      casejump(184) =  15
      casejump(185) =  16
      casejump(186) =  31
      casejump(187) =  32
      casejump(188) =  33
      casejump(189) =  18
      casejump(190) =  34
      casejump(191) =  20
      casejump(192) =  15
      casejump(193) =  16
      casejump(194) =  11
      casejump(195) =  35
      casejump(196) =   1
      casejump(197) =  36
      casejump(198) =   2
      casejump(199) =   3
      casejump(200) =   4
      casejump(201) =  16
      casejump(202) =  57
      casejump(203) =   5
      casejump(204) =  55
      casejump(205) =   7
      casejump(206) =  47
      casejump(207) =   9
      casejump(208) =  10
      casejump(209) =  37
      casejump(210) =  11
      casejump(211) =  61
      casejump(212) =  59
      casejump(213) =  12
      casejump(214) =  13
      casejump(215) =  38
      casejump(216) =  15
      casejump(217) =  16
      casejump(218) =  57
      casejump(219) =  61
      casejump(220) =  76
      casejump(221) =  18
      casejump(222) =  19
      casejump(223) =  20
      casejump(224) =  63
      casejump(225) =  40
      casejump(226) =  22
      casejump(227) =  23
      casejump(228) =   1
      casejump(229) =  18
      casejump(230) =  13
      casejump(231) =  41
      casejump(232) =   4
      casejump(233) =  24
      casejump(234) =  48
      casejump(235) =   5
      casejump(236) =  26
      casejump(237) =  18
      casejump(238) =  19
      casejump(239) =   9
      casejump(240) =  10
      casejump(241) =  28
      casejump(242) =  11
      casejump(243) =  42
      casejump(244) =  30
      casejump(245) =  43
      casejump(246) =  13
      casejump(247) =  20
      casejump(248) =  15
      casejump(249) =  16
      casejump(250) =  31
      casejump(251) =  32
      casejump(252) =  64
      casejump(253) =  18
      casejump(254) =  19
      casejump(255) =  20
      casejump(256) =  63
      casejump(257) =  54
      casejump(258) =  22
      casejump(259) =  23
      casejump(260) =   1
      casejump(261) =  51
      casejump(262) =   2
      casejump(263) =   3
      casejump(264) =   4
      casejump(265) =  24
      casejump(266) =  22
      casejump(267) =   5
      casejump(268) =   6
      casejump(269) =   7
      casejump(270) =   8
      casejump(271) =   9
      casejump(272) =  10
      casejump(273) =  28
      casejump(274) =  11
      casejump(275) =  23
      casejump(276) =   1
      casejump(277) =  12
      casejump(278) =  13
      casejump(279) =  14
      casejump(280) =  15
      casejump(281) =  16
      casejump(282) =  17
      casejump(283) =   5
      casejump(284) =   1
      casejump(285) =  18
      casejump(286) =  19
      casejump(287) =  20
      casejump(288) =  21
      casejump(289) =  54
      casejump(290) =  22
      casejump(291) =  23
      casejump(292) =   1
      casejump(293) =  77
      casejump(294) =  65
      casejump(295) =  66
      casejump(296) =   4
      casejump(297) =  24
      casejump(298) =  25
      casejump(299) =   5
      casejump(300) =  26
      casejump(301) =  67
      casejump(302) =  65
      casejump(303) =   9
      casejump(304) =  27
      casejump(305) =  28
      casejump(306) =  11
      casejump(307) =  29
      casejump(308) =  30
      casejump(309) =  68
      casejump(310) =  13
      casejump(311) =  66
      casejump(312) =   4
      casejump(313) =  16
      casejump(314) =  31
      casejump(315) =  32
      casejump(316) =  33
      casejump(317) =  18
      casejump(318) =  34
      casejump(319) =   9
      casejump(320) =   4
      casejump(321) =  24
      casejump(322) =  22
      casejump(323) =  35
      casejump(324) =   1
      casejump(325) =  36
      casejump(326) =   2
      casejump(327) =   3
      casejump(328) =   4
      casejump(329) =  24
      casejump(330) =  58
      casejump(331) =   5
      casejump(332) =  26
      casejump(333) =   7
      casejump(334) =  56
      casejump(335) =   9
      casejump(336) =  10
      casejump(337) =  37
      casejump(338) =  11
      casejump(339) =   5
      casejump(340) =   1
      casejump(341) =  12
      casejump(342) =  13
      casejump(343) =  38
      casejump(344) =  15
      casejump(345) =  16
      casejump(346) =  31
      casejump(347) =   5
      casejump(348) =  26
      casejump(349) =  18
      casejump(350) =  19
      casejump(351) =  20
      casejump(352) =  39
      casejump(353) =  40
      casejump(354) =  22
      casejump(355) =  23
      casejump(356) =   1
      casejump(357) =  67
      casejump(358) =  65
      casejump(359) =  41
      casejump(360) =   4
      casejump(361) =  24
      casejump(362) =  58
      casejump(363) =   5
      casejump(364) =  26
      casejump(365) =  67
      casejump(366) =  78
      casejump(367) =   9
      casejump(368) =  69
      casejump(369) =  28
      casejump(370) =  11
      casejump(371) =  42
      casejump(372) =  30
      casejump(373) =  43
      casejump(374) =  13
      casejump(375) =   9
      casejump(376) =   4
      casejump(377) =  16
      casejump(378) =  31
      casejump(379) =  32
      casejump(380) =  44
      casejump(381) =  18
      casejump(382) =  70
      casejump(383) =   9
      casejump(384) =  69
      casejump(385) =  28
      casejump(386) =  11
      casejump(387) =  23
      casejump(388) =   1
      casejump(389) =  45
      casejump(390) =   2
      casejump(391) =   3
      casejump(392) =   4
      casejump(393) =  16
      casejump(394) =  11
      casejump(395) =   5
      casejump(396) =   6
      casejump(397) =   7
      casejump(398) =   8
      casejump(399) =   9
      casejump(400) =  10
      casejump(401) =  28
      casejump(402) =  11
      casejump(403) =  62
      casejump(404) =  30
      casejump(405) =  12
      casejump(406) =  13
      casejump(407) =  60
      casejump(408) =  15
      casejump(409) =  16
      casejump(410) =  17
      casejump(411) =  32
      casejump(412) =  30
      casejump(413) =  18
      casejump(414) =  19
      casejump(415) =  20
      casejump(416) =  21
      casejump(417) =  46
      casejump(418) =  22
      casejump(419) =  23
      casejump(420) =   1
      casejump(421) =  68
      casejump(422) =  13
      casejump(423) =  66
      casejump(424) =   4
      casejump(425) =  24
      casejump(426) =  25
      casejump(427) =   5
      casejump(428) =  26
      casejump(429) =  18
      casejump(430) =  13
      casejump(431) =   9
      casejump(432) =  27
      casejump(433) =  28
      casejump(434) =  11
      casejump(435) =  62
      casejump(436) =  30
      casejump(437) =  68
      casejump(438) =  13
      casejump(439) =  79
      casejump(440) =  71
      casejump(441) =  16
      casejump(442) =  31
      casejump(443) =  32
      casejump(444) =  33
      casejump(445) =  18
      casejump(446) =  34
      casejump(447) =  72
      casejump(448) =  71
      casejump(449) =  16
      casejump(450) =  11
      casejump(451) =  35
      casejump(452) =   1
      casejump(453) =  36
      casejump(454) =   2
      casejump(455) =   3
      casejump(456) =   4
      casejump(457) =  16
      casejump(458) =  31
      casejump(459) =   5
      casejump(460) =  26
      casejump(461) =   7
      casejump(462) =  47
      casejump(463) =   9
      casejump(464) =  10
      casejump(465) =  37
      casejump(466) =  11
      casejump(467) =  32
      casejump(468) =  30
      casejump(469) =  12
      casejump(470) =  13
      casejump(471) =  38
      casejump(472) =  15
      casejump(473) =  16
      casejump(474) =  31
      casejump(475) =  32
      casejump(476) =  64
      casejump(477) =  18
      casejump(478) =  19
      casejump(479) =  20
      casejump(480) =  63
      casejump(481) =  40
      casejump(482) =  22
      casejump(483) =  23
      casejump(484) =   1
      casejump(485) =  18
      casejump(486) =  13
      casejump(487) =  41
      casejump(488) =   4
      casejump(489) =  24
      casejump(490) =  48
      casejump(491) =   5
      casejump(492) =  26
      casejump(493) =  18
      casejump(494) =  70
      casejump(495) =   9
      casejump(496) =  69
      casejump(497) =  28
      casejump(498) =  11
      casejump(499) =  42
      casejump(500) =  30
      casejump(501) =  43
      casejump(502) =  13
      casejump(503) =  72
      casejump(504) =  71
      casejump(505) =  16
      casejump(506) =  31
      casejump(507) =  32
      casejump(508) =  64
      casejump(509) =  18
      casejump(510) =  70
      casejump(511) =  72
      casejump(512) =  80

end subroutine getcjp512
