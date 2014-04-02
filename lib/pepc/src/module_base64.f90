! This file is part of PEPC - The Pretty Efficient Parallel Coulomb Solver.
! 
! Copyright (C) 2002-2014 Juelich Supercomputing Centre, 
!                         Forschungszentrum Juelich GmbH,
!                         Germany
! 
! PEPC is free software: you can redistribute it and/or modify
! it under the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
! 
! PEPC is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
! 
! You should have received a copy of the GNU Lesser General Public License
! along with PEPC.  If not, see <http://www.gnu.org/licenses/>.
!

module module_base64
  implicit none

    character, parameter :: chars(0:63) = ['A','B','C','D','E','F','G','H', &
                                              'I','J','K','L','M','N','O','P', &
                                              'Q','R','S','T','U','V','W','X', &
                                              'Y','Z','a','b','c','d','e','f', &
                                              'g','h','i','j','k','l','m','n', &
                                              'o','p','q','r','s','t','u','v', &
                                              'w','x','y','z','0','1','2','3', &
                                              '4','5','6','7','8','9','+','/']
    character, parameter :: finalchar = '='

    type base64_encoder
      private
        integer*8 :: buffer  = 0
        integer   :: bits    = 0
        integer   :: istream = 6
        logical   :: bigendian = .true.

      contains
           procedure :: getnextbyte    => base64_encoder_getnextbyte
           procedure :: flushbuffer  => base64_encoder_flushbuffer
           generic   :: encode       => encode_real8, encode_real4, encode_int8, encode_int4
           procedure :: finish       => base64_encoder_finish
           procedure :: start        => base64_encoder_start
           procedure :: encode_real8 => base64_encoder_encode_real8
           procedure :: encode_real4 => base64_encoder_encode_real4
           procedure :: encode_int8  => base64_encoder_encode_int8
           procedure :: encode_int4  => base64_encoder_encode_int4
    end type base64_encoder

    contains


      function base64_encoder_getnextbyte(base64)
        implicit none
        class(base64_encoder) :: base64
        integer*1 :: base64_encoder_getnextbyte
        base64_encoder_getnextbyte   = int(iand(ibits(base64%buffer, 58, 6), Z'3F'), kind(base64_encoder_getnextbyte))
        base64%buffer                = ishft(base64%buffer, 6)
      end function base64_encoder_getnextbyte


      subroutine base64_encoder_flushbuffer(base64)
        implicit none
        class(base64_encoder) :: base64
        integer :: i
        character  :: reschars(4)
        do while (base64%bits >= 24)
          do i=1,4
            reschars(i) = chars(base64%getnextbyte())
          end do
          base64%bits = base64%bits - 24
          write(base64%istream, '(4a)', advance='no') reschars(1:4)
        end do
      end subroutine base64_encoder_flushbuffer


      subroutine base64_encoder_finish(base64)
        implicit none
        class(base64_encoder) :: base64
        integer :: i
        integer :: leftbytes
        character  :: reschars(4)
        leftbytes = (base64%bits+5)/6

        if (leftbytes > 0) then
          do i=1,leftbytes
            reschars(i) = chars(base64%getnextbyte())
          end do
          do i=leftbytes+1,4
            reschars(i) = finalchar
          end do
          write(base64%istream, '(4a)', advance='no') reschars(1:4)
        endif

        base64%bits   = 0
        base64%buffer = 0_8
      end subroutine base64_encoder_finish


      subroutine base64_encoder_start(base64, istream_, bigendian_)
        implicit none
        class(base64_encoder) :: base64
        integer :: istream_
        logical :: bigendian_
        base64%bits      = 0
        base64%buffer    = 0_8
        base64%istream   = istream_
        base64%bigendian = bigendian_
      end subroutine base64_encoder_start


      subroutine base64_encoder_encode_real8(base64, data)
        implicit none
        class(base64_encoder) :: base64
        real*8, intent(in), value :: data
        integer*8 :: tmp
        tmp = transfer(data, tmp)
        call base64%encode_int8(tmp)
      end subroutine base64_encoder_encode_real8


      subroutine base64_encoder_encode_real4(base64, data)
        implicit none
        class(base64_encoder) :: base64
        real*4, intent(in), value :: data
        integer*4 :: tmp
        tmp = transfer(data, tmp)
        call base64%encode_int4(tmp)
      end subroutine base64_encoder_encode_real4


      subroutine base64_encoder_encode_int8(base64, data)
        implicit none
        class(base64_encoder) :: base64
        integer*8, intent(in), value :: data
        integer*4 :: tmp1, tmp2
        tmp1 = int(ibits(data,  0, 32), kind(tmp1))
        tmp2 = int(ibits(data, 32, 32), kind(tmp2))
        call base64%encode_int4(tmp2)
        call base64%encode_int4(tmp1)
      end subroutine base64_encoder_encode_int8


      subroutine base64_encoder_encode_int4(base64, data)
        implicit none
        class(base64_encoder) :: base64
        integer*4, value, intent(in) :: data
        integer*8 :: tmp
        tmp = transfer(data, tmp)
        if (.not. base64%bigendian) tmp = ishft(tmp, -32)
        tmp = ishft(iand(tmp, Z'FFFFFFFF'), 32-base64%bits)
        base64%buffer = ior(base64%buffer, tmp)
        base64%bits = base64%bits + 32
        call base64%flushbuffer()
     end subroutine base64_encoder_encode_int4

end module module_base64

