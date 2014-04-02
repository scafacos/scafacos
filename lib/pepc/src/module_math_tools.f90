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

!>
!> Several math tools, primarily matrix manipulations,
!> Legendre polynomials and related stuff. Additionally
!> functions to calculate the biggest power in given interval
!>
module module_math_tools
      implicit none
      save
      private

      public inverse3
      public bpi
      public cross_product

      interface cross_product
        module procedure cross_product3
      end interface

      contains

        !>
        !> Vector (cross) product of two 3D-vectors
        !>
        function cross_product3(a, b)
          implicit none
          real*8, dimension(1:3) :: cross_product3
          real*8, dimension(1:3), intent(in) :: a, b

          cross_product3(1) = a(2)*b(3) - a(3)*b(2)
          cross_product3(2) = a(3)*b(1) - a(1)*b(3)
          cross_product3(3) = a(1)*b(2) - a(2)*b(1)

        end function cross_product3


        !>
        !> Determinant of a real*8 2x2 matrix
        !>
        real*8 function det2(mat)
          implicit none
          real*8, intent(in) :: mat(1:2,1:2)

          det2 = mat(1,1)*mat(2,2) - mat(1,2)*mat(2,1)

        end function det2


        !>
        !> Determinant of a real*8 2x2 matrix that is saved as a vector (columns first)
        !>
        real*8 function det2f(mat)
          implicit none
          real*8, intent(in) :: mat(1:4)

          det2f = mat(1)*mat(4)-mat(2)*mat(3)

        end function det2f


        !>
        !> Determinant of a real*8 3x3 matrix
        !>
        real*8 function det3(mat)
          implicit none
          real*8, intent(in) :: mat(1:3,1:3)

          det3 = mat(1,1)*mat(2,2)*mat(3,3) + mat(1,2)*mat(2,3)*mat(3,1) + mat(1,3)*mat(2,1)*mat(3,2) &
               - mat(1,3)*mat(2,2)*mat(3,1) - mat(1,2)*mat(2,1)*mat(3,3) - mat(1,1)*mat(2,3)*mat(3,2)

        end function det3


        !>
        !>
        !>
        real*8 function cofact(mat, i, j)
          implicit none
          real*8, intent(in) :: mat(1:3,1:3)
          integer, intent(in) :: i, j

          real*8 :: cf(1:2,1:2)

          cf(1:i-1,1:j-1) = mat(1  :i-1,1  :j-1)
          cf(i:   ,1:j-1) = mat(i+1:   ,1  :j-1)
          cf(1:i-1,j:   ) = mat(1  :i-1,j+1:   )
          cf(i:   ,j:   ) = mat(i+1:   ,j+1:   )

          cofact = det2(cf)

          write(*,*) 'cofact(A, i, j) = ', i, j, cofact

        end function cofact


        !>
        !>  Simple 3x3-matrix inversion using Cramers rule
        !>
        function inverse3(m)
          implicit none
          real*8, intent(in) :: m(1:3,1:3)

          real*8 :: inverse3(1:3,1:3)
          real*8 :: test(1:3,1:3)
          real*8 :: det
          integer :: i,j
          logical :: failed

          det = det3(m)
          inverse3 = 0

          if (det == 0) then
            write (*,*) 'You are trying to invert a singular matrix - this is really evil'
            stop
          end if

          inverse3 = reshape([ & ! first column
                               det2f([m(2,2),m(3,2),m(2,3),m(3,3)]), &
                               det2f([m(2,3),m(3,3),m(2,1),m(3,1)]), &
                               det2f([m(2,1),m(3,1),m(2,2),m(3,2)]), &
                               ! 2nd column
                               det2f([m(1,3),m(3,3),m(1,2),m(3,2)]), &
                               det2f([m(1,1),m(3,1),m(1,3),m(3,3)]), &
                               det2f([m(1,2),m(3,2),m(1,1),m(3,1)]), &
                               ! 3rd column
                               det2f([m(1,2),m(2,2),m(1,3),m(2,3)]), &
                               det2f([m(1,3),m(2,3),m(1,1),m(2,1)]), &
                               det2f([m(1,1),m(2,1),m(1,2),m(2,2)]) ], [3,3])

          inverse3 = inverse3 / det

          ! just for testing purposes
          test = matmul(inverse3,m)
          failed = .false.

          do i=1,3
            do j=1,3
              if (i==j) then
                failed = failed .or. (abs(test(i,j)-1.D0) > 1.D-15)
              else
                failed = failed .or. (abs(test(i,j))      > 1.D-15)
              end if
            end do
          end do

          if (failed) then
            write (*,*) 'inverse3 failed for the following input matrix:', m, 'output: ', inverse3, 'test: ', test
            stop
          end if

        end function inverse3


        !>
        !> Calculates the biggest power to \a base in a given interval with the limits \a a and \a b. In this
        !> routine it is done by a more direct expression. Be careful, that \f$a<b\f$!
        !>
        !> @param[in] a First Limit of interval.
        !> @param[in] b Second Limit of interval.
        !>
        function bpi(a, b)
          use module_pepc_types
          use treevars, only: idim
          use module_spacefilling
          use module_debug
          implicit none
          include 'mpif.h'

          integer(kind_key), intent(in) :: a, b
          integer(kind_key) :: bpi

          integer(kind_key) :: axorb, mask
          integer(kind_level) :: bpilevel

          DEBUG_ASSERT_MSG(a < b ,'("Error: b < a in math_tools::bpi(a = ",o22,"(oct), ",i22,"(dec), b = ",o22,"(oct), ",i22,"(dec))")', a,a,b,b)

          axorb             = ieor(a,b)
          bpilevel          = level_from_key(axorb)
          mask              = not(2_8**(idim*bpilevel) - 1)
          bpi               = iand(b,mask) ! extract highest bits from b (which has to be the larger of both parameters)

        end function bpi
end module module_math_tools
