! This file is part of PEPC - The Pretty Efficient Parallel Coulomb Solver.
! 
! Copyright (C) 2002-2013 Juelich Supercomputing Centre, 
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!> Contains mapper functions for space-filling curves
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module module_spacefilling
      implicit none

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  public variable declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      integer, public :: curve_type = 1 !(0: Morton, 1: Hilbert) 

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  public subroutine declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      interface coord_to_key
        module procedure coord_to_key_lastlevel, coord_to_key_level, &
          veccoord_to_key_lastlevel, veccoord_to_key_level
      end interface coord_to_key

      interface key_to_coord
        module procedure key_to_coord, key_to_veccoord
      end interface key_to_coord

      interface is_ancestor_of_particle
        module procedure is_ancestor_of_particle_nolevel, is_ancestor_of_particle_withlevel
      end interface is_ancestor_of_particle

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  private variable declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  subroutine-implementation  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      contains

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> calculates level from key by finding position of placeholder bit
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        elemental function level_from_key(key)
          use treevars, only: idim
          implicit none
          integer*8, intent(in) :: key
          integer :: level_from_key

          ! using log_{2**idim}(key):
          ! level_from_key = int( log(1._8*key) / log((2._8)**idim))
          ! counting leading zeros (faster):
          level_from_key = int((bit_size(key) - leadz(key) - 1) / idim)

        end function


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> returns the key of a key's parent
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        elemental function parent_key_from_key(key)
          implicit none
          integer*8, intent(in) :: key
          integer*8 :: parent_key_from_key

          parent_key_from_key = shift_key_by_level(key, -1)
        end function parent_key_from_key


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> returns the child number of `key` with respect to its parent
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        elemental function child_number_from_key(key)
          use treevars, only: idim
          implicit none

          integer*8, intent(in) :: key
          integer :: child_number_from_key

          child_number_from_key = int(ibits(key, 0, idim))
        end function child_number_from_key


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> returns `key` shifted up (negative argument) or down by a number
        !> of levels
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        pure function shift_key_by_level(key, lvl)
          use treevars, only: idim
          implicit none

          integer*8, intent(in) :: key
          integer, intent(in) :: lvl
          integer*8 :: shift_key_by_level

          shift_key_by_level = ishft(key, idim * lvl)
        end function shift_key_by_level


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> returns the key for child `n` of a node with key `key` 
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        pure function child_key_from_parent_key(key, n)
          implicit none

          integer*8, intent(in) :: key
          integer, intent(in) :: n
          integer*8 :: child_key_from_parent_key

          child_key_from_parent_key = shift_key_by_level(key, 1) + n
        end function child_key_from_parent_key


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> checks whether key_a is an ancestor of key_c (which must be at highest tree level, i.e. a particle key) 
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        pure function is_ancestor_of_particle_nolevel(key_c,key_a)
          use treevars, only: nlev
          implicit none
          logical :: is_ancestor_of_particle_nolevel
          integer*8, intent(in) :: key_a, key_c
 
          is_ancestor_of_particle_nolevel = &
            (shift_key_by_level(key_c, level_from_key(key_a) - nlev) == key_a)

        end function


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> checks whether key_a is an ancestor of key_c (which must be at highest tree level, i.e. a particle key) 
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        pure function is_ancestor_of_particle_withlevel(key_c,key_a,level_a)
          use treevars, only: nlev
          implicit none
          logical :: is_ancestor_of_particle_withlevel
          integer*8, intent(in) :: key_a, key_c
          integer, intent(in) :: level_a
 
          is_ancestor_of_particle_withlevel = (shift_key_by_level(key_c, level_a - nlev) == key_a)

        end function

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> calculates keys from local particles (faster than per-particle call to coord_to_key())
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine compute_particle_keys(particles)
          use treevars
          use module_debug
          implicit none
          type(t_particle), intent(inout) :: particles(1:npp)
          integer*8, dimension(:,:), allocatable :: intcoord
          real*8 :: s(idim)
          integer :: j

          allocate(intcoord(idim, npp))

          s=boxsize(1:idim)/2_8**nlev       ! refinement length

          ! (xmin, ymin, zmin) is the translation vector from the tree box to the simulation region (in 1st octant)
          do j = 1,npp
            intcoord(:,j) = int(( particles(j)%x(1:idim) - boxmin(1:idim) )/s, kind = 8) ! partial keys
          end do

          ! construct particle keys
          select case (curve_type)
            case (0) ! Z-curve
              do j = 1,npp
                particles(j)%key = intcoord_to_key_morton(intcoord(:,j))
              end do

            case (1) ! Hilbert curve (original pattern)
              select case (idim)
                case (1) ! fallback to Z-curve
                  do j=1,npp
                    particles(j)%key = intcoord_to_key_morton(intcoord(:,j))
                  end do
                case (2) ! 2D hilbert curve
                  do j=1,npp
                    particles(j)%key = intcoord_to_key_hilbert2D(intcoord(:,j))
                  end do
                case (3) ! 3D hilbert curve
                  do j=1,npp
                    particles(j)%key = intcoord_to_key_hilbert3D(intcoord(:,j))
                  end do
              end select
          end select

          deallocate(intcoord)

        end subroutine compute_particle_keys


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> calculates key from particle coordinate on top level
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        function coord_to_key_lastlevel(x, y, z)
          implicit none

          integer*8 :: coord_to_key_lastlevel
          real*8, intent(in) :: x, y, z

          coord_to_key_lastlevel = veccoord_to_key_lastlevel([x, y, z])
        end function coord_to_key_lastlevel


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> calculates key from particle coordinate as vector on top level
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        function veccoord_to_key_lastlevel(x)
          use treevars, only : idim, nlev, boxsize, boxmin
          implicit none
          integer*8 :: veccoord_to_key_lastlevel
          real*8, intent(in) :: x(3)
          integer*8 :: ic(idim)
          real*8 :: s(idim)

          s=boxsize(1:idim)/2_8**nlev       ! refinement length

          ! (xmin, ymin, zmin) is the translation vector from the tree box to the simulation region (in 1st octant)
          ic = int((x(1:idim) - boxmin(1:idim))/s, kind = 8)           ! partial keys

          ! construct particle keys
          select case (curve_type)
            case (0) ! Z-curve
              veccoord_to_key_lastlevel = intcoord_to_key_morton(ic)
            case (1) ! Hilbert curve (original pattern)
              select case (idim)
                case (1) ! fallback to Z-curve
                  veccoord_to_key_lastlevel = intcoord_to_key_morton(ic)
                case (2) ! 2D hilbert curve
                  veccoord_to_key_lastlevel = intcoord_to_key_hilbert2D(ic)
                case (3) ! 3D hilbert curve
                  veccoord_to_key_lastlevel = intcoord_to_key_hilbert3D(ic)
                case default
                  veccoord_to_key_lastlevel = 0
              end select
            case default
              veccoord_to_key_lastlevel = 0
          end select

        end function veccoord_to_key_lastlevel


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> calculates key from particle coordinate on certain tree level
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        function coord_to_key_level(x, y, z, level)
          implicit none

          integer*8 :: coord_to_key_level
          real*8, intent(in) :: x, y, z
          integer, intent(in) :: level

          coord_to_key_level = veccoord_to_key_level([x, y, z], level)
        end function coord_to_key_level


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> calculates key from particle coordinate as vector on certain tree level
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        function veccoord_to_key_level(x, level)
          use treevars, only : nlev
          implicit none
          integer*8 :: veccoord_to_key_level
          real*8, intent(in) :: x(3)
          integer, intent(in) :: level

          veccoord_to_key_level = shift_key_by_level(veccoord_to_key_lastlevel(x), level-nlev)
        end function veccoord_to_key_level


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> calculates particle coordinate from key
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine key_to_coord(key, x, y, z)
          implicit none

          integer*8, intent(in) :: key
          real*8, intent(out) :: x, y, z

          real*8 :: xv(3)

          call key_to_veccoord(key, xv)
          x = xv(1)
          y = xv(2)
          z = xv(3)
        end subroutine key_to_coord


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> calculates particle coordinate as vector from key
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine key_to_veccoord(key, x)
          use treevars, only : idim, nlev, boxsize, boxmin
          implicit none
          integer*8, intent(in) :: key
          real*8, intent(inout) :: x(3)
          integer*8 :: ic(idim)
          real*8 :: s(idim)

          ! construct particle coordiantes
          select case (curve_type)
            case (0) ! Z-curve
              call key_to_vecintcoord_morton(key, ic)
            case (1) ! Hilbert curve (original pattern)
              call key_to_vecintcoord_hilbert(key, ic)
            case default
              ic = 0
          end select

          s=boxsize(1:idim)/2_8**nlev       ! refinement length

          ! (xmin, ymin, zmin) is the translation vector from the tree box to the simulation region (in 1st octant)
          x(1:idim) = (real(ic,kind(1._8)) + 0.5_8) * s + boxmin(1:idim)

        end subroutine key_to_veccoord


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> (Morton-)Z-curve
        !> construct keys by interleaving coord bits and add placeholder bit
        !> note use of 64-bit constants to ensure correct arithmetic
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        function intcoord_to_key_morton(ic)
          use treevars, only : idim, nlev
          implicit none
          integer*8, intent(in) :: ic(idim)
          integer*8 :: intcoord_to_key_morton
          integer :: i, j

            ! set placeholder bit
            intcoord_to_key_morton = 1_8

            ! key generation
            do i=nlev-1,0,-1
              do j=idim,1,-1
                intcoord_to_key_morton = ior(ishft(intcoord_to_key_morton, 1), ibits(ic(j), i, 1_8))
              end do
            end do

        end function intcoord_to_key_morton


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> (Morton-)Z-curve
        !> keys were constructed by interleaving coord bits and adding placeholder bit
        !> input key must be right-adjusted and must contain a placeholder bit
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine key_to_intcoord_morton(key, ix, iy, iz)
          implicit none

          integer*8, intent(in) :: key
          integer*8, intent(out) :: ix, iy, iz

          integer*8 :: ic(3)

          call key_to_vecintcoord_morton(key, ic)
          ix = ic(1)
          iy = ic(2)
          iz = ic(3)
          
        end subroutine key_to_intcoord_morton


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> (Morton-)Z-curve
        !> keys were constructed by interleaving coord bits and adding placeholder bit
        !> input key must be right-adjusted and must contain a placeholder bit
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine key_to_vecintcoord_morton(key, ic)
          use treevars, only : idim, nlev
          implicit none
          integer*8, intent(out) :: ic(idim)
          integer*8, intent(in) :: key
          integer :: i, j, lev

          lev = level_from_key(key)

          ic = 0

          do i=0,lev-1
            do j=1,idim
              ic(j) = ior(ic(j), ishft(ibits(key, idim*i + j - 1, 1), nlev-lev+i))
            end do
          end do

        end subroutine key_to_vecintcoord_morton


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> Hilbert-curve,
        !>
        !> algorithm from
        !>
        !> Kamata, S.-I.; Eason, R.O.; Bandou, Y.; ,
        !> "A new algorithm for N-dimensional Hilbert scanning",
        !> Image Processing, IEEE Transactions on , vol.8, no.7, pp.964-973, Jul 1999
        !> doi: 10.1109/83.772242
        !> http://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=772242&isnumber=16772
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        function intcoord_to_key_hilbert2D(ic)
          use treevars
          implicit none
          integer*8, intent(in) :: ic(2)
          integer*8 :: intcoord_to_key_hilbert2D
          integer :: i, j

          integer*8, parameter :: CI2(0:3)    = [0,3,1,2] ! 2D - inverse hilbert cell
          integer*8, parameter :: G2(0:3,0:1) = reshape([3,0,0,3,0,0,0,3],shape(G2))
          integer*8 :: horder           ! order of the hilbert cell C
          integer*8 :: exchange, reverse
          integer*8 :: itemp(2), change
          integer*8 :: cval

          ! copy, because construction alters original values
          itemp=ic

          ! set placeholder bit
          intcoord_to_key_hilbert2D = 1

          ! key generation
          do i=nlev-1,0,-1

            cval = 0_8

            do j=2,1,-1
              cval = ior(ishft(cval, 1), ibits(itemp(j), i, 1_8))
              itemp(j) = ibclr(itemp(j), i)
            end do

            ! get H-order (gathering upper bits as in z-mapping and combine to binary number)
            horder = CI2(cval)
            exchange = G2(horder,0)
            reverse = G2(horder,1)

            ! appending H-order to hkey
            intcoord_to_key_hilbert2D = ior(ishft(intcoord_to_key_hilbert2D, 2), horder)

            ! transform partial curve with the gene rules for the next level step
            ! exchange
            select case (exchange)
              case (3)
                change = itemp(1)
                itemp(1) = itemp(2)
                itemp(2) = change
            end select
                
            ! reverse
            do j=1,2
              if (btest(reverse, j - 1)) itemp(j) = iand(not(itemp(j)), 2_8**(i) - 1)
            end do
          end do

        end function intcoord_to_key_hilbert2D


        function intcoord_to_key_hilbert3D(ic)
          use treevars
          implicit none
          integer*8, intent(in) :: ic(3)
          integer*8 :: intcoord_to_key_hilbert3D
          integer :: i, j

          integer*8, parameter :: CI3(0:7)    = [0,1,3,2,7,6,4,5] ! 3D - inverse hilbert cell
          integer*8, parameter :: G3(0:7,0:1) = reshape([5,6,0,5,5,0,6,5,0,0,0,5,0,0,6,5],shape(G3))     ! 3D - hilbert gene
          integer*8 :: horder           ! order of the hilbert cell C
          integer*8 :: exchange, reverse
          integer*8 :: itemp(3), change
          integer*8 :: cval

          ! copy, because construction alters original values
          itemp=ic

          ! set placeholder bit
          intcoord_to_key_hilbert3D = 1

          ! key generation
          do i=nlev-1,0,-1

            cval = 0_8

            do j=3,1,-1
              cval = ior(ishft(cval, 1), ibits(itemp(j), i, 1_8))
              itemp(j) = ibclr(itemp(j), i)
            end do

            ! get H-order (gathering upper bits as in z-mapping and combine to binary number)
            horder = CI3(cval)
            exchange = G3(horder,0)
            reverse = G3(horder,1)

            ! appending H-order to hkey
            intcoord_to_key_hilbert3D = ior(ishft(intcoord_to_key_hilbert3D, 3), horder)

            ! transform partial curve with the gene rules for the next level step
            ! exchange
            select case (exchange)
              case (5)
                change = itemp(1)
                itemp(1) = itemp(3)
                itemp(3) = change
              case (6)
                change = itemp(2)
                itemp(2) = itemp(3)
                itemp(3) = change
            end select

            ! reverse
            do j=1,3
              if (btest(reverse, j - 1)) itemp(j) = iand(not(itemp(j)), 2_8**(i) - 1)
            end do
          end do

        end function intcoord_to_key_hilbert3D


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> inverse Hilbert mapping
        !> input key must be right-adjusted and must contain a placeholder bit
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine key_to_intcoord_hilbert(key, ix, iy, iz)
          implicit none
          integer*8, intent(out) :: ix, iy, iz
          integer*8, intent(in) :: key

          integer*8 :: ic(3)

          call key_to_vecintcoord_hilbert(key, ic)
          ix = ic(1)
          iy = ic(2)
          iz = ic(3)

        end subroutine key_to_intcoord_hilbert


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> inverse Hilbert mapping
        !> input key as vector must be right-adjusted and must contain a placeholder bit
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine key_to_vecintcoord_hilbert(key, ic)
          use treevars
          implicit none
          integer*8, intent(out) :: ic(idim)
          integer*8, intent(in) :: key

          integer*8 :: change, horder, cval
          integer*8 :: exchange, reverse
          integer*8 :: i, j, k
          integer :: lev

          integer*8, parameter :: C2(0:3) = [0,2,3,1] ! 2D - hilbert cell
          integer*8, parameter :: G2(0:3,0:1) = reshape([3,0,0,3,0,0,0,3],shape(G2)) ! 2D - hilbert gene
          integer*8, parameter :: C3(0:7)    = [0,1,3,2,6,7,5,4] ! 3D - hilbert cell
          integer*8, parameter :: G3(0:7,0:1) = reshape([5,6,0,5,5,0,6,5,0,0,0,5,0,0,6,5],shape(G3)) ! 3D - hilbert gene
          !integer*8, parameter :: G(0:7,0:1) = reshape([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],shape(G)) ! 3D - hilbert gene

          lev = level_from_key(key)

          ic = 0_8

          do i=0,lev-1
            horder = ibits(key,idim*i,idim)

            if (idim == 2) then
              cval = C2(horder)
              reverse = G2(horder,1)
              exchange = G2(horder,0)
            else
              cval = C3(horder)
              reverse = G3(horder,1)
              exchange = G3(horder,0)
            end if

            if (i>0) then
              ! reverse
              do j=1,idim
                if (btest(reverse, j - 1)) ic(j) = iand(not(ic(j)), 2_8**(i) - 1)
              end do

              ! exchange
              do j=1,idim-1
                if (.not. btest(exchange, j - 1)) cycle

                do k=j+1,idim
                  if (btest(exchange, k - 1)) then
                    change = ic(j)
                    ic(j) = ic(k)
                    ic(k) = change
                  end if
                end do
              end do
            end if

            do j=1,idim
              ic(j) = ior(ishft(ibits(cval, j-1, 1), i), ic(j))
            end do
          end do

          do j=1,idim
            ic(j) = ishft(ic(j), nlev-lev)
          end do

        end subroutine key_to_vecintcoord_hilbert

end module module_spacefilling
