! This file is part of PEPC - The Pretty Efficient Parallel Coulomb Solver.
! 
! Copyright (C) 2002-2012 Juelich Supercomputing Centre, 
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
        module procedure coord_to_key_lastlevel, coord_to_key_level
      end interface coord_to_key

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
        pure elemental function level_from_key(key)
          implicit none
          integer*8, intent(in) :: key
          integer :: level_from_key

          ! using log_8(key):
          ! level_from_key = int( log(1._8*key) / log(8._8))
          ! counting leading zeros (faster):
          level_from_key = int((bit_size(key) - leadz(key) - 1) / 3)

        end function


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> checks whether key_a is an ancestor of key_c (which must be at highest tree level, i.e. a particle key) 
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        pure function is_ancestor_of_particle(key_c,key_a)
          use treevars
          implicit none
          logical :: is_ancestor_of_particle
          integer*8, intent(in) :: key_a, key_c
 
          is_ancestor_of_particle = (ishft(key_c,3*(level_from_key(key_a)-nlev)) == key_a)

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
          integer*8, dimension(3,npp) :: intcoord
          real*8 :: s(3)
          integer :: j

          s=boxsize/2**nlev       ! refinement length

          ! (xmin, ymin, zmin) is the translation vector from the tree box to the simulation region (in 1st octant)
          do j = 1,npp
            intcoord(:,j) = int(( particles(j)%x - boxmin )/s) ! partial keys
          end do

          ! construct particle keys
          select case (curve_type)
            case (0) ! Z-curve
              do j = 1,npp
                 particles(j)%key = intcoord_to_key_morton(intcoord(:,j))
              end do

            case (1) ! Hilbert curve (original pattern)
             do j=1,npp
                 particles(j)%key = intcoord_to_key_hilbert(intcoord(:,j))
             end do

          end select

        end subroutine compute_particle_keys


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> calculates key from particle coordinate on top level
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        function coord_to_key_lastlevel(x, y, z)
          use treevars, only : nlev, boxsize, boxmin
          implicit none
          integer*8 :: coord_to_key_lastlevel
          real*8, intent(in) :: x, y, z
          integer*8 :: ix(3)
          real*8 :: s(3)

          s=boxsize/2**nlev       ! refinement length

          ! (xmin, ymin, zmin) is the translation vector from the tree box to the simulation region (in 1st octant)
          ix = int(( [x, y, z] - boxmin )/s)           ! partial keys

          ! construct particle keys
          select case (curve_type)
            case (0) ! Z-curve
              coord_to_key_lastlevel = intcoord_to_key_morton(ix)
            case (1) ! Hilbert curve (original pattern)
              coord_to_key_lastlevel = intcoord_to_key_hilbert(ix)
            case default
              coord_to_key_lastlevel = 0
          end select

        end function coord_to_key_lastlevel


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> calculates key from particle coordinate on certain tree level
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        function coord_to_key_level(x, y, z, level)
          use treevars, only : nlev
          implicit none
          integer*8 :: coord_to_key_level
          real*8, intent(in) :: x, y, z
          integer, intent(in) :: level

          coord_to_key_level = ishft(coord_to_key_lastlevel(x, y, z),3*(level-nlev))

        end function coord_to_key_level


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> calculates particle coordinate from key
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine key_to_coord(key, x, y, z)
          use treevars, only : nlev, boxsize, boxmin
          implicit none
          integer*8, intent(in) :: key
          real*8, intent(out) :: x, y, z
          integer*8 :: ix, iy, iz
          real*8 :: s(3)

          ! construct particle coordiantes
          select case (curve_type)
            case (0) ! Z-curve
              call key_to_intcoord_morton(key, ix, iy, iz)
            case (1) ! Hilbert curve (original pattern)
              call key_to_intcoord_hilbert(key, ix, iy, iz)
            case default
              ix = 0
              iy = 0
              iz = 0
          end select

          s=boxsize/2**nlev       ! refinement length

          ! (xmin, ymin, zmin) is the translation vector from the tree box to the simulation region (in 1st octant)
          x = (real(ix,kind(1._8)) + 0.5_8) * s(1) + boxmin(1)
          y = (real(iy,kind(1._8)) + 0.5_8) * s(2) + boxmin(2)
          z = (real(iz,kind(1._8)) + 0.5_8) * s(3) + boxmin(3)

        end subroutine key_to_coord



        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> calculates particle coordinate from key  and respects given dimension
        !> \@todo This function should use real 1D, 2D, 3D-spacefilling curves instead
        !> of simply restricting the 3D-curve to lower dimensions
        !> using the default_coordinates parameter, the default value for
        !> coordinates that lie outside the idim range can be given
        !> e.g. default_coordinates=[0.,0.,5.],  idim=2 --> particles are restricted to a xy-plane with z=5.
        !>      default_coordinates=coordinates, idim=1 --> only x-coordinate of coordinates is changed according to spacefilling curve
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine key_to_coord_dim(key, coordinates, idim, default_coordinates)
          use treevars, only : nlev, boxsize, boxmin
          implicit none
          integer*8, intent(in) :: key
          real*8, intent(out) :: coordinates(3)
          integer, intent(in) :: idim !< dimension of the system (1,2,3)
          real*8, intent(in) :: default_coordinates(3) !< values for coordinates that are not affected by spacefilling curve (due to idim restriction)
          integer*8 :: ix, iy, iz
          real*8 :: s(3), x(3)

          ! construct particle coordiantes
          select case (curve_type)
            case (0) ! Z-curve
              call key_to_intcoord_morton(key, ix, iy, iz)
            case (1) ! Hilbert curve (original pattern)
              call key_to_intcoord_hilbert(key, ix, iy, iz)
            case default
              ix = 0
              iy = 0
              iz = 0
          end select

          s=boxsize/2**nlev       ! refinement length

          ! (xmin, ymin, zmin) is the translation vector from the tree box to the simulation region (in 1st octant)
          x = [ (real(ix,kind(1._8)) + 0.5_8) * s(1) + boxmin(1), &
                (real(iy,kind(1._8)) + 0.5_8) * s(2) + boxmin(2), &
                (real(iz,kind(1._8)) + 0.5_8) * s(3) + boxmin(3)  ]

          coordinates(1:idim)   = x(1:idim)
          coordinates(idim+1:3) = default_coordinates(idim+1:3)

        end subroutine key_to_coord_dim


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> (Morton-)Z-curve
        !> construct keys by interleaving coord bits and add placeholder bit
        !> note use of 64-bit constants to ensure correct arithmetic
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        function intcoord_to_key_morton(ic)
          use treevars, only : nlev
          implicit none
          integer*8, intent(in) :: ic(3)
          integer*8 :: intcoord_to_key_morton
          integer :: i
          integer*8 :: cval

            ! set placeholder bit
            intcoord_to_key_morton = 1

            ! key generation
            do i=nlev-1,0,-1
               cval = 4_8*ibits( ic(3),i, 1_8 ) + 2_8*ibits( ic(2),i, 1_8 ) + 1_8*ibits( ic(1),i, 1_8 )
               ! appending bit triple to key
               intcoord_to_key_morton = ior(ishft(intcoord_to_key_morton, 3), cval)
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
          use treevars, only : nlev
          implicit none
          integer*8, intent(out) :: ix, iy, iz
          integer*8, intent(in) :: key
          integer :: i, lev

          lev = level_from_key(key)

          ix = 0
          iy = 0
          iz = 0

          do i=0,lev-1
            ix = ior(ix, ishft(ibits(key, 3*i + 0, 1), nlev-lev+i))
            iy = ior(iy, ishft(ibits(key, 3*i + 1, 1), nlev-lev+i))
            iz = ior(iz, ishft(ibits(key, 3*i + 2, 1), nlev-lev+i))
          end do

        end subroutine key_to_intcoord_morton


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
        function intcoord_to_key_hilbert(ic)
          use treevars
          implicit none
          integer*8, intent(in) :: ic(3)
          integer*8 :: intcoord_to_key_hilbert
          integer :: i

		  integer*8, parameter :: CI(0:7)    = [0,1,3,2,7,6,4,5] ! 3D - inverse hilbert cell
          integer*8, parameter :: G(0:7,0:1) = reshape([5,6,0,5,5,0,6,5,0,0,0,5,0,0,6,5],shape(G))     ! 3D - hilbert gene
		  integer*8 :: horder           ! order of the hilbert cell C
		  integer*8 :: xtemp,ytemp,ztemp,change
		  integer*8 :: cval

	        ! copy, because construction alters original values
	        xtemp=ic(1)
	        ytemp=ic(2)
	        ztemp=ic(3)

	        ! set placeholder bit
	        intcoord_to_key_hilbert = 1

	        ! key generation
	        do i=nlev-1,0,-1

	          cval = 4_8*ibits( ztemp,i, 1_8 ) + 2_8*ibits( ytemp,i, 1_8 ) + 1_8*ibits( xtemp,i, 1_8 )
              xtemp = ibclr(xtemp, i)
              ytemp = ibclr(ytemp, i)
              ztemp = ibclr(ztemp, i)

	           ! get H-order (gathering upper bits as in z-mapping and combine to binary number)
	           horder = CI( cval )

	           ! appending H-order to hkey
	           intcoord_to_key_hilbert = ior(ishft(intcoord_to_key_hilbert, 3), horder)

	           ! transform partial curve with the gene rules for the next level step
	           ! exchange
	           select case (G(horder,0))
	             case (5) ! (= 101[zyx]) --> change z and x
	                change = ztemp
	                ztemp  = xtemp
	                xtemp  = change
                 case (6) ! (= 110[zyx]) --> change z and y
                    change = ztemp
                    ztemp  = ytemp
                    ytemp  = change
	           end select

	           ! reverse
	           select case (G(horder,1))
	             case (5) ! (= 101[zyx]) --> reverse z and x
	               ztemp = iand(not(ztemp),2_8**(i)-1)
	               xtemp = iand(not(xtemp),2_8**(i)-1)
                 case (6) ! (= 110[zyx]) --> reverse z and y
                   ztemp = iand(not(ztemp),2_8**(i)-1)
                   ytemp = iand(not(ytemp),2_8**(i)-1)
	           end select
	        end do

        end function intcoord_to_key_hilbert



        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> inverse Hilbert mapping
        !> input key must be right-adjusted and must contain a placeholder bit
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine key_to_intcoord_hilbert(key, ix, iy, iz)
          use treevars
          implicit none
          integer*8, intent(out) :: ix, iy, iz
          integer*8, intent(in) :: key

		  integer*8 :: change, horder, cval
		  integer*8 :: i
		  integer :: lev

          integer*8, parameter :: C(0:7)    = [0,1,3,2,6,7,5,4] ! 3D - hilbert cell
          integer*8, parameter :: G(0:7,0:1) = reshape([5,6,0,5,5,0,6,5,0,0,0,5,0,0,6,5],shape(G)) ! 3D - hilbert gene
          !integer*8, parameter :: G(0:7,0:1) = reshape([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],shape(G)) ! 3D - hilbert gene

          lev = level_from_key(key)

		  iz = 0
		  iy = 0
		  ix = 0

		  do i=0,lev-1
		     horder = ibits(key,3*i,3)
		     cval   = C(horder)

             if (i>0) then
               select case(G(horder,1))
               case(5)
                  ix=iand(not(ix),2_8**(i)-1)
                  iz=iand(not(iz),2_8**(i)-1)
               case(6)
                  iy=iand(not(iy),2_8**(i)-1)
                  iz=iand(not(iz),2_8**(i)-1)
               end select

               select case(G(horder,0))
               case(5)
                  change = ix
                  ix     = iz
                  iz     = change
               case(6)
                  change = iy
                  iy     = iz
                  iz     = change
               end select
             end if

             iz=ior(ishft(ibits(cval,2,1),i),iz)
             iy=ior(ishft(ibits(cval,1,1),i),iy)
             ix=ior(ishft(ibits(cval,0,1),i),ix)
		  end do

          ix = ishft(ix, nlev-lev)
          iy = ishft(iy, nlev-lev)
          iz = ishft(iz, nlev-lev)

        end subroutine key_to_intcoord_hilbert



end module module_spacefilling
