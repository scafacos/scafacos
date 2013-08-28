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

!>
!>  Encapsulates functions for periodic boundary conditions
!>
module module_mirror_boxes
    implicit none
    private

      ! lattice basis vectors
      real*8, public :: t_lattice_1(3) = [1, 0, 0] !< 1st vector of lattice basis
      real*8, public :: t_lattice_2(3) = [0, 1, 0] !< 2nd vector of lattice basis
      real*8, public :: t_lattice_3(3) = [0, 0, 1] !< 3rd vector of lattice basis
      real*8, public :: Lattice(3,3) !< holds the lattice transformation matrix
      real*8, public :: LatticeInv(3,3) !< holds the inverse lattice transformation matrix
      real*8, public :: LatticeCenter(3) !< holds the central point of the lattice box
      real*8, public :: LatticeOrigin(3) = [0., 0., 0.] !< holds the lattice origin - can be modified by calling program
      logical, public :: periodicity(3) = [.false., .false., .false.]  !< boolean switches for determining periodicity directions
      integer, public :: mirror_box_layers = 1 !< size of near-field layer (nmber of shells)
      real*8, public :: unit_box_volume = 1.0 !< volume box that is spanned by t_lattice_1..3, is initialized in calc_neighour_boxes()
      !> variables that should not be written to
      integer, public :: num_neighbour_boxes = 1
      integer, dimension(:,:), allocatable, public :: neighbour_boxes !dimensions in 3D case (3,27)
      !> is set to OR(periodicity), to be able to exit from all procedures as fast as possible
      !> if nothing is to be done
      logical, public :: do_periodic
      ! number of boxes to include into each direction
      integer, public :: periodicity_switches(3)
      #ifndef NO_SPATIAL_INTERACTION_CUTOFF
      !> all nodes, where any(abs(coc(1:3)-particle_position(1:3)) > spatial_interaction_cutoff(1:3) are
      !> ignored when calculating interactions (for convenient minimum image method, where only half
      !> of the mirror boxes is included)
      real*8 , public :: spatial_interaction_cutoff(3) = huge(0._8) * [1., 1., 1.]
      #endif

      public calc_neighbour_boxes
      public init_movement_constraint
      public constrain_periodic
      public lattice_vect
      public lattice_indices
      public system_is_unit_cube
      public system_uses_principal_axes
      public check_lattice_boundaries

      ! true for lattices with axes parallel to cartesian system
      logical :: simplelattice

    contains

        !>
        !> Calculates vector with respect to lattice base vectors
        !> @param[in] ijk lattice indices
        !>
        function lattice_vect(ijk)
          implicit none
          integer, intent(in) :: ijk(3)
          real*8 :: lattice_vect(3)

          lattice_vect = ijk(1)*t_lattice_1 + ijk(2)*t_lattice_2 + ijk(3)*t_lattice_3
        end function lattice_vect


        !>
        !> Calculates the lattice indices of the unit cell into which the
        !> argument vector points
        !> @param[in] xyz vector
        !>
        function lattice_indices(xyz)
          implicit none
          real*8, intent(in) :: xyz(3)
          integer :: lattice_indices(3)

          lattice_indices = floor(matmul(xyz, LatticeInv))
        end function lattice_indices


        !>
        !> Prepares the list of neighbour boxes within ws
        !> stores their number in num_neighbour_boxes and their logical
        !> indices/coordinates in neighbour_boxes
        !>
        subroutine calc_neighbour_boxes()
        use module_math_tools, only : cross_product
        implicit none
        integer :: i,j,k,idx

         LatticeCenter = 0.5*(t_lattice_1 + t_lattice_2 + t_lattice_3) + LatticeOrigin

          if (allocated(neighbour_boxes)) deallocate(neighbour_boxes)

          where (periodicity(:))
            periodicity_switches = mirror_box_layers
          end where

          idx = 0

          num_neighbour_boxes = product(2*periodicity_switches+1)
          allocate(neighbour_boxes(3,num_neighbour_boxes))

          ! the boxes of the shells that surround the central box
          ! are ordered in a way, that a box and its central-symmetric counterpart are grouped together
          do i = -periodicity_switches(1),-1
            do j = -periodicity_switches(2),periodicity_switches(2)
              do k = -periodicity_switches(3),periodicity_switches(3)
                  idx = idx + 1
                  neighbour_boxes(:,idx) = [ i,  j ,  k]
                  idx = idx + 1
                  neighbour_boxes(:,idx) = [-i, -j , -k]
              end do
            end do
          end do

          i = 0
          do j = -periodicity_switches(2),-1
            do k = -periodicity_switches(3),periodicity_switches(3)
                idx = idx + 1
                neighbour_boxes(:,idx) = [ i,  j ,  k]
                idx = idx + 1
                neighbour_boxes(:,idx) = [-i, -j , -k]
            end do
          end do

          i = 0
          j = 0
          do k = -periodicity_switches(3),-1
            idx = idx + 1
            neighbour_boxes(:,idx) = [ i,  j ,  k]
            idx = idx + 1
            neighbour_boxes(:,idx) = [-i, -j , -k]
          end do

          neighbour_boxes(:,idx+1) = [0, 0, 0] ! center box is put to the back of the boxlist for easier iteration
          
          unit_box_volume = abs(dot_product(cross_product(t_lattice_1, t_lattice_2),t_lattice_3))

          call init_movement_constraint()
        end subroutine calc_neighbour_boxes


        !>
        !> Initializes transformation matrices between cartesian system and lattice basis
        !>
        subroutine init_movement_constraint
          use module_math_tools
          implicit none

          Lattice(1,:) = t_lattice_1
          Lattice(2,:) = t_lattice_2
          Lattice(3,:) = t_lattice_3
          LatticeInv = Inverse3(Lattice)

          ! simplify the movement constraint if the lattice is really simple
          simplelattice = system_uses_principal_axes()
        end subroutine init_movement_constraint


        !>
        !> Shifts all particles that left the original box
        !> back into it
        !>
        subroutine constrain_periodic(particles)
            use module_math_tools
            use module_pepc_types

            implicit none

            type(t_particle), intent(inout) :: particles(:)
            
            integer(kind_particle) :: np_local, p
            real*8 :: lattice_coord(3), real_coord(3), latticewalls(3)

            np_local = size(particles, kind=kind_particle)
            
            if (simplelattice) then
                latticewalls(1:3) = [t_lattice_1(1), t_lattice_2(2), t_lattice_3(3)]

                do p = 1,np_local
                  where (periodicity) particles(p)%x = modulo(particles(p)%x-LatticeOrigin, latticewalls) + LatticeOrigin
                end do

            else
              ! the lattice axes are not parallel to the outer cartesian axes
                do p = 1,np_local

                    lattice_coord = matmul(particles(p)%x-LatticeOrigin, LatticeInv)

                    if (any(lattice_coord > 1.D+0) .or. any(lattice_coord < 0.D+0)) then

                      lattice_coord = modulo(lattice_coord, 1.D+0)

                      real_coord = matmul(lattice_coord, Lattice) + LatticeOrigin

                      where (periodicity) particles(p)%x = real_coord

                    end if
                end do
            end if
        end subroutine


        !>
        !>  returns |r|^2
        !>
        real*8 function normsq(r)
          implicit none
          real*8, intent(in) :: r(3)

          normsq = dot_product(r, r)
        end function


        !>
        !>  returns .true. if t_lattice_1,2,3 are aligned to the cartesian axes x,y,z
        !>
        logical function system_uses_principal_axes()
          implicit none

          system_uses_principal_axes = (t_lattice_1(1) >  0.) .and. &
                                       (t_lattice_1(2) == 0.) .and. &
                                       (t_lattice_1(3) == 0.) .and. &
                                       (t_lattice_2(1) == 0.) .and. &
                                       (t_lattice_2(2) >  0.) .and. &
                                       (t_lattice_2(3) == 0.) .and. &
                                       (t_lattice_3(1) == 0.) .and. &
                                       (t_lattice_3(2) == 0.) .and. &
                                       (t_lattice_3(3) >  0.)
        end function


        !>
        !>  returns .true. if t_lattice_1,2,3 span a unit box
        !>
        logical function system_is_unit_cube()
          implicit none

          system_is_unit_cube = system_uses_principal_axes() .and.  &
                                (normsq(t_lattice_1) == 1.0) .and.  &
                                (normsq(t_lattice_2) == 1.0) .and.  &
                                (normsq(t_lattice_3) == 1.0)
        end function
        

        !>
        !>  Verifies, if all particles lie inside the central box,
        !>  otherweise LatticeOrigin or LatticeCenter are invalid
        !>
        logical function check_lattice_boundaries(particles)
          use module_pepc_types
          use module_debug
          implicit none
          type (t_particle), dimension(:), intent(in) :: particles
          integer(kind_particle) :: p
          real*8 :: lattice_coord(3)
          
          check_lattice_boundaries = .true.
          
          do p = 1,size(particles)

            lattice_coord = matmul(particles(p)%x-LatticeOrigin, LatticeInv)

            if (any(((lattice_coord > 1.D+0) .or. (lattice_coord < 0.D+0)) .and. periodicity)) then
              DEBUG_WARNING_ALL('("Particle ", I0, " at coords", 3(X,G12.3), " lies outside the central box. Maybe, LatticeCenter or LatticeOrigin is wrong." )', p, particles(p)%x)
            end if
          end do
        end function
end module module_mirror_boxes
