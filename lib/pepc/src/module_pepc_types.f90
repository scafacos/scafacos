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
!> Contains all lpepc-specific types and routines for registering them to MPI
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module module_pepc_types
  use module_interaction_specific_types
  implicit none

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  public variable declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      integer :: MPI_TYPE_particle_data,   &
                 MPI_TYPE_tree_node_interaction_data,   &
                 MPI_TYPE_particle_results, &
                 MPI_TYPE_particle,         &
                 MPI_TYPE_tree_node_transport_package

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  public subroutine declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !> Data structure for shipping single particles
      integer, private, parameter :: nprops_particle = 8 ! # particle properties to ship
      type t_particle
         real*8 :: x(1:3)      !< coordinates
         real*8 :: work        !< work load from force sum
         integer*8 :: key      !< particle key, i.e. key on highgest tree level
         integer*8 :: key_leaf !< key of corresponding leaf (tree node)
         integer :: label      !< particle label (only for diagnostic purposes, can be used freely by the frontend
         integer :: pid        !< particle owner
         type(t_particle_data) :: data       !< real physics (charge, etc.)
         type(t_particle_results) :: results !< results of calc_force_etc and companions
      end type t_particle

      ! Data structure for shipping multiple moments of child nodes
      integer, private, parameter :: nprops_tree_node_transport_package = 5 ! Number of multipole properties to ship
      type t_tree_node_transport_package
         integer*8 :: key     ! key
         integer   :: byte    ! byte code
         integer   :: leaves  ! # leaves contained
         integer   :: owner   ! owner where multipole resides
         type(t_tree_node_interaction_data) :: m ! real physics
      end type t_tree_node_transport_package


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
      !> Creates and registers lpepc-MPI types
      !>
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine register_lpepc_mpi_types()
        use module_interaction_specific_types
        implicit none
        include 'mpif.h'

        integer, parameter :: max_props = nprops_particle + nprops_tree_node_transport_package

        integer :: ierr
        ! address calculation
        integer, dimension(1:max_props) :: blocklengths, displacements, types
        integer(KIND=MPI_ADDRESS_KIND), dimension(0:max_props) :: address
        ! dummies for address calculation
        type(t_particle)  :: dummy_particle
        type(t_tree_node_transport_package) :: dummy_tree_node

        ! first register the interaction-specific MPI types since they are embedded into the lpepc-datatypes
        call register_interaction_specific_mpi_types(MPI_TYPE_particle_data, MPI_TYPE_tree_node_interaction_data, MPI_TYPE_particle_results)

        ! register particle type
        blocklengths(1:nprops_particle)  = [3, 1, 1, 1, 1, 1, 1, 1]
        types(1:nprops_particle)         = [MPI_REAL8, MPI_REAL8, MPI_INTEGER8, MPI_INTEGER8, MPI_INTEGER, MPI_INTEGER, MPI_TYPE_particle_data, MPI_TYPE_particle_results]
        call MPI_GET_ADDRESS( dummy_particle,          address(0), ierr )
        call MPI_GET_ADDRESS( dummy_particle%x,        address(1), ierr )
        call MPI_GET_ADDRESS( dummy_particle%work,     address(2), ierr )
        call MPI_GET_ADDRESS( dummy_particle%key,      address(3), ierr )
        call MPI_GET_ADDRESS( dummy_particle%key_leaf, address(4), ierr )
        call MPI_GET_ADDRESS( dummy_particle%label,    address(5), ierr )
        call MPI_GET_ADDRESS( dummy_particle%pid,      address(6), ierr )
        call MPI_GET_ADDRESS( dummy_particle%data,     address(7), ierr )
        call MPI_GET_ADDRESS( dummy_particle%results,  address(8), ierr )
        displacements(1:nprops_particle) = int(address(1:nprops_particle) - address(0))
        call MPI_TYPE_STRUCT( nprops_particle, blocklengths, displacements, types, MPI_TYPE_particle, ierr )
        call MPI_TYPE_COMMIT( MPI_TYPE_particle, ierr)

        ! register tree_node type
        blocklengths(1:nprops_tree_node_transport_package)  = [1, 1, 1, 1, 1]
        types(1:nprops_tree_node_transport_package)         = [MPI_INTEGER8, MPI_INTEGER, MPI_INTEGER, MPI_INTEGER, MPI_TYPE_tree_node_interaction_data]
        call MPI_GET_ADDRESS( dummy_tree_node,        address(0), ierr )
        call MPI_GET_ADDRESS( dummy_tree_node%key,    address(1), ierr )
        call MPI_GET_ADDRESS( dummy_tree_node%byte,   address(2), ierr )
        call MPI_GET_ADDRESS( dummy_tree_node%leaves, address(3), ierr )
        call MPI_GET_ADDRESS( dummy_tree_node%owner,  address(4), ierr )
        call MPI_GET_ADDRESS( dummy_tree_node%m    ,  address(5), ierr )
        displacements(1:nprops_tree_node_transport_package) = int(address(1:nprops_tree_node_transport_package) - address(0))
        call MPI_TYPE_STRUCT( nprops_tree_node_transport_package, blocklengths, displacements, types, MPI_TYPE_tree_node_transport_package, ierr )
        call MPI_TYPE_COMMIT( MPI_TYPE_tree_node_transport_package, ierr)

    end subroutine register_lpepc_mpi_types


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !>
      !> Deregisters lpepc- and interaction-specific MPI types
      !>
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine free_lpepc_mpi_types()
        implicit none
        include 'mpif.h'
        integer :: ierr

        call MPI_TYPE_FREE( MPI_TYPE_tree_node_transport_package,        ierr)
        call MPI_TYPE_FREE( MPI_TYPE_particle,         ierr)
        call MPI_TYPE_FREE( MPI_TYPE_particle_results, ierr)
        call MPI_TYPE_FREE( MPI_TYPE_tree_node_interaction_data,   ierr)
        call MPI_TYPE_FREE( MPI_TYPE_particle_data,    ierr)

      end subroutine free_lpepc_mpi_types

end module module_pepc_types
