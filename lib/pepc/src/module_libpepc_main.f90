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
!> main internal pepc routines
!>
module module_libpepc_main
    use module_debug, only : debug_level
    use treevars, only : np_mult, interaction_list_length_factor, num_threads, idim
    use module_spacefilling, only : curve_type
    use module_domains, only: weighted
    use module_box, only: force_cubic_domain
    use module_mirror_boxes, only: mirror_box_layers, periodicity
    use module_pepc_types

    implicit none
    private

    public libpepc_restore_particles
    public libpepc_traverse_tree
    public libpepc_grow_tree
    public libpepc_timber_tree
    public libpepc_read_parameters
    public libpepc_write_parameters

    namelist /libpepc/ debug_level, periodicity, np_mult, curve_type, force_cubic_domain, weighted, interaction_list_length_factor, mirror_box_layers, num_threads, idim

    contains


    !>
    !> reads libpepc-specific parameters from file
    !>
    subroutine libpepc_read_parameters(filehandle)
        use module_debug, only: pepc_status
        implicit none
        integer, intent(in) :: filehandle

        call pepc_status("READ PARAMETERS, section libpepc")
        read(filehandle,NML=libpepc)
    end subroutine libpepc_read_parameters


    !>
    !> writes libpepc-specific parameters to file
    !>
    subroutine libpepc_write_parameters(filehandle)
        implicit none
        integer, intent(in) :: filehandle

        write(filehandle, NML=libpepc)
    end subroutine


    !>
    !> Builds the tree from the given particles, redistributes particles
    !> to other MPI ranks if necessary (i.e. reallocates particles changing size(p))
    !>
    subroutine libpepc_grow_tree(t, p)
      use module_pepc_types, only: t_particle
      use module_tree, only: t_tree
      use module_tree_grow, only: tree_grow
      use module_tree_communicator, only: tree_communicator_start
      use module_debug, only: pepc_status
      use module_interaction_specific, only : calc_force_after_grow
      implicit none

      type(t_tree), intent(inout) :: t
      type(t_particle), allocatable, intent(inout) :: p(:) !< input particle data, initializes %x, %data, %work appropriately (and optionally set %label) before calling this function

      call tree_grow(t, p)
      call tree_communicator_start(t)

      call pepc_status('AFTER GROW: CALC FORCE')
      call calc_force_after_grow(p)
      ! call pepc_status('AFTER GROW: TRAVERSE')
      call pepc_status('AFTER GROW DONE')
    end subroutine


    !>
    !> Frees all resources allocated for the tree `t`
    !>
    subroutine libpepc_timber_tree(t)
      use module_timings
      use module_debug, only : pepc_status
      use module_tree, only: t_tree, tree_destroy
      use module_tree_communicator, only: tree_communicator_stop
      implicit none

      type(t_tree), intent(inout) :: t

      call pepc_status('TIMBER TREE')

      ! deallocate particle and result arrays
      call tree_communicator_stop(t)
      call timer_start(t_deallocate)
      call tree_destroy(t)
      call timer_stop(t_deallocate)
      call timer_stop(t_all)

      call pepc_status('TREE HAS FALLEN')
    end subroutine libpepc_timber_tree


    !>
    !> Traverses the complete tree for the given particles, i.e. computes
    !> the field values at their positions. Although missing information
    !> is automatically requested from remote MPI ranks, it is important
    !> that the particle coordinates fit to the local MPI ranks domain
    !> to avoid excessive communication
    !> If field values on some regular grid are needed, they can be
    !> generated using pepc_prepare_local_grid()
    !> Otherwise, it makes sense to provide the same particles as given/returned
    !> from to pepc_grow_tree()
    !>
    subroutine libpepc_traverse_tree(t, p)
        use module_pepc_types
        use treevars
        use module_walk
        use module_mirror_boxes
        use module_timings
        use module_interaction_specific
        use module_debug
        use module_tree, only: t_tree, tree_allocated
        implicit none

        type(t_tree), target, intent(inout) :: t
        type(t_particle), target, intent(inout) :: p(:) !< input particle data, initializes %x, %data, %work appropriately (and optionally set %label) before calling this function

        integer :: ibox

        call pepc_status('TRAVERSE TREE')
        DEBUG_ASSERT(tree_allocated(t))
        call timer_start(t_walk)

        call tree_walk_init(t, p, num_threads)

        call timer_start(t_fields_passes)
        do ibox = 1,num_neighbour_boxes ! sum over all boxes within ws=1
            ! tree walk finds interaction partners and calls interaction routine for particles on short list
            call tree_walk_run(lattice_vect(neighbour_boxes(:,ibox)))
        end do ! ibox = 1,num_neighbour_boxes
        call timer_stop(t_fields_passes)

        call tree_walk_uninit(t, p)

        ! add lattice contribution
        call timer_start(t_lattice)
        ! add lattice contribution and other per-particle-forces
        ! TODO: do not call calc_force_per_particle here!
        call calc_force_per_particle(p)
        call timer_stop(t_lattice)
        call timer_stop(t_walk)
        call pepc_status('TRAVERSAL DONE')
    end subroutine libpepc_traverse_tree


    !>
    !> Restores the initial particle distribution (before calling pepc_grow_tree() ).
    !>
    subroutine libpepc_restore_particles(t, particles)
        use module_pepc_types, only: t_particle
        use module_timings
        use module_debug, only : pepc_status
        use module_domains, only : domain_restore
        use module_tree, only: t_tree
        implicit none

        type(t_tree), intent(inout) :: t
        type(t_particle), allocatable, intent(inout) :: particles(:) !< input particle data on local MPI rank - is replaced by original particle data that was given before calling pepc_grow_tree()

        call pepc_status('RESTORE PARTICLES')

        ! restore initial particle order specified by calling routine to reassign computed forces
        call timer_start(t_restore)
        call domain_restore(t%decomposition, particles)
        call timer_stop(t_restore)

        call pepc_status('RESTORATION DONE')
    end subroutine libpepc_restore_particles
end module module_libpepc_main
