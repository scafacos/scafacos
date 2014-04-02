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
!>  Interaction specific wrappers to support frontends with old interfaces
!> nothing is obligatory for the treecode here
!>
module module_pepc_wrappers
     use module_pepc_types, only : t_particle
     use module_interaction_specific_types, only : t_particle_data, EMPTY_PARTICLE_RESULTS
     implicit none
     private

      type(t_particle), public, dimension(:), allocatable :: particles

      public pepc_fields_coulomb_wrapper
      public pepc_grid_fields_coulomb_wrapper

      contains

    !>
    !>   Calculate fields and potential for supplied particle coordinates p_x, p_y, p_z and charges p_q
    !>
    !>   Returns fields Ex, Ey, Ez and potential pot excluding external terms
    !>   @param[in] np_local local number of particles
    !>   @param[in] npart_total total particle number
    !>   @param[in] p_x dimension(1:np_local) - x-component of particle coordinates
    !>   @param[in] p_y dimension(1:np_local) - y-component of particle coordinates
    !>   @param[in] p_z dimension(1:np_local) - z-component of particle coordinates
    !>   @param[in] p_q dimension(1:np_local) - particle charge
    !>   @param[inout] p_w dimension(1:np_local) - particle workload from previous iteration (should be set to 1.0 for the first timestep)
    !>   @param[in] p_label dimension(1:np_local) - particle label (may any number except zero)
    !>   @param[out] p_Ex dimension(1:np_local) - x-component of electric field
    !>   @param[out] p_Ey dimension(1:np_local) - y-component of electric field
    !>   @param[out] p_Ez dimension(1:np_local) - z-component of electric field
    !>   @param[out] p_pot dimension(1:np_local) - electric potential
    !>   @param[in] itime current simulation timestep number
    !>   @param[in] weighted selector for load balancing
    !>   @param[in] curve_type selector for type of space filling curve
    !>   @param[in] num_neighbours number of neighbour boxes to be considered during tree walk
    !>   @param[in] neighbours shift vectors to neighbour boxes
    !>   @param[in] no_dealloc if set to .true., deallocation of tree-structures is prevented to allow for front-end triggered diagnostics
    !>
    !>  TODO: update function documentation
    subroutine pepc_fields_coulomb_wrapper(np_local,p_x, p_y, p_z, p_q, p_w, p_label, &
                    p_Ex, p_Ey, p_Ez, p_pot, itime, no_dealloc, force_const)
        use treevars
        use module_pepc
        implicit none
        integer, intent(inout) :: np_local  ! # particles on this CPU
        integer, intent(in) :: itime  ! timestep
        real*8, intent(in), dimension(np_local) :: p_x, p_y, p_z  ! coords and velocities: x1,x2,x3, y1,y2,y3, etc
        real*8, intent(in), dimension(np_local) :: p_q ! charges, masses
        integer, intent(in), dimension(np_local) :: p_label  ! particle label
        real*8, intent(out), dimension(np_local) :: p_ex, p_ey, p_ez, p_pot  ! fields and potential to return
        real*8, dimension(np_local) :: p_w ! work loads
        logical, intent(in) :: no_dealloc
        real*8, intent(in) :: force_const
        
        integer :: i
        
        if (allocated(particles))        deallocate(particles)
        allocate(particles(1:np_local))

        do i=1,np_local
            particles(i) = t_particle( [p_x(i), p_y(i), p_z(i)],       &  ! position
                                              max(p_w(i), 1._8),       &  ! workload from last step
                                                           -1_8,       &  ! key - will be assigned later
                                                           -1_8,       &  ! leaf key - will be assigned later
                                                     p_label(i),       &  ! particle label for tracking purposes
                                       t_particle_data( p_q(i)),       &  ! charge etc
                                       EMPTY_PARTICLE_RESULTS )
        end do

        call pepc_particleresults_clear(particles)
        call pepc_grow_and_traverse(particles, itime=itime, no_dealloc=no_dealloc, no_restore=.false.)
        
        np_local = size(particles)

        ! read data from particle_coordinates, particle_results, particle_properties
        do i=1,np_local
          p_ex(i)  = force_const*particles(i)%results%e(1)
          p_ey(i)  = force_const*particles(i)%results%e(2)
          p_ez(i)  = force_const*particles(i)%results%e(3)
          p_pot(i) = force_const*particles(i)%results%pot
          p_w(i)   = particles(i)%work
        end do

    end subroutine



    subroutine pepc_grid_fields_coulomb_wrapper(ngp,p_x, p_y, p_z, p_label, p_Ex, p_Ey, p_Ez, p_pot, &
                              num_neighbour_boxes, neighbour_boxes, force_const)
      use treevars, only : me
      use module_pepc
      implicit none
      integer, intent(in) :: ngp
      real*8, intent(in) :: p_x(ngp), p_y(ngp), p_z(ngp)
      integer, intent(in) :: p_label(ngp)
      real*8, intent(out) :: p_Ex(ngp), p_Ey(ngp), p_Ez(ngp), p_pot(ngp)
      integer, intent(in) :: num_neighbour_boxes !< number of shift vectors in neighbours list (must be at least 1 since [0, 0, 0] has to be inside the list)
      integer, intent(in) :: neighbour_boxes(3, num_neighbour_boxes) ! list with shift vectors to neighbour boxes that shall be included in interaction calculation, at least [0, 0, 0] should be inside this list
      real*8, intent(in) :: force_const

      type(t_particle), dimension(:), allocatable :: grid_particles

      integer :: i

      allocate(grid_particles(ngp))

      do i=1,ngp
        grid_particles(i) = t_particle( [p_x(i), p_y(i), p_z(i)],       &  ! position
                                                              1.,       &  ! workload from last step
                                                            -1_8,       &  ! key - will be assigned later
                                                            -1_8,       &  ! leaf key - will be assigned later
                                                      p_label(i),       &  ! particle label for tracking purposes
                                          t_particle_data( 0.0 ),       &  ! charge etc
                                       EMPTY_PARTICLE_RESULTS )
      end do

      call pepc_particleresults_clear(grid_particles)

      call pepc_traverse_tree(grid_particles)

      do i=1,ngp
        p_ex(i)  = force_const*grid_particles(i)%results%e(1)
        p_ey(i)  = force_const*grid_particles(i)%results%e(2)
        p_ez(i)  = force_const*grid_particles(i)%results%e(3)
        p_pot(i) = force_const*grid_particles(i)%results%pot
      end do

      deallocate(grid_particles)
    end subroutine
end module module_pepc_wrappers
