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


! ==============================================================
!
!
!                  PEPC-S
!
!    Parallel Efficient Parallel Coulomb-solver: Single Call Version
!
!  ==============================================================
#include "fcs_fconfig.h"

module module_pepcs

  contains

    subroutine pepc(local_particles, local_max_particles, &
                          positions, charges,  &
                          field, potentials,   &
                          virial_tensor,       &
                          lat_x, lat_y, lat_z, &
                          lat_period,          &
                          lat_corr,            &
                          eps, theta, db_level) bind(c,name='pepc')
        use module_pepc_types
        use module_fmm_framework
        use module_mirror_boxes
        use module_pepc
        use treevars
        use module_interaction_specific, only : theta2, eps2, mac_select, force_law
        use module_debug, only : debug_level
        use, intrinsic :: iso_c_binding
        implicit none
        include 'mpif.h'

        fcs_integer, intent(inout) :: local_particles, local_max_particles
        fcs_real, intent(in) :: positions(3,local_max_particles), charges(local_max_particles)
        fcs_real, intent(out) :: field(3,local_max_particles), potentials(local_max_particles)
        fcs_real, intent(out), dimension(3,3) :: virial_tensor
        fcs_real, intent(in), dimension(3) :: lat_x, lat_y, lat_z
        fcs_integer, intent(in), dimension(3) :: lat_period
        fcs_integer, intent(in) :: lat_corr
        fcs_real, intent(in) :: theta, eps
        fcs_integer, intent(in) :: db_level

        integer, parameter :: itime = 1
        real, parameter :: np_mult_ = -45

        type(t_particle), allocatable :: particles(:)
        integer :: i

        np_mult = np_mult_

        allocate(particles(local_particles))

        ! copy coordinates and charges to internal data structures
        do i=1,local_particles
            particles(i)%x(1:3) = positions(1:3,i)
            particles(i)%work   = 1._8
            particles(i)%label  = i ! this is just for debugging purposes
            particles(i)%data%q = charges(i)
        end do

        ! initialize calc force params
        theta2      = theta**2
        mac_select  = 0 ! Barnes-Hut MAC, everything else leads to an N^2 code
        eps2        = eps**2
        force_law   = 3 ! 3D-Coulomb

        ! =============================================================
        ! TODO: do this in some scafacos_init function instead
        ! of in every timestep
        ! initialize framework for lattice contributions (is automatically ignored if periodicity = [false, false, false]
        t_lattice_1 = lat_x
        t_lattice_2 = lat_y
        t_lattice_3 = lat_z
        periodicity(1:3) = (lat_period(1:3) == 1)
        fmm_extrinsic_correction = lat_corr
        debug_level = db_level
        call pepc_prepare(3_kind_dim)
        ! =============================================================

        call pepc_particleresults_clear(particles)
        call pepc_grow_and_traverse(particles, itime=itime)

        ! read fields and potentials from internal data structures
        do i=1,local_particles
            field(1:3,i)  = particles(i)%results%e
            potentials(i) = particles(i)%results%pot
        end do

        virial_tensor = 0. !TODO
        if (me==0) write(*,*) "TODO: Virial unsupported in this version of pepc"

        deallocate(particles)

    end subroutine

end module
