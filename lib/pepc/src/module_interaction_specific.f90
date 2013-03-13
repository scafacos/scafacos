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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!> Encapsulates anything that is directly involved in force calculation
!> and multipole expansion manipulation
!> i.e. shifting along the tree, computing forces between particles and cluster, etc.
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module module_interaction_specific
     use module_pepc_types
     use module_interaction_specific_types
     implicit none
     save
     private

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  public variable declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      real*8, parameter :: WORKLOAD_PENALTY_MAC  = 1._8 !< TODO: currently unused
      real*8, parameter :: WORKLOAD_PENALTY_INTERACTION = 30._8

      integer, public :: force_law    = 3      !< 3 = 3D-Coulomb, 2 = 2D-Coulomb
      integer, public :: mac_select   = 0      !< selector for multipole acceptance criterion, mac_select==0: Barnes-Hut
      logical, public :: include_far_field_if_periodic = .true. !< if set to false, the far-field contribution to periodic boundaries is ignored (aka 'minimum-image-mode')
      real*8, public  :: theta2       = 0.36  !< square of multipole opening angle
      real*8, public  :: eps2         = 0.0    !< square of short-distance cutoff parameter for plummer potential (0.0 corresponds to classical Coulomb)
      real*8, public  :: kelbg_invsqrttemp = 0.0 !< inverse square root of temperature for kelbg potential

! CS DEBUG KRAM FÜR INTERACTION PARTNER
      integer*8, allocatable,public :: interaction_keylist(:,:)
      integer, allocatable,public :: no_interaction_partners(:)
      real*8, allocatable,public :: interaction_vbox(:,:,:)
! ENDE CS
      namelist /calc_force_coulomb/ force_law, mac_select, include_far_field_if_periodic, theta2, eps2, kelbg_invsqrttemp


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  public subroutine declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! currently, all public functions in module_interaction_specific are obligatory
      public multipole_from_particle
      public shift_multipoles_up
      public results_add
      public calc_force_per_interaction
      public calc_force_per_particle
      public mac
      public particleresults_clear
      public calc_force_read_parameters
      public calc_force_write_parameters
      public calc_force_finalize
      public calc_force_prepare
      public calc_force_after_grow
      public get_number_of_interactions_per_particle

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
      !> Computes multipole properties of a single particle
      !>
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine multipole_from_particle(particle_pos, particle, multipole)
        implicit none
        real*8, intent(in) :: particle_pos(3)
        type(t_particle_data), intent(in) :: particle
        type(t_tree_node_interaction_data), intent(out) :: multipole

        multipole = t_tree_node_interaction_data(particle_pos, &
                                     particle%q,   &
                                 abs(particle%q),  &
                                     [0., 0., 0.], &
                                     [0., 0., 0.], &
                                      0., 0., 0., 0. )
      end subroutine


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !>
      !> Accumulates multipole properties of child nodes to parent node
      !>
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine shift_multipoles_up(parent, children)
        implicit none
        type(t_tree_node_interaction_data), intent(out) :: parent
        type(t_tree_node_interaction_data), intent(in) :: children(:)

        integer :: nchild, j

        real*8 :: shift(1:3)

        nchild = size(children)

        parent%charge     = SUM( children(1:nchild)%charge )
        parent%abs_charge = SUM( children(1:nchild)%abs_charge )

        ! centre of charge
        parent%coc        = [0., 0., 0.]

        if (parent%abs_charge .ne. 0.) then
          ! use center-of-charge because we may divide by abs_charge
          do j=1,nchild
            parent%coc(1:3) = parent%coc(1:3) + ( children(j)%coc(1:3) * children(j)%abs_charge )
          end do

          parent%coc(1:3) = parent%coc(1:3) / parent%abs_charge
        else
          ! use geometric center
          do j=1,nchild
            parent%coc(1:3) = parent%coc(1:3) +   children(j)%coc(1:3)
          end do

         parent%coc(1:3) = parent%coc(1:3) / nchild
        endif

        ! multipole properties
        parent%dip    = [0., 0., 0.]
        parent%quad   = [0., 0., 0.]
        parent%xyquad = 0.
        parent%yzquad = 0.
        parent%zxquad = 0.

        do j=1,nchild
          shift(1:3) = parent%coc(1:3) - children(j)%coc

          ! dipole moment
          parent%dip = parent%dip + children(j)%dip - children(j)%charge*shift(1:3)

          ! quadrupole moment
          parent%quad(1:3) = parent%quad(1:3) + children(j)%quad(1:3) - 2*children(j)%dip(1:3)*shift(1:3) + children(j)%charge*shift(1:3)**2

          parent%xyquad = parent%xyquad + children(j)%xyquad - children(j)%dip(1)*shift(2) - children(j)%dip(2)*shift(1) + children(j)%charge*shift(1)*shift(2)
          parent%yzquad = parent%yzquad + children(j)%yzquad - children(j)%dip(2)*shift(3) - children(j)%dip(3)*shift(2) + children(j)%charge*shift(2)*shift(3)
          parent%zxquad = parent%zxquad + children(j)%zxquad - children(j)%dip(3)*shift(1) - children(j)%dip(1)*shift(3) + children(j)%charge*shift(3)*shift(1)
        end do

        parent%bmax = maxval(sqrt((parent%coc(1)-children(1:nchild)%coc(1))**2+(parent%coc(2)-children(1:nchild)%coc(2))**2+(parent%coc(3)-children(1:nchild)%coc(3))**2) + children(1:nchild)%bmax)

      end subroutine


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !>
      !> adds res2 to res1
      !>
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine results_add(res1, res2)
        implicit none
        type(t_particle_results), intent(inout) :: res1
        type(t_particle_results), intent(in) :: res2

        res1%e    = res1%e    + res2%e
        res1%pot  = res1%pot  + res2%pot
      end subroutine


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !>
      !> reads interaction-specific parameters from file
      !>
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine calc_force_read_parameters(filehandle)
        use module_debug, only: pepc_status
        implicit none
        integer, intent(in) :: filehandle

        call pepc_status("READ PARAMETERS, section calc_force_coulomb")
        read(filehandle, NML=calc_force_coulomb)

      end subroutine

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !>
      !> writes interaction-specific parameters to file
      !>
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine calc_force_write_parameters(filehandle)
        use module_debug, only: pepc_status
        implicit none
        integer, intent(in) :: filehandle

        write(filehandle, NML=calc_force_coulomb)

      end subroutine


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !>
      !> computes derived parameters for calc force module
      !>
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine calc_force_prepare()
        use treevars, only : me, MPI_COMM_lpepc
        use module_fmm_framework, only : fmm_framework_init
        implicit none

        call fmm_framework_init(me, MPI_COMM_lpepc)

      end subroutine


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !>
      !> initializes static variables of calc force module that depend 
      !> on particle data and might be reused on subsequent traversals
      !>
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine calc_force_after_grow(particles, nparticles)
        use module_pepc_types
        use module_fmm_framework, only : fmm_framework_timestep
        use module_mirror_boxes, only : do_periodic
        implicit none
        type(t_particle), dimension(:), intent(in) :: particles
        integer, intent(in) :: nparticles

        ! calculate spherical multipole expansion of central box
        ! this cannot be done in calc_force_per_particle() since there, possibly
        ! other particles are used than we need for the multipoles
        ! e.g. in the case of a second traverse for test/grid particles
        if ((do_periodic) .and. (include_far_field_if_periodic)) then
          call fmm_framework_timestep(particles, nparticles)
        end if

      end subroutine      


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !>
      !> subroutine must return the estimated number of iteractions per particle
      !> for the current mac and/or parameters and the supplied total number of particles
      !>
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine get_number_of_interactions_per_particle(npart_total, nintmax)
        implicit none
        integer, intent(in) :: npart_total !< total number of particles
        integer, intent(out) :: nintmax !< maximum number of interactions per particle

        real*8 :: invnintmax !< inverse of nintmax to avoid division by zero for theta == 0.0

        ! Estimate of interaction list length - Hernquist expression
        ! applies for BH-MAC
        invnintmax = max(theta2 / (35.*log(1.*npart_total)) , 1._8/npart_total)
        nintmax    = int(1._8/invnintmax)

      end subroutine



      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !>
      !> finalizes the calc force module at end of simulation
      !>
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine calc_force_finalize()
        implicit none
        ! nothing to do here
      end subroutine calc_force_finalize


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !>
      !> generic Multipole Acceptance Criterion
      !>
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      function mac(particle, node, dist2, boxlength2)

        use treevars, only : tree_nodes
        implicit none

        logical :: mac
        integer, intent(in) :: node
        type(t_particle), intent(in) :: particle
        real*8, intent(in) :: dist2
        real*8, intent(in) :: boxlength2

        select case (mac_select)
            case (0)
              ! Barnes-Hut-MAC
              mac = (theta2 * dist2 > boxlength2)
            case (1)
               ! Bmax-MAC
              mac = (theta2 * dist2 > min(tree_nodes(node)%bmax**2,3.0*boxlength2)) !TODO: Can we put the min into bmax itself? And **2?
            case default
              ! N^2 code
              mac = .false.
        end select

      end function

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !>
      !> clears result in t_particle datatype - usually, this function does not need to be touched
      !> due to dependency on module_pepc_types and(!) on module_interaction_specific, the
      !> function cannot reside in module_interaction_specific that may not include module_pepc_types
      !>
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine particleresults_clear(particles, nparticles)
        use module_pepc_types
        implicit none
        type(t_particle), intent(inout) :: particles(nparticles)
        integer, intent(in) :: nparticles

        particles(1:nparticles)%results = EMPTY_PARTICLE_RESULTS

      end subroutine


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> Force calculation wrapper.
        !> This function is thought for pre- and postprocessing of
        !> calculated fields, and for being able to call several
        !> (different) force calculation routines
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine calc_force_per_interaction(particle, node, key, delta, dist2, vbox, node_is_leaf)
          use module_pepc_types
          use treevars
          implicit none

          type(t_tree_node_interaction_data), intent(in) :: node
          integer*8, intent(in) :: key
          type(t_particle), intent(inout) :: particle
          logical, intent(in) :: node_is_leaf
          real*8, intent(in) :: vbox(3), delta(3), dist2


          real*8 :: exyz(3), phic

          select case (force_law)
            case (2)  !  compute 2D-Coulomb fields and potential of particle p from its interaction list

                if (node_is_leaf) then
                    ! It's a leaf, do direct summation
                    call calc_force_coulomb_2D_direct(node, delta(1:2), dot_product(delta(1:2), delta(1:2)), exyz(1), exyz(2),phic)
                else
                    ! It's a twig, do ME
                    call calc_force_coulomb_2D(node, delta(1:2), dot_product(delta(1:2), delta(1:2)), exyz(1), exyz(2),phic)
                end if
                exyz(3) = 0.

            case (3)  !  compute 3D-Coulomb fields and potential of particle p from its interaction list

                if (node_is_leaf) then
                    ! It's a leaf, do direct summation
                    call calc_force_coulomb_3D_direct(node, delta, dist2, exyz(1), exyz(2), exyz(3), phic)
                else
                    ! It's a twig, do ME
                    call calc_force_coulomb_3D(node, delta, dist2, exyz(1), exyz(2), exyz(3), phic)
                end if

            case (4)  ! LJ potential for quiet start
                call calc_force_LJ(node, delta, dist2, exyz(1), exyz(2), exyz(3), phic)
                exyz(3) = 0.

            case (5)  !  compute 3D-Coulomb fields and potential for particle-cluster interaction
                      !  and Kelbg for particle-particle interaction

                if (node_is_leaf) then
                    ! It's a leaf, do direct summation with kelbg
                    call calc_force_kelbg_3D_direct(particle, node, delta, dist2, exyz(1), exyz(2), exyz(3), phic)
                else
                    ! It's a twig, do ME with coulomb
                    call calc_force_coulomb_3D(node, delta, dist2, exyz(1), exyz(2), exyz(3), phic)
                end if

! START CHRISTAN SALMAGNE; ADDED FOR DEBUGGING

            case (6)  !  used to save interaction partners
                no_interaction_partners(particle%label)=no_interaction_partners(particle%label)+1
                interaction_keylist(particle%label,no_interaction_partners(particle%label))=key
                interaction_vbox(particle%label,no_interaction_partners(particle%label),1:3)=vbox(1:3)


! END CS


            case default
              exyz = 0.
              phic = 0.
          end select

          particle%results%e         = particle%results%e    + exyz
          particle%results%pot       = particle%results%pot  + phic
          particle%work              = particle%work         + WORKLOAD_PENALTY_INTERACTION

        end subroutine calc_force_per_interaction

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> Force calculation wrapper for contributions that only have
        !> to be added once per particle
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine calc_force_per_particle(particles, nparticles)
          use treevars, only: num_threads
          use module_debug, only : pepc_status
          use module_pepc_types
          use treevars, only : me
          use module_fmm_framework
          use module_mirror_boxes
          implicit none

          integer, intent(in) :: nparticles
          type(t_particle), intent(inout) :: particles(:)
          real*8 :: e_lattice(3), phi_lattice
          integer :: p

          call pepc_status('CALC FORCE PER PARTICLE')

          potfarfield  = 0.
          potnearfield = 0.

          if ((do_periodic) .and. (include_far_field_if_periodic)) then

             if ((me==0) .and. (force_law .ne. 3)) write(*,*) "Warning: far-field lattice contribution is currently only supported for force_law==3"
          !$ call omp_set_num_threads(num_threads)
          !$OMP  PARALLEL DO DEFAULT(PRIVATE) SHARED(particles) SCHEDULE(RUNTIME) REDUCTION(+:potfarfield,potnearfield)
             do p=1,nparticles
                call fmm_sum_lattice_force(particles(p)%x, e_lattice, phi_lattice)

                potfarfield  = potfarfield  + phi_lattice               * particles(p)%data%q
                potnearfield = potnearfield + particles(p)%results%pot  * particles(p)%data%q

                particles(p)%results%e     = particles(p)%results%e     + e_lattice
                particles(p)%results%pot   = particles(p)%results%pot   +  phi_lattice
             end do
          !$OMP  END PARALLEL DO
          !$ call omp_set_num_threads(1)

          end if

          call pepc_status('CALC FORCE PER PARTICLE DONE')

        end subroutine calc_force_per_particle


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> Calculates 3D Coulomb interaction of particle p with tree node inode
        !> that is shifted by the lattice vector vbox
        !> results are returned in eps, sumfx, sumfy, sumfz, sumphi
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine calc_force_coulomb_3D(t, d, dist2, sumfx, sumfy, sumfz, sumphi)
          use module_pepc_types
          use treevars
          implicit none

          type(t_tree_node_interaction_data), intent(in) :: t !< index of particle to interact with
          real*8, intent(in) :: d(3), dist2 !< separation vector and magnitude**2 precomputed in walk_single_particle
          real*8, intent(out) ::  sumfx,sumfy,sumfz,sumphi

          real*8 :: rd,dx,dy,dz,r,dx2,dy2,dz2
          real*8 :: dx3,dy3,dz3,rd3,rd5,rd7,fd1,fd2,fd3,fd4,fd5,fd6

             sumfx  = 0.
             sumfy  = 0.
             sumfz  = 0.
             sumphi = 0.

             !  preprocess distances
             dx = d(1)
             dy = d(2)
             dz = d(3)


             r = sqrt(dist2+eps2)
             rd = 1./r
             rd3 = rd**3
             rd5 = rd**5
             rd7 = rd**7

             dx2 = dx**2
             dy2 = dy**2
             dz2 = dz**2
             dx3 = dx**3
             dy3 = dy**3
             dz3 = dz**3

             fd1 = 3.*dx2*rd5 - rd3
             fd2 = 3.*dy2*rd5 - rd3
             fd3 = 3.*dz2*rd5 - rd3
             fd4 = 3.*dx*dy*rd5
             fd5 = 3.*dy*dz*rd5
             fd6 = 3.*dx*dz*rd5

             ! potential

             sumphi = sumphi + t%charge*rd    &                           !  monopole term
                                        !
                  + (dx*t%dip(1) + dy*t%dip(2) + dz*t%dip(3))*rd3  &    !  dipole
                                        !     Dx             Dy            Dz
                  + 0.5*fd1*t%quad(1) + 0.5*fd2*t%quad(2) + 0.5*fd3*t%quad(3)  &  !  quadrupole
                                        !           Qxx                 Qyy                 Qzz
                  + fd4*t%xyquad + fd5*t%yzquad + fd6*t%zxquad
             !   Qxy            Qyz             Qzx

             !  forces

             sumfx = sumfx + t%charge*dx*rd3 &      ! monopole term
                                        !
                  + fd1*t%dip(1) + fd4*t%dip(2) + fd6*t%dip(3)   &   !  dipole term
                                        !
                  + (15.*dx3*rd7 - 9.*dx*rd5 )*0.5*t%quad(1) &     !
                  + ( 15.*dy*dx2*rd7 - 3.*dy*rd5 )*t%xyquad &     !
                  + ( 15.*dz*dx2*rd7 - 3.*dz*rd5 )*t%zxquad &     !   quadrupole term
                  + ( 15*dx*dy*dz*rd7 )*t%yzquad &                !
                  + ( 15.*dx*dy2*rd7 - 3.*dx*rd5 )*0.5*t%quad(2) & !
                  + ( 15.*dx*dz2*rd7 - 3.*dx*rd5 )*0.5*t%quad(3)   !

             sumfy = sumfy + t%charge*dy*rd3 &
                  + fd2*t%dip(2) + fd4*t%dip(1) + fd5*t%dip(3)  &
                  + ( 15.*dy3*rd7 - 9.*dy*rd5 )*0.5*t%quad(2) &
                  + ( 15.*dx*dy2*rd7 - 3.*dx*rd5 )*t%xyquad &
                  + ( 15.*dz*dy2*rd7 - 3.*dz*rd5 )*t%yzquad &
                  + ( 15.*dx*dy*dz*rd7 )*t%zxquad &
                  + ( 15.*dy*dx2*rd7 - 3.*dy*rd5 )*0.5*t%quad(1) &
                  + ( 15.*dy*dz2*rd7 - 3.*dy*rd5 )*0.5*t%quad(3)

             sumfz = sumfz + t%charge*dz*rd3 &
                  + fd3*t%dip(3) + fd5*t%dip(2) + fd6*t%dip(1)  &
                  + ( 15.*dz3*rd7 - 9.*dz*rd5 )*0.5*t%quad(3) &
                  + ( 15.*dx*dz2*rd7 - 3.*dx*rd5 )*t%zxquad &
                  + ( 15.*dy*dz2*rd7 - 3.*dy*rd5 )*t%yzquad &
                  + ( 15.*dx*dy*dz*rd7 )*t%xyquad &
                  + ( 15.*dz*dy2*rd7 - 3.*dz*rd5 )*0.5*t%quad(2) &
                  + ( 15.*dz*dx2*rd7 - 3.*dz*rd5 )*0.5*t%quad(1)

        end subroutine calc_force_coulomb_3D


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> Calculates 2D Coulomb interaction of particle p with tree node inode
        !> that is shifted by the lattice vector vbox
        !> results are returned in eps, sumfx, sumfy, sumphi 
        !> Unregularized force law is: 
        !>   Phi = -2q log R 
        !>   Ex = -dPhi/dx = 2 q x/R^2 etc 
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine calc_force_coulomb_2D(t, d, dist2, sumfx, sumfy, sumphi)
          use module_pepc_types
          use treevars
          implicit none

          type(t_tree_node_interaction_data), intent(in) :: t !< index of particle to interact with
          real*8, intent(in) :: d(2), dist2 !< separation vector and magnitude**2 precomputed in walk_single_particle
          real*8, intent(out) ::  sumfx,sumfy,sumphi

          real*8 :: dx,dy,d2,rd2,rd4,rd6,dx2,dy2,dx3,dy3

          sumfx  = 0.
          sumfy  = 0.
          sumphi = 0.

          !  preprocess distances and reciprocals
          dx = d(1)
          dy = d(2)

          d2  = dist2+eps2
          rd2 = 1./d2 
          rd4 = rd2**2 
          rd6 = rd2**3 
          dx2 = dx**2 
          dy2 = dy**2 
          dx3 = dx**3 
          dy3 = dy**3 
	  
          sumphi = sumphi - 0.5*t%charge*log(d2)    &                           !  monopole term 
               ! 
               + (dx*t%dip(1) + dy*t%dip(2) )*rd2  &    !  dipole
               !                               
               + 0.5*t%quad(1)*(dx2*rd4 - rd2) + 0.5*t%quad(2)*(dy2*rd4 - rd2) + t%xyquad*dx*dy*rd4  !  quadrupole
          
          sumfx = sumfx + t%charge*dx*rd2  &   ! monopole 
               ! 
               + t%dip(1)*(2*dx2*rd4 - rd2) + t%dip(2)*2*dx*dy*rd4  &  ! dipole
               ! 
               + 0.5*t%quad(1)*(8*dx3*rd6 - 6*dx*rd4) &                    ! quadrupole
               + 0.5*t%quad(2)*(8*dx*dy**2*rd6 - 2*dx*rd4) &
               +     t%xyquad*(8*dx2*dy*rd6 - 2*dy*rd4) 
          
          sumfy = sumfy + t%charge*dy*rd2  &   ! monopole 
               ! 
               + t%dip(2)*(2*dy2*rd4 - rd2) + t%dip(1)*2*dx*dy*rd4  &  ! dipole
               ! 
               + 0.5*t%quad(2)*(8*dy3*rd6 - 6*dy*rd4) &                    ! quadrupole
               + 0.5*t%quad(1)*(8*dy*dx**2*rd6 - 2*dy*rd4) &
               +     t%xyquad*(8*dy2*dx*rd6 - 2*dx*rd4) 

        end subroutine calc_force_coulomb_2D
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> CALC_FORCE_LJ
        !>
        !> Calculates 3D Lennard-Jones interaction of particle p with tree node inode
        !> shifted by the lattice vector vbox
        !> results are returned sumfx, sumfy, sumfz, sumphi
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine calc_force_LJ(t, d, dist2, sumfx, sumfy, sumfz, sumphi)
          use module_pepc_types
          use treevars
          implicit none

          type(t_tree_node_interaction_data), intent(in) :: t !< index of particle to interact with
          real*8, intent(in) :: d(3), dist2 !< separation vector and magnitude**2 precomputed in walk_single_particle
          real*8, intent(out) ::  sumfx,sumfy,sumfz,sumphi
          real*8 :: dx,dy,dz,r2
          real*8 :: flj, epsc2, plj, aii2, aii2_r2, r

          sumfx  = 0.
          sumfy  = 0.
          sumfz  = 0.
          sumphi = 0.

          !  preprocess distances
          dx  = d(1)
          dy  = d(2)
          dz  = d(3)
          r2 = dist2

          !    epsc should be > a_ii to get evenly spaced ions
          aii2  = eps2
          epsc2 = 0.8*aii2
          plj   = 0.

          ! Force is repulsive up to and just beyond aii
          if (r2 > epsc2) then
              aii2_r2 = aii2/r2
          else
              aii2_r2 = aii2/epsc2
          endif

          flj = 2.*(aii2_r2)**4 - 1.*(aii2_r2  )**2

          ! potential
          sumphi = sumphi + plj

          !  forces
          r     = sqrt(r2)
          sumfx = sumfx + dx/r*flj
          sumfy = sumfy + dy/r*flj
          !    	  sumfz = sumfz + dz/r*flj
          sumfz=0.

      end subroutine calc_force_LJ

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !>
      !> Calculates 3D Coulomb interaction of particle p with particle inode
      !> that is shifted by the lattice vector vbox
      !> results are returned in eps, sumfx, sumfy, sumfz, sumphi
      !>
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine calc_force_coulomb_3D_direct(t, d, dist2, sumfx, sumfy, sumfz, sumphi)
          use module_pepc_types
          use treevars
          implicit none

          type(t_tree_node_interaction_data), intent(in) :: t !< index of particle to interact with
          real*8, intent(in) :: d(3), dist2 !< separation vector and magnitude**2 precomputed in walk_single_particle
          real*8, intent(out) ::  sumfx,sumfy,sumfz,sumphi

          real*8 :: rd,dx,dy,dz,r,charge, rd3

          !  preprocess distances
          dx = d(1)
          dy = d(2)
          dz = d(3)

          r = sqrt(dist2+eps2)
          rd = 1./r
          rd3 = rd**3

          charge = t%charge

          ! potential
          sumphi = charge*rd

          !  forces

          sumfx = charge*dx*rd3

          sumfy = charge*dy*rd3

          sumfz = charge*dz*rd3

      end subroutine calc_force_coulomb_3D_direct


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !>
      !> Calculates 3D Kelbg interaction of particle p with particle inode
      !> that is shifted by the lattice vector vbox
      !> results are returned in eps, sumfx, sumfy, sumfz, sumphi
      !>
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine calc_force_kelbg_3D_direct(particle, t, d, dist2, sumfx, sumfy, sumfz, sumphi)
          use module_pepc_types
          use treevars
          implicit none

          type(t_particle), intent(inout) :: particle
          type(t_tree_node_interaction_data), intent(in) :: t !< index of particle to interact with
          real*8, intent(in) :: d(3), dist2 !< separation vector and magnitude**2 precomputed in walk_single_particle
          real*8, intent(out) ::  sumfx,sumfy,sumfz,sumphi
          real*8 :: rd,r,rd3
          real*8, parameter :: sqrtpi = sqrt(acos(-1.0_8))
          real*8 :: ome, rol, lambda, q, fprefac

          q           = t%charge

          ! TODO: lambda must be adjusted depending on mass and temperature of interacting partners - currently it is fixed for electron-proton interactions
          if (particle%data%q * q < 0.) then
            ! e-i or i-e interaction
            lambda = 1.00027227_8 * kelbg_invsqrttemp
          else
            if ( q > 0. ) then
              ! i-i interaction
              lambda = 0.03300355_8 * kelbg_invsqrttemp
            else
              ! e-e interaction
              lambda = 1.41421356_8 * kelbg_invsqrttemp
            endif
          endif

          r   = sqrt(dist2)
          rd  = 1. / r
          rd3 = rd**3
          rol = r  / lambda        !< "r over lambda"
          ome = 1  - exp(-rol*rol) !< "one minus exp(stuff)"

          ! potential
          sumphi  = q * rd  * (ome + sqrtpi*rol*(1-erf(rol)))
          !  forces
          fprefac = q * rd3 * ome

          sumfx = fprefac*d(1)
          sumfy = fprefac*d(2)
          sumfz = fprefac*d(3)

      end subroutine calc_force_kelbg_3D_direct


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !>
      !> Calculates 2D Coulomb interaction of particle p with tree node inode
      !> that is shifted by the lattice vector vbox
      !> results are returned in eps, sumfx, sumfy, sumphi
      !> Unregularized force law is:
      !>   Phi = -2q log R
      !>   Ex = -dPhi/dx = 2 q x/R^2 etc
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine calc_force_coulomb_2D_direct(t, d, dist2, sumfx, sumfy, sumphi)
          use module_pepc_types
          use treevars
          implicit none

          type(t_tree_node_interaction_data), intent(in) :: t !< index of particle to interact with
          real*8, intent(in) :: d(2), dist2 !< separation vector and magnitude**2 precomputed in walk_single_particle
          real*8, intent(out) ::  sumfx,sumfy,sumphi

          real*8 :: dx,dy,d2,rd2,charge

          !  preprocess distances and reciprocals
          dx = d(1)
          dy = d(2)

          d2  = dist2+eps2
          rd2 = 1./d2


          charge = t%charge

          sumphi = - 0.5*charge*log(d2)

          sumfx = charge*dx*rd2

          sumfy = charge*dy*rd2

      end subroutine calc_force_coulomb_2D_direct


  end module module_interaction_specific
