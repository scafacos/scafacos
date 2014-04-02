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
!> Defines a derived type that represents boxes and associated procedures.
!>
module module_box
  implicit none
  private

  logical, public :: force_cubic_domain = .false. !< if set to .true., pepc uses an overall cubic enclosure of the particle cloud instead of the cuboid (closer) one

  !>
  !> Derived type that represents a rectangular cuboid.
  !>
  type, public :: t_box
    real*8 :: boxmin(3) !< front lower left corner of the box
    real*8 :: boxmax(3) !< back upper right corner of the box
    real*8 :: boxsize(3) !< the box diagonal connecting `boxmin` to `boxmax`
  end type t_box

  public :: box_create

  contains

  !>
  !> Determines the bounding box `b` that contains the coordinates of all
  !> particles in `p` on all ranks in the communication environment `c`.
  !>
  subroutine box_create(b, p, c)
    use module_comm_env, only: t_comm_env
    use module_pepc_types
    use module_debug
    implicit none
    include 'mpif.h'

    type(t_box), intent(out) :: b !< the bounding box that contains all `p`
    type(t_particle), intent(in) :: p(:) !< particles to embed in the bounding box
    type(t_comm_env), intent(in) :: c !< communication environment

    integer(kind_default) :: ierr
    real*8 :: min_local(3), max_local(3)

    ! Find limits of local simulation region
    min_local(1) = minval(p(:)%x(1))
    max_local(1) = maxval(p(:)%x(1))
    min_local(2) = minval(p(:)%x(2))
    max_local(2) = maxval(p(:)%x(2))
    min_local(3) = minval(p(:)%x(3))
    max_local(3) = maxval(p(:)%x(3))

    if (force_cubic_domain) then
      min_local(1:3) = minval(min_local(1:3))
      max_local(1:3) = maxval(max_local(1:3))
    endif

    ! Find global limits
    call MPI_ALLREDUCE(min_local, b%boxmin, 3, MPI_REAL8, MPI_MIN, c%comm, ierr)
    call MPI_ALLREDUCE(max_local, b%boxmax, 3, MPI_REAL8, MPI_MAX, c%comm, ierr)

    ! Safety margin - put buffer region around particles
    b%boxsize = b%boxmax - b%boxmin
    b%boxmin  = b%boxmin - b%boxsize / 10000.0
    b%boxmax  = b%boxmax + b%boxsize / 10000.0
    b%boxsize = b%boxmax - b%boxmin

    if (dbg(DBG_DOMAIN)) then
      DEBUG_WARNING('(6(a15,f12.4/),a15,3f12.4)',
      'xmin = ',b%boxmin(1),'xmax = ',b%boxmax(1),
      'ymin = ',b%boxmin(2),'ymax = ',b%boxmax(2),
      'zmin = ',b%boxmin(3),'zmax = ',b%boxmax(3),
      'boxsize = ',b%boxsize )
    end if
  end subroutine box_create

end module module_box

