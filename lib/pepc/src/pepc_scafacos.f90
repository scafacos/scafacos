#include <config.h>

subroutine pepc_scafacos_initialize(communicator) bind(c)

  use module_pepc

  implicit none

  integer, intent(inout) :: communicator
  integer :: my_rank, n_ranks

  call pepc_initialize("pepc-scafacos", my_rank, n_ranks, .false., comm=communicator)

end subroutine pepc_scafacos_initialize

subroutine pepc_scafacos_finalize(communicator) bind(c)

  use module_pepc

  implicit none

  integer, intent(inout) :: communicator
  integer :: my_rank, n_ranks

  call pepc_finalize(communicator)

end subroutine pepc_scafacos_finalize

subroutine pepc_scafacos_run(nlocal, ntotal, positions, charges, &
  efield, potentials, work, virial, box_a, box_b, box_c, periodicity_in, &
  lattice_corr, eps, theta, db_level, nwt, npm) bind(c)

  use iso_c_binding

  use module_pepc
  use module_walk, only : max_particles_per_thread
  use module_pepc_types
  use module_interaction_specific, only : theta2, eps2
  use module_mirror_boxes, only : t_lattice_1, t_lattice_2, t_lattice_3, periodicity
  use module_fmm_framework, only : fmm_extrinsic_correction
  use module_debug, only : debug_level
  use treevars, only : np_mult, num_threads

  implicit none

  !!! passed variables
  integer(kind = fcs_integer_kind_isoc), intent(inout) :: nlocal, ntotal
  real(kind = fcs_real_kind_isoc),       intent(in)    :: positions(3*nlocal), charges(nlocal)
  real(kind = fcs_real_kind_isoc),       intent(inout) :: efield(3*nlocal), potentials(nlocal), work(nlocal)
  real(kind = fcs_real_kind_isoc),       intent(inout) :: virial(9)
  real(kind = fcs_real_kind_isoc),       intent(in)    :: box_a(3), box_b(3), box_c(3)
  integer(kind = fcs_integer_kind_isoc), intent(in)    :: periodicity_in(3), lattice_corr
  real(kind = fcs_real_kind_isoc),       intent(in)    :: eps, theta, npm
  integer(kind = fcs_integer_kind_isoc), intent(in)    :: db_level, nwt

  !!! pepc internal variables
  type(t_particle), allocatable   :: particles(:)
  real(kind = fcs_real_kind_isoc) :: box_scale = 1.0

  !!! loop and tmp variables
  integer(kind = fcs_integer_kind_isoc) :: ip
  integer                               :: itime = 0
  integer                               :: rc
  integer                               :: pepc_nlocal, pepc_ntotal

  !!! debug output, may be removed ...
  if(db_level .gt. 4) then
     write(*,*) "inside pepc fortran: nlocal", nlocal
     write(*,*) "inside pepc fortran: box_c", box_c
     write(*,*) "inside pepc fortran: theta", theta
     write(*,*) "inside pepc fortran: db_level", db_level
     write(*,*) "inside pepc fortran: num_threads", nwt
     
     if(nlocal .ge. 3) then
        write(*,*) "inside pepc fortran: positions(4:9)", positions(4:9)
     end if
  end if

  !!! set pepc interaction (Coulomb) parameter
  theta2 = theta*theta
  eps2   = eps*eps
  num_threads              = nwt
  max_particles_per_thread = 100
  np_mult                  = npm
  if (db_level > 0) debug_level = ibset(db_level,0)

  !!! setup periodic domain
  fmm_extrinsic_correction = lattice_corr
  periodicity              = periodicity_in > 0

  !!! set scaling length -> all coordinates are in [0,1]^3
  if(any(periodicity .eqv. .true.)) then
    !! assume first component of first box vector contains the scaling length
    box_scale = box_a(1)
  end if

  !! scale epsilon and lattice vectors
  eps2        = eps2  / box_scale**2
  t_lattice_1 = box_a / box_scale
  t_lattice_2 = box_b / box_scale
  t_lattice_3 = box_c / box_scale
  
  !!! process the changed parameter
  call pepc_prepare(3_kind_dim)

  !!! copy input values (pos, q) into pepc particle structure, clear results
  allocate(particles(nlocal), stat=rc)
  if(rc .ne. 0) then
     write(*,*) "[PEPC] could not allocate memory for particle structure"
  end if

  !!! copy and scale positions
  do ip=1, nlocal
     particles(ip)%x           = positions(3*ip-2 : 3*ip) / box_scale
     particles(ip)%work        = work(ip)
     particles(ip)%data%q      = charges(ip)
     particles(ip)%results%e   = 0_8
     particles(ip)%results%pot = 0_8
  end do

  !!! call pepc routines
  pepc_nlocal = INT(nlocal, KIND(pepc_nlocal))
  pepc_ntotal = INT(ntotal, KIND(pepc_ntotal))
  call pepc_grow_and_traverse(particles, itime, .false., .false.)
  nlocal = INT(pepc_nlocal, KIND(nlocal))

  !!! copy result values (efield, pot), including scaling, into scafacos buffers
  do ip=1, nlocal
     efield(3*ip-2 : 3*ip) = particles(ip)%results%e   / (box_scale**2)
     potentials(ip)        = particles(ip)%results%pot / box_scale
     work(ip)              = particles(ip)%work
  end do

  !!! free particle structure
  deallocate(particles)

end subroutine pepc_scafacos_run
