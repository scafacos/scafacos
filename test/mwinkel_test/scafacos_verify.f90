!*************************************************************
! *
! * ScaFaCoS-Test programme
! *
! * by Lukas Arnold & Mathias Winkel
! *
! *************************************************************

#include <fcs_fconfig.h>
#include "setup.h"

#define NUMMETHODS 8
#define NUMPARTICLES_DEFAULT 100
#define METHOD_DEFAULT 0

! maximum number of particles to print out
#define MAX_PART_PRINT 10

module globals
  use fcs_module
  character(len=*), parameter, dimension(0:NUMMETHODS-1) :: fcs_method = ["DIRECT", &
                 "FMM   ", "MEMD  ", "NFFT  ", "P3M   ", "PEPC  ", "PP3MG ", "VMG   "]
  character(len=*), parameter, dimension(0:NUMSETUPS-1) :: test_setup = ["TRIVIAL", &
                 "NACL   ", "SIO2   "]
  fcs_real :: dummy = 1.
  integer, parameter :: rk = kind(dummy)

  ! default values taken from interface test
#if FCS_ENABLE_PP3MG
    integer :: PP3MG_dims(0:2)
    fcs_integer :: PP3MG_cells_x = 64! 64
    fcs_integer :: PP3MG_cells_y = 64! 64
    fcs_integer :: PP3MG_cells_z = 64! 64
    fcs_integer :: PP3MG_ghost_cells = 8! 8
    fcs_integer :: PP3MG_pol_degree = 4! 4
    fcs_real    :: PP3MG_err_bound = 1e-7_rk
    fcs_integer :: PP3MG_max_iterations = 10! 4
    fcs_integer :: PP3MG_max_particles = 1000! 10000
    logical :: PP3MG_periodic = .false.
#endif
#if FCS_ENABLE_P3M
    fcs_real    :: p3m_r_cut = 1.001;
    fcs_integer :: p3m_mesh  = 64;
    fcs_integer :: p3m_cao   = 7;
    fcs_real    :: p3m_alpha = 2.5;
#endif
#if FCS_ENABLE_PEPC
    fcs_real    :: PEPC_epsilon = 0.0_rk
    fcs_real    :: PEPC_theta = 0.3_rk
    fcs_integer :: PEPC_debuglevel = -1
#endif
#if FCS_ENABLE_FMM
    fcs_integer :: FMM_absrel = FCS_FMM_CUSTOM_RELATIVE
    fcs_integer :: FMM_dipole_correction = FCS_FMM_STANDARD_DIPOLE_CORRECTION
    fcs_real    :: FMM_deltaE = 1e-2
#endif

end module

subroutine checkres(res, my_rank)
   use fcs_module
   implicit none
   type(FCSResult), intent(in) ::  res
   fcs_integer, intent(in) :: my_rank
   fcs_integer :: ret

   call fcsresult_get_returncode(res,ret)
   if (ret /= FCS_SUCCESS) then
     if (my_rank == 0) call fcsresult_printResult(res)
     stop
   endif

end subroutine



subroutine setup_methodspecific(handle, my_rank, num_ranks)
  use fcs_module
  use globals
  implicit none
  type(FCSHandle) :: handle
  fcs_integer :: my_rank, num_ranks
  type(FCSResult) :: fcs_res
  fcs_integer :: method
  integer :: ierr

  call fcshandle_get_method(handle, method)

  select case (method)
#if FCS_ENABLE_PP3MG
    case (FCS_PP3MG)
        if (num_ranks > 1) then
            PP3MG_dims(0:2) = 0
            call MPI_Dims_create(num_ranks,3,PP3MG_dims, ierr)
        else
            PP3MG_dims(0:2) = 1
        endif

        call fcs_PP3MG_setup(handle, PP3MG_dims, PP3MG_cells_x, PP3MG_cells_y, PP3MG_cells_z, &
                              PP3MG_ghost_cells, PP3MG_max_iterations, PP3MG_max_particles, &
                              PP3MG_periodic, PP3MG_pol_degree, PP3MG_err_bound, fcs_res) ! TODO: in C heisst die Fkt: fcs_setup_PP3MG
        call checkres(fcs_res, my_rank)
#endif
#if FCS_ENABLE_P3M
    case (FCS_P3M) ! TODO: Fortran interface does not seem to be complete
        !fcs_P3M_set_r_cut(handle, p3m_r_cut, fcs_res)
        !call checkres(fcs_res, my_rank)
        !fcs_P3M_set_mesh(handle, p3m_mesh, fcs_res)
        !call checkres(fcs_res, my_rank)
        !fcs_P3M_set_cao(handle, p3m_cao, fcs_res)
        !call checkres(fcs_res, my_rank)
        !fcs_P3M_set_alpha(handle, p3m_alpha, fcs_res)
        !call checkres(fcs_res, my_rank)
#endif
#if FCS_ENABLE_PEPC
    case (FCS_PEPC)
        call fcs_PEPC_setup(handle, PEPC_epsilon, PEPC_theta, PEPC_debuglevel, fcs_res) !TODO: in C heisst die Funktion fcs_setup_PEPC
        call checkres(fcs_res, my_rank)
#endif
#if FCS_ENABLE_VMG
    case (FCS_VMG)
        call fcs_VMG_setup(handle, 6,.false.,6,25,3,1,1.0e-6_8, fcs_res)! TODO: in C heisst die Fkt: fcs_setup_VMG
        call checkres(fcs_res, my_rank)
#endif
#if FCS_ENABLE_FMM
    case (FCS_FMM)
        call fcs_FMM_setup(handle, FMM_absrel, FMM_deltaE, FMM_dipole_correction, fcs_res)! TODO: in C heisst die Fkt: fcs_setup_FMM
        call checkres(fcs_res, my_rank)
#endif
  end select
end subroutine



subroutine readparams(method, nparticles, setup_id, per)
    implicit none
    fcs_integer, intent(inout) :: method, nparticles, setup_id
    logical, intent(inout) :: per(3)

    integer :: argc, l,i
    character(50) :: value
    l = 50

    argc = command_argument_count()

    if (argc > 0) then
      call get_command_argument(1, value, l)
      value = trim(value)
      read(value,*) method
    endif

    if (argc > 1) then
      call get_command_argument(2, value, l)
      value = trim(value)
      read(value,*) nparticles
    endif

    if (argc > 2) then
      call get_command_argument(3, value, l)
      value = trim(value)
      read(value,*) setup_id
    endif

    if (argc > 3) then
      call get_command_argument(4, value, l)
      value = trim(value)
      per(1:3) = .false.

      if (len(value) >= 3) then
        do i=1,3
          if (value(i:i) == "1") per(i) = .true.
        end do
      endif
    endif

    if ((method > NUMMETHODS-1) .or. (method < 0)) method = 0
    if ((setup_id > NUMSETUPS-1) .or. (setup_id < 0)) setup_id = 0
end subroutine


subroutine printparams(num_ranks, my_rank, method, ntotal, setup_id, per)
  use globals
  implicit none
  fcs_integer, intent(in) :: num_ranks, my_rank, method, ntotal, setup_id
  logical, intent(in) ::  per(3)
  character(100) :: cmd
  integer :: l, i

  if (my_rank /= 0) return

  call get_command(cmd, l)
  l = scan(cmd, " ")
  cmd = cmd(1:l)

  write(*,'(/," ***** ScaFaCoS Test ***** ",/)')
  write(*,'("Call with ''", a," [METHOD [NUMPARTICLES [SETUP_ID [PERIODICITY]]]]''")') trim(cmd)
  write(*,'("     e.g. ''", a," 5 10 011''")') trim(cmd)
  write(*,'("Parameters:")')
  write(*,'("     NUMPARTICLES - total particle number (default ", I0, ")")') NUMPARTICLES_DEFAULT
  write(*,'("     SETUP_ID     - one of the following integer values (default ", I0, ", current highlighted)")') SETUP_DEFAULT
  do i=0,NUMSETUPS-1
    if (i==setup_id) then
      write(*,'("     ==>")', advance='no')
    else
      write(*,'("        ")', advance='no')
    endif

    write(*,'("  ", I2, "  - ", a8)') i, trim(test_setup(i))
  end do
  write(*,'("     PERIODICITY  - bitmask for periodicity in three spatial dimensions (default ''000'')")')
  write(*,'("     METHOD       - one of the following integer values (default ", I0," current highlighted)")') METHOD_DEFAULT

  do i=0,NUMMETHODS-1
    if (i==method) then
      write(*,'("     -->")', advance='no')
    else
      write(*,'("        ")', advance='no')
    endif

    write(*,'("  ", I2, "  - ", a8)') i, trim(fcs_method(i))
  end do

  write(*,'(/,"Running on ",I0, " MPI ranks with in total ", I0, " particles",/)') num_ranks, ntotal
end subroutine



program main
  use fcs_module
  use globals
  use iso_fortran_env
  implicit none

  fcs_integer :: num_ranks, my_rank
  fcs_integer :: method = METHOD_DEFAULT
  fcs_integer :: ntotal = NUMPARTICLES_DEFAULT
  logical :: per(3) = [.false.,.false.,.false.]
  fcs_real, dimension(3) ::  boxa, boxb, boxc;
  fcs_real :: energy_local  = 0.
  fcs_real :: energy_global = 0.
  fcs_real :: trvirial      = 0.
  fcs_real :: trvirial2l    = 0.
  fcs_real :: trvirial2     = 0.
  fcs_integer:: setup_id    = SETUP_DEFAULT
  integer :: comm

  logical :: comp_flag = .false.
  fcs_integer :: nparts, max_nparts
  fcs_real, allocatable, dimension(:) ::  pos, qs, es, pot, virial

   integer :: ierr, i

   type(FCSOutput)     ::  fcs_out

   type(FCSResult)     ::  fcs_res
   type(FCSHandle)     ::  fcs_handle

   call MPI_Init(ierr)

   call MPI_Comm_size(MPI_COMM_WORLD, num_ranks, ierr)
   call MPI_Comm_rank(MPI_COMM_WORLD, my_rank, ierr)

   call readparams(method, ntotal, setup_id, per)
   ! Particle setup
   max_nparts = MAXPARTS(ntotal/num_ranks + 1);

   allocate(pos(0:max_nparts*3-1))
   allocate( qs(0:max_nparts*1-1))
   allocate( es(0:max_nparts*3-1))
   allocate(pot(0:max_nparts*1-1))
   allocate(virial(9))

   call setup_particles(my_rank, num_ranks, setup_id, ntotal, nparts, max_nparts, pos, qs, es, pot, boxa, boxb, boxc, comm);

   call printparams(num_ranks, my_rank, method, ntotal, setup_id, per)

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! ScaFaCoS-Initialization
   call fcs_init(fcs_handle, trim(fcs_method(method)), MPI_COMM_WORLD, fcs_res) !TODO: Reihenfolge der Parameter wie in der C-Variante
   call checkres(fcs_res, my_rank)
   !! ScaFaCoS: Generic Parameter Setup
   call fcs_set_common(fcs_handle, .true., boxa, boxb, boxc, ntotal, ntotal, per, fcs_res)
   call checkres(fcs_res, my_rank)
   !! ScaFaCoS: Method-specific Parameter Setup
   call setup_methodspecific(fcs_handle, my_rank, num_ranks)
   !! ScaFaCoS: Method-specific Internal Tuning
   call fcs_tune(fcs_handle, nparts, max_nparts, pos, qs, comp_flag, fcs_res)
   call checkres(fcs_res, my_rank)
   !! ScaFaCoS: Output of internal Parameters
   if (my_rank ==0) call fcshandle_print(fcs_handle, output_unit) !:TODO: Die Funktion heisst im C-Interface fcs_output() und hat einen Parameter weniger
   !! ScaFaCoS: Creation of Output Handle
   call fcsoutput_create(fcs_out, nparts, max_nparts, fcs_res) ! TODO: wozu nparts & max_nparts ?
   call checkres(fcs_res, my_rank)
   !! ScaFaCoS: Invocation of Coulomb Solver
   call fcs_run(fcs_handle,nparts,max_nparts,pos,qs,comp_flag,fcs_out,fcs_res) !TODO: wo sind es und pot?
   call checkres(fcs_res, my_rank)
   !! ScaFaCoS: Obtain result data
   ! TODO: die folgenden drei Zeilen sind ist in der C-Variante so nicht noetig
   call fcsoutput_get_potentials(fcs_out, pot) ! TODO: diese Funktionen kopieren nicht die Daten in mein pot-Feld, sondern behandeln pot, es und virial wie einen Zeiger, der danach irgndwohin zeigt, wo die Felddaten liegen
   call fcsoutput_get_field(fcs_out, es) ! TODO: d.h auch, dass die Feldindizes 1:nparts laufen statt wie nutzerspezifisch alloziiert
   call fcsoutput_get_virial(fcs_out, virial) ! TODO: Funktion heisst anders als in der C-Variante (Unterstrich)
   !! ScaFaCoS: Free any allocated objects
   call fcs_destroy(fcs_handle,fcs_res)
   call checkres(fcs_res, my_rank)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !! Evaluation of results
   do i=0,nparts-1
     energy_local = energy_local + pot(i+1)*qs(i)/2.
     trvirial2l   = trvirial2l   + dot_product(pos(3*i+0:3*i+2), es(3*i+1:3*i+3)) * qs(i)
   end do

   energy_local = energy_local / (1.*ntotal)

   call MPI_Allreduce(energy_local, energy_global, 1, FCS_MPI_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)

   do i=0,2
     trvirial = trvirial + virial(4*i+1)
   end do

   trvirial   = trvirial   / (1.*ntotal)
   trvirial2l = trvirial2l / (1.*ntotal)

   call MPI_Allreduce(trvirial2l, trvirial2, 1, FCS_MPI_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)

   !! Output of results
   if (my_rank ==0) then
      write(*,'(/,/,"******* debug output: content of output handle ********",/)')
      call fcsOutput_print(fcs_out) ! TODO: Die C-Funktion hat noch den MPI-Rank als Argument

      write(*,'(/,/,"******* debug output: energy and trace of virial ******",/)')
      write(*,'("energy/particle:          ", g14.8)') energy_global
      write(*,'("trace of virial/particle: ", g14.8)') trvirial
      write(*,'("dito (from forces):       ", g14.8)') trvirial2

      write(*,'(/,/,"******* debug output: particles and fields ************",/)')
      write(*,'("[Rank] Part. | (      x     ,     y      ,     z      )       q       |", &
                              "  (     Ex     ,     Ey     ,     Ez     )     pot     ")')
      write(*,'("-------------+--------------------------------------------------------+", &
                              "-------------------------------------------------------")')
      flush(output_unit)
   endif

   call MPI_Barrier(MPI_COMM_WORLD, ierr)

   do i=0,min(nparts,MAX_PART_PRINT/num_ranks)-1
      write(*,'("[",i4.4,"] ",i5.5" | (", g12.3,",",g12.3,",",g12.3,") ",g12.3,   &
                                 "  |  (", g12.3,",",g12.3,",",g12.3,") ",g12.3)') &
                 my_rank, i, pos(3*i+0), pos(3*i+1), pos(3*i+2), qs(i),            &
                 es(3*i+0+1), es(3*i+1+1), es(3*i+2+1), pot(i+1) !TODO: Feldindizes laufen 1:nparts statt 0:nparts-1
   end do

   call MPI_Finalize()

end program
