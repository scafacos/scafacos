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
!> Encapsulates calculation of the lattice contribution by means
!> of the FMM-approach to the lattice
!>
module module_fmm_framework
      use module_pepc_types
      use module_debug
      use module_mirror_boxes
      use module_interaction_specific_types, only : t_tree_node_interaction_data
      implicit none
      include 'mpif.h'
      private

      ! TODO: set dipole- and low-order stuff in MLattice to zero

      !> far- and near-field contribution to potential energy (has to be calculated in fields.p90)
      real*8, public :: potfarfield, potnearfield

      integer, public, parameter :: FMM_EXTRINSIC_CORRECTION_NONE        =  0 !< no extrinsic-to-intrinsic correction
      integer, public, parameter :: FMM_EXTRINSIC_CORRECTION_REDLACK     =  1 !< correction expression as given by Redlack and Grindlay (only for cubic boxes)
                                                                              !< (see [J.Chem.Phys. 107, 10131, eqn.(19,20)] for details, inside this publication, the volume factor is missing;
                                                                              !<      [J. Chem. Phys. 101, 5024, eqn (5)] contains it) -- this is the method that was used if do_extrinsic_correction=.true. before
      integer, public, parameter :: FMM_EXTRINSIC_CORRECTION_FICTCHARGE  =  2 !< fictitious charges as given by Kudin (should work for all unit chell shapes) [Kudin 1998, eq. (2.8)]
      integer, public, parameter :: FMM_EXTRINSIC_CORRECTION_MEASUREMENT =  3 !< measurement of correction value [Kudin 1998, eq. (2.6,, 2.7)], FIXME: currently not implemented
      !> type of dipole correction, see [J.Chem.Phys. 107, 10131, eq. (19,20)]
      integer, public :: fmm_extrinsic_correction = FMM_EXTRINSIC_CORRECTION_FICTCHARGE

      public fmm_framework_init
      public fmm_framework_timestep
      public fmm_sum_lattice_force
      public lattice_vect
      public fmm_framework_param_dump

      ! general stuff
      integer(kind_pe) :: myrank
      integer :: MPI_COMM_fmm
      ! precision flags
      integer, parameter :: kfp                  = 8 ! numeric precision (kind value)
      integer, parameter :: MPI_REAL_fmm         = MPI_REAL8
      logical, parameter  :: chop_arrays         = .false.
      real(kfp), parameter :: prec = 1.E-16_kfp
      ! shortcut notations
      real(kfp), parameter :: zero = 0._kfp
      real(kfp), parameter :: one  = 1._kfp
      real(kfp), parameter :: two  = 2._kfp
      real(kfp), parameter :: three= 3._kfp
      ! FMM-PARAMETERS
      integer, parameter :: Lmax_multipole = 20
      integer, parameter :: Lmax_taylor    = Lmax_multipole * 2
      integer, parameter :: MaxIter        = 32
      integer :: ws = 1
      logical, parameter :: use_pretabulated_lattice_coefficients = .false.
      ! FMM-VARIABLES
      integer, parameter :: fmm_array_length_multipole = Lmax_multipole*(Lmax_multipole+1)/2+Lmax_multipole+1
      integer, parameter :: fmm_array_length_taylor    = Lmax_taylor   *(Lmax_taylor   +1)/2+Lmax_taylor   +1
      ! internally calculated FMM variables
      complex(kfp) :: mu_cent(1:fmm_array_length_taylor)
      complex(kfp) :: omega_tilde(1:fmm_array_length_multipole)
      complex(kfp) :: MLattice(1:fmm_array_length_taylor)
      !> variables for extrinsic to intrinsic correction
      real(kfp) :: box_dipole(3) = zero
      real(kfp) :: quad_trace    = zero
      !> fictitious charges and their position: fictcharge(0,i)=q_i, fictcharge(1:3,i)=r_i
      type(t_tree_node_interaction_data) :: fictcharge(1:4)
      integer :: nfictcharge = 0
      
      contains

        !>
        !> Module Initialization, should be called on program startup
        !> after setting up all geometric parameters etc.
        !>  - well-separation criterion is automatically set to module_mirror_boxes::mirror_box_layers
        !>  - module_mirror_boxes::calc_neighbour_boxes() must have been called before calling this routine
        !>
        !> @param[in] mpi_rank MPI rank of process for controlling debug output
        !> @param[in] mpi_comm MPI communicator to be used
        !>
        subroutine fmm_framework_init(mpi_rank, mpi_comm)
          use module_debug
          use module_mirror_boxes, only : mirror_box_layers
          implicit none
          integer(kind_pe), intent(in) :: mpi_rank
          integer, intent(in) :: mpi_comm

          myrank       = mpi_rank
          MPI_COMM_fmm = mpi_comm
          ws           = mirror_box_layers

          do_periodic = any(periodicity(1:3))

           ! anything above has to be done in any case
          if (do_periodic) then
            MLattice = 0

            if (use_pretabulated_lattice_coefficients .or. system_is_unit_cube()) then
              call load_lattice_coefficients(MLattice)
              if (use_pretabulated_lattice_coefficients) then
                DEBUG_WARNING('(a)', "Using pretabulated lattice coefficients. These are only valid for 3D-periodic unit-box simulations regions.")
              endif
            else
              call calc_lattice_coefficients(MLattice)
            endif

            if ((myrank == 0) .and. dbg(DBG_PERIODIC)) then
              call WriteTableToFile('MLattice.tab', MLattice, Lmax_taylor)
            end if
          end if
        end subroutine fmm_framework_init


        !>
        !> Refreshes Multipole information and Taylor coefficients,
        !> has to be called every timestep with particles that were used in tree buildup
        !>
        subroutine fmm_framework_timestep(particles)
          use module_pepc_types
          use module_mirror_boxes
          use module_debug
          implicit none
          type(t_particle), intent(in) :: particles(:)

          if (do_periodic) then
            if (.not. check_lattice_boundaries(particles)) then
              DEBUG_ERROR(*, 'Lattice contribution will be wrong. Aborting.')
            endif
            
            if (      (fmm_extrinsic_correction == FMM_EXTRINSIC_CORRECTION_REDLACK) &
                 .or. (fmm_extrinsic_correction == FMM_EXTRINSIC_CORRECTION_FICTCHARGE)) then
              call calc_box_dipole(particles)
            endif
            call calc_omega_tilde(particles)
            call calc_mu_cent(omega_tilde, mu_cent)
          endif
        end subroutine fmm_framework_timestep


        !>
        !> Calculates the lattice coefficients for computing mu_cent
        !>
        subroutine calc_lattice_coefficients(ML)
          use module_debug
          implicit none

          complex(kfp), intent(out) :: ML(1:fmm_array_length_taylor)
          ! Variables for lattice coefficient calculation
          complex(kfp) :: Mstar(1:fmm_array_length_multipole)
          complex(kfp) :: Lstar(1:fmm_array_length_taylor)

          integer :: l, m, iter
          character(20) :: fn

          call pepc_status('LATTICE COEFFICIENTS: Starting calculation')

          Mstar = 0
          Lstar = 0

          ! pretabulation of necessary values
          do l = 0,Lmax_multipole
            do m = 0,l
              Mstar( tblinv(l, m, Lmax_multipole) ) = MstarFunc(l, m)
            end do
          end do

          do l = 0,Lmax_taylor
            do m = 0,l
              Lstar( tblinv(l, m, Lmax_taylor) )    = LstarFunc(l, m)
            end do
          end do

          call chop(Mstar)
          call chop(Lstar)

          if ((myrank == 0) .and. dbg(DBG_PERIODIC)) then
            call WriteTableToFile('Mstar.tab', Mstar, Lmax_multipole)
            call WriteTableToFile('Lstar.tab', Lstar, Lmax_taylor)
          end if

          ! zeroth step of iteration
          ML = Lstar
          
          ! FIXME untested : call zero_terms_taylor(ML)

          do iter = 1,MaxIter

            ML = M2L( UL( ML ) , MStar ) + Lstar

            ! FIXME untested : call zero_terms_taylor(ML)

            if ((myrank == 0) .and. dbg(DBG_PERIODIC)) then
              write(fn,'("MLattice.", I2.2, ".tab")') iter
              call WriteTableToFile(trim(fn), ML, Lmax_taylor)
            endif
          end do

          call chop(ML)

          ! ML(1:tblinv(3,3))=0
          call pepc_status('LATTICE COEFFICIENTS: finished calculation')
        end subroutine calc_lattice_coefficients


        !>
        !> Sets those elements in a Taylor expansion to zero that vanish
        !> anyway, see J. Chem. Phys 121, 2886, Section V
        !>
        subroutine zero_terms_taylor(M)
          use module_mirror_boxes, only : periodicity
          implicit none
          complex(kfp), intent(inout) :: M(1:fmm_array_length_taylor)
          
          ! 1D, 2D, and 3D periodicity
          M(tblinv(0, 0, Lmax_taylor)) = (zero, zero)
          
          if (count(periodicity) > 2) then
            ! only for 3D-periodic case
            M(tblinv(1, 0, Lmax_taylor)) = (zero, zero)
            M(tblinv(1, 1, Lmax_taylor)) = (zero, zero)
          endif
        end subroutine


        !>
        !> Sets those elements in a Multipole expansion to zero that vanish
        !> anyway, see J. Chem. Phys 121, 2886, Section V
        !>
        subroutine zero_terms_multipole(O)
          use module_mirror_boxes, only : periodicity
          implicit none
          complex(kfp), intent(inout) :: O(1:fmm_array_length_multipole)
          
          ! 1D, 2D, and 3D periodicity
          O(tblinv(0, 0, Lmax_multipole)) = (zero, zero)
          
          if (count(periodicity) > 2) then
            ! only for 3D-periodic case
            O(tblinv(1, 0, Lmax_multipole)) = (zero, zero)
            O(tblinv(1, 1, Lmax_multipole)) = (zero, zero)
          endif
        end subroutine
         

        !>
        !> Sets the lattice coefficients for computing mu_cent
        !> data computed with this code (Lamx_multipole=50, MaxIter=32)
        !> compare with
        !>  - [Challacombe, White, Head-Gordon: J. Chem. Phys. 107, 10131]
        !>  - PhD thesis of Ivo Kabadshow
        !>
        subroutine load_lattice_coefficients(M)
          use module_debug
          implicit none

          complex(kfp), intent(inout) :: M(1:fmm_array_length_taylor)

          call pepc_status('LATTICE COEFFICIENTS: Loading')

          if (Lmax_taylor >= 0) then
            M(tblinv(0, 0, Lmax_taylor)) =   0.66196708399502771246378770654822400E+33_kfp
          endif
          if (Lmax_taylor >= 2) then
            M(tblinv(2,  0, Lmax_taylor)) =   0.29501109556337144811408208229350297E-13_kfp
          endif
          if (Lmax_taylor >= 4) then
            M(tblinv(4,  0, Lmax_taylor)) =   0.28119304871888664010270986182149500E+01_kfp
            M(tblinv(4,  4, Lmax_taylor)) =   0.14059652435944322235172876389697194E+02_kfp
          endif
          if (Lmax_taylor >= 6) then
            M(tblinv(6,  0, Lmax_taylor)) =   0.54795908739322285452288952001254074E+00_kfp
            M(tblinv(6,  4, Lmax_taylor)) =   -0.38357136117525469920508385257562622E+01_kfp
          endif
          if (Lmax_taylor >= 8) then
            M(tblinv(8,  0, Lmax_taylor)) =   0.12156157302097911099281191127374768E+03_kfp
            M(tblinv(8,  4, Lmax_taylor)) =   0.12156157302097928152306849369779229E+03_kfp
            M(tblinv(8,  8, Lmax_taylor)) =   0.79015022463636523752938956022262573E+04_kfp
          endif
          if (Lmax_taylor >= 10) then
            M(tblinv(10,  0, Lmax_taylor)) =   0.31179916736109191788273165002465248E+03_kfp
            M(tblinv(10,  4, Lmax_taylor)) =   -0.68595816819440403833141317591071129E+03_kfp
            M(tblinv(10,  8, Lmax_taylor)) =   -0.11661288859304897414403967559337616E+05_kfp
          endif
          if (Lmax_taylor >= 12) then
            M(tblinv(12,  0, Lmax_taylor)) =   0.24245612747359092463739216327667236E+06_kfp
            M(tblinv(12,  4, Lmax_taylor)) =   0.20375858264140240498818457126617432E+06_kfp
            M(tblinv(12,  8, Lmax_taylor)) =   0.70682666545985033735632896423339844E+06_kfp
            M(tblinv(12, 12, Lmax_taylor)) =   0.23702435984527084231376647949218750E+09_kfp
          endif
          if (Lmax_taylor >= 14) then
            M(tblinv(14,  0, Lmax_taylor)) =   0.20954087119885545689612627029418945E+07_kfp
            M(tblinv(14,  4, Lmax_taylor)) =   -0.26940969154138588346540927886962891E+07_kfp
            M(tblinv(14,  8, Lmax_taylor)) =   -0.17062613797621075063943862915039063E+08_kfp
            M(tblinv(14, 12, Lmax_taylor)) =   -0.65406686224214184284210205078125000E+09_kfp
          endif
          if (Lmax_taylor >= 16) then
            M(tblinv(16,  0, Lmax_taylor)) =   0.54279858299650156497955322265625000E+09_kfp
            M(tblinv(16,  4, Lmax_taylor)) =   0.22841041529105401039123535156250000E+09_kfp
            M(tblinv(16,  8, Lmax_taylor)) =   0.12973301854895758628845214843750000E+10_kfp
            M(tblinv(16, 12, Lmax_taylor)) =   0.25882484900055580139160156250000000E+11_kfp
            M(tblinv(16, 16, Lmax_taylor)) =   0.69973653547984179687500000000000000E+13_kfp
          endif
          if (Lmax_taylor >= 18) then
            M(tblinv(18,  0, Lmax_taylor)) =   0.14686049951258449554443359375000000E+11_kfp
            M(tblinv(18,  4, Lmax_taylor)) =   -0.15376346487994709014892578125000000E+11_kfp
            M(tblinv(18,  8, Lmax_taylor)) =   -0.24226921558569625854492187500000000E+11_kfp
            M(tblinv(18, 12, Lmax_taylor)) =   -0.10692416604738659667968750000000000E+13_kfp
            M(tblinv(18, 16, Lmax_taylor)) =   -0.39585194668444328125000000000000000E+14_kfp
          endif
          if (Lmax_taylor >= 20) then
            M(tblinv(20,  0, Lmax_taylor)) =   0.29414124910043237304687500000000000E+13_kfp
            M(tblinv(20,  4, Lmax_taylor)) =   0.50799363324581787109375000000000000E+12_kfp
            M(tblinv(20,  8, Lmax_taylor)) =   0.43319375525806137695312500000000000E+13_kfp
            M(tblinv(20, 12, Lmax_taylor)) =   0.47785845726839648437500000000000000E+14_kfp
            M(tblinv(20, 16, Lmax_taylor)) =   0.16103883836731942500000000000000000E+16_kfp
            M(tblinv(20, 20, Lmax_taylor)) =   0.50103139602723200000000000000000000E+18_kfp
          endif
          if (Lmax_taylor >= 22) then
            M(tblinv(22,  0, Lmax_taylor)) =   0.15704492101503415625000000000000000E+15_kfp
            M(tblinv(22,  4, Lmax_taylor)) =   -0.13167755327760200000000000000000000E+15_kfp
            M(tblinv(22,  8, Lmax_taylor)) =   -0.19393166587859906250000000000000000E+15_kfp
            M(tblinv(22, 12, Lmax_taylor)) =   -0.27630226717294555000000000000000000E+16_kfp
            M(tblinv(22, 16, Lmax_taylor)) =   -0.70704971434484640000000000000000000E+17_kfp
            M(tblinv(22, 20, Lmax_taylor)) =   -0.47508773151765964800000000000000000E+19_kfp
          endif
          if (Lmax_taylor >= 24) then
            M(tblinv(24,  0, Lmax_taylor)) =   0.49787311095229424000000000000000000E+17_kfp
            M(tblinv(24,  4, Lmax_taylor)) =   0.17493021689018656000000000000000000E+17_kfp
            M(tblinv(24,  8, Lmax_taylor)) =   0.47119084871278656000000000000000000E+17_kfp
            M(tblinv(24, 12, Lmax_taylor)) =   0.27781944545203244800000000000000000E+18_kfp
            M(tblinv(24, 16, Lmax_taylor)) =   0.39862401037549358080000000000000000E+19_kfp
            M(tblinv(24, 20, Lmax_taylor)) =   0.20240160016231468236800000000000000E+21_kfp
            M(tblinv(24, 24, Lmax_taylor)) =   0.14575287885847205406310400000000000E+24_kfp
          endif
          if (Lmax_taylor >= 26) then
            M(tblinv(26,  0, Lmax_taylor)) =   0.38825844842441052160000000000000000E+19_kfp
            M(tblinv(26,  4, Lmax_taylor)) =   -0.23376165598727029760000000000000000E+19_kfp
            M(tblinv(26,  8, Lmax_taylor)) =   -0.63197488917431695360000000000000000E+19_kfp
            M(tblinv(26, 12, Lmax_taylor)) =   -0.31129980912121520128000000000000000E+20_kfp
            M(tblinv(26, 16, Lmax_taylor)) =   -0.33888386262957850624000000000000000E+21_kfp
            M(tblinv(26, 20, Lmax_taylor)) =   -0.10355256572352317620224000000000000E+23_kfp
            M(tblinv(26, 24, Lmax_taylor)) =   -0.16526385068628147039109120000000000E+25_kfp
          endif
          if (Lmax_taylor >= 28) then
            M(tblinv(28, 0, Lmax_taylor)) =   0.15660190510982260326400000000000000E+22_kfp
            M(tblinv(28, 4, Lmax_taylor)) =   0.52389541759563156684800000000000000E+21_kfp
            M(tblinv(28, 8, Lmax_taylor)) =   0.10189364356236629770240000000000000E+22_kfp
            M(tblinv(28, 12, Lmax_taylor)) =   0.60516507919279487713280000000000000E+22_kfp
            M(tblinv(28, 16, Lmax_taylor)) =   0.44738315975593851617280000000000000E+23_kfp
            M(tblinv(28, 20, Lmax_taylor)) =   0.78725443506303414147481600000000000E+24_kfp
            M(tblinv(28, 24, Lmax_taylor)) =   0.74143442402174531771826176000000000E+26_kfp
            M(tblinv(28, 28, Lmax_taylor)) =   0.69657271111833057487556182016000000E+29_kfp
          endif
          if (Lmax_taylor >= 30) then
            M(tblinv(30,   0, Lmax_taylor)) =   0.17565213196687985606656000000000000E+24_kfp
            M(tblinv(30,  4, Lmax_taylor)) =   -0.90284881313281234436096000000000000E+23_kfp
            M(tblinv(30,  8, Lmax_taylor)) =   -0.22493547031100229523865600000000000E+24_kfp
            M(tblinv(30, 12, Lmax_taylor)) =   -0.83171156039264345522176000000000000E+24_kfp
            M(tblinv(30, 16, Lmax_taylor)) =   -0.67252085460554353821614080000000000E+25_kfp
            M(tblinv(30, 20, Lmax_taylor)) =   -0.93070427828784211593003008000000000E+26_kfp
            M(tblinv(30, 24, Lmax_taylor)) =   -0.43512348030325508745734389760000000E+28_kfp
            M(tblinv(30, 28, Lmax_taylor)) =   -0.94353504240103762217068081971200000E+30_kfp
          endif
          if (Lmax_taylor >= 32) then
            M(tblinv(32,  0, Lmax_taylor)) =   0.78458227133695181778321408000000000E+26_kfp
            M(tblinv(32,  4, Lmax_taylor)) =   0.19800659318770194003787776000000000E+26_kfp
            M(tblinv(32,  8, Lmax_taylor)) =   0.40707502807452171343233024000000000E+26_kfp
            M(tblinv(32, 12, Lmax_taylor)) =   0.19735657821534461612995379200000000E+27_kfp
            M(tblinv(32, 16, Lmax_taylor)) =   0.12450909829876972207619440640000000E+28_kfp
            M(tblinv(32, 20, Lmax_taylor)) =   0.14182362308984682793708552192000000E+29_kfp
            M(tblinv(32, 24, Lmax_taylor)) =   0.40345092874684789236447995494400000E+30_kfp
            M(tblinv(32, 28, Lmax_taylor)) =   0.45444444205682558857562241892352000E+32_kfp
            M(tblinv(32, 32, Lmax_taylor)) =   0.50653072148998691431842768333307904E+35_kfp
          endif
          if (Lmax_taylor >= 34) then
            M(tblinv(34,   0, Lmax_taylor)) =   0.12278316702141080761339478016000000E+29_kfp
            M(tblinv(34,  4, Lmax_taylor)) =   -0.59526574797376592337887559680000000E+28_kfp
            M(tblinv(34,  8, Lmax_taylor)) =   -0.11302400297437305963909480448000000E+29_kfp
            M(tblinv(34, 12, Lmax_taylor)) =   -0.38755995239471066353042980864000000E+29_kfp
            M(tblinv(34, 16, Lmax_taylor)) =   -0.23039296235918866510876613017600000E+30_kfp
            M(tblinv(34, 20, Lmax_taylor)) =   -0.24872491441724127447545492275200000E+31_kfp
            M(tblinv(34, 24, Lmax_taylor)) =   -0.54905190027187592671741710696448000E+32_kfp
            M(tblinv(34, 28, Lmax_taylor)) =   -0.31190929147945806824910150623559680E+34_kfp
            M(tblinv(34, 32, Lmax_taylor)) =   -0.81660979692832048923194261235983974E+36_kfp
          endif
          if (Lmax_taylor >= 36) then
            M(tblinv(36,  0, Lmax_taylor)) =   0.68042610132272997473266475991040000E+31_kfp
            M(tblinv(36,  4, Lmax_taylor)) =   0.16899043765958400366236600893440000E+31_kfp
            M(tblinv(36,  8, Lmax_taylor)) =   0.34324695002843274371624075264000000E+31_kfp
            M(tblinv(36, 12, Lmax_taylor)) =   0.10891961584200954183772928999424000E+32_kfp
            M(tblinv(36, 16, Lmax_taylor)) =   0.58629472355230521588913774002176000E+32_kfp
            M(tblinv(36, 20, Lmax_taylor)) =   0.49998188032652907531537992594227200E+33_kfp
            M(tblinv(36, 24, Lmax_taylor)) =   0.91906783013472935407104169837854720E+34_kfp
            M(tblinv(36, 28, Lmax_taylor)) =   0.34340831446128471752777543275302093E+36_kfp
            M(tblinv(36, 32, Lmax_taylor)) =   0.43228808611207350488871067804720169E+38_kfp
            M(tblinv(36, 36, Lmax_taylor)) =   0.68507719909013291704154043991494769E+41_kfp
          endif
          if (Lmax_taylor >= 38) then
            M(tblinv(38,  0, Lmax_taylor)) =    0.13898229629288749146802627217981440E+34_kfp
            M(tblinv(38,  4, Lmax_taylor)) =   -0.60270926668898774338105280928153600E+33_kfp
            M(tblinv(38,  8, Lmax_taylor)) =   -0.10775625047577013971523006848040960E+34_kfp
            M(tblinv(38, 12, Lmax_taylor)) =   -0.34766229802879845205641022637342720E+34_kfp
            M(tblinv(38, 16, Lmax_taylor)) =   -0.14823209389154711189182603646205952E+35_kfp
            M(tblinv(38, 20, Lmax_taylor)) =   -0.11752399281934163140926597781625242E+36_kfp
            M(tblinv(38, 24, Lmax_taylor)) =   -0.17828636580470227178141304832036700E+37_kfp
            M(tblinv(38, 28, Lmax_taylor)) =   -0.52302605673609923514633048016233169E+38_kfp
            M(tblinv(38, 32, Lmax_taylor)) =   -0.35118870041597762879848563642386361E+40_kfp
            M(tblinv(38, 36, Lmax_taylor)) =   -0.12463088506162654718653256758114972E+43_kfp
          endif
          if (Lmax_taylor >= 40) then
            M(tblinv(40,  0, Lmax_taylor)) =   0.93977144756954277273979168067944448E+36_kfp
            M(tblinv(40,  4, Lmax_taylor)) =   0.23714015731391362877925506034119475E+36_kfp
            M(tblinv(40,  8, Lmax_taylor)) =   0.42743502595313863625949324121630310E+36_kfp
            M(tblinv(40, 12, Lmax_taylor)) =   0.11495283403223308187177106969334907E+37_kfp
            M(tblinv(40, 16, Lmax_taylor)) =   0.52276986864725787950760775438271775E+37_kfp
            M(tblinv(40, 20, Lmax_taylor)) =   0.33727766809055668839809403273932177E+38_kfp
            M(tblinv(40, 24, Lmax_taylor)) =   0.40969785140672314541232018394122394E+39_kfp
            M(tblinv(40, 28, Lmax_taylor)) =   0.96152785881824123788054152991004573E+40_kfp
            M(tblinv(40, 32, Lmax_taylor)) =   0.44874401799468575928849907858428924E+42_kfp
            M(tblinv(40, 36, Lmax_taylor)) =   0.72067379047599356261543581978928705E+44_kfp
            M(tblinv(40, 40, Lmax_taylor)) =   0.14800934278189952376201098541279335E+48_kfp
          endif
          if (Lmax_taylor >= 42) then
            M(tblinv(42,  0, Lmax_taylor)) =    0.24101076028312310966209393434722042E+39_kfp
            M(tblinv(42,  4, Lmax_taylor)) =   -0.93281129758571302744991251804222128E+38_kfp
            M(tblinv(42,  8, Lmax_taylor)) =   -0.16863712622021360123396968923689019E+39_kfp
            M(tblinv(42, 12, Lmax_taylor)) =   -0.46853042206730476424367486564414836E+39_kfp
            M(tblinv(42, 16, Lmax_taylor)) =   -0.17641892567348816407157155646613016E+40_kfp
            M(tblinv(42, 20, Lmax_taylor)) =   -0.10797655474740776079434477868503980E+41_kfp
            M(tblinv(42, 24, Lmax_taylor)) =   -0.11246877236296608016299482503228058E+42_kfp
            M(tblinv(42, 28, Lmax_taylor)) =   -0.20969526317812940926745387753660932E+43_kfp
            M(tblinv(42, 32, Lmax_taylor)) =   -0.75724665322139793932309232476273515E+44_kfp
            M(tblinv(42, 36, Lmax_taylor)) =   -0.66768070679697744883455851575286634E+46_kfp
            M(tblinv(42, 40, Lmax_taylor)) =   -0.29755158895080900545179527774736656E+49_kfp
          endif
          if (Lmax_taylor >= 44) then
            M(tblinv(44,  0, Lmax_taylor)) =   0.18802895551242673459138461852530244E+42_kfp
            M(tblinv(44,  4, Lmax_taylor)) =   0.43311029073565793748509578220261652E+41_kfp
            M(tblinv(44,  8, Lmax_taylor)) =   0.73183638551002990062234712155205585E+41_kfp
            M(tblinv(44, 12, Lmax_taylor)) =   0.19166034402485976735982481534810537E+42_kfp
            M(tblinv(44, 16, Lmax_taylor)) =   0.73232212910142631303132939886618574E+42_kfp
            M(tblinv(44, 20, Lmax_taylor)) =   0.41012643572849966900066752483471331E+43_kfp
            M(tblinv(44, 24, Lmax_taylor)) =   0.36428462773419336183953301342681366E+44_kfp
            M(tblinv(44, 28, Lmax_taylor)) =   0.55060498466948680359498813385646596E+45_kfp
            M(tblinv(44, 32, Lmax_taylor)) =   0.15457075645609170865398960598067162E+47_kfp
            M(tblinv(44, 36, Lmax_taylor)) =   0.94835735621465907545197009937397149E+48_kfp
            M(tblinv(44, 40, Lmax_taylor)) =   0.18713883520362201373510777806841623E+51_kfp
            M(tblinv(44, 44, Lmax_taylor)) =   0.45107670269567111253033410002734787E+54_kfp
          endif
          if (Lmax_taylor >= 46) then
            M(tblinv(46,  0, Lmax_taylor)) =    0.59784640964094988095060157708745782E+44_kfp
            M(tblinv(46,  4, Lmax_taylor)) =   -0.21718732499012689597368020389604494E+44_kfp
            M(tblinv(46,  8, Lmax_taylor)) =   -0.36782514431632691178151233050748931E+44_kfp
            M(tblinv(46, 12, Lmax_taylor)) =   -0.88873556144958931198450002034029892E+44_kfp
            M(tblinv(46, 16, Lmax_taylor)) =   -0.31603183651814223604291954933692593E+45_kfp
            M(tblinv(46, 20, Lmax_taylor)) =   -0.16091883309137843277165046271290753E+46_kfp
            M(tblinv(46, 24, Lmax_taylor)) =   -0.13287100221864350079915775281789533E+47_kfp
            M(tblinv(46, 28, Lmax_taylor)) =   -0.17200569376295918617028117347212389E+48_kfp
            M(tblinv(46, 32, Lmax_taylor)) =   -0.38146081679001297318366488451434897E+49_kfp
            M(tblinv(46, 36, Lmax_taylor)) =   -0.17463078833474928188216353940559520E+51_kfp
            M(tblinv(46, 40, Lmax_taylor)) =   -0.19348009546647245429748014166800000E+53_kfp
            M(tblinv(46, 44, Lmax_taylor)) =   -0.10021492491483168007613903982234254E+56_kfp
          endif
          if (Lmax_taylor >= 48) then
            M(tblinv(48,  0, Lmax_taylor)) =   0.54091590625627984312694298022521638E+47_kfp
            M(tblinv(48,  4, Lmax_taylor)) =   0.11663289142002507285726138558407866E+47_kfp
            M(tblinv(48,  8, Lmax_taylor)) =   0.19297849313251926818093115443991633E+47_kfp
            M(tblinv(48, 12, Lmax_taylor)) =   0.46489785748709634074809721449710401E+47_kfp
            M(tblinv(48, 16, Lmax_taylor)) =   0.15135360353073524771531329454490740E+48_kfp
            M(tblinv(48, 20, Lmax_taylor)) =   0.74875422886287882154503340453273852E+48_kfp
            M(tblinv(48, 24, Lmax_taylor)) =   0.54529024679443529183710127100291682E+49_kfp
            M(tblinv(48, 28, Lmax_taylor)) =   0.61848643203375079716184234543977896E+50_kfp
            M(tblinv(48, 32, Lmax_taylor)) =   0.11344260412736385826599342396561789E+52_kfp
            M(tblinv(48, 36, Lmax_taylor)) =   0.39355623229229834367370697954204835E+53_kfp
            M(tblinv(48, 40, Lmax_taylor)) =   0.30004515203011408745525386533152145E+55_kfp
            M(tblinv(48, 44, Lmax_taylor)) =   0.68834213748584856263825356451645833E+57_kfp
            M(tblinv(48, 48, Lmax_taylor)) =   0.20001957436429808127412432998753191E+61_kfp
          endif
          if (Lmax_taylor >= 50) then
            M(tblinv(50,  0, Lmax_taylor)) =    0.20864059126219213229101827182756468E+50_kfp
            M(tblinv(50,  4, Lmax_taylor)) =   -0.71371462465877191246460694120821817E+49_kfp
            M(tblinv(50,  8, Lmax_taylor)) =   -0.11390675800486528205469614339919759E+50_kfp
            M(tblinv(50, 12, Lmax_taylor)) =   -0.25698108938627681912829253813380236E+50_kfp
            M(tblinv(50, 16, Lmax_taylor)) =   -0.82682295714696801394294513691512288E+50_kfp
            M(tblinv(50, 20, Lmax_taylor)) =   -0.36119437965220962693352871174498619E+51_kfp
            M(tblinv(50, 24, Lmax_taylor)) =   -0.24469060499208437710499965690675910E+52_kfp
            M(tblinv(50, 28, Lmax_taylor)) =   -0.24985831554220241753398146184375358E+53_kfp
            M(tblinv(50, 32, Lmax_taylor)) =   -0.39487363485332672405775590478159064E+54_kfp
            M(tblinv(50, 36, Lmax_taylor)) =   -0.10818721669463859115111182196995084E+56_kfp
            M(tblinv(50, 40, Lmax_taylor)) =   -0.60006392362363866029775125016277996E+57_kfp
            M(tblinv(50, 44, Lmax_taylor)) =   -0.78921187860113696431951889503631055E+59_kfp
            M(tblinv(50, 48, Lmax_taylor)) =   -0.48554519650961009793033053599894219E+62_kfp
          endif

          if (Lmax_taylor > 50) then
            DEBUG_WARNING(*, 'load_lattice_coefficients(): Lmax_taylor > 50')
          endif

          call pepc_status('LATTICE COEFFICIENTS: finished')
        end subroutine load_lattice_coefficients


        !>
        !> Calculates the overall multipole expansion of the whole
        !> central box
        !>
        subroutine calc_omega_tilde(particles)
          use treevars, only: num_threads
          use module_pepc_types
          use module_debug
          implicit none

          type(t_particle), intent(in) :: particles(:)

          integer(kind_particle) :: p
          integer :: ierr
          
          omega_tilde = zero

          ! calculate multipole contributions of all local particles

          !$ call omp_set_num_threads(num_threads)
          !$OMP  PARALLEL DO DEFAULT(PRIVATE) SHARED(particles,LatticeCenter) SCHEDULE(RUNTIME) REDUCTION(+:omega_tilde)
          do p=1,size(particles, kind=kind(p))
            call addparticle(omega_tilde, particles(p)%x, particles(p)%data%q)
          end do
          !$OMP  END PARALLEL DO          
          !$ call omp_set_num_threads(1)
          
          ! extrinsic correction via fictitious charges according to [Kudin & Scuseria, ChemPhysLet 283, 61 (1998)] on rank 0
          if (fmm_extrinsic_correction == FMM_EXTRINSIC_CORRECTION_FICTCHARGE) then
            ! the fictitious charges are needed on all ranks for central-cell interaction
            ! but will only be added to the cell multipole expansion on rank 0
            
            nfictcharge = 0
            
            do p=1,3
              if (periodicity(p)) then
                nfictcharge = nfictcharge + 1
                fictcharge(nfictcharge)%coc(1:3) =   LatticeOrigin + Lattice(p, :)                                ! position
                fictcharge(nfictcharge)%charge   = - box_dipole(p) / sqrt(dot_product(Lattice(p,:),Lattice(p,:))) ! charge
              end if
            end do
            
            nfictcharge = nfictcharge + 1
            fictcharge(nfictcharge)%coc(1:3)     =  LatticeOrigin                                                 ! position
            fictcharge(nfictcharge)%charge       = -sum(fictcharge(1:nfictcharge-1)%charge)                       ! charge
            
            if (myrank==0) then
              do p=1,nfictcharge
                 call addparticle(omega_tilde, fictcharge(p)%coc, fictcharge(p)%charge)
              end do
            end if
          endif

          call chop(omega_tilde)

          ! sum multipole contributions from all processors - treat complex as two real numbers since e.g. complex*32 is not supported by mpi
          call MPI_ALLREDUCE(MPI_IN_PLACE, omega_tilde, 2*fmm_array_length_multipole, MPI_REAL_fmm, MPI_SUM, MPI_COMM_fmm, ierr)

          if ((myrank == 0) .and. dbg(DBG_PERIODIC)) then
            call WriteTableToFile('omega_tilde.tab', omega_tilde, Lmax_multipole)
          end if

          if (abs(omega_tilde( tblinv(0, 0, Lmax_multipole))) > 0.) then
            DEBUG_WARNING(*, 'The central box is not charge-neutral: Q_total=omega_tilde( tblinv(0, 0))=', omega_tilde( tblinv(0, 0, Lmax_multipole)), ' Setting to zero, resulting potentials might be wrong.' )
            omega_tilde( tblinv(0, 0, Lmax_multipole)) = (zero, zero) ! FIXME this line should be removed if zero_terms_..() functions are used
          end if
          
          ! FIXME untested : call zero_terms_multipole(omega_tilde)
          
          ! sum contributions from all processors

        contains
          subroutine addparticle(om, R, q)
            implicit none
            complex(kfp), intent(inout) :: om(1:fmm_array_length_multipole)
            real*8, intent(in) :: R(3)
            real*8, intent(in) :: q
            real(kfp) :: S(3)
            integer :: ll, mm
             
            S   = cartesian_to_spherical(R - LatticeCenter)

            do ll=0,Lmax_multipole
              do mm=0,ll
                om( tblinv(ll, mm, Lmax_multipole) ) = om( tblinv(ll, mm, Lmax_multipole) ) + omega(ll, mm, S, q)
              end do
            end do
          end subroutine
        end subroutine calc_omega_tilde


        !>
        !> Calculates the (cartesian) overall dipole moment
        !> \f$\sum_p q(p){\vec r}_p\f$ and the
        !> trace of the quadrupole matrix
        !> \f$\frac{2\pi}{3}\sum_p q(p){\vec r}_p\cdot{\vec r}_p\f$
        !> for performing the extrinsic-to-intrinsic correction
        !>
        subroutine calc_box_dipole(particles)
          use module_debug
          use module_pepc_types
          implicit none

          type(t_particle), intent(in) :: particles(:)

          real*8 :: r(3)

          integer(kind_particle) :: p
          integer :: ierr

          box_dipole = zero
          quad_trace = zero

          ! calculate dipole contributions of all local particles
          do p=1,size(particles, kind=kind(p))
            r = particles(p)%x - LatticeCenter
            box_dipole = box_dipole + particles(p)%data%q * r
            quad_trace = quad_trace + particles(p)%data%q * dot_product(r, r)
          end do
          
          ! sum contributions from all processors
          call MPI_ALLREDUCE(MPI_IN_PLACE, box_dipole, 3, MPI_REAL_fmm, MPI_SUM, MPI_COMM_fmm, ierr)
          call MPI_ALLREDUCE(MPI_IN_PLACE, quad_trace, 1, MPI_REAL_fmm, MPI_SUM, MPI_COMM_fmm, ierr)
        end subroutine calc_box_dipole

         
        !>
        !> Calculates the lattice contribution with respect to the
        !> centre of the original box
        !>
        subroutine calc_mu_cent(omega, mu)
          implicit none
          complex(kfp), intent(in) :: omega(1:fmm_array_length_multipole)
          complex(kfp), intent(out) :: mu(1:fmm_array_length_taylor)

          ! contribution of all outer lattice cells, with regards to the centre of the original box
          mu = M2L(MLattice, omega)

          call chop(mu)

          if ((myrank == 0) .and. dbg(DBG_PERIODIC)) then
            call WriteTableToFile('mu_cent.tab', mu_cent, Lmax_taylor)
          end if
        end subroutine calc_mu_cent


        !>
        !> Calculates P^~ as given by [Challacombe et al, eq. (7)]
        !> This is consistent with [Kudin, eq(2)],
        !> but contradicts with [White, eq. (2a)]
        !>
        !> uses negative order relation as given by
        !>  http://mathworld.wolfram.com/AssociatedLegendrePolynomial.html
        !>
        recursive real(kfp) function Ptilda(l, m, x) result(Pt)
          implicit none
          integer, intent(in) :: l, m
          real(kfp), intent(in) :: x

          if (m >= 0) then
            Pt = LegendreP(l, m, x) / factorial(l + m)
          else
            Pt = (-one)**m * Ptilda(l, abs(m), x)
          endif
        end function Ptilda


        !>
        !> Calculates P^~ as given by [Challocombe et al, eq. (8)]
        !> This is (almost) consistent with [Kudin, eq(3)],
        !> but contradicts with [White, eq. (2b)]
        !>
        !> uses negative order relation as given by
        !>  http://mathworld.wolfram.com/AssociatedLegendrePolynomial.html
        !>
        recursive real(kfp) function P2tilda(l, m, x) result(P2t)
          implicit none
          integer, intent(in) :: l, m
          real(kfp), intent(in) :: x

          if (m >= 0) then
            P2t = LegendreP(l, m, x) * factorial(l - m)
          else
            P2t = (-one)**abs(m) * P2tilda(l, abs(m), x)
          endif
        end function P2tilda


        !>
        !> Calculates the chargeless moments of the multipole expansion
        !> for a certain particle (position given in spherical coordinates
        !>
        complex(kfp) function OMultipole(l, m, s)
          implicit none
          integer, intent(in) :: l, m
          real(kfp), intent(in) :: s(3)

          complex(kfp), parameter :: i = (zero,one)

          if ((l<0) .or. (m<-l) .or. (m>l)) then
              OMultipole = 0 ! these are zero per definition
          else
              if (s(1) > 0 .or. l>0) then
                OMultipole = s(1)**real(l, kind=kfp) * &
                             Ptilda(l, m, s(2)) * exp(-i*real(m,kind=kfp)*s(3) )
              else
                ! avoid having to compute 0^0 = 1 since some runtimes do not like that
                OMultipole = Ptilda(l, m, s(2)) * exp(-i*real(m,kind=kfp)*s(3) )
              endif

          endif
        end function OMultipole


        !>
        !> Calculates the chargeless moments of the Taylor expansion
        !> for a certain particle (coordinates given in spherical system)
        !>
        complex(kfp) function MTaylor(l, m, s)
          implicit none
          integer, intent(in) :: l, m
          real(kfp), intent(in) :: s(3)

          complex(kfp), parameter :: i = (zero,one)

          if ((l<0) .or. (m<-l) .or. (m>l)) then
              MTaylor = 0 ! these are zero per definition
          else
              MTaylor = one/(s(1)**real(l+1,kind=kfp)) * &
                        P2tilda(l, m, s(2)) * exp( i*real(m,kind=kfp)*s(3) )
          endif
        end function MTaylor


        !>
        !> Calculates the charged moments of the multipole expansion
        !> for a certain particle (position in spherical coordinates)
        !>
        complex(kfp) function omega(l, m, s, q)
          implicit none
          integer, intent(in) :: l, m
          real(kfp), intent(in) :: s(3)
          real*8, intent(in) :: q

          omega = real(q, kind=kfp) * OMultipole(l, m, s)
        end function omega


        !>
        !> Calculates the charged moments of the Taylor expansion
        !> for a certain particle (coordinates in spherical system)
        !>
        complex(kfp) function mu(l, m, s, q)
          implicit none
          integer, intent(in) :: l, m
          real(kfp), intent(in) :: s(3)
          real*8, intent(in) :: q

          mu = real(q, kind=kfp) * MTaylor(l, m, s)
        end function mu


        !>
        !> Calculates force at individual position that results
        !> from mirror boxes beyond the near field region,
        !> i.e. the lattice contribution
        !>
        subroutine fmm_sum_lattice_force(pos, e_lattice, phi_lattice)
          use module_mirror_boxes, only : num_neighbour_boxes, lattice_vect, neighbour_boxes
          use module_coulomb_kernels, only : calc_force_coulomb_3D_direct
          implicit none

          real*8, intent(in) :: pos(3)
          real*8, intent(out) ::  e_lattice(3), phi_lattice

          complex(kfp) :: mu_shift(1:fmm_array_length_taylor), O_R(1:fmm_array_length_multipole)
          integer :: l, m
          real*8 :: R(3)
          real(kfp) :: S(3)
          real(kfp) :: prefact
          integer :: p, ibox
          real(kfp), parameter :: pi=acos(-one)
          real*8 :: etmp(3), phitmp, delta(3)
          
          prefact = two*pi/(three*unit_box_volume)

          if (.not. do_periodic) then
              e_lattice   = 0
              phi_lattice = 0
            else
              ! shift mu_cent to the position of our particle
              R        = pos - LatticeCenter
              S        = cartesian_to_spherical(R)

              do l = 0,Lmax_multipole
                do m = 0,l
                  O_R(tblinv(l, m, Lmax_multipole)) = OMultipole(l, m, S)
                end do
              end do

              mu_shift = L2L(O_R, mu_cent, 1)

              ! E = -grad(Phi)
              e_lattice(1) = - real(tbl(mu_shift, 1, 1, Lmax_taylor))
              e_lattice(2) = -aimag(tbl(mu_shift, 1, 1, Lmax_taylor))
              e_lattice(3) = - real(tbl(mu_shift, 1, 0, Lmax_taylor))                            
              phi_lattice  =   real(tbl(mu_shift, 0, 0, Lmax_taylor))

              select case (fmm_extrinsic_correction)
                case (FMM_EXTRINSIC_CORRECTION_NONE)
                  ! nothing to do here
                case (FMM_EXTRINSIC_CORRECTION_REDLACK)
                  e_lattice   = e_lattice   + two*prefact * box_dipole
                  phi_lattice = phi_lattice - two*prefact * dot_product(R, box_dipole) + prefact * quad_trace
                case (FMM_EXTRINSIC_CORRECTION_FICTCHARGE)
                  do p=1,nfictcharge
                    ! interact with fictcharge(p)
                    ! we loop over all vbox-vectors, in fact we are only interested in the surface charges since the others cancel anyway, but
                    ! the exception for this is too complicated for now - FIXME: correct this
                    do ibox = 1,num_neighbour_boxes ! sum over all boxes within ws=1
                      delta = pos - lattice_vect(neighbour_boxes(:,ibox)) - fictcharge(p)%coc
                      call calc_force_coulomb_3D_direct(fictcharge(p), delta, dot_product(delta, delta), etmp, phitmp)
                      e_lattice   = e_lattice   + etmp
                      phi_lattice = phi_lattice + phitmp
                    end do
                  end do
                case (FMM_EXTRINSIC_CORRECTION_MEASUREMENT)
                  DEBUG_ERROR('("fmm_extrinsic_correction == FMM_EXTRINSIC_CORRECTION_MEASUREMENT currently not supported")')
                case default
                  DEBUG_ERROR('("fmm_extrinsic_correction == ", I0, " not supported")', fmm_extrinsic_correction)
              end select

          end if
        end subroutine fmm_sum_lattice_force


        !>
        !> Table access function for giving arbitrary l and m
        !>
        !> The function cares for picking the right indices and
        !> respects symmetry
        !>
        !> @param[in] l multipole order
        !> @param[in] m
        !> @param[in] A table
        !>
        complex(kfp) function tbl(A, l, m, Lmax)
          implicit none
          integer, intent(in) :: l, m, Lmax
          complex(kfp), intent(in) :: A(*)

          if ((l<0)) then
            DEBUG_ERROR('("tbl(A,l,m) - invalid arguments. l=", I0, " m=", I0)', l, m)
          endif

          if ((l>Lmax) .or. (abs(m)>l)) then
            tbl = 0
          else
            tbl = A( tblinv(l, abs(m), Lmax) )

            if (m<0) tbl = (-one)**m * conjg(tbl)
          end if
        end function tbl


        !>
        !> Table access function for giving arbitrary l and m
        !>
        !> Calculates the flat index, where an entry for l and m
        !> has to be stored in a one-dimensional array
        !>
        !> @param[in] l
        !> @param[in] m
        !>
        integer function tblinv(l,m,Lmax)
          use module_debug
          implicit none
          integer, intent(in) :: l, m,Lmax

          if ((l<0) .or. (m<0) .or. (m>l) .or. (l>Lmax)) then
            DEBUG_ERROR('("tblinv(l,m) - invalid arguments. l=", I0, " m=", I0)', l, m)
          endif

          tblinv = l*(l+1)/2 + 1 + m
        end function tblinv


        !>
        !> M2M-Operator (denoted with \f$\triangleleft\f$ )
        !> eq. (6) in [Kudin]
        !>
        function M2M(O_a, O_b)
          implicit none
          complex(kfp), intent(in) :: O_a(1:fmm_array_length_multipole)
          complex(kfp), intent(in) :: O_b(1:fmm_array_length_multipole)
          complex(kfp), dimension(1:fmm_array_length_multipole) :: M2M

          complex(kfp) :: t
          integer :: l, m, j, k

          do l = 0,Lmax_multipole
            do m = 0,l

              t = zero

              do j = 0,l ! Attention, this sum only runs to l instead of infty|Lmax
                do k = -j,j
                  t = t + tbl(O_b, l-j, m-k, Lmax_multipole) * tbl(O_a, j, k, Lmax_multipole)
                end do
              end do

              M2M( tblinv(l, m, Lmax_multipole) ) = t
            end do
          end do

          call chop(M2M)
        end function M2M


        !>
        !> M2L-Operator (denoted with \f$\otimes\f$ )
        !> eq. (7) in [Kudin]
        !>
        function M2L(M_b, O_a)
          implicit none
          complex(kfp), intent(in) :: M_b(1:fmm_array_length_taylor)
          complex(kfp), intent(in) :: O_a(1:fmm_array_length_multipole)
          complex(kfp), dimension(1:fmm_array_length_taylor) :: M2L

          complex(kfp) :: t
          integer :: l, m, j, k

          do l = 0,Lmax_taylor
            do m = 0,l

              t = 0_kfp

              do j = 0,Lmax_taylor ! should be infinity, but see last page of [Kudin]: there they state, that (p-l) is enough
                do k = -j,j
                  t = t + (-1)**j * tbl(M_b, j+l, k+m, Lmax_taylor) * tbl(O_a, j, k, Lmax_multipole) !TODO: evtl. (-1)**j hinzufuegen (vgl. S. 82 bei Ivo) ??
                end do
              end do

              M2L( tblinv(l, m, Lmax_taylor) ) = t
            end do
          end do

          call chop(M2L)
        end function M2L


        !>
        !> L2L-Operator (denoted with \f$\triangleright\f$ )
        !> eq. (8) in [Kudin]
        !>
        function L2L(O_b, M_r, max_l)
          implicit none
          complex(kfp), intent(in) :: O_b(1:fmm_array_length_multipole)
          complex(kfp), intent(in) :: M_r(1:fmm_array_length_taylor)
          complex(kfp), dimension(1:fmm_array_length_taylor) :: L2L
          integer, intent(in), optional :: max_l

          complex(kfp) :: t
          integer :: l, m, j, k
          integer :: maxl

          if (present(max_l)) then
            maxl = max_l
          else
            maxl = Lmax_taylor
          endif

          do l = 0,maxl
            do m = 0,l

              t = 0_kfp

              do j = l,Lmax_multipole ! Attention, this sum starts at l since O_b(j-l<0,..) = 0 anyway
                do k = -j,j
                  t = t + tbl(O_b, j-l, k-m, Lmax_multipole) * tbl(M_r, j, k, Lmax_taylor)
                end do
              end do

              L2L( tblinv(l, m, Lmax_taylor) ) = t
            end do
          end do

          call chop(L2L)
        end function L2L


        !>
        !> Scaling Operator \f$\mathcal{U}_L\f$ for Taylor coefficients
        !> @param[in] L table with Taylor coefficients
        !>
        function UL(L)
          implicit none
          complex(kfp), intent(in) :: L(1:fmm_array_length_taylor)
          complex(kfp), dimension(1:fmm_array_length_taylor) :: UL
          real(kfp) :: scalefac

          integer :: ll, mm

          scalefac = real(2*ws+1, kind=kfp)

          do ll = 0,Lmax_taylor
            do mm = 0,ll
              UL( tblinv(ll, mm, Lmax_taylor) ) = L( tblinv(ll, mm, Lmax_taylor) ) / scalefac**real(ll+1, kind=kfp)
            end do
          end do

          call chop(UL)
        end function UL


        !>
        !> Formal summation of \f$M\f$ over NF, ie all (27) neighbouring boxes
        !> with some overhead to avoid numerical elimination of small values
        !>
        complex(kfp) function MstarFunc(l,m)
          implicit none
          integer, intent(in) :: l, m
          real(kfp) :: rpart(num_neighbour_boxes), ipart(num_neighbour_boxes),rp,ip, s(3)
          complex(kfp) :: tmp
          complex(kfp), parameter :: ic = (zero,one)

          integer :: i

          do i = 1,num_neighbour_boxes
             s = cartesian_to_spherical(-lattice_vect(neighbour_boxes(:,i)))

             tmp      = OMultipole(l, m, s)
              ! we store the summands and order them before performing the sum
              ! to avoid numeric elimination
             rpart(i) =       real(tmp,  kind=kfp)
             ipart(i) = real(aimag(tmp), kind=kfp)
          end do

          call sort_abs(rpart)
          call sort_abs(ipart)

          rp = zero
          ip = zero

          do i = 1, num_neighbour_boxes,1 ! we sum up all contributions, starting from smallest
            rp = rp + rpart(i)
            ip = ip + ipart(i)
          end do

          MstarFunc = rp + ic*ip ! do not use cmplx()-function since it yields wrong results with complex*32 types
        end function MstarFunc


        !>
        !> Formal summation of \f$M\f$ over FF`, ie a lot of boxes
        !> with some overhead to avoid numerical elimination of small values
        !>
        complex(kfp) function LstarFunc(l,m)
          implicit none
          integer, intent(in) :: l, m
          real(kfp) :: rpart(num_neighbour_boxes*(num_neighbour_boxes-1)), ipart(num_neighbour_boxes*(num_neighbour_boxes-1)),rp,ip, s(3)
          complex(kfp) :: tmp
          complex(kfp), parameter :: ic = (zero,one)
          integer :: i, ii, k

          k = 0

          ! sum over all boxes within FF' (cells in the far field of the central cell but in the near field of the central supercell 3x3x3 that embeds cell {0,0,0} in the center)
          do i = 1,num_neighbour_boxes-1 ! central box is being omitted in this loop
            do ii = 1,num_neighbour_boxes
              s   = cartesian_to_spherical((2*ws+1)*lattice_vect(neighbour_boxes(:,i)) + lattice_vect(neighbour_boxes(:,ii)))
              tmp = MTaylor(l, m, s)
              k   = k+1
              ! we store the summands and order them before performing the sum
              ! to avoid numeric elimination
              rpart(k) =       real(tmp,  kind=kfp)
              ipart(k) = real(aimag(tmp), kind=kfp)
            end do
          end do

          call sort_abs(rpart)
          call sort_abs(ipart)

          rp = zero
          ip = zero

          do i = 1,k ! we sum up all contributions, starting from smallest
            rp = rp + rpart(i)
            ip = ip + ipart(i)
          end do

          LStarFunc = rp + ic*ip ! do not use cmplx()-function since it yields wrong results with complex*32 types
        end function LstarFunc


        !>
        !> Converts cartesian coordinates to spherical system
        !> @param[in]  cartesian  cartesian vector [x, y, z]
        !>
        function cartesian_to_spherical(cartesian)
          implicit none
          real*8, intent(in)  :: cartesian(3)
          real(kfp), dimension(3) :: cartesian_to_spherical, c

          c = real(cartesian, kind=kfp)

          cartesian_to_spherical(1) = sqrt(dot_product(c, c))

          if (cartesian_to_spherical(1) == 0) then
            cartesian_to_spherical(2) = one
          else
            cartesian_to_spherical(2) = c(3) / cartesian_to_spherical(1)
          end if

          if ((c(1) == 0) .and. (c(2) == 0)) then
            cartesian_to_spherical(3) = zero
          else
            cartesian_to_spherical(3) = atan2(c(2), c(1))
          end if
        end function cartesian_to_spherical


        !>
        !> Writes contents of table T to a output stream s in a structured way
        !>
        subroutine PrintTable(s, T, Lmax)
          implicit none
          complex(kfp), intent(in) :: T(:)
          integer, intent(in) :: Lmax
          integer, intent(in) :: s

          integer :: ll,mm,idx

          idx = 0

          do ll = 0,Lmax
            do mm = 0,ll
              idx = idx + 1

              write(s,'(I6, I6, I6, D50.35, D50.35)') idx, ll, mm, tbl(T, ll, mm, Lmax)
            end do
          end do
        end subroutine PrintTable


        !>
        !> Writes contents of table T to file s
        !>
        subroutine WriteTableToFile(s, T, Lmax)
          implicit none
          complex(kfp), intent(in) :: T(:)
          integer, intent(in) :: Lmax
          character(len=*) :: s
          integer, parameter :: temp_file = 60

          open(temp_file, file=trim(s))

          write(temp_file,*) "# filename    = ", trim(s)
          write(temp_file,*) "# t_lattice_1 = ", t_lattice_1
          write(temp_file,*) "# t_lattice_2 = ", t_lattice_2
          write(temp_file,*) "# t_lattice_3 = ", t_lattice_3
          write(temp_file,*) "# periodicity_switches = ", periodicity_switches
          write(temp_file,*) "# Lmax_multipole = ", Lmax_multipole
          write(temp_file,*) "# Lmax_taylor    = ", Lmax_taylor

          write(temp_file,*) "# Lmax    = ", Lmax
          write(temp_file,*) "# MaxIter = ", MaxIter
          write(temp_file,*) "##########################"
          write(temp_file,*) "  idx     l     m           real-part                                         imaginary-part"

          call PrintTable(temp_file, T, Lmax)

          close(temp_file)
        end subroutine WriteTableToFile


        !>
        !> Dumps all parameters to the stream ifile
        !>
        subroutine fmm_framework_param_dump(ifile)
          implicit none
          integer, intent(in) :: ifile

          if (myrank .ne. 0) return

          write(ifile,*) "LATTICE: ------------- Lattice fmm-framework switches ----------------"
          write(ifile,*) "LATTICE: Lmax_multipole = ", Lmax_multipole
          write(ifile,*) "LATTICE: Lmax_taylor    = ", Lmax_taylor
          write(ifile,*) "LATTICE: MaxIter        = ", MaxIter
          write(ifile,*) "LATTICE: ws             = ", ws
          write(ifile,*) "LATTICE: t_lattice_1    = ", t_lattice_1
          write(ifile,*) "LATTICE: t_lattice_2    = ", t_lattice_2
          write(ifile,*) "LATTICE: t_lattice_3    = ", t_lattice_3
          write(ifile,*) "LATTICE: periodicity    = ", periodicity
          write(ifile,*) "LATTICE: # neighbours    = ", num_neighbour_boxes
          write(ifile,*)
        end subroutine fmm_framework_param_dump


        !>
        !> Computes the associated Legendre polynomial \f$P_{lm}(x)\f$.
        !> Here m and l are integers satisfying  \f$0 \leq m \leq l\f$,
        !> while x lies in the range \f$-1 \leq x \leq 1\f$.
        !>
        !> Code fragment for \f$P_l^m(x)\f$ taken from
        !>
        !> Numerical Recipes in Fortran 77: The Art of Scientific Computing
        !>              (ISBN 0-521-43064-X)
        !> pg. 246ff
        !>
        !> and modified to give \f$P_{lm}(x)\f$:
        !> \f$ P_{lm}(x) = (-1)^m P_l^m (x) \f$, see
        !>
        !> Abramowitz and Stegun: Handbook of Mathematical Functions
        !> Section 8. Legendre Functions (pg. 332)
        !>
        real(kfp) function LegendreP(l,m,x)
          use module_debug
          implicit none
          integer, intent(in) :: l, m
          real(kfp) ::x

          integer :: i,ll
          real(kfp) :: fact,pll,pmm,pmmp1,somx2

          pll = zero

          if ( (m < 0) .or. (m > l) .or. (abs(x) > 1) .or.  (l<0)) then
            DEBUG_ERROR(*,'Invalid arguments for LegendreP(',l,m,x,')')
          endif

          pmm = one     ! Compute P_m^m

          if (m > 0) then
            somx2 = sqrt((one-x)*(one+x))
            fact  = one

            do i = 1,m
               pmm  = -pmm * fact * somx2
               fact = fact+two
            enddo
          endif

          if (l == m) then
            LegendreP = pmm
          else
            pmmp1 = x*(two*m+one)*pmm  ! Compute P_m+1^m

            if (l == m+1) then
              LegendreP = pmmp1
            else                  ! Compute P_l^m , l > m + 1
              do ll = m+2,l
                pll   = (x*(two*ll-one)*pmmp1-(ll+m-one)*pmm)/real(ll-m,kind=kfp)
                pmm   = pmmp1
                pmmp1 = pll
              enddo

              LegendreP = pll
            endif
          endif

          LegendreP = (-one)**m * LegendreP
        end function LegendreP


        !>
        !> Calculates the factorial of the argument
        !>
        real(kfp) function factorial(n)
            use module_debug
            implicit none
            integer, intent(in) :: n

            if (n<0) then
              DEBUG_ERROR(*,"Tried to calculate factorial of negative argument.")
            end if

            if (n>170) then
              DEBUG_ERROR(*,"Tried to calculate factorial with n>170. This would lead to numeric overflow.")
            end if

            select case (n)
              case ( 0)
                factorial =                            one
              case ( 1)
                factorial =                            one
              case ( 2)
                factorial =                            two
              case ( 3)
                factorial =                            6._kfp
              case ( 4)
                factorial =                           24._kfp
              case ( 5)
                factorial =                          120._kfp
              case ( 6)
                factorial =                          720._kfp
              case ( 7)
                factorial =                         5040._kfp
              case ( 8)
                factorial =                        40320._kfp
              case ( 9)
                factorial =                       362880._kfp
              case (10)
                factorial =                      3628800._kfp
              case (11)
                factorial =                     39916800._kfp
              case (12)
                factorial =                    479001600._kfp
              case (13)
                factorial =                   6227020800._kfp
              case (14)
                factorial =                  87178291200._kfp
              case (15)
                factorial =                1307674368000._kfp
              case (16)
                factorial =               20922789888000._kfp
              case (17)
                factorial =              355687428096000._kfp
              case (18)
                factorial =             6402373705728000._kfp
              case (19)
                factorial =           121645100408832000._kfp
              case (20)
                factorial =          2432902008176640000._kfp
              case (21)
                factorial =         51090942171709440000._kfp
              case (22)
                factorial =       1124000727777607680000._kfp
              case (23)
                factorial =      25852016738884976640000._kfp
              case (24)
                factorial =     620448401733239439360000._kfp
              case default
                factorial = gamma(real(n+1, kfp))
            end select
        end function factorial


        !>
        !> Sorts the given values with a heap sort approach
        !> in order of ther absolute value
        !> compare (Numerical Recipes f90, p1171)
        !>
          subroutine sort_abs(arr)
            implicit none
            real(kfp), intent(inout) :: arr(:)
            integer :: i,n

            n = size(arr)

            do i=n/2,1,-1                      ! Left range - hiring phase (heap creation)
               call sift_down(i,n)
            end do

            do i=n,2,-1                        ! Right range of sift-down is decr. from n-1 ->1
               ! during retirement/promotion (heap selection) phase.
               call swap( arr(1),arr(i) )      ! Clear space at end of array and retire top of heap into it
               call sift_down( 1,i-1)
            end do

          contains
            subroutine sift_down(l,r)        ! Carry out the sift-down on element arr(l) to maintain
              integer, intent(in) :: l,r     ! the heap structure
              integer :: j,jold    ! index
              real(kfp) :: a

              a    = arr(l)
              jold = l
              j    = l + l
              do                   ! do while j <= r
                 if (j > r) exit
                 if (j < r) then
                   if (abs(arr(j)) < abs(arr(j+1))) j = j+1
                 endif
                 if (abs(a) >= abs(arr(j))) exit       ! Found a`s level, so terminate sift-down

                 arr(jold) = arr(j)
                 jold      = j                    ! Demote a and continue
                 j         = j+j
              end do
              arr(jold) = a                  ! Put a into its slot
            end subroutine sift_down

            subroutine swap(p,q)
                real(kfp) :: p,q, dum
                dum = p
                p   = q
                q   = dum
            end subroutine swap
          end subroutine sort_abs


        !>
        !>  Sets all matrix entries that are smaller than 1.e-16 to 0.
        !> (separately for real and imaginary part)
        !> This is the same as Mathematicas Chop[]-function
        !>
        subroutine chop(a)
          implicit none
          complex(kfp), intent(inout) :: a(:)
          integer :: i
          complex(kfp), parameter :: ic = (zero,one)
          real(kfp) :: re, im

          if (.not. chop_arrays) return

          DEBUG_WARNING(*, 'chopping some array')

          do i=1,size(a)
            re =       real(a(i),  kind=kfp)
            im = real(aimag(a(i)), kind=kfp)

            if (abs(re) < prec) re = zero
            if (abs(im) < prec) im = zero

            a(i) = re + ic*im ! do not use cmplx()-function here (see above)
          end do
        end subroutine chop
end module module_fmm_framework





