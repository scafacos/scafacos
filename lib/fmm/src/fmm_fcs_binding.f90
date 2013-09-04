      module fmm_fcs_binding
       use fmmkinds
       use mwigner
       type FMM_internal_params_t
!         sequence ! not possible, since c_ptr is involved (itself not sequenced)
! status error analysis
         integer(kind=fmm_integer) :: serroranalysis
         integer(kind=fmm_integer) :: nerroranalysis
! fmm internal parameter
         integer(kind=fmm_integer) :: pgd
         integer(kind=fmm_integer) :: ws
         integer(kind=fmm_integer) :: depth
         integer(kind=fmm_integer) :: nmultipoles
         integer(kind=fmm_integer) :: ncsar
         integer(kind=fmm_integer) :: parabola
! nearfield potential
         integer(kind=fmm_integer) :: ilinearpotentialsv
         integer(kind=fmm_integer) :: iplummerpotentialsv
! loadbalancing
         integer(kind=fmm_integer) :: iboxloadlength
         integer(kind=fmm_integer) :: percentageofimbalance
         integer(kind=fmm_integer) :: minofimbalance
         integer(kind=fmm_integer) :: maxofimbalance
         integer(kind=fmm_integer) :: doload
         logical(kind=fmm_logical) :: firsterroranalysis
         logical(kind=fmm_logical),dimension(0:100) :: hugep
         real(kind=fmm_real) :: fracdepth
         real(kind=fmm_real) :: shmonopole
         real(kind=fmm_real),dimension(3) :: linearodistancesv
         real(kind=fmm_real) :: aoplummersv
         real(kind=fmm_real),dimension(100) :: hugef
         type(twignerd) :: wignerd
         type(c_ptr) :: iboxloadcptr
         
         integer(kind=fmm_integer) :: presorted
         integer(kind=fmm_integer) :: resort
         type(c_ptr) :: resort_ptr
       end type FMM_internal_params_t
      end module fmm_fcs_binding

