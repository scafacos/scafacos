
      subroutine fmm_init(FMM_internal_params)
      use fmmkinds
      use fmm_fcs_binding
      implicit none
      type(FMM_internal_params_t), intent(inout) :: FMM_internal_params
      
        FMM_internal_params%firsterroranalysis = .true.
        FMM_internal_params%depth = 0
        !FMM_internal_params%nerroranalysis = huge(1_8)
      end subroutine fmm_init
