      module fmm_hooks_mod

      use iso_c_binding, only : c_ptr, c_null_ptr
      implicit none
      type(c_ptr) :: fmm_hooks = c_null_ptr

        
      interface

        type(c_ptr) function fcs_fmm_hooks_create() bind(c, Name='fcs_fmm_hooks_create')
        use iso_c_binding, only : c_ptr
        implicit none
        end function fcs_fmm_hooks_create

        subroutine fcs_fmm_hooks_destroy(hooks) bind(c, Name='fcs_fmm_hooks_destroy')
        use iso_c_binding, only : c_ptr
        implicit none
        type(c_ptr), value :: hooks
        end subroutine fcs_fmm_hooks_destroy

        subroutine fcs_fmm_hooks_near_start(hooks) bind(c, Name='fcs_fmm_hooks_near_start')
        use iso_c_binding, only : c_ptr
        implicit none
        type(c_ptr), value :: hooks
        end subroutine fcs_fmm_hooks_near_start

        subroutine fcs_fmm_hooks_near_stop(hooks) bind(c, Name='fcs_fmm_hooks_near_stop')
        use iso_c_binding, only : c_ptr
        implicit none
        type(c_ptr), value :: hooks
        end subroutine fcs_fmm_hooks_near_stop

        subroutine fcs_fmm_hooks_far_start(hooks) bind(c, Name='fcs_fmm_hooks_far_start')
        use iso_c_binding, only : c_ptr
        implicit none
        type(c_ptr), value :: hooks
        end subroutine fcs_fmm_hooks_far_start

        subroutine fcs_fmm_hooks_far_stop(hooks) bind(c, Name='fcs_fmm_hooks_far_stop')
        use iso_c_binding, only : c_ptr
        implicit none
        type(c_ptr), value :: hooks
        end subroutine fcs_fmm_hooks_far_stop

        subroutine fcs_fmm_hooks_print(hooks) bind(c, Name='fcs_fmm_hooks_print')
        use iso_c_binding, only : c_ptr
        implicit none
        type(c_ptr), value :: hooks
        end subroutine fcs_fmm_hooks_print

      end interface

      end module fmm_hooks_mod
