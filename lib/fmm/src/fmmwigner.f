      module mwigner
       use fmmkinds
       implicit none
       type twignerd
        sequence
        real(kind=fmm_real), pointer:: wignerd(:,:,:,:,:)
       end type twignerd
      end module mwigner
