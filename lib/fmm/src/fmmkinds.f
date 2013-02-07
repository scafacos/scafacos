#include "fmm.h"
#include "fmmkinds.h"

      module fmmkinds
#ifdef FMM_ISO_C_BINDING
       use iso_c_binding
#endif
       implicit none
       integer, parameter:: fmm_real=FMM_REAL
       integer, parameter:: fmm_real_extended=FMM_REAL_EXTENDED
       integer, parameter:: fmm_real_itor=FMM_REAL_ITOR
       integer, parameter:: fmm_integer=FMM_INTEGER
#ifdef FMM_ISO_C_BINDING
       integer, parameter:: fmm_c_integer=FMM_C_INTEGER
       integer, parameter:: fmm_c_real=FMM_C_REAL
#endif
#ifdef FMM_COMPRESSION
#ifdef FMM_EXTREMECOMPRESSION
#ifdef FMM_EXTREMEEXTREMECOMPRESSION
       integer, parameter:: fmm_xyz_to_integer=FMM_XYZ_TO_INTEGER
#endif
#endif
#endif
       integer, parameter:: fmm_testalloc_integer=FMM_TESTALLOC_INTEGER
       integer, parameter:: fmm_logical=FMM_LOGICAL
       integer, parameter:: fmm_pointer=FMM_POINTER
#ifdef FMM_ISO_C_BINDING
#ifdef FMM_ALLOCALIGNED
       integer, parameter:: fmm_alignment=FMM_ALIGNMENT
#endif
#endif
#ifdef FMM_PARALLEL
       integer, parameter:: mp_real_max=FMM_MP_REAL_MAX
       integer, parameter::mp_integer_processes=FMM_MP_INTEGER_PROCESSES
#endif
       type tsndibox
        integer(kind=fmm_integer), pointer:: sndibox(:)
       end type tsndibox

       type tsndomegatree
        real(kind=fmm_real), pointer:: sndomegatree(:,:)
       end type tsndomegatree

       type tsndq
        real(kind=fmm_real), pointer:: sndq(:)
       end type tsndq

       type tsndxyz
        real(kind=fmm_real), pointer:: sndxyz(:,:)
       end type tsndxyz

#ifdef FMM_LOADSORT
       type tsndiboxload
        real(kind=fmm_real), pointer:: sndiboxload(:)
       end type tsndiboxload
#endif
#ifdef FMM_UNIFORMGRID
       type tsnduniformgrid
        real(kind=fmm_real), pointer:: snduniformgrid(:,:)
       end type tsnduniformgrid
#endif
       type telem
        integer(kind=fmm_integer) pos,val
        type(telem), pointer:: prev,next
       end type telem
      end module fmmkinds
