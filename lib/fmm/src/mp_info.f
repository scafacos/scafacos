#include "fmm.h"

#ifdef FMM_PARALLEL
      module mp_info
       use fmmkinds
       use mp_wrapper
       implicit none
       integer(kind=fmm_integer) me,nnodes,gbinfo3,gbinfo4,gbinfo5,
     . gbinfo7,gbinfo8,kind_integer_processes
#ifdef FMM_COMPRESSION
       integer(kind=1), allocatable:: gbinfo(:,:,:)
#else
       integer(kind=fmm_integer), allocatable:: gbinfo(:,:)
#endif
       integer(kind=mp_integer_processes), allocatable:: idp(:)
       logical(kind=fmm_logical) mp_setup
       type(c_ptr), allocatable, target:: gbpt(:),gbsndrcvol(:),
     . gbsndrcvbuffer(:),gbsndibox(:),gbsndromegatree(:),
     . gbsndiomegatree(:),gbsndq(:),gbsndxyz(:)
#ifdef FMM_LOADSORT
       type(c_ptr), allocatable, target:: gbsndiboxload(:)
#endif
#ifdef FMM_UNIFORMGRID
       type(c_ptr) gbptrcvuniformgrid
       type(c_ptr), allocatable, target:: gbsnduniformgrid(:)
#endif
       type(c_ptr) gbptsndibox,gbptsndromegatree,gbptsndiomegatree
       type(c_ptr) gbptsndrmutree,gbptsndimutree,gbptsndmutree
       type(c_ptr) gbptsndxyz,gbptsndpot,gbptsndgrad,gbptrcvpot,
     . gbptrcvgrad,gbptgbinfo7
       type(c_ptr) ptscr1,ptscr2
#ifdef FMM_LOADSORT
       type(c_ptr) gbptsndiboxload,gbptrcviboxload,gbptscr
       type(c_ptr) gbpt5sndiboxload,gbpt5scr
#endif
#ifdef FMM_UNIFORMGRID
       type(c_ptr) gbptsnduniformgrid
#endif
       type(tsndibox), allocatable, volatile:: psndibox(:)
       type(tsndomegatree), allocatable, volatile:: psndromegatree(:),
     . psndiomegatree(:)
       type(tsndq), allocatable, volatile:: psndq(:)
       type(tsndxyz), allocatable, volatile:: psndxyz(:)
#ifdef FMM_LOADSORT
       type(tsndiboxload), allocatable, volatile:: psndiboxload(:)
#endif
#ifdef FMM_UNIFORMGRID
       type(tsnduniformgrid), allocatable, volatile:: psnduniformgrid(:)
#endif
#ifdef FMM_DEBUG
       integer(kind=fmm_integer),allocatable:: notifycount(:),
     . waitcount(:)
#endif
      end module mp_info
#endif