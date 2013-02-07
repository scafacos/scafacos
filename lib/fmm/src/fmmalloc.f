#include "fmm.h"
#include "fmmkinds.h"

      module fmmalloc
c
       use fmmkinds
c
       implicit none
c
       private
c
       integer(kind=fmm_integer), public:: maxallocation,corememory,
     . sortmemory
c
       integer(kind=fmm_integer), public:: nalloc,maxnalloc,nallocr
#ifdef FMM_PARALLEL
       integer(kind=fmm_integer), public:: nmp_alloc,maxnmp_alloc
#endif
       integer(kind=fmm_integer), public:: rtob,rextendedtob,ritortob,
     . itob,ltob
#ifdef FMM_PARALLEL
       integer(kind=fmm_integer), public:: iptob
#endif
#ifdef FMM_PARALLEL
       integer(kind=fmm_integer), public:: gbpttob,typeinttob,
     . typerealtob,typeqtob,typexyztob
#ifdef FMM_LOADSORT
       integer(kind=fmm_integer), public:: typeint2tob
#endif
#ifdef FMM_UNIFORMGRID
       integer(kind=fmm_integer), public:: typeuniformgridtob
#endif
       integer(kind=fmm_integer), public:: typeelemtob
#endif
c
       public:: fmmallocate,fmmdeallocate,fmmallocatept,fmmdeallocatept
#ifdef FMM_PARALLEL
       public:: mp_fmmallocate,mp_fmmdeallocate
#endif
c
#ifdef FMM_ISO_C_BINDING
#ifdef FMM_ALLOCALIGNED
       public:: fmmallocate_aligned,fmmdeallocate_aligned
#endif
#endif
       public:: indallocate,indallocatept,inddeallocate,inddeallocatept
       public:: dbpallocate,dbpdeallocate
       public:: srtallocate,srtdeallocate
c
       interface fmmallocate
        module procedure fmmallocdbl1d
#if FMM_REAL != FMM_REAL_EXTENDED
        module procedure fmmallocdblextended1d
#endif
#if FMM_REAL != FMM_REAL_ITOR && FMM_REAL_EXTENDED != FMM_REAL_ITOR
        module procedure fmmallocdblitor1d
#endif
        module procedure fmmallocdbl2d,fmmallocdbl3d
        module procedure fmmallocdbl4d,fmmallocdbl5d
        module procedure fmmallocint1d,fmmallocint2d
        module procedure fmmallocint3d,fmmallocint4d
        module procedure fmmalloclog1d,fmmalloclog2d
        module procedure fmmallocbyte,fmmallocbyte2d,fmmallocbyte3d
#ifdef FMM_PARALLEL
#if FMM_INTEGER != FMM_MP_INTEGER_PROCESSES
        module procedure fmmallocintp
#endif
        module procedure fmmallocgpt,fmmalloctypeint,fmmalloctypereal
        module procedure fmmalloctypeq,fmmalloctypexyz
#ifdef FMM_LOADSORT
        module procedure fmmalloctypeint2
#endif
#ifdef FMM_UNIFORMGRID
        module procedure fmmalloctypeuniformgrid
#endif
        module procedure fmmalloctypeelem
#endif
       end interface fmmallocate
c
       interface fmmdeallocate
        module procedure fmmdeallocdbl1d
#if FMM_REAL != FMM_REAL_EXTENDED
        module procedure fmmdeallocdblextended1d
#endif
#if FMM_REAL != FMM_REAL_ITOR && FMM_REAL_EXTENDED != FMM_REAL_ITOR
        module procedure fmmdeallocdblitor1d
#endif
        module procedure fmmdeallocdbl2d,fmmdeallocdbl3d
        module procedure fmmdeallocdbl4d,fmmdeallocdbl5d
        module procedure fmmdeallocint1d,fmmdeallocint2d
        module procedure fmmdeallocint3d,fmmdeallocint4d
        module procedure fmmdealloclog1d,fmmdealloclog2d
        module procedure fmmdeallocbyte,fmmdeallocbyte2d,
     .  fmmdeallocbyte3d
#ifdef FMM_PARALLEL
#if FMM_INTEGER != FMM_MP_INTEGER_PROCESSES
        module procedure fmmdeallocintp
#endif
        module procedure fmmdeallocgpt,fmmdealloctypeint,
     .  fmmdealloctypereal
        module procedure fmmdealloctypeq,fmmdealloctypexyz
#ifdef FMM_LOADSORT
        module procedure fmmdealloctypeint2
#endif
#ifdef FMM_UNIFORMGRID
        module procedure fmmdealloctypeuniformgrid
#endif
        module procedure fmmdealloctypeelem
#endif
       end interface fmmdeallocate
c
       interface fmmallocatept
        module procedure fmmallocdpt1d
        module procedure fmmallocdpt5d
        module procedure fmmallocipt1d
        module procedure fmmalloclpt1d
       end interface fmmallocatept
c
       interface fmmdeallocatept
        module procedure fmmdeallocdpt1d
        module procedure fmmdeallocdpt5d
        module procedure fmmdeallocipt1d
        module procedure fmmdealloclpt1d
       end interface fmmdeallocatept
c
#ifdef FMM_PARALLEL
       interface mp_fmmallocate
        module procedure mp_fmmalloc
       end interface mp_fmmallocate
c
       interface mp_fmmdeallocate
        module procedure mp_fmmdealloc
       end interface mp_fmmdeallocate
#endif
c
#ifdef FMM_ISO_C_BINDING
#ifdef FMM_ALLOCALIGNED
       interface fmmallocate_aligned
        module procedure fmmalloc_aligned1d
        module procedure fmmalloc_aligned2d
       end interface fmmallocate_aligned
c
       interface fmmdeallocate_aligned
        module procedure fmmdealloc_aligned
       end interface fmmdeallocate_aligned
#endif
#endif
c
       interface indallocate
        module procedure indallocint1d
       end interface indallocate
c
       interface inddeallocate
        module procedure inddeallocint1d
       end interface inddeallocate
c
       interface indallocatept
        module procedure indallocipt1d
       end interface indallocatept
c
       interface inddeallocatept
        module procedure inddeallocipt1d
       end interface inddeallocatept
c
       interface dbpallocate
        module procedure dbpallocdbl1d
        module procedure dbpallocdbl3d
       end interface dbpallocate
c
       interface dbpdeallocate
        module procedure dbpdeallocdbl1d
        module procedure dbpdeallocdbl3d
       end interface dbpdeallocate
c
       interface srtallocate
        module procedure srtallocbyte1d
        module procedure srtallocdbl2d
       end interface srtallocate
c
       interface srtdeallocate
        module procedure srtdeallocbyte1d
        module procedure srtdeallocdbl2d
       end interface srtdeallocate
c
       contains
        subroutine fmmallocdbl1d(a,ia,ib,istat)
         use fmmkinds
         implicit none
         real(kind=fmm_real), allocatable:: a(:)
         integer(kind=fmm_integer), intent(in):: ia,ib
         integer(kind=fmm_integer), intent(out):: istat
         integer(kind=fmm_integer) i
         i = nalloc+rtob*(ib-ia+1)
         if(maxallocation.gt.0) then
          if(maxallocation.ge.i) then
           allocate(a(ia:ib),stat = istat)
          else
           istat = 1
          endif
         else
          allocate(a(ia:ib),stat = istat)
         endif
         if(istat.eq.0) then
          nalloc = i
          maxnalloc = max(maxnalloc,nalloc)
          sortmemory = max(sortmemory,maxnalloc)
         else
#ifdef FMM_INFO
          write(6,*) ' allocated memory:',nalloc,'requested memory:',i
#endif
         endif
        end subroutine fmmallocdbl1d
c
#if FMM_REAL != FMM_REAL_EXTENDED
        subroutine fmmallocdblextended1d(a,ia,ib,istat)
         use fmmkinds
         implicit none
         real(kind=fmm_real_extended), allocatable:: a(:)
         integer(kind=fmm_integer), intent(in):: ia,ib
         integer(kind=fmm_integer), intent(out):: istat
         integer(kind=fmm_integer) i
         i = nalloc+rextendedtob*(ib-ia+1)
         if(maxallocation.gt.0) then
          if(maxallocation.ge.i) then
           allocate(a(ia:ib),stat = istat)
          else
           istat = 1
          endif
         else
          allocate(a(ia:ib),stat = istat)
         endif
         if(istat.eq.0) then
          nalloc = i
          maxnalloc = max(maxnalloc,nalloc)
          sortmemory = max(sortmemory,maxnalloc)
         else
#ifdef FMM_INFO
          write(6,*) ' allocated memory:',nalloc,'requested memory:',i
#endif
         endif
        end subroutine fmmallocdblextended1d
#endif
c
#if FMM_REAL != FMM_REAL_ITOR && FMM_REAL_EXTENDED != FMM_REAL_ITOR
        subroutine fmmallocdblitor1d(a,ia,ib,istat)
         use fmmkinds
         implicit none
         real(kind=fmm_real_itor), allocatable:: a(:)
         integer(kind=fmm_integer), intent(in):: ia,ib
         integer(kind=fmm_integer), intent(out):: istat
         integer(kind=fmm_integer) i
         i = nalloc+ritortob*(ib-ia+1)
         if(maxallocation.gt.0) then
          if(maxallocation.ge.i) then
           allocate(a(ia:ib),stat = istat)
          else
           istat = 1
          endif
         else
          allocate(a(ia:ib),stat = istat)
         endif
         if(istat.eq.0) then
          nalloc = i
          maxnalloc = max(maxnalloc,nalloc)
          sortmemory = max(sortmemory,maxnalloc)
         else
#ifdef FMM_INFO
          write(6,*) ' allocated memory:',nalloc,'requested memory:',i
#endif
         endif
        end subroutine fmmallocdblitor1d
#endif
c
        subroutine fmmallocdbl2d(a,ia,ib,ja,jb,istat)
         use fmmkinds
         implicit none
         real(kind=fmm_real), allocatable:: a(:,:)
         integer(kind=fmm_integer), intent(in):: ia,ib,ja,jb
         integer(kind=fmm_integer), intent(out):: istat
         integer(kind=fmm_integer) i
         i = nalloc+rtob*(ib-ia+1)*(jb-ja+1)
         if(maxallocation.gt.0) then
          if(maxallocation.ge.i) then
           allocate(a(ia:ib,ja:jb),stat = istat)
          else
           istat = 1
          endif
         else
          allocate(a(ia:ib,ja:jb),stat = istat)
         endif
         if(istat.eq.0) then
          nalloc = i
          maxnalloc = max(maxnalloc,nalloc)
          sortmemory = max(sortmemory,maxnalloc)
         else
#ifdef FMM_INFO
          write(6,*) ' allocated memory:',nalloc,'requested memory:',i
#endif
         endif
        end subroutine fmmallocdbl2d
c
        subroutine fmmallocdbl3d(a,ia,ib,ja,jb,ka,kb,istat)
         use fmmkinds
         implicit none
         real(kind=fmm_real), allocatable:: a(:,:,:)
         integer(kind=fmm_integer), intent(in):: ia,ib,ja,jb,ka,kb
         integer(kind=fmm_integer), intent(out):: istat
         integer(kind=fmm_integer) i
         i = nalloc+rtob*(ib-ia+1)*(jb-ja+1)*(kb-ka+1)
         if(maxallocation.gt.0) then
          if(maxallocation.ge.i) then
           allocate(a(ia:ib,ja:jb,ka:kb),stat = istat)
          else
           istat = 1
          endif
         else
          allocate(a(ia:ib,ja:jb,ka:kb),stat = istat)
         endif
         if(istat.eq.0) then
          nalloc = i
          maxnalloc = max(maxnalloc,nalloc)
          sortmemory = max(sortmemory,maxnalloc)
         else
#ifdef FMM_INFO
          write(6,*) ' allocated memory:',nalloc,'requested memory:',i
#endif
         endif
        end subroutine fmmallocdbl3d
c
        subroutine fmmallocdbl4d(a,ia,ib,ja,jb,ka,kb,la,lb,istat)
         use fmmkinds
         implicit none
         real(kind=fmm_real), allocatable:: a(:,:,:,:)
         integer(kind=fmm_integer), intent(in):: ia,ib,ja,jb,ka,kb,la,lb
         integer(kind=fmm_integer), intent(out):: istat
         integer(kind=fmm_integer) i
         i = nalloc+rtob*(ib-ia+1)*(jb-ja+1)*(kb-ka+1)*(lb-la+1)
         if(maxallocation.gt.0) then
          if(maxallocation.ge.i) then
           allocate(a(ia:ib,ja:jb,ka:kb,la:lb),stat = istat)
          else
           istat = 1
          endif
         else
          allocate(a(ia:ib,ja:jb,ka:kb,la:lb),stat = istat)
         endif
         if(istat.eq.0) then
          nalloc = i
          maxnalloc = max(maxnalloc,nalloc)
          sortmemory = max(sortmemory,maxnalloc)
         else
#ifdef FMM_INFO
          write(6,*) ' allocated memory:',nalloc,'requested memory:',i
#endif
         endif
        end subroutine fmmallocdbl4d
c
        subroutine fmmallocdbl5d(a,ia,ib,ja,jb,ka,kb,la,lb,ma,mb,istat)
         use fmmkinds
         implicit none
         real(kind=fmm_real), allocatable:: a(:,:,:,:,:)
         integer(kind=fmm_integer), intent(in):: ia,ib,ja,jb,ka,kb,la,
     .   lb,ma,mb
         integer(kind=fmm_integer), intent(out):: istat
         integer(kind=fmm_integer) i
         i=nalloc+rtob*(ib-ia+1)*(jb-ja+1)*(kb-ka+1)*(lb-la+1)*(mb-ma+1)
         if(maxallocation.gt.0) then
          if(maxallocation.ge.i) then
           allocate(a(ia:ib,ja:jb,ka:kb,la:lb,ma:mb),stat = istat)
          else
           istat = 1
          endif
         else
          allocate(a(ia:ib,ja:jb,ka:kb,la:lb,ma:mb),stat = istat)
         endif
         if(istat.eq.0) then
          nalloc = i
          maxnalloc = max(maxnalloc,nalloc)
          sortmemory = max(sortmemory,maxnalloc)
         else
#ifdef FMM_INFO
          write(6,*) ' allocated memory:',nalloc,'requested memory:',i
#endif
         endif
        end subroutine fmmallocdbl5d
c
        subroutine fmmallocint1d(a,ia,ib,istat)
         use fmmkinds
         implicit none
         integer(kind=fmm_integer), allocatable:: a(:)
         integer(kind=fmm_integer), intent(in):: ia,ib
         integer(kind=fmm_integer), intent(out):: istat
         integer(kind=fmm_integer) i
         i = nalloc+itob*(ib-ia+1)
         if(maxallocation.gt.0) then
          if(maxallocation.ge.i) then
           allocate(a(ia:ib),stat = istat)
          else
           istat = 1
          endif
         else
          allocate(a(ia:ib),stat = istat)
         endif
         if(istat.eq.0) then
          nalloc = i
          maxnalloc = max(maxnalloc,nalloc)
          sortmemory = max(sortmemory,maxnalloc)
         else
#ifdef FMM_INFO
          write(6,*) ' allocated memory:',nalloc,'requested memory:',i
#endif
         endif
        end subroutine fmmallocint1d
c
#ifdef FMM_PARALLEL
#if FMM_INTEGER != FMM_MP_INTEGER_PROCESSES
        subroutine fmmallocintp(a,ia,ib,istat)
         use fmmkinds
         implicit none
         integer(kind=mp_integer_processes), allocatable:: a(:)
         integer(kind=fmm_integer), intent(in):: ia,ib
         integer(kind=fmm_integer), intent(out):: istat
         integer(kind=fmm_integer) i
         i = nalloc+iptob*(ib-ia+1)
         if(maxallocation.gt.0) then
          if(maxallocation.ge.i) then
           allocate(a(ia:ib),stat = istat)
          else
           istat = 1
          endif
         else
          allocate(a(ia:ib),stat = istat)
         endif
         if(istat.eq.0) then
          nalloc = i
          maxnalloc = max(maxnalloc,nalloc)
          sortmemory = max(sortmemory,maxnalloc)
         else
#ifdef FMM_INFO
          write(6,*) ' allocated memory:',nalloc,'requested memory:',i
#endif
         endif
        end subroutine fmmallocintp
#endif
#endif
c
        subroutine fmmallocint2d(a,ia,ib,ja,jb,istat)
         use fmmkinds
         implicit none
         integer(kind=fmm_integer), allocatable:: a(:,:)
         integer(kind=fmm_integer), intent(in):: ia,ib,ja,jb
         integer(kind=fmm_integer), intent(out):: istat
         integer(kind=fmm_integer) i
         i = nalloc+itob*(ib-ia+1)*(jb-ja+1)
         if(maxallocation.gt.0) then
          if(maxallocation.ge.i) then
           allocate(a(ia:ib,ja:jb),stat = istat)
          else
           istat = 1
          endif
         else
          allocate(a(ia:ib,ja:jb),stat = istat)
         endif
         if(istat.eq.0) then
          nalloc = i
          maxnalloc = max(maxnalloc,nalloc)
          sortmemory = max(sortmemory,maxnalloc)
         else
#ifdef FMM_INFO
          write(6,*) ' allocated memory:',nalloc,'requested memory:',i
#endif
         endif
        end subroutine fmmallocint2d
c
        subroutine fmmallocint3d(a,ia,ib,ja,jb,ka,kb,istat)
         use fmmkinds
         implicit none
         integer(kind=fmm_integer), allocatable:: a(:,:,:)
         integer(kind=fmm_integer), intent(in):: ia,ib,ja,jb,ka,kb
         integer(kind=fmm_integer), intent(out):: istat
         integer(kind=fmm_integer) i
         i = nalloc+itob*(ib-ia+1)*(jb-ja+1)*(kb-ka+1)
         if(maxallocation.gt.0) then
          if(maxallocation.ge.i) then
           allocate(a(ia:ib,ja:jb,ka:kb),stat = istat)
          else
           istat = 1
          endif
         else
          allocate(a(ia:ib,ja:jb,ka:kb),stat = istat)
         endif
         if(istat.eq.0) then
          nalloc = i
          maxnalloc = max(maxnalloc,nalloc)
          sortmemory = max(sortmemory,maxnalloc)
         else
#ifdef FMM_INFO
          write(6,*) ' allocated memory:',nalloc,'requested memory:',i
#endif
         endif
        end subroutine fmmallocint3d
c
        subroutine fmmallocint4d(a,ia,ib,ja,jb,ka,kb,la,lb,istat)
         use fmmkinds
         implicit none
         integer(kind=fmm_integer), allocatable:: a(:,:,:,:)
         integer(kind=fmm_integer), intent(in):: ia,ib,ja,jb,ka,kb,la,lb
         integer(kind=fmm_integer), intent(out):: istat
         integer(kind=fmm_integer) i
         i = nalloc+itob*(ib-ia+1)*(jb-ja+1)*(kb-ka+1)*(lb-la+1)
         if(maxallocation.gt.0) then
          if(maxallocation.ge.i) then
           allocate(a(ia:ib,ja:jb,ka:kb,la:lb),stat = istat)
          else
           istat = 1
          endif
         else
          allocate(a(ia:ib,ja:jb,ka:kb,la:lb),stat = istat)
         endif
         if(istat.eq.0) then
          nalloc = i
          maxnalloc = max(maxnalloc,nalloc)
          sortmemory = max(sortmemory,maxnalloc)
         else
#ifdef FMM_INFO
          write(6,*) ' allocated memory:',nalloc,'requested memory:',i
#endif
         endif
        end subroutine fmmallocint4d
c
        subroutine fmmalloclog1d(a,ia,ib,istat)
         use fmmkinds
         implicit none
         logical(kind=fmm_logical), allocatable:: a(:)
         integer(kind=fmm_integer), intent(in):: ia,ib
         integer(kind=fmm_integer), intent(out):: istat
         integer(kind=fmm_integer) i
         i = nalloc+ltob*(ib-ia+1)
         if(maxallocation.gt.0) then
          if(maxallocation.ge.i) then
           allocate(a(ia:ib),stat = istat)
          else
           istat = 1
          endif
         else
          allocate(a(ia:ib),stat = istat)
         endif
         if(istat.eq.0) then
          nalloc = i
          maxnalloc = max(maxnalloc,nalloc)
          sortmemory = max(sortmemory,maxnalloc)
         else
#ifdef FMM_INFO
          write(6,*) ' allocated memory:',nalloc,'requested memory:',i
#endif
         endif
        end subroutine fmmalloclog1d
c
        subroutine fmmalloclog2d(a,ia,ib,ja,jb,istat)
         use fmmkinds
         implicit none
         logical(kind=fmm_logical), allocatable:: a(:,:)
         integer(kind=fmm_integer), intent(in):: ia,ib,ja,jb
         integer(kind=fmm_integer), intent(out):: istat
         integer(kind=fmm_integer) i
         i = nalloc+ltob*(ib-ia+1)*(jb-ja+1)
         if(maxallocation.gt.0) then
          if(maxallocation.ge.i) then
           allocate(a(ia:ib,ja:jb),stat = istat)
          else
           istat = 1
          endif
         else
          allocate(a(ia:ib,ja:jb),stat = istat)
         endif
         if(istat.eq.0) then
          nalloc = i
          maxnalloc = max(maxnalloc,nalloc)
          sortmemory = max(sortmemory,maxnalloc)
         else
#ifdef FMM_INFO
          write(6,*) ' allocated memory:',nalloc,'requested memory:',i
#endif
         endif
        end subroutine fmmalloclog2d
c
        subroutine fmmallocbyte(a,ia,ib,istat)
         use fmmkinds
         implicit none
         integer(kind=1), allocatable:: a(:)
         integer(kind=fmm_integer), intent(in):: ia,ib
         integer(kind=fmm_integer), intent(out):: istat
         integer(kind=fmm_integer) i
         i = nalloc+(ib-ia+1)
         if(maxallocation.gt.0) then
          if(maxallocation.ge.i) then
           allocate(a(ia:ib),stat = istat)
          else
           istat = 1
          endif
         else
          allocate(a(ia:ib),stat = istat)
         endif
         if(istat.eq.0) then
          nalloc = i
          maxnalloc = max(maxnalloc,nalloc)
          sortmemory = max(sortmemory,maxnalloc)
         else
#ifdef FMM_INFO
          write(6,*) ' allocated memory:',nalloc,'requested memory:',i
#endif
         endif
        end subroutine fmmallocbyte
c
        subroutine fmmallocbyte2d(a,ia,ib,ja,jb,istat)
         use fmmkinds
         implicit none
         integer(kind=1), allocatable:: a(:,:)
         integer(kind=fmm_integer), intent(in):: ia,ib,ja,jb
         integer(kind=fmm_integer), intent(out):: istat
         integer(kind=fmm_integer) i
         i = nalloc+(ib-ia+1)*(jb-ja+1)
         if(maxallocation.gt.0) then
          if(maxallocation.ge.i) then
           allocate(a(ia:ib,ja:jb),stat = istat)
          else
           istat = 1
          endif
         else
          allocate(a(ia:ib,ja:jb),stat = istat)
         endif
         if(istat.eq.0) then
          nalloc = i
          maxnalloc = max(maxnalloc,nalloc)
          sortmemory = max(sortmemory,maxnalloc)
         else
#ifdef FMM_INFO
          write(6,*) ' allocated memory:',nalloc,'requested memory:',i
#endif
         endif
        end subroutine fmmallocbyte2d
c
        subroutine fmmallocbyte3d(a,ia,ib,ja,jb,ka,kb,istat)
         use fmmkinds
         implicit none
         integer(kind=1), allocatable:: a(:,:,:)
         integer(kind=fmm_integer), intent(in):: ia,ib,ja,jb,ka,kb
         integer(kind=fmm_integer), intent(out):: istat
         integer(kind=fmm_integer) i
         i = nalloc+(ib-ia+1)*(jb-ja+1)*(kb-ka+1)
         if(maxallocation.gt.0) then
          if(maxallocation.ge.i) then
           allocate(a(ia:ib,ja:jb,ka:kb),stat = istat)
          else
           istat = 1
          endif
         else
          allocate(a(ia:ib,ja:jb,ka:kb),stat = istat)
         endif
         if(istat.eq.0) then
          nalloc = i
          maxnalloc = max(maxnalloc,nalloc)
          sortmemory = max(sortmemory,maxnalloc)
         else
#ifdef FMM_INFO
          write(6,*) ' allocated memory:',nalloc,'requested memory:',i
#endif
         endif
        end subroutine fmmallocbyte3d
c
#ifdef FMM_PARALLEL
        subroutine fmmallocgpt(gpt,ia,ib,istat)
         use fmmkinds
         implicit none
         type(c_ptr), allocatable, target:: gpt(:)
         integer(kind=fmm_integer), intent(in):: ia,ib
         integer(kind=fmm_integer), intent(out):: istat
         integer(kind=fmm_integer) i
         i = nalloc+gbpttob*(ib-ia+1)
         if(maxallocation.gt.0) then
          if(maxallocation.ge.i) then
           allocate(gpt(ia:ib),stat = istat)
          else
           istat = 1
          endif
         else
          allocate(gpt(ia:ib),stat = istat)
         endif
         if(istat.eq.0) then
          nalloc = i
          maxnalloc = max(maxnalloc,nalloc)
          sortmemory = max(sortmemory,maxnalloc)
         else
#ifdef FMM_INFO
          write(6,*) ' allocated memory:',nalloc,'requested memory:',i
#endif
         endif
        end subroutine fmmallocgpt
c
        subroutine fmmalloctypeint(typeint,ia,ib,istat)
         use fmmkinds
         implicit none
         type(tsndibox), allocatable:: typeint(:)
         integer(kind=fmm_integer), intent(in):: ia,ib
         integer(kind=fmm_integer), intent(out):: istat
         integer(kind=fmm_integer) i
         i = nalloc+typeinttob*(ib-ia+1)
         if(maxallocation.gt.0) then
          if(maxallocation.ge.i) then
           allocate(typeint(ia:ib),stat = istat)
          else
           istat = 1
          endif
         else
          allocate(typeint(ia:ib),stat = istat)
         endif
         if(istat.eq.0) then
          nalloc = i
          maxnalloc = max(maxnalloc,nalloc)
          sortmemory = max(sortmemory,maxnalloc)
         else
#ifdef FMM_INFO
          write(6,*) ' allocated memory:',nalloc,'requested memory:',i
#endif
         endif
        end subroutine fmmalloctypeint
c
#ifdef FMM_LOADSORT
        subroutine fmmalloctypeint2(typeint2,ia,ib,istat)
         use fmmkinds
         implicit none
         type(tsndiboxload), allocatable:: typeint2(:)
         integer(kind=fmm_integer), intent(in):: ia,ib
         integer(kind=fmm_integer), intent(out):: istat
         integer(kind=fmm_integer) i
         i = nalloc+typeint2tob*(ib-ia+1)
         if(maxallocation.gt.0) then
          if(maxallocation.ge.i) then
           allocate(typeint2(ia:ib),stat = istat)
          else
           istat = 1
          endif
         else
          allocate(typeint2(ia:ib),stat = istat)
         endif
         if(istat.eq.0) then
          nalloc = i
          maxnalloc = max(maxnalloc,nalloc)
          sortmemory = max(sortmemory,maxnalloc)
         else
#ifdef FMM_INFO
          write(6,*) ' allocated memory:',nalloc,'requested memory:',i
#endif
         endif
        end subroutine fmmalloctypeint2
#endif
c
#ifdef FMM_UNIFORMGRID
        subroutine fmmalloctypeuniformgrid(typeuniformgrid,ia,ib,istat)
         use fmmkinds
         implicit none
         type(tsnduniformgrid), allocatable:: typeuniformgrid(:)
         integer(kind=fmm_integer), intent(in):: ia,ib
         integer(kind=fmm_integer), intent(out):: istat
         integer(kind=fmm_integer) i
         i = nalloc+typeuniformgridtob*(ib-ia+1)
         if(maxallocation.gt.0) then
          if(maxallocation.ge.i) then
           allocate(typeuniformgrid(ia:ib),stat = istat)
          else
           istat = 1
          endif
         else
          allocate(typeuniformgrid(ia:ib),stat = istat)
         endif
         if(istat.eq.0) then
          nalloc = i
          maxnalloc = max(maxnalloc,nalloc)
          sortmemory = max(sortmemory,maxnalloc)
         else
#ifdef FMM_INFO
          write(6,*) ' allocated memory:',nalloc,'requested memory:',i
#endif
         endif
        end subroutine fmmalloctypeuniformgrid
#endif
c
        subroutine fmmalloctypeelem(typeelem,istat)
         use fmmkinds
         implicit none
         type(telem), pointer:: typeelem
         integer(kind=fmm_integer), intent(out):: istat
         integer(kind=fmm_integer) i
         i = nalloc+typeelemtob
         if(maxallocation.gt.0) then
          if(maxallocation.ge.i) then
           allocate(typeelem,stat = istat)
          else
           istat = 1
          endif
         else
          allocate(typeelem,stat = istat)
         endif
         if(istat.eq.0) then
          nalloc = i
          maxnalloc = max(maxnalloc,nalloc)
          sortmemory = max(sortmemory,maxnalloc)
         else
#ifdef FMM_INFO
          write(6,*) ' allocated memory:',nalloc,'requested memory:',i
#endif
         endif
        end subroutine fmmalloctypeelem
c
        subroutine fmmalloctypereal(typereal,ia,ib,istat)
         use fmmkinds
         implicit none
         type(tsndomegatree), allocatable:: typereal(:)
         integer(kind=fmm_integer), intent(in):: ia,ib
         integer(kind=fmm_integer), intent(out):: istat
         integer(kind=fmm_integer) i
         i = nalloc+typerealtob*(ib-ia+1)
         if(maxallocation.gt.0) then
          if(maxallocation.ge.i) then
           allocate(typereal(ia:ib),stat = istat)
          else
           istat = 1
          endif
         else
          allocate(typereal(ia:ib),stat = istat)
         endif
         if(istat.eq.0) then
          nalloc = i
          maxnalloc = max(maxnalloc,nalloc)
          sortmemory = max(sortmemory,maxnalloc)
         else
#ifdef FMM_INFO
          write(6,*) ' allocated memory:',nalloc,'requested memory:',i
#endif
         endif
        end subroutine fmmalloctypereal
c
        subroutine fmmalloctypeq(typeq,ia,ib,istat)
         use fmmkinds
         implicit none
         type(tsndq), allocatable:: typeq(:)
         integer(kind=fmm_integer), intent(in):: ia,ib
         integer(kind=fmm_integer), intent(out):: istat
         integer(kind=fmm_integer) i
         i = nalloc+typeqtob*(ib-ia+1)
         if(maxallocation.gt.0) then
          if(maxallocation.ge.i) then
           allocate(typeq(ia:ib),stat = istat)
          else
           istat = 1
          endif
         else
          allocate(typeq(ia:ib),stat = istat)
         endif
         if(istat.eq.0) then
          nalloc = i
          maxnalloc = max(maxnalloc,nalloc)
          sortmemory = max(sortmemory,maxnalloc)
         else
#ifdef FMM_INFO
          write(6,*) ' allocated memory:',nalloc,'requested memory:',i
#endif
         endif
        end subroutine fmmalloctypeq
c
        subroutine fmmalloctypexyz(typexyz,ia,ib,istat)
         use fmmkinds
         implicit none
         type(tsndxyz), allocatable:: typexyz(:)
         integer(kind=fmm_integer), intent(in):: ia,ib
         integer(kind=fmm_integer), intent(out):: istat
         integer(kind=fmm_integer) i
         i = nalloc+typexyztob*(ib-ia+1)
         if(maxallocation.gt.0) then
          if(maxallocation.ge.i) then
           allocate(typexyz(ia:ib),stat = istat)
          else
           istat = 1
          endif
         else
          allocate(typexyz(ia:ib),stat = istat)
         endif
         if(istat.eq.0) then
          nalloc = i
          maxnalloc = max(maxnalloc,nalloc)
          sortmemory = max(sortmemory,maxnalloc)
         else
#ifdef FMM_INFO
          write(6,*) ' allocated memory:',nalloc,'requested memory:',i
#endif
         endif
        end subroutine fmmalloctypexyz
#endif
c
        subroutine fmmdeallocdbl1d(a,istat)
         use fmmkinds
         implicit none
         real(kind=fmm_real), allocatable:: a(:)
         integer(kind=fmm_integer), intent(out):: istat
         integer(kind=fmm_integer) i
         i = size(a)
         deallocate(a,stat = istat)
         if(istat.eq.0) nalloc = nalloc-rtob*i
        end subroutine fmmdeallocdbl1d
c
#if FMM_REAL != FMM_REAL_EXTENDED
        subroutine fmmdeallocdblextended1d(a,istat)
         use fmmkinds
         implicit none
         real(kind=fmm_real_extended), allocatable:: a(:)
         integer(kind=fmm_integer), intent(out):: istat
         integer(kind=fmm_integer) i
         i = size(a)
         deallocate(a,stat = istat)
         if(istat.eq.0) nalloc = nalloc-rextendedtob*i
        end subroutine fmmdeallocdblextended1d
#endif
c
#if FMM_REAL != FMM_REAL_ITOR && FMM_REAL_EXTENDED != FMM_REAL_ITOR
        subroutine fmmdeallocdblitor1d(a,istat)
         use fmmkinds
         implicit none
         real(kind=fmm_real_itor), allocatable:: a(:)
         integer(kind=fmm_integer), intent(out):: istat
         integer(kind=fmm_integer) i
         i = size(a)
         deallocate(a,stat = istat)
         if(istat.eq.0) nalloc = nalloc-ritortob*i
        end subroutine fmmdeallocdblitor1d
#endif
c
        subroutine fmmdeallocdbl2d(a,istat)
         use fmmkinds
         implicit none
         real(kind=fmm_real), allocatable:: a(:,:)
         integer(kind=fmm_integer), intent(out):: istat
         integer(kind=fmm_integer) i
         i = size(a)
         deallocate(a,stat = istat)
         if(istat.eq.0) nalloc = nalloc-rtob*i
        end subroutine fmmdeallocdbl2d
c
        subroutine fmmdeallocdbl3d(a,istat)
         use fmmkinds
         implicit none
         real(kind=fmm_real), allocatable:: a(:,:,:)
         integer(kind=fmm_integer), intent(out):: istat
         integer(kind=fmm_integer) i
         i = size(a)
         deallocate(a,stat = istat)
         if(istat.eq.0) nalloc = nalloc-rtob*i
        end subroutine fmmdeallocdbl3d
c
        subroutine fmmdeallocdbl4d(a,istat)
         use fmmkinds
         implicit none
         real(kind=fmm_real), allocatable:: a(:,:,:,:)
         integer(kind=fmm_integer), intent(out):: istat
         integer(kind=fmm_integer) i
         i = size(a)
         deallocate(a,stat = istat)
         if(istat.eq.0) nalloc = nalloc-rtob*i
        end subroutine fmmdeallocdbl4d
c
        subroutine fmmdeallocdbl5d(a,istat)
         use fmmkinds
         implicit none
         real(kind=fmm_real), allocatable:: a(:,:,:,:,:)
         integer(kind=fmm_integer), intent(out):: istat
         integer(kind=fmm_integer) i
         i = size(a)
         deallocate(a,stat = istat)
         if(istat.eq.0) nalloc = nalloc-rtob*i
        end subroutine fmmdeallocdbl5d
c
        subroutine fmmdeallocint1d(a,istat)
         use fmmkinds
         implicit none
         integer(kind=fmm_integer), allocatable:: a(:)
         integer(kind=fmm_integer), intent(out):: istat
         integer(kind=fmm_integer) i
         i = size(a)
         deallocate(a,stat = istat)
         if(istat.eq.0) nalloc = nalloc-itob*i
        end subroutine fmmdeallocint1d
c
#ifdef FMM_PARALLEL
#if FMM_INTEGER != FMM_MP_INTEGER_PROCESSES
        subroutine fmmdeallocintp(a,istat)
         use fmmkinds
         implicit none
         integer(kind=mp_integer_processes), allocatable:: a(:)
         integer(kind=fmm_integer), intent(out):: istat
         integer(kind=fmm_integer) i
         i = size(a)
         deallocate(a,stat = istat)
         if(istat.eq.0) nalloc = nalloc-iptob*i
        end subroutine fmmdeallocintp
#endif
#endif
c
        subroutine fmmdeallocint2d(a,istat)
         use fmmkinds
         implicit none
         integer(kind=fmm_integer), allocatable:: a(:,:)
         integer(kind=fmm_integer), intent(out):: istat
         integer(kind=fmm_integer) i
         i = size(a)
         deallocate(a,stat = istat)
         if(istat.eq.0) nalloc = nalloc-itob*i
        end subroutine fmmdeallocint2d
c
        subroutine fmmdeallocint3d(a,istat)
         use fmmkinds
         implicit none
         integer(kind=fmm_integer), allocatable:: a(:,:,:)
         integer(kind=fmm_integer), intent(out):: istat
         integer(kind=fmm_integer) i
         i = size(a)
         deallocate(a,stat = istat)
         if(istat.eq.0) nalloc = nalloc-itob*i
        end subroutine fmmdeallocint3d
c
        subroutine fmmdeallocint4d(a,istat)
         use fmmkinds
         implicit none
         integer(kind=fmm_integer), allocatable:: a(:,:,:,:)
         integer(kind=fmm_integer), intent(out):: istat
         integer(kind=fmm_integer) i
         i = size(a)
         deallocate(a,stat = istat)
         if(istat.eq.0) nalloc = nalloc-itob*i
        end subroutine fmmdeallocint4d
c
        subroutine fmmdealloclog1d(a,istat)
         use fmmkinds
         implicit none
         logical(kind=fmm_logical), allocatable:: a(:)
         integer(kind=fmm_integer), intent(out):: istat
         integer(kind=fmm_integer) i
         i = size(a)
         deallocate(a,stat = istat)
         if(istat.eq.0) nalloc = nalloc-ltob*i
        end subroutine fmmdealloclog1d
c
        subroutine fmmdealloclog2d(a,istat)
         use fmmkinds
         implicit none
         logical(kind=fmm_logical), allocatable:: a(:,:)
         integer(kind=fmm_integer), intent(out):: istat
         integer(kind=fmm_integer) i
         i = size(a)
         deallocate(a,stat = istat)
         if(istat.eq.0) nalloc = nalloc-ltob*i
        end subroutine fmmdealloclog2d
c
        subroutine fmmdeallocbyte(a,istat)
         use fmmkinds
         implicit none
         integer(kind=1), allocatable:: a(:)
         integer(kind=fmm_integer), intent(out):: istat
         integer(kind=fmm_integer) i
         i = size(a)
         deallocate(a,stat = istat)
         if(istat.eq.0) nalloc = nalloc-i
        end subroutine fmmdeallocbyte
c
        subroutine fmmdeallocbyte2d(a,istat)
         use fmmkinds
         implicit none
         integer(kind=1), allocatable:: a(:,:)
         integer(kind=fmm_integer), intent(out):: istat
         integer(kind=fmm_integer) i
         i = size(a)
         deallocate(a,stat = istat)
         if(istat.eq.0) nalloc = nalloc-i
        end subroutine fmmdeallocbyte2d
c
        subroutine fmmdeallocbyte3d(a,istat)
         use fmmkinds
         implicit none
         integer(kind=1), allocatable:: a(:,:,:)
         integer(kind=fmm_integer), intent(out):: istat
         integer(kind=fmm_integer) i
         i = size(a)
         deallocate(a,stat = istat)
         if(istat.eq.0) nalloc = nalloc-i
        end subroutine fmmdeallocbyte3d
c
#ifdef FMM_PARALLEL
        subroutine fmmdeallocgpt(gpt,istat)
         use fmmkinds
         implicit none
         type(c_ptr), allocatable, target:: gpt(:)
         integer(kind=fmm_integer), intent(out):: istat
         integer(kind=fmm_integer) i
         i = size(gpt,1)
         deallocate(gpt,stat = istat)
         if(istat.eq.0) nalloc = nalloc-gbpttob*i
        end subroutine fmmdeallocgpt
c
        subroutine fmmdealloctypeint(typeint,istat)
         use fmmkinds
         implicit none
         type(tsndibox), allocatable:: typeint(:)
         integer(kind=fmm_integer), intent(out):: istat
         integer(kind=fmm_integer) i
         i = size(typeint,1)
         deallocate(typeint,stat = istat)
         if(istat.eq.0) nalloc = nalloc-typeinttob*i
        end subroutine fmmdealloctypeint
c
#ifdef FMM_LOADSORT
        subroutine fmmdealloctypeint2(typeint2,istat)
         use fmmkinds
         implicit none
         type(tsndiboxload), allocatable:: typeint2(:)
         integer(kind=fmm_integer), intent(out):: istat
         integer(kind=fmm_integer) i
         i = size(typeint2,1)
         deallocate(typeint2,stat = istat)
         if(istat.eq.0) nalloc = nalloc-typeint2tob*i
        end subroutine fmmdealloctypeint2
#endif
c
#ifdef FMM_UNIFORMGRID
        subroutine fmmdealloctypeuniformgrid(typeuniformgrid,istat)
         use fmmkinds
         implicit none
         type(tsnduniformgrid), allocatable:: typeuniformgrid(:)
         integer(kind=fmm_integer), intent(out):: istat
         integer(kind=fmm_integer) i
         i = size(typeuniformgrid)
         deallocate(typeuniformgrid,stat = istat)
         if(istat.eq.0) nalloc = nalloc-typeuniformgridtob*i
        end subroutine fmmdealloctypeuniformgrid
#endif
c
        subroutine fmmdealloctypeelem(typeelem,istat)
         use fmmkinds
         implicit none
         type(telem), pointer:: typeelem
         integer(kind=fmm_integer), intent(out):: istat
         deallocate(typeelem,stat = istat)
         if(istat.eq.0) nalloc = nalloc-typeelemtob
        end subroutine fmmdealloctypeelem
c
        subroutine fmmdealloctypereal(typereal,istat)
         use fmmkinds
         implicit none
         type(tsndomegatree), allocatable:: typereal(:)
         integer(kind=fmm_integer), intent(out):: istat
         integer(kind=fmm_integer) i
         i = size(typereal,1)
         deallocate(typereal,stat = istat)
         if(istat.eq.0) nalloc = nalloc-typerealtob*i
        end subroutine fmmdealloctypereal
c
        subroutine fmmdealloctypeq(typeq,istat)
         use fmmkinds
         implicit none
         type(tsndq), allocatable:: typeq(:)
         integer(kind=fmm_integer), intent(out):: istat
         integer(kind=fmm_integer) i
         i = size(typeq,1)
         deallocate(typeq,stat = istat)
         if(istat.eq.0) nalloc = nalloc-typeqtob*i
        end subroutine fmmdealloctypeq
c
        subroutine fmmdealloctypexyz(typexyz,istat)
         use fmmkinds
         implicit none
         type(tsndxyz), allocatable:: typexyz(:)
         integer(kind=fmm_integer), intent(out):: istat
         integer(kind=fmm_integer) i
         i = size(typexyz,1)
         deallocate(typexyz,stat = istat)
         if(istat.eq.0) nalloc = nalloc-typexyztob*i
        end subroutine fmmdealloctypexyz
#endif
c
        subroutine fmmallocdpt1d(a,ia,ib,istat)
         use fmmkinds
         implicit none
         real(kind=fmm_real), pointer:: a(:)
         integer(kind=fmm_integer), intent(in):: ia,ib
         integer(kind=fmm_integer), intent(out):: istat
         integer(kind=fmm_integer) i
         i = nalloc+rtob*(ib-ia+1)
         if(maxallocation.gt.0) then
          if(maxallocation.ge.i) then
           allocate(a(ia:ib),stat = istat)
          else
           istat = 1
          endif
         else
          allocate(a(ia:ib),stat = istat)
         endif
         if(istat.eq.0) then
          nalloc = i
          maxnalloc = max(maxnalloc,nalloc)
          sortmemory = max(sortmemory,maxnalloc)
         else
#ifdef FMM_INFO
          write(6,*) ' allocated memory:',nalloc,'requested memory:',i
#endif
         endif
        end subroutine fmmallocdpt1d
c
        subroutine fmmallocdpt5d(a,ia,ib,ja,jb,ka,kb,la,lb,ma,mb,istat)
         use fmmkinds
         implicit none
         real(kind=fmm_real), pointer:: a(:,:,:,:,:)
         integer(kind=fmm_integer), intent(in):: ia,ib,ja,jb,ka,kb,la,
     .   lb,ma,mb
         integer(kind=fmm_integer), intent(out):: istat
         integer(kind=fmm_integer) i
         i=nalloc+rtob*(ib-ia+1)*(jb-ja+1)*(kb-ka+1)*(lb-la+1)*(mb-ma+1)
         if(maxallocation.gt.0) then
          if(maxallocation.ge.i) then
           allocate(a(ia:ib,ja:jb,ka:kb,la:lb,ma:mb),stat = istat)
          else
           istat = 1
          endif
         else
          allocate(a(ia:ib,ja:jb,ka:kb,la:lb,ma:mb),stat = istat)
         endif
         if(istat.eq.0) then
          nalloc = i
          maxnalloc = max(maxnalloc,nalloc)
          sortmemory = max(sortmemory,maxnalloc)
         else
#ifdef FMM_INFO
          write(6,*) ' allocated memory:',nalloc,'requested memory:',i
#endif
         endif
        end subroutine fmmallocdpt5d
c
        subroutine fmmallocipt1d(a,ia,ib,istat)
         use fmmkinds
         implicit none
         integer(kind=fmm_integer), pointer:: a(:)
         integer(kind=fmm_integer), intent(in):: ia,ib
         integer(kind=fmm_integer), intent(out):: istat
         integer(kind=fmm_integer) i
         i = nalloc+itob*(ib-ia+1)
         if(maxallocation.gt.0) then
          if(maxallocation.ge.i) then
           allocate(a(ia:ib),stat = istat)
          else
           istat = 1
          endif
         else
          allocate(a(ia:ib),stat = istat)
         endif
         if(istat.eq.0) then
          nalloc = i
          maxnalloc = max(maxnalloc,nalloc)
          sortmemory = max(sortmemory,maxnalloc)
         else
#ifdef FMM_INFO
          write(6,*) ' allocated memory:',nalloc,'requested memory:',i
#endif
         endif
        end subroutine fmmallocipt1d
c
        subroutine fmmalloclpt1d(a,ia,ib,istat)
         use fmmkinds
         implicit none
         logical(kind=fmm_logical), pointer:: a(:)
         integer(kind=fmm_integer), intent(in):: ia,ib
         integer(kind=fmm_integer), intent(out):: istat
         integer(kind=fmm_integer) i
         i = nalloc+ltob*(ib-ia+1)
         if(maxallocation.gt.0) then
          if(maxallocation.ge.i) then
           allocate(a(ia:ib),stat = istat)
          else
           istat = 1
          endif
         else
          allocate(a(ia:ib),stat = istat)
         endif
         if(istat.eq.0) then
          nalloc = i
          maxnalloc = max(maxnalloc,nalloc)
          sortmemory = max(sortmemory,maxnalloc)
         else
#ifdef FMM_INFO
          write(6,*) ' allocated memory:',nalloc,'requested memory:',i
#endif
         endif
        end subroutine fmmalloclpt1d
c
        subroutine fmmdeallocdpt1d(a,istat)
         use fmmkinds
         implicit none
         real(kind=fmm_real), pointer:: a(:)
         integer(kind=fmm_integer), intent(out):: istat
         integer(kind=fmm_integer) i
         i = size(a)
         deallocate(a,stat = istat)
         a => null()
         if(istat.eq.0) nalloc = nalloc-rtob*i
        end subroutine fmmdeallocdpt1d
c
        subroutine fmmdeallocdpt5d(a,istat)
         use fmmkinds
         implicit none
         real(kind=fmm_real), pointer:: a(:,:,:,:,:)
         integer(kind=fmm_integer), intent(out):: istat
         integer(kind=fmm_integer) i
         i = size(a)
         deallocate(a,stat = istat)
         a => null()
         if(istat.eq.0) nalloc = nalloc-rtob*i
        end subroutine fmmdeallocdpt5d
c
        subroutine fmmdeallocipt1d(a,istat)
         use fmmkinds
         implicit none
         integer(kind=fmm_integer), pointer:: a(:)
         integer(kind=fmm_integer), intent(out):: istat
         integer(kind=fmm_integer) i
         i = size(a)
         deallocate(a,stat = istat)
         a => null()
         if(istat.eq.0) nalloc = nalloc-itob*i
        end subroutine fmmdeallocipt1d
c
        subroutine fmmdealloclpt1d(a,istat)
         use fmmkinds
         implicit none
         logical(kind=fmm_logical), pointer:: a(:)
         integer(kind=fmm_integer), intent(out):: istat
         integer(kind=fmm_integer) i
         i = size(a)
         deallocate(a,stat = istat)
         a => null()
         if(istat.eq.0) nalloc = nalloc-ltob*i
        end subroutine fmmdealloclpt1d
c
#ifdef FMM_PARALLEL
        subroutine mp_fmmalloc(gpt,n)
         use fmmkinds
         use mp_info, only: nnodes,mp_allocate,mp_error
         implicit none
         type(c_ptr), allocatable, target:: gpt(:)
         integer(kind=fmm_integer) n,i,j,k,istat
         if(n.ge.0) then
          i = nalloc+n
          j = nmp_alloc+n
          if(maxallocation.gt.0) then
           if(maxallocation.ge.i) then
            if(n.gt.0) then
             if(nnodes.gt.1) then
              k = n
             else
              k = -n
             endif
            else
             k = n
            endif
            call mp_allocate(gpt,k,istat)
           else
            istat = 1
           endif
          else
           if(n.gt.0) then
            if(nnodes.gt.1) then
             k = n
            else
             k = -n
            endif
           else
            k = n
           endif
           call mp_allocate(gpt,k,istat)
          endif
          if(istat.eq.0) then
           nalloc = i
           maxnalloc = max(maxnalloc,nalloc)
           sortmemory = max(sortmemory,maxnalloc)
           nmp_alloc = j
           maxnmp_alloc = max(maxnmp_alloc,nmp_alloc)
          else
#ifdef FMM_INFO
           write(6,*) ' allocated memory:',nalloc,'requested memory:',i
#endif
           call mp_error(istat)
          endif
         else
          call mp_error(n)
         endif
        end subroutine mp_fmmalloc
#endif
c
#ifdef FMM_PARALLEL
        subroutine mp_fmmdealloc(gpt,n)
         use fmmkinds
         use mp_info, only: nnodes,mp_deallocate,mp_error
         implicit none
         type(c_ptr) gpt
         integer(kind=fmm_integer) n,istat
         if(n.ge.0) then
          if(n.gt.0) then
           if(nnodes.gt.1) then
            istat = n
           else
            istat = -n
           endif
          else
           istat = n
          endif
          call mp_deallocate(gpt,istat)
          if(istat.eq.0) then
           nalloc = nalloc-n
           nmp_alloc = nmp_alloc-n
          else
#ifdef FMM_INFO
           write(6,*) ' mp_deallocate failed, istat = ',istat
#endif
           call mp_error(istat)
          endif
         else
          call mp_error(n)
         endif
        end subroutine mp_fmmdealloc
#endif
c
#ifdef FMM_ISO_C_BINDING
#ifdef FMM_ALLOCALIGNED
        subroutine fmmalloc_aligned1d(cpt,ia,ib,istat)
         use fmmkinds
         use allocation
         implicit none
         type(c_ptr) cpt
         integer(kind=fmm_integer), intent(in):: ia,ib
         integer(kind=fmm_integer), intent(out):: istat
         integer(kind=fmm_integer) i
         i = nalloc+rtob*(ib-ia+1)
         if(maxallocation.gt.0) then
          if(maxallocation.ge.i) then
           call allocate_aligned(cpt,ia,ib,istat)
          else
           istat = 1
          endif
         else
          call allocate_aligned(cpt,ia,ib,istat)
         endif
         if(istat.eq.0) then
          nalloc = i
          maxnalloc = max(maxnalloc,nalloc)
          sortmemory = max(sortmemory,maxnalloc)
         else
#ifdef FMM_INFO
          write(6,*) ' allocated memory:',nalloc,'requested memory:',i
#endif
         endif
        end subroutine fmmalloc_aligned1d
#endif
#endif
c
#ifdef FMM_ISO_C_BINDING
#ifdef FMM_ALLOCALIGNED
        subroutine fmmalloc_aligned2d(cpt,ia,ib,ja,jb,istat)
         use fmmkinds
         use allocation
         implicit none
         type(c_ptr) cpt
         integer(kind=fmm_integer), intent(in):: ia,ib,ja,jb
         integer(kind=fmm_integer), intent(out):: istat
         integer(kind=fmm_integer) i
         i = nalloc+rtob*(ib-ia+1)*(jb-ja+1)
         if(maxallocation.gt.0) then
          if(maxallocation.ge.i) then
           call allocate_aligned(cpt,ia,ib,ja,jb,istat)
          else
           istat = 1
          endif
         else
          call allocate_aligned(cpt,ia,ib,ja,jb,istat)
         endif
         if(istat.eq.0) then
          nalloc = i
          maxnalloc = max(maxnalloc,nalloc)
          sortmemory = max(sortmemory,maxnalloc)
         else
#ifdef FMM_INFO
          write(6,*) ' allocated memory:',nalloc,'requested memory:',i
#endif
         endif
        end subroutine fmmalloc_aligned2d
#endif
#endif
c
#ifdef FMM_ISO_C_BINDING
#ifdef FMM_ALLOCALIGNED
        subroutine fmmdealloc_aligned(cpt,n)
         use fmmkinds
         use allocation
         implicit none
         type(c_ptr) cpt
         integer(kind=fmm_integer) n,i
         i = rtob*n
         call deallocate_aligned(cpt,n)
         if(n.eq.0) then
          cpt = c_null_ptr
          nalloc = nalloc-i
         endif
        end subroutine fmmdealloc_aligned
#endif
#endif
c
        subroutine indallocint1d(a,ia,ib,istat)
         use fmmkinds
         implicit none
         integer(kind=fmm_integer), allocatable:: a(:)
         integer(kind=fmm_integer), intent(in):: ia,ib
         integer(kind=fmm_integer), intent(out):: istat
         integer(kind=fmm_integer) i
         i = nalloc+itob*(ib-ia+1)
         if(corememory.gt.0) then
          if(corememory.ge.i) then
           allocate(a(ia:ib),stat = istat)
          else
           istat = 1
          endif
         else
          allocate(a(ia:ib),stat = istat)
         endif
         if(istat.eq.0) then
          nalloc = i
          maxnalloc = max(maxnalloc,nalloc)
          sortmemory = max(sortmemory,maxnalloc)
         endif
        end subroutine indallocint1d
c
        subroutine inddeallocint1d(a,istat)
         use fmmkinds
         implicit none
         integer(kind=fmm_integer), allocatable:: a(:)
         integer(kind=fmm_integer), intent(out):: istat
         integer(kind=fmm_integer) i
         i = size(a)
         deallocate(a,stat = istat)
         if(istat.eq.0) nalloc = nalloc-itob*i
        end subroutine inddeallocint1d
c
        subroutine indallocipt1d(a,ia,ib,istat)
         use fmmkinds
         implicit none
         integer(kind=fmm_integer), pointer:: a(:)
         integer(kind=fmm_integer), intent(in):: ia,ib
         integer(kind=fmm_integer), intent(out):: istat
         integer(kind=fmm_integer) i
         i = nalloc+itob*(ib-ia+1)
         if(corememory.gt.0) then
          if(corememory.ge.i) then
           allocate(a(ia:ib),stat = istat)
          else
           istat = 1
          endif
         else
          allocate(a(ia:ib),stat = istat)
         endif
         if(istat.eq.0) then
          nalloc = i
          maxnalloc = max(maxnalloc,nalloc)
          sortmemory = max(sortmemory,maxnalloc)
         endif
        end subroutine indallocipt1d
c
        subroutine inddeallocipt1d(a,istat)
         use fmmkinds
         implicit none
         integer(kind=fmm_integer), pointer:: a(:)
         integer(kind=fmm_integer), intent(out):: istat
         integer(kind=fmm_integer) i
         i = size(a)
         deallocate(a,stat = istat)
         a => null()
         if(istat.eq.0) nalloc = nalloc-itob*i
        end subroutine inddeallocipt1d
c
        subroutine srtallocbyte1d(a,ia,ib,istat)
         use fmmkinds
         implicit none
         integer(kind=1), allocatable:: a(:)
         integer(kind=fmm_integer), intent(in):: ia,ib
         integer(kind=fmm_integer), intent(out):: istat
         integer(kind=fmm_integer) i
         i = nalloc+(ib-ia+1)
         if(sortmemory.gt.0) then
          if(sortmemory.ge.i) then
           allocate(a(ia:ib),stat = istat)
          else
           istat = 1
          endif
         else
          allocate(a(ia:ib),stat = istat)
         endif
         if(istat.eq.0) then
          nalloc = i
          maxnalloc = max(maxnalloc,nalloc)
          sortmemory = max(sortmemory,maxnalloc)
         endif
        end subroutine srtallocbyte1d
c
        subroutine srtallocdbl2d(a,ia,ib,ja,jb,istat)
         use fmmkinds
         implicit none
         real(kind=fmm_real), allocatable:: a(:,:)
         integer(kind=fmm_integer), intent(in):: ia,ib,ja,jb
         integer(kind=fmm_integer), intent(out):: istat
         integer(kind=fmm_integer) i
         i = nalloc+(ib-ia+1)*(jb-ja+1)
         if(sortmemory.gt.0) then
          if(sortmemory.ge.i) then
           allocate(a(ia:ib,ja:jb),stat = istat)
          else
           istat = 1
          endif
         else
          allocate(a(ia:ib,ja:jb),stat = istat)
         endif
         if(istat.eq.0) then
          nalloc = i
          maxnalloc = max(maxnalloc,nalloc)
          sortmemory = max(sortmemory,maxnalloc)
         endif
        end subroutine srtallocdbl2d
c
        subroutine srtdeallocbyte1d(a,istat)
         use fmmkinds
         implicit none
         integer(kind=1), allocatable:: a(:)
         integer(kind=fmm_integer), intent(out):: istat
         integer(kind=fmm_integer) i
         i = size(a)
         deallocate(a,stat = istat)
         if(istat.eq.0) nalloc = nalloc-i
        end subroutine srtdeallocbyte1d
c
        subroutine srtdeallocdbl2d(a,istat)
         use fmmkinds
         implicit none
         real(kind=fmm_real), allocatable:: a(:,:)
         integer(kind=fmm_integer), intent(out):: istat
         integer(kind=fmm_integer) i
         i = size(a)
         deallocate(a,stat = istat)
         if(istat.eq.0) nalloc = nalloc-i
        end subroutine srtdeallocdbl2d
c
        subroutine dbpallocdbl1d(a,ia,ib,istat)
         use fmmkinds
         implicit none
         real(kind=fmm_real_extended), allocatable:: a(:)
         integer(kind=fmm_integer), intent(in):: ia,ib
         integer(kind=fmm_integer), intent(out):: istat
         integer(kind=fmm_integer) i
         i = nalloc+rtob*(ib-ia+1)
         if(maxallocation.gt.0) then
          if(maxallocation.ge.i) then
           allocate(a(ia:ib),stat = istat)
          else
           istat = 1
          endif
         else
          allocate(a(ia:ib),stat = istat)
         endif
         if(istat.eq.0) then
          nalloc = i
          maxnalloc = max(maxnalloc,nalloc)
          sortmemory = max(sortmemory,maxnalloc)
         else
#ifdef FMM_INFO
          write(6,*) ' allocated memory:',nalloc,'requested memory:',i
#endif
         endif
        end subroutine dbpallocdbl1d
c
        subroutine dbpallocdbl3d(a,ia,ib,ja,jb,ka,kb,istat)
         use fmmkinds
         implicit none
         real(kind=fmm_real_extended), allocatable:: a(:,:,:)
         integer(kind=fmm_integer), intent(in):: ia,ib,ja,jb,ka,kb
         integer(kind=fmm_integer), intent(out):: istat
         integer(kind=fmm_integer) i
         i = nalloc+rtob*(ib-ia+1)*(jb-ja+1)*(kb-ka+1)
         if(maxallocation.gt.0) then
          if(maxallocation.ge.i) then
           allocate(a(ia:ib,ja:jb,ka:kb),stat = istat)
          else
           istat = 1
          endif
         else
          allocate(a(ia:ib,ja:jb,ka:kb),stat = istat)
         endif
         if(istat.eq.0) then
          nalloc = i
          maxnalloc = max(maxnalloc,nalloc)
          sortmemory = max(sortmemory,maxnalloc)
         else
#ifdef FMM_INFO
          write(6,*) ' allocated memory:',nalloc,'requested memory:',i
#endif
         endif
        end subroutine dbpallocdbl3d
c
        subroutine dbpdeallocdbl1d(a,istat)
         use fmmkinds
         implicit none
         real(kind=fmm_real_extended), allocatable:: a(:)
         integer(kind=fmm_integer), intent(out):: istat
         integer(kind=fmm_integer) i
         i = size(a)
         deallocate(a,stat = istat)
         if(istat.eq.0) nalloc = nalloc-rtob*i
        end subroutine dbpdeallocdbl1d
c
        subroutine dbpdeallocdbl3d(a,istat)
         use fmmkinds
         implicit none
         real(kind=fmm_real_extended), allocatable:: a(:,:,:)
         integer(kind=fmm_integer), intent(out):: istat
         integer(kind=fmm_integer) i
         i = size(a)
         deallocate(a,stat = istat)
         if(istat.eq.0) nalloc = nalloc-rtob*i
        end subroutine dbpdeallocdbl3d
      end module fmmalloc

