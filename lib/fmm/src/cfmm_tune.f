#include "fmm.h"
#include "fmmkinds.h"

      subroutine cfmm_tune(ncharges,q,xyzin,ierr,de,der,
     .energy,
     .fmmpot,
     .fmmgrad,
     .periodic,
     .periodica,
     .periodlength,
     .dipolecorrection,
     .ilinearpotential,
     .linearodistance,
     .iplummerpotential,
     .aoplummer,
     .snewerroranalysis,
     .homogen,
     .maxdepth,
     .unrolled2,
     .balance_load,
     .FMM_internal_params)
c
      use fmmkinds
      use fmmint34
      use fmmjmp
      use fmmhybrid
      use fmmnsqrndiv
      use fmmalloc
      use getneighbors_vars
      use qinfo
      use msort
      use fmmicharge1icharge2
      use fmmicharge5icharge6
      use fmmjcharge1jcharge2
      use mcoordinates
      use mem_info
      use mplummer
      use mwigner
#ifdef FMM_UNIFORMGRID
      use muniformgrid
#endif
      use fmm_fcs_binding 
#ifdef FMM_COMPRESSION
      use compression
#endif
#ifdef FMM_TREETOGRAD
      use mtreetograd
#endif
#ifdef FMM_DAMPING
      use mdamping
#endif
#ifdef FMM_PARALLEL
      use mp_info
#ifdef FMM_LOADSORT
      use mp_load
#endif
#endif
c
      implicit none
c
      integer(kind=fmm_integer) maxnmultipoles
      parameter(maxnmultipoles=FMM_MAXNMULTIPOLES)
c
      integer(kind=fmm_integer) maxws
      parameter(maxws=1)
c
      integer(kind=fmm_integer) maxncsar
      parameter(maxncsar=7*maxws*maxws*maxws+18*maxws*maxws+14*maxws+5)
c
      integer(kind=fmm_integer) maxncar
      parameter(maxncar=16*maxws*maxws+24*maxws+8)
c
      integer(kind=fmm_integer) maxnrar
      parameter(maxnrar=7*maxws*maxws*maxws+21*maxws*maxws+21*maxws+7)
c
      integer(kind=fmm_integer) maxn2multipoles
      parameter(maxn2multipoles=maxnmultipoles+maxnmultipoles)
c
      integer(kind=fmm_integer) maxwsd
      parameter(maxwsd=2*maxws+1)
c
      integer(kind=fmm_integer) mmaxwsd
      parameter(mmaxwsd=-maxwsd)
c
      integer(kind=fmm_integer) maxwsd2
      parameter(maxwsd2=2*maxwsd*maxwsd)
c
      integer(kind=fmm_integer) maxwsd3
      parameter(maxwsd3=3*maxwsd*maxwsd)
c
      real(kind=fmm_real), allocatable:: dbl(:)
c
#ifdef FMM_COMPRESSION
      integer(kind=1), allocatable:: iboxsrt(:,:)
#else
      integer(kind=fmm_integer), allocatable, target:: iboxsrt(:)
#endif
      integer(kind=fmm_integer), allocatable:: iboxjmp(:)
#ifdef FMM_IBOXSCR
      integer(kind=fmm_integer), allocatable, target:: iboxscr(:)
#endif
      integer(kind=fmm_integer), pointer:: piboxscr(:),piboxsrt(:)
      integer(kind=fmm_integer), allocatable:: int3x(:),int3y(:),
     .int3z(:),int3p(:),int3q(:),nbofmb(:)
      integer(kind=fmm_integer) ipo(3),jpo(3),mask(3)
c
      real(kind=fmm_real) fmmdist(mmaxwsd:maxwsd,mmaxwsd:maxwsd,
     .mmaxwsd:maxwsd)
      integer(kind=fmm_integer) gb(2,maxwsd*maxwsd*maxwsd)
      real(kind=fmm_real) gbsh(3,maxwsd*maxwsd*maxwsd)
c
      integer(kind=fmm_integer) ncsar,icsar(0:maxwsd,0:maxwsd2),
     .jcsar(maxncsar)
      integer(kind=fmm_integer) ncar,icar(mmaxwsd:maxwsd,mmaxwsd:maxwsd)
      integer(kind=fmm_integer) isar(mmaxwsd:maxwsd,0:maxwsd)
      integer(kind=fmm_integer) nrar
      integer(kind=fmm_integer), allocatable:: irar(:,:)
      real(kind=fmm_real), allocatable:: grar(:,:)
      logical(kind=fmm_logical), allocatable:: sgrar(:)
c
      integer(kind=fmm_integer) nfmmcos(maxwsd3),fmmcos(2,maxncsar)
c
      integer(kind=fmm_integer) mi
      parameter(mi=-1)
c
      integer(kind=fmm_integer) jcar(mi:1,mi:1)
      real(kind=fmm_real) hcar(0:maxnmultipoles,4),
     .hsar(0:maxnmultipoles,4)
c
      integer(kind=fmm_integer) ierr
c
      integer(kind=fmm_integer) ldf
      parameter(ldf=200)
      integer(kind=fmm_integer) inf
      parameter(inf=ldf)
      integer(kind=fmm_integer) ldff
      parameter(ldff=ldf+ldf)
      real(kind=fmm_real) fmmerr(0:maxnmultipoles,maxws),
     .pfmmerr(0:maxnmultipoles),
     .merr(0:maxnmultipoles,maxws)
c
      integer(kind=fmm_integer) ncharges
      real(kind=fmm_real), target:: xyzin(3,ncharges)
      real(kind=fmm_real), allocatable, target:: xyz(:,:)
      real(kind=fmm_real), pointer:: xyzt(:,:)
      real(kind=fmm_real) q(*)
      real(kind=fmm_real) coul,shx,shy,shz,sf,sh,fac(0:170),rfac(0:170)
      real(kind=fmm_real) pow(0:maxnmultipoles),sg(0:maxn2multipoles)
      real(kind=fmm_real) fr(0:maxn2multipoles)
      real(kind=fmm_real), allocatable:: coeff1(:,:),coeff2(:,:),
     .coeff3(:,:,:),
     .coeff4(:,:),coeff5(:,:,:),coeff6(:,:)
      real(kind=fmm_real) energy
#ifdef FMM_NOPOT
      real(kind=fmm_real), target:: fmmgrad(3*ncharges)
      real(kind=fmm_real) fmmpot(*)
#else
#ifdef FMM_TREETOGRAD
      real(kind=fmm_real), target:: fmmgrad(3*ncharges)
#else
      real(kind=fmm_real) fmmgrad(*)
#endif
      real(kind=fmm_real), target:: fmmpot(ncharges)
#endif
      real(kind=fmm_real), pointer:: pfmmpot(:)
      real(kind=fmm_real), allocatable:: d2(:,:,:),d3(:,:,:),d2f(:,:,:),
     .d3f(:,:,:)
      real(kind=fmm_real), allocatable:: bfgn(:)
      real(kind=fmm_real) rl(0:maxn2multipoles),cmphi(0:maxn2multipoles)
      real(kind=fmm_real) smphi(0:maxn2multipoles),
     .cmphipi(0:maxn2multipoles)
      real(kind=fmm_real) smphipi(0:maxn2multipoles)
      real(kind=fmm_real), allocatable, target:: omegatree(:)
      real(kind=fmm_real), allocatable, target:: mutree(:)
      real(kind=fmm_real), pointer:: romegatree(:),iomegatree(:),
     .rmutree(:),imutree(:)
      integer(kind=fmm_integer), allocatable:: taylor(:)
      real(kind=fmm_real) delta,de,der,dem,periodlength,
     .linearodistance(*),aoplummer,lineardistance(0:3),efarfield,e1per,
     .enearfield,enfinbox,enfbibj,virialtensor(3,3)
      real(kind=fmm_real_extended) efarfieldpot
      real(kind=fmm_real) ctheta,stheta
      real(kind=fmm_real) calpha,salpha,cbeta,sbeta
c
      integer(kind=fmm_integer), allocatable:: indscr(:)
c
      integer(kind=fmm_integer) buflen
      parameter(buflen=131072)
      integer(kind=fmm_integer) bfglen
      parameter(bfglen=4*buflen)
      real(kind=fmm_real) bfg(buflen,4),csar(buflen),car(buflen),
     .sar(buflen),
     .rar(buflen)
      integer(kind=fmm_integer) jibfglen
      parameter(jibfglen=6*buflen)
      integer(kind=fmm_integer) jibfg(buflen,6),isrt(buflen),
     .kbxyzar(buflen),
     .indar(buflen),kboxxyzar(buflen),kboxindar(buflen),kbar(buflen)
      equivalence(bfg,csar)
      equivalence(bfg(1,2),car)
      equivalence(bfg(1,3),sar)
      equivalence(bfg(1,4),rar)
      equivalence(jibfg,isrt)
      equivalence(jibfg(1,2),kbxyzar)
      equivalence(jibfg(1,3),indar)
      equivalence(jibfg(1,4),kboxxyzar)
      equivalence(jibfg(1,5),kboxindar)
      equivalence(jibfg(1,6),kbar)
c
      real(kind=fmm_real) fracdepth,shmonopole
c
      real(kind=fmm_real), allocatable:: flvlar(:),powsq(:)
c
      integer(kind=fmm_integer) parabola,ilevelmn,bfgnlen
      real(kind=fmm_real) cx,cy,cz
c
#ifdef FMM_TREETOGRAD
      integer(kind=fmm_integer) startbox,endbox,pagejump,pageshift,
     .pageshiftg,pagemask,pageaddr,indsize,pagepossize,indskpjump
c
      logical(kind=fmm_logical) pages
#endif
c
      integer(kind=fmm_integer) maxdepth,maxdepthp,mmaxdepth,nbytes,
     .nbits,maxint,maxmint,i,depth,nmultipoles,mnmultipoles,
     .n2multipoles,nsqmultipoles,j,jnbi,nbfg,ntree
      integer(kind=fmm_integer), allocatable:: bitpos(:),mbitpos(:),
     .nboxesinlevel(:),nboxeslevel(:)
c
#if defined(FMM_IBOXSCR) && !defined(FMM_COMPRESSION)
      integer(kind=fmm_integer), allocatable:: ibox(:)
#else
      integer(kind=fmm_integer), allocatable, target:: ibox(:)
#endif
c
      logical(kind=fmm_logical) dfmmmerr(maxws)
c
      integer(kind=fmm_integer) k,l,m,igtaylor,mgtaylor,ntaylor,ws
      real(kind=fmm_real) enearfieldpot,energypot
c
      integer(kind=fmm_integer) ishx,maskx,ishy,masky,mishx,mishy,
     .maskxy
c
      real(kind=fmm_real) qqq,corrsx,corrsy,corrsz,corrs,corrsh
      real(kind=fmm_real) gp,gsq
#ifdef FMM_DEBUG
      real(kind=fmm_real) sgradx,sgrady,sgradz
      real(kind=fmm_real) sagradx,sagrady,sagradz
#endif
c
      integer(kind=fmm_integer) periodic,periodica,dipolecorrection,
     .ilinearpotential,iplummerpotential,homogen
      logical(kind=fmm_logical) compute
      integer(kind=fmm_integer) negpos
      logical(kind=fmm_logical) nothing,linearpotential,copyxyz,sh4,sh3,
     .changepos,shmp,cachopt,g2db,precomputeallds,withaop,withbop,
     .withcop,withtaylor,associatedwignerd,gtaylor,doit
c
      logical(kind=fmm_logical) hugep(0:100)
      real(kind=fmm_real) hugef(100)
c
      real(kind=fmm_real) zero
      parameter(zero=0.e0_fmm_real)
      real(kind=fmm_real) one
      parameter(one=1.e0_fmm_real)
      real(kind=fmm_real) two
      parameter(two=2.e0_fmm_real)
      real(kind=fmm_real) three
      parameter(three=3.e0_fmm_real)
      real(kind=fmm_real) half
      parameter(half=one/two)
      real(kind=fmm_real_extended) zero_extended
      parameter(zero_extended=0.e0_fmm_real_extended)
c
      integer(kind=fmm_integer) n
      integer(kind=fmm_integer) unrolled2,pgd
c
      integer(kind=fmm_integer) snewerroranalysis
      integer(kind=fmm_integer) serroranalysis,nerroranalysis
c      logical(kind=fmm_logical), save:: firsterroranalysis = .true.
      logical(kind=fmm_logical) firsterroranalysis
      logical(kind=fmm_logical) erroranalysis
c      save pgd,depth,fracdepth,nmultipoles,shmonopole,parabola
      real(kind=fmm_real) linearodistancesv(3),aoplummersv
      integer(kind=fmm_integer) ilinearpotentialsv,iplummerpotentialsv
c
      type(twignerd) wignerd
c
      integer(kind=fmm_integer) nallocst
c-ik         
      type(FMM_internal_params_t):: FMM_internal_params
      integer(kind=fmm_integer) :: balance_load
c-ik
#ifdef FMM_PARALLEL
      call mp_rank(MP_ALLNODES,me)
      call mp_nnodes(MP_ALLNODES,nnodes)
#ifdef FMM_INFO
      if(me.eq.0) write(6,*) ' nnodes: ',nnodes
#endif
#endif
c-ik 
#ifdef FMM_LOADSORT 
      if(balance_load.gt.0) then
c         doload = .true.  
         doload = .false.
      else
         doload = .false.
      endif
#endif
      serroranalysis = FMM_internal_params%serroranalysis
      nerroranalysis = FMM_internal_params%nerroranalysis
      pgd = FMM_internal_params%pgd
      depth = FMM_internal_params%depth
      nmultipoles = FMM_internal_params%nmultipoles
      parabola = FMM_internal_params%parabola
      ilinearpotentialsv = FMM_internal_params%ilinearpotentialsv
      iplummerpotentialsv = FMM_internal_params%iplummerpotentialsv
      firsterroranalysis = FMM_internal_params%firsterroranalysis
      hugep = FMM_internal_params%hugep
      fracdepth = FMM_internal_params%fracdepth
      shmonopole = FMM_internal_params%shmonopole
      linearodistancesv = FMM_internal_params%linearodistancesv
      aoplummersv = FMM_internal_params%aoplummersv
      hugef = FMM_internal_params%hugef
      wignerd%wignerd => FMM_internal_params%wignerd%wignerd

c-ik
c      call seticharge1icharge2(maxncharges)
      icharge1 = 1
      icharge2 = ncharges
      icharge5 = icharge1
      icharge6 = icharge2
      icharges = ncharges

#ifdef FMM_PARALLEL
#ifdef FMM_LOADSORT
      if(doload) then 
        stop "error doload & tune"
      else
        ichargesout = icharges
        micharge1 = icharge1
        micharge2 = icharge2
        micharges = icharges
      endif
#else
      ichargesout = icharges 
      micharge1 = icharge1
      micharge2 = icharge2
      micharges = icharges
#endif
#else
      ichargesout = icharges
      micharge1 = icharge1
      micharge2 = icharge2 
      micharges = icharges
#endif
c
      call intbyt(nbytes,nbits,maxint,maxmint)
c
      call setmdmaxallocation(maxallocation)
      call setmdcorememorysortmemory(corememory,sortmemory)
      call chkmdmaxallocationcorememorysortmemory(maxallocation,
     .corememory,sortmemory)
c
      call setmdfmmalloc(nbytes)
c
      call setnalloctomaxnalloc(nalloc,maxnalloc)
c
#ifdef FMM_STATISTICS
      call fmmstatistics(ncharges,.true.)
#endif
c
#ifdef FMM_PARALLEL
#if defined(FMM_CPUTIME) || defined(FMM_WALLTIME)
      call c_cputime_walltime(.true.)
#endif
#endif
c
      call stmdfmmalloc(nalloc,nallocst)
      call prtmdfmmalloc(nalloc,maxnalloc,'start of cfmm_tune')
c
#ifdef FMM_PARALLEL
      mp_setup = .false.
#endif
c
c      call seticharge1icharge2(ncharges)
c
#ifdef FMM_PARALLEL
      if(nnodes.gt.0) then
         i = nnodes-1
         call seticharge3(i)
      else
         call bummer('fmm: error, nnodes = ',nnodes)
      endif
#else
      icharge3 = icharge1
#endif
c
#ifdef FMM_NOPOT
      pfmmpot => fmmgrad
#else
      pfmmpot => fmmpot
#endif
c
      call setfmmi0(ierr,de)
c
      if(de.gt.zero) then
       if(ldf.gt.maxnmultipoles) then
        call calnegpos(periodic,q,energy,pfmmpot,fmmgrad,negpos,nothing)
        if(nothing) then
         pfmmpot => null()
#ifdef FMM_PARALLEL
#ifdef FMM_LOADSORT
         if(doload) iboxload => null()
#endif
#endif
         call edmdfmmalloc(nalloc,nallocst,'cfmm_tune')
         call prtmdfmmalloc(nalloc,maxnalloc,'  end of cfmm_tune')
#ifdef FMM_PARALLEL
#if defined(FMM_CPUTIME) || defined(FMM_WALLTIME)
         call c_cputime_walltime(.false.)
#endif
#endif
#ifdef FMM_STATISTICS
         call fmmstatistics(ncharges,.false.)
#endif
         return
        endif
       else
        call bummer('fmm: (ldf-maxnmultipoles) = ',(ldf-maxnmultipoles))
       endif
      else
       call bummer('fmm: error, ierr = ',ierr)
      endif
c
c      call intbyt(nbytes,nbits,maxint,maxmint)
c
c      call intmaxdepth(maxdepth,nbits)
c
      maxdepthp = maxdepth+1
      mmaxdepth = -maxdepth
c
#ifdef FMM_PARALLEL
      call calkind_integer_processes(nbits,maxmint,maxint)
#endif
c
      if(firsterroranalysis) then
         if(snewerroranalysis.gt.0) then
            serroranalysis = snewerroranalysis
            nerroranalysis = 1
            erroranalysis = .true.
            ilinearpotentialsv = ilinearpotential
            linearodistancesv(1) = linearodistance(1)
            linearodistancesv(2) = linearodistance(2)
            linearodistancesv(3) = linearodistance(3)
            iplummerpotentialsv = iplummerpotential
            aoplummersv = aoplummer
            firsterroranalysis = .false.
         else
            call bummer('fmm: error, snewerroranalysis = ',
     .      snewerroranalysis)
         endif
      elseif(serroranalysis.gt.0) then
         if(serroranalysis.gt.nerroranalysis) then
            erroranalysis = .false.
            nerroranalysis = nerroranalysis+1
         elseif(serroranalysis.eq.nerroranalysis) then
            nerroranalysis = 1
            erroranalysis = .true.
         else
            call bummer('fmm: (serroranalysis-nerroranalysis) = ',
     .      (serroranalysis-nerroranalysis))
         endif
      else
         call bummer('fmm: error, serroranalysis = ',serroranalysis)
      endif
c
#ifdef FMM_NOPOT
      if(erroranalysis.and.(homogen.eq.0)) then
         call bummer('fmm: not yet implemented, homogen = ',homogen)
      endif
#endif
c
      if(erroranalysis) then
         if(maxdepth.ge.0) then
c            call fmmallocate(hugep,0,maxdepthp,i)
c            if(i.eq.0) then
c               call fmmallocate(hugef,1,maxdepthp,i)
c               if(i.eq.0) then
c                  i = ltob*(maxdepthp+1)+rtob*(maxdepthp)
c                  nallocr = i
c                  nallocst = nallocst+i
                  hugep(0) = .false.
                  do 337 i = 1,maxdepthp
                     hugep(i) = .false.
                     hugef(i) = one
 337              continue
c               else
c                  call bummer('fmm: error, i = ',i)
c               endif
c            else
c               call bummer('fmm: error, i = ',i)
c            endif
         else
            call bummer('fmm: error, maxdepth = ',maxdepth)
         endif
      endif
c
#ifdef FMM_UNIFORMGRID
      if(periodic.eq.0) then
       if(iplummerpotential.eq.0) then
        if(ilinearpotential.eq.0) then
         plummer_potential = .false.
         linearpotential = .false.
         if(uniformgridox.gt.zero) then
          if(uniformgridoy.gt.zero) then
           if(uniformgridoz.gt.zero) then
            uniformgridm = max(uniformgridox,uniformgridoy,
     .      uniformgridoz)
            uniformgridx = uniformgridox
            uniformgridy = uniformgridoy
            uniformgridz = uniformgridoz
            uniformgridalpha = uniformgridoalpha
           else
            call bummer('fmm: error, uniformgridoz <= ',0)
           endif
          else
           call bummer('fmm: error, uniformgridoy <= ',0)
          endif
         else
          call bummer('fmm: error, uniformgridox <= ',0)
         endif
        else
         call bummer('fmm: error, ilinearpotential = ',ilinearpotential)
        endif
       else
        call bummer('fmm: error, iplummerpotential =',iplummerpotential)
       endif
      else
       call bummer(
     . 'fmm: FMM_UNIFORMGRID & periodic not yet implemented, periodic=',
     . periodic)
      endif
#else
      if(iplummerpotential.eq.1) then
       if(ilinearpotential.eq.0) then
        if(aoplummer.gt.zero) then
         plummer_potential = .true.
         linearpotential = .false.
         a_plummer = aoplummer
        elseif(aoplummer.eq.zero) then
         plummer_potential = .false.
         linearpotential = .false.
         a_plummer = aoplummer
        else
         call bummer('fmm: error, aoplummer < ',0)
        endif
       else
        call bummer('fmm: error, ilinearpotential = ',ilinearpotential)
       endif
      elseif(iplummerpotential.eq.0) then
       plummer_potential = .false.
       if(ilinearpotential.eq.0) then
        linearpotential = .false.
       elseif((ilinearpotential.gt.0).and.(ilinearpotential.lt.8)) then
        linearpotential = .true.
c
        if(ilinearpotential.eq.1) then
         if(linearodistance(3).gt.zero) then
          lineardistance(0) = linearodistance(3)
          lineardistance(1) = zero
          lineardistance(2) = zero
          lineardistance(3) = linearodistance(3)
         else
          call bummer('fmm: error, ilinearpotential =',ilinearpotential)
         endif
        elseif(ilinearpotential.eq.2) then
         if(linearodistance(2).gt.zero) then
          lineardistance(0) = linearodistance(2)
          lineardistance(1) = zero
          lineardistance(2) = linearodistance(2)
          lineardistance(3) = zero
         else
          call bummer('fmm: error, ilinearpotential =',ilinearpotential)
         endif
        elseif(ilinearpotential.eq.3) then
         if(linearodistance(2).gt.zero) then
          if(linearodistance(3).gt.zero) then
           lineardistance(0)=max(linearodistance(2),linearodistance(3))
           lineardistance(1) = zero
           lineardistance(2) = linearodistance(2)
           lineardistance(3) = linearodistance(3)
          else
           call bummer('fmm: error, ilinearpotential=',ilinearpotential)
          endif
         else
          call bummer('fmm: error, ilinearpotential =',ilinearpotential)
         endif
        elseif(ilinearpotential.eq.4) then
         if(linearodistance(1).gt.zero) then
          lineardistance(0) = linearodistance(1)
          lineardistance(1) = linearodistance(1)
          lineardistance(2) = zero
          lineardistance(3) = zero
         else
          call bummer('fmm: error, ilinearpotential =',ilinearpotential)
         endif
        elseif(ilinearpotential.eq.5) then
         if(linearodistance(1).gt.zero) then
          if(linearodistance(3).gt.zero) then
           lineardistance(0)=max(linearodistance(1),linearodistance(3))
           lineardistance(1) = linearodistance(1)
           lineardistance(2) = zero
           lineardistance(3) = linearodistance(3)
          else
           call bummer('fmm: error, ilinearpotential=',ilinearpotential)
          endif
         else
          call bummer('fmm: error, ilinearpotential =',ilinearpotential)
         endif
        elseif(ilinearpotential.eq.6) then
         if(linearodistance(1).gt.zero) then
          if(linearodistance(2).gt.zero) then
           lineardistance(0)=max(linearodistance(1),linearodistance(2))
           lineardistance(1) = linearodistance(1)
           lineardistance(2) = linearodistance(2)
           lineardistance(3) = zero
          else
           call bummer('fmm: error, ilinearpotential=',ilinearpotential)
          endif
         else
          call bummer('fmm: error, ilinearpotential =',ilinearpotential)
         endif
        elseif(linearodistance(1).gt.zero) then
         if(linearodistance(2).gt.zero) then
          if(linearodistance(3).gt.zero) then
           lineardistance(0)=max(linearodistance(1),linearodistance(2),
     .     linearodistance(3))
           lineardistance(1) = linearodistance(1)
           lineardistance(2) = linearodistance(2)
           lineardistance(3) = linearodistance(3)
          else
           call bummer('fmm: error, ilinearpotential=',ilinearpotential)
          endif
         else
          call bummer('fmm: error, ilinearpotential =',ilinearpotential)
         endif
        else
         call bummer('fmm: error, ilinearpotential = ',ilinearpotential)
        endif
       else
        call bummer('fmm: error, ilinearpotential = ',ilinearpotential)
       endif
      else
       call bummer('fmm: error, iplummerpotential = ',iplummerpotential)
      endif
#endif
c
      call factorial(170,fac)
      call rfactorial(170,rfac)
c
      call fmmallocate(flvlar,1,maxdepthp,i)
      if(i.eq.0) then
         call fmmlvl(maxdepthp,flvlar)
      else
         call bummer('fmm: error, i = ',i)
      endif
c
      call fmmallocate(powsq,mmaxdepth,maxdepthp,i)
      if(i.eq.0) then
         call calpowsq(maxdepth,mmaxdepth,powsq)
      else
         call bummer('fmm: error, i = ',i)
      endif
c
      if(nbits.gt.0) then
         i = nbits-1
         call fmmallocate(bitpos,0,i,j)
         if(j.eq.0) then
            call fmmallocate(mbitpos,0,i,j)
            if(j.eq.0) then
               call calbitpmbp(i,maxmint,bitpos,mbitpos)
            else
               call bummer('fmm: error, j = ',j)
            endif
         else
            call bummer('fmm: error, j = ',j)
         endif
      else
         call bummer('fmm: error, nbits = ',nbits)
      endif
c
      ws = 1
c
      if(ws.le.0) call bummer('fmm: error, ws = ',ws)
      if(ws.gt.maxws) call bummer('fmm: (ws-maxws) = ',(ws-maxws))
c
#ifdef FMM_INFO
#ifdef FMM_PARALLEL
      if(me.eq.0) then
#endif
         write(6,'(''          ws = '',i13)') ws
#ifdef FMM_PARALLEL
      endif
#endif
#endif
c
      mem_s = 5
      mmem_s = -mem_s
      mem_l = bitpos(mem_s)-1
      if(mem_s.gt.1) then
         mem_m = bitpos(nbits+mmem_s)-1
      else
         mem_m = maxint
      endif
      mem_s2 = 4
      mmem_s2 = -mem_s2
      mem_l2 = bitpos(mem_s2)-1
      if(mem_s.gt.1) then
         mem_m2 = bitpos(nbits+mmem_s2)-1
      else
         mem_m2 = maxint
      endif
c
      copyxyz = .true.
      sh4 = .true.
      sh3 = .false.
      useqinsh = .true.
      changepos = .true.
      shmp = .false.
      cachopt = .true.
      g2db = .true.
c setup from scafaocs interface
      unrolled2 = 0
      precomputeallds = .false.
#ifdef FMM_PARALLEL
      withaop = .true.
      withbop = .true.
      withcop = .true.
      withtaylor = .false.
      jmp = .false.
      hybrid = .false.
#ifdef FMM_LOADSORT
c-ik      doload = .false.
#endif
#else
      withaop = .true.
      withbop = .true.
      withcop = .true.
      withtaylor = .true.
      jmp = .true.
      hybrid = .true.
#endif
      nsqr = 19
      ndiv = 7
c      nsqr = 31
c      ndiv = 4
      int3xyzd = maxdepth
      int3xyzd1 = int3xyzd+1
      if(int3xyzd.ge.0) then
         int3xyz = bitpos(int3xyzd)-1
      else
         int3xyz = int3xyzd
      endif
      if(erroranalysis) then
         pgd = 9
         call setpgd(bitpos,ws,maxdepth,pgd)
      endif
c
#ifdef FMM_DAMPING
      sh4 = .true.
      useqinsh = .true.
#else
      if(sh4.or.(.not.qscratch)) useqinsh = .true.
#endif
      if(.not.changepos) shmp = .false.
      if(.not.withbop) hybrid = .true.
c
      if(jmp) then
         call caljmps(ncharges,nbits,erroranalysis,homogen,maxmint,
     .   periodic,changepos,bitpos,mbitpos)
      else
         call caljmplj(jmpg,jmph,jmpp,jmpupd,jmpj,jmpn,jmpb,jmpjp)
      endif
c
      sf = one
      sh = one
      gp = one
      gsq = one
      gsq = -gsq
c
      efarfield = zero
      efarfieldpot = zero
      e1per = zero
      enearfield = zero
      enearfieldpot = zero
      enfinbox = zero
      enfbibj = zero
      coul = zero
c
      do 778 i = 1,maxws
         dfmmmerr(i) = .true.
 778  continue
c
c      call coulomb(ncharges,q,xyzin,coul)
#ifdef FMM_INFO
#ifdef FMM_PARALLEL
      if(me.eq.0) then
#endif
c         write(6,'(1x,q50.40,1x,z32)') coul,coul
#ifdef FMM_PARALLEL
      endif
#endif
#endif
c      call coulgrad(ncharges,q,xyzin)
c      call coulpot(ncharges,q,xyzin)
c
      call calishmsk(nbits,ishx,maskx,ishy,masky,mishx,mishy,maskxy)
c
      ipo(1) = 0
      ipo(2) = ishy
      ipo(3) = ishx
c
      jpo(1) = 0
      jpo(2) = mishy
      jpo(3) = mishx
c
      mask(1) = maskxy
      mask(2) = masky
      mask(3) = maskx
c
      call fmmg(maxwsd,maxncsar,ws,ncsar,icsar,fmmcos,.false.)
c
      if(maxnmultipoles.gt.0) then
         i = -maxnmultipoles
      else
         i = 0
      endif
c
      if(ncsar.gt.0) then
         j = ncsar+2
      else
         call bummer('fmm: error, ncsar = ',ncsar)
      endif
c
      associatedwignerd = associated(wignerd%wignerd)
c
      if(precomputeallds) then
         if(.not.associatedwignerd) then
            call fmmallocatept(wignerd%wignerd,0,maxnmultipoles,i,
     .      maxnmultipoles,0,maxnmultipoles,1,4,1,j,k)
            if(k.ne.0) call bummer('fmm: error, k = ',k)
         endif
c
         call calallds(ws,maxnmultipoles,ncsar,maxwsd,icsar,fmmcos,
     .   wignerd)
      endif
c
      call fmmh(maxnmultipoles,maxwsd,mmaxwsd,mi,maxncar,ws,ncar,icar,
     .jcar,hcar,hsar,.false.)
c
#ifdef FMM_INFO
#ifdef FMM_PARALLEL
      if(me.eq.0) then
#endif
         write(6,'(''    ncharges = '',i13)') ncharges
#ifdef FMM_PARALLEL
      endif
#endif
#endif
c
#ifdef FMM_COMPRESSION
#ifdef FMM_PARALLEL
#ifdef FMM_LOADSORT
      if(doload) call bummer('fmm: compression & loadsort, me = ',me)
#endif
#endif
#ifdef FMM_EXTREMECOMPRESSION
#ifdef FMM_EXTREMEEXTREMECOMPRESSION
      if(nchcompression.lt.0) then
         call bummer('fmm: error, nchcompression = ',nchcompression)
      endif
#endif
#endif
#endif
c
#ifdef FMM_COMPRESSION
      if(homogen.eq.1) then
       if(jmp) then
        call bummer('fmm: compression and jmp requested, homogen = ',
     .  homogen)
       elseif(nbytes.gt.0) then
        i = nbytes-1
        call fmmallocate(int1m,0,i,j)
        if(j.ne.0) call bummer('fmm: error, j = ',j)
        call fmmallocate(mint1,0,i,j)
        if(j.ne.0) call bummer('fmm: error, j = ',j)
        call fmmallocate(int1scr,0,i,j)
        if(j.ne.0) call bummer('fmm: error, j = ',j)
#ifdef FMM_EXTREMECOMPRESSION
        maxboxnumber = 0
        maxnp = 0
        call calmindp(ws,periodic,nbits,bitpos,mindp)
        if(mindp.ge.0) then
         if(mindp.le.maxdepth) then
          i = nbits-1
          j = 3*mindp
          if(j.lt.i) then
           mindpl = bitpos(j)
           mindpm = mindpl-1
           call fmmallocate(pinb,0,mindpm,i)
           if(i.eq.0) then
            do 776 i = 0,mindpm
             pinb(i) = 0
 776        continue
           else
            call bummer('fmm: error, i = ',i)
           endif
          else
           call bummer('fmm: (j-i) = ',(j-i))
          endif
         else
          call bummer('fmm: (mindp-maxdepth) = ',(mindp-maxdepth))
         endif
        else
         call bummer('fmm: error, mindp = ',mindp)
        endif
#endif
       else
        call bummer('fmm: error, nbytes = ',nbytes)
       endif
      else
       call bummer('fmm: compression requested, homogen = ',homogen)
      endif
#endif
c
#ifdef FMM_PARALLEL
      if(icharges.gt.0) then
#else
      if(ncharges.gt.0) then
#endif
         if(copyxyz) then
#ifdef FMM_PARALLEL
            call fmmallocate(xyz,1,3,micharge1,micharge2,i)
#else
            call fmmallocate(xyz,1,3,1,ncharges,i)
#endif
            if(i.eq.0) then
               xyzt => xyz
            else
               call bummer('fmm: error, i = ',i)
            endif
         else
            xyzt => xyzin
         endif
      else
#ifdef FMM_PARALLEL
         call bummer('fmm: error, icharges = ',icharges)
#else
         call bummer('fmm: error, ncharges = ',ncharges)
#endif
      endif
c
#ifdef FMM_DAMPING
      enfdba = .false.
#endif
c
      if(homogen.eq.1) then
       if(jmp) then
        call bummer('fmm: compression and jmp requested, homogen = ',
     .  homogen)
       else
#ifdef FMM_COMPRESSION
        call sccoorda(ws,nbits,bitpos,xyzin,xyzt,sf,periodic,periodica,
     .  periodlength,changepos)
c
        if(periodic.gt.0) then
         if(nj.gt.0) then
          i = nbits-1
          call fmmallocate(bitposm1,0,i,j)
          if(j.ne.0) call bummer('fmm: error, j = ',j)
          call calbitposm1(nbits,maxint,bitpos)
          call fmmallocate(mskj,1,nbytes,i)
          if(i.ne.0) call bummer('fmm: error, i = ',i)
          call calmskj(nbytes)
          call calm4(nbytes,nbits,bitpos,periodic)
          call fmmallocate(folder,1,nfolder,i)
          if(i.ne.0) call bummer('fmm: error, i = ',i)
         elseif(nj.lt.0) then
          call bummer('fmm: error, nj = ',nj)
         endif
        elseif(periodic.ne.0) then
         call bummer('fmm: error, periodic = ',periodic)
        endif
c
        call sccoordb(ws,nbits,bitpos,mbitpos,bfg,bfglen,q,xyzt,sf,
     .  periodic,linearpotential,lineardistance,changepos,negpos,
     .  dipolecorrection,qqq,corrsx,corrsy,corrsz,corrs,corrsh)
#endif
        if(erroranalysis) then
#ifdef FMM_COMPRESSION
         doit = .true.
#else
         doit = .false.
#endif
         if(doit) then
          call pass2dehomogen(maxdepth,ws,ncharges,maxnmultipoles,
     .    nmultipoles,inf,ldf,ldff,dfmmmerr,maxwsd,ncsar,icsar,
     .    nfmmcos,fmmcos,fmmerr,pfmmerr,merr,bitpos,nbits,maxint,
     .    depth,fracdepth,parabola,withbop,hugep,hugef,negpos,
     .    periodic,pgd,changepos,ierr,de,shmonopole,sf,
     .    linearpotential,lineardistance)
         endif
        endif
c
#ifdef FMM_COMPRESSION
#ifdef FMM_EXTREMECOMPRESSION
        call sccoord2(periodic,nbits,depth,fracdepth,shmonopole,sh,xyzt,
     .  linearpotential,lineardistance)
#endif
c
#ifdef FMM_EXTREMECOMPRESSION
        if(depth.ge.mindp) then
         mindps = 3*(mindp-depth)
        elseif(depth.ge.0) then
         mindps = 0
         mindp = depth
         i = nbits-1
         j = 3*mindp
         if(j.lt.i) then
          mindpl = bitpos(j)
          mindpm = mindpl-1
          do 777 i = 0,mindpm
           pinb(i) = 0
 777      continue
         else
          call bummer('fmm: (j-i) = ',(j-i))
         endif
        else
         call bummer('fmm: error, depth = ',depth)
        endif
#endif
#endif
       endif
#ifdef FMM_COMPRESSION
      else
       call bummer('fmm: compression requested, homogen = ',homogen)
#else
      elseif(homogen.ne.0) then
       call bummer('fmm: error, homogen = ',homogen)
#endif
      endif
c
#ifdef FMM_PARALLEL
      if(icharges.gt.0) then
#else
      if(ncharges.gt.0) then
#endif
#ifdef FMM_PARALLEL
#ifdef FMM_IBOXSCR
         call fmmallocate(iboxscr,micharge1,micharge2,i)
         if(i.ne.0) call bummer('fmm: error, i = ',i)
         piboxscr => iboxscr
         call fmmallocate(ibox,micharge1,micharge2,i)
         if(i.ne.0) call bummer('fmm: error, i = ',i)
#else
         call fmmallocate(ibox,micharge1,micharge2,i)
         if(i.ne.0) call bummer('fmm: error, i = ',i)
         piboxscr => ibox
#endif
#else
#ifdef FMM_IBOXSCR
         call fmmallocate(iboxscr,1,ncharges,i)
         if(i.ne.0) call bummer('fmm: error, i = ',i)
         piboxscr => iboxscr
         call fmmallocate(ibox,1,ncharges,i)
         if(i.ne.0) call bummer('fmm: error, i = ',i)
#else
         call fmmallocate(ibox,1,ncharges,i)
         if(i.ne.0) call bummer('fmm: error, i = ',i)
         piboxscr => ibox
#endif
#endif
      else
#ifdef FMM_PARALLEL
         call bummer('fmm: error, icharges = ',icharges)
#else
         call bummer('fmm: error, ncharges = ',ncharges)
#endif
      endif
c
#ifdef FMM_COMPRESSION
#ifdef FMM_EXTREMECOMPRESSION
#ifdef FMM_IBOXSCR
      call iboxscinfo(ishx,ishy,mmaxdepth,depth,bitpos,powsq,xyzt,ibox,
     .iboxsrt,iboxscr)
#else
      call iboxinfo(mmaxdepth,depth,bitpos,powsq,xyzt,ibox,ibox)
#endif
c
      if(depth.eq.0) then
       maxboxnumber = 1
       maxnp = ncharges
      elseif(depth.lt.0) then
       call bummer('fmm: error, depth = ',depth)
      endif
#endif
c
      call calnb(ncharges,nbytes,nbits,bitpos,mbitpos,depth,maxint,
     .maxmint)
c
#ifdef FMM_EXTREMECOMPRESSION
#ifdef FMM_EXTREMEEXTREMECOMPRESSION
      call scompression(nbytes,nbits,maxint,bitpos,mbitpos)
#endif
#endif
c
      call calint1mmint1(nbytes)
c
      if(nint1.ge.-1) then
       if(nint1.lt.nbytes) then
        addr_desc(1) = nint1+1
       else
        call bummer('fmm: (nint1-nbytes) = ',(nint1-nbytes))
       endif
      else
       call bummer('fmm: error, nint1 = ',nint1)
      endif
      if(nbib.gt.1) then
       if(nbib.le.nbits) then
        addr_desc(2) = nbits-nbib
       else
        call bummer('fmm: (nbib-nbits) = ',(nbib-nbits))
       endif
      else
       call bummer('fmm: error, nbib = ',nbib)
      endif
#else
      addr_desc(1) = nbytes
      addr_desc(2) = 0
#endif
c
#ifdef FMM_PARALLEL
      if(nnodes.gt.0) then
         i = nnodes-1
#ifdef FMM_COMPRESSION
         if(nintgb1.ge.0) then
            if(nintgb1.lt.nbytes) then
               call fmmallocate(gbinfo,0,nintgb1,5,6,0,i,j)
               if(j.ne.0) call bummer('fmm: error, j = ',j)
            else
               call bummer('fmm: (nintgb1-nbytes) = ',(nintgb1-nbytes))
            endif
         else
            call bummer('fmm: error, nintgb1 = ',nintgb1)
         endif
#else
         call fmmallocate(gbinfo,5,6,0,i,j)
         if(j.ne.0) call bummer('fmm: error, j = ',j)
#endif
      else
         call bummer('fmm: error, nnodes = ',nnodes)
      endif
#endif
c
#ifdef FMM_PARALLEL
      if(icharges.gt.0) then
#else
      if(ncharges.gt.0) then
#endif
#ifdef FMM_PARALLEL
#ifdef FMM_COMPRESSION
         if(twoint) then
            if(nint1.ge.0) then
               call fmmallocate(iboxsrt,0,nint1,icharge1,icharge2,i)
               if(i.ne.0) call bummer('fmm: error, i = ',i)
            else
               call bummer('fmm: error, nint1 = ',nint1)
            endif
         elseif(nint1.ne.-1) then
            call bummer('fmm: error, nint1 = ',nint1)
         endif
#else
         call fmmallocate(iboxsrt,micharge1,micharge2,i)
         if(i.ne.0) call bummer('fmm: error, i = ',i)
#endif
#else
#ifdef FMM_COMPRESSION
         if(twoint) then
            if(nint1.ge.0) then
               call fmmallocate(iboxsrt,0,nint1,1,ncharges,i)
               if(i.ne.0) call bummer('fmm: error, i = ',i)
            else
               call bummer('fmm: error, nint1 = ',nint1)
            endif
         elseif(nint1.ne.-1) then
            call bummer('fmm: error, nint1 = ',nint1)
         endif
#else
         call fmmallocate(iboxsrt,1,ncharges,i)
         if(i.ne.0) call bummer('fmm: error, i = ',i)
#endif
#endif
      else
#ifdef FMM_PARALLEL
         call bummer('fmm: error, icharges = ',icharges)
#else
         call bummer('fmm: error, ncharges = ',ncharges)
#endif
      endif
c
#ifdef FMM_COMPRESSION
#ifdef FMM_EXTREMECOMPRESSION
#ifdef FMM_EXTREMEEXTREMECOMPRESSION
#ifdef FMM_SIGNEXPONENT
      call s6tocoordinates(xyzt)
#else
      call s3tocoordinates(xyzt)
#endif
#endif
      call iboxsrttoibox(ibox,iboxsrt)
#endif
#else
      call setiboxsrt(iboxsrt)
#endif
c
      if(homogen.eq.1) then
         int3xyzd = depth
         int3xyzd1 = int3xyzd+1
         if(int3xyzd.ge.0) then
            int3xyz = bitpos(int3xyzd)-1
         else
            int3xyz = int3xyzd
         endif
      elseif(homogen.ne.0) then
         call bummer('fmm: error, homogen = ',homogen)
      endif
c
      if(int3xyzd.ge.0) then
         call fmmallocate(int3x,0,int3xyz,i)
         if(i.ne.0) call bummer('fmm: error, i = ',i)
         call fmmallocate(int3y,0,int3xyz,i)
         if(i.ne.0) call bummer('fmm: error, i = ',i)
         call fmmallocate(int3z,0,int3xyz,i)
         if(i.ne.0) call bummer('fmm: error, i = ',i)
         call fmmallocate(int3p,0,int3xyz,i)
         if(i.ne.0) call bummer('fmm: error, i = ',i)
      endif
c
      call getneighborsinit(nbits,bitpos,mbitpos)
      call getcjpinit(i)
      call fmmallocate(int3q,0,i,j)
      if(j.eq.0) then
         int3q(0) = i-1
         call getcjp(int3q)
      else
         call bummer('fmm: error, j = ',j)
      endif
c
      i = maxdepth+1
      call fmmallocate(nbofmb,1,i,j)
      if(j.eq.0) then
         do 338 j = 1,i
            nbofmb(j) = 0
 338     continue
      else
         call bummer('fmm: error, j = ',j)
      endif
c
      dem = de
c
      if(jmp) then
        if(ncharges.gt.0) then
          call fmmallocate(iboxjmp,1,ncharges,i)
          if(i.eq.0) then
            call iboxjmpz(ncharges,iboxjmp)
          else
            call caljmpl(jmp,jmpg,jmph,jmpp,jmpupd,jmpj,jmpn,jmpb,jmpjp)
          endif
        else
          call bummer('fmm: error, ncharges = ',ncharges)
        endif
      endif
c
#ifndef FMM_COMPRESSION
      call sccoorda(ws,nbits,bitpos,xyzin,xyzt,sf,periodic,periodica,
     .periodlength,changepos)
c
      if(periodic.gt.0) then
       if(nj.gt.0) then
        i = nbits-1
        call fmmallocate(bitposm1,0,i,j)
        if(j.ne.0) call bummer('fmm: error, j = ',j)
        call calbitposm1(nbits,maxint,bitpos)
        call fmmallocate(mskj,1,nbytes,i)
        if(i.ne.0) call bummer('fmm: error, i = ',i)
        call calmskj(nbytes)
        call calm4(nbytes,nbits,bitpos,periodic)
        call fmmallocate(folder,1,nfolder,i)
        if(i.ne.0) call bummer('fmm: error, i = ',i)
       elseif(nj.lt.0) then
        call bummer('fmm: error, nj = ',nj)
       endif
      elseif(periodic.ne.0) then
       call bummer('fmm: error, periodic = ',periodic)
      endif
c
      call sccoordb(ws,nbits,bitpos,mbitpos,bfg,bfglen,q,xyzt,sf,
     .periodic,linearpotential,lineardistance,changepos,negpos,
     .dipolecorrection,qqq,corrsx,corrsy,corrsz,corrs,corrsh)
#endif
c
      call scalecoordinates(ncharges,ws,buflen,maxnmultipoles,
     .nmultipoles,maxdepth,mmaxdepth,depth,nbytes,nbits,maxint,maxmint,
     .bitpos,mbitpos,q,pfmmpot,fmmgrad,fmmgrad(icharges+1),
     .fmmgrad(2*icharges+1),xyzt,ibox,iboxsrt,bfg,bfglen,ierr,de,inf,
     .ldf,ldff,dfmmmerr,maxwsd,ncsar,icsar,nfmmcos,fmmcos,fmmerr,
     .pfmmerr,merr,sf,fracdepth,sh,powsq,parabola,piboxscr,iboxjmp,
     .mmaxwsd,fmmdist,ipo,jpo,mask,ishx,ishy,mishx,mishy,maskxy,
     .shmonopole,enearfield,enfinbox,enfbibj,gb,gbsh,int3x,int3y,int3z,
     .int3p,int3q,withbop,hugep,hugef,gp,gsq,negpos,dipolecorrection,
     .qqq,corrsx,corrsy,corrsz,corrs,corrsh,periodic,periodica,
     .periodlength,pgd,nbofmb,sh4,changepos,shmp,linearpotential,
     .ilinearpotential,lineardistance,erroranalysis,homogen)
c
      if(nmultipoles.le.maxnmultipoles) then
#ifdef FMM_INFO
#ifdef FMM_PARALLEL
         if(me.eq.0) then
#endif
            write(6,*) ' parabola = ',parabola
            write(6,*) ' nmultipoles = ',nmultipoles
#ifdef FMM_PARALLEL
         endif
#endif
#endif
      else
         call bummer('fmm: (nmultipoles-maxnmultipoles) = ',
     .   (nmultipoles-maxnmultipoles))
      endif
c
      if(depth.ge.0) then
         i = depth+1
c
         call fmmallocate(nboxesinlevel,1,i,j)
         if(j.ne.0) call bummer('fmm: error, j = ',j)
         call fmmallocate(nboxeslevel,1,i,j)
         if(j.ne.0) call bummer('fmm: error, j = ',j)
      else
         call bummer('fmm: error, depth = ',depth)
      endif
c
      call sortcharges(ncharges,mmaxdepth,depth,q,pfmmpot,fmmgrad,
     .fmmgrad(icharges+1),fmmgrad(2*icharges+1),xyzt,sh,powsq,ibox,
     .piboxscr,iboxsrt,iboxjmp,jibfg,jibfglen,nboxesinlevel,nboxeslevel,
     .ws,maxnmultipoles,nmultipoles,maxint,maxmint,nbits,bitpos,mbitpos,
     .ishx,ishy,mishx,mishy,maskxy,fracdepth,sf,parabola,ierr,dem,de,
     .shmonopole,merr(0,ws),bfg,bfglen,enearfield,enfinbox,enfbibj,gb,
     .gbsh,int3x,int3y,int3z,int3p,int3q,withbop,hugep,hugef,periodic,
     .gp,gsq,negpos,corrsh,pgd,nbofmb,sh4,changepos,linearpotential,
     .ilinearpotential,lineardistance,erroranalysis,homogen)
c
      if(jmp) then
         call fmmdeallocate(iboxjmp,i)
         if(i.ne.0) call bummer('fmm: error, i = ',i)
      endif
c
      if(depth.lt.0) call bummer('fmm: error, depth = ',depth)
      if(depth.gt.maxdepth) call bummer('fmm: (depth-maxdepth) = ',
     .(depth-maxdepth))
c
      if(nmultipoles.ge.0) then
         if(nmultipoles.gt.0) then
            mnmultipoles = -nmultipoles
         else
            mnmultipoles = 0
         endif
c
         if(associatedwignerd.or.precomputeallds) then
            compute = .false.
#ifdef FMM_INFO
#ifdef FMM_PARALLEL
            if(me.eq.0) then
#endif
               write(6,*) ' rotation matrices are pre-computed.'
               write(6,*) ' associatedwignerd = ',associatedwignerd
               write(6,*) '   precomputeallds = ',precomputeallds
               write(6,*) '           compute = ',compute
#ifdef FMM_PARALLEL
            endif
#endif
#endif
         else
            compute = .true.
#ifdef FMM_INFO
#ifdef FMM_PARALLEL
            if(me.eq.0) then
#endif
               write(6,*) ' rotation matrices are not pre-computed.'
               write(6,*) ' associatedwignerd = ',associatedwignerd
               write(6,*) '   precomputeallds = ',precomputeallds
               write(6,*) '           compute = ',compute
#ifdef FMM_PARALLEL
            endif
#endif
#endif
         endif
c
         call fmmallocate(d2,0,nmultipoles,mnmultipoles,nmultipoles,0,
     .   nmultipoles,i)
         if(i.ne.0) call bummer('fmm: error, i = ',i)
         call fmmallocate(d3,0,nmultipoles,mnmultipoles,nmultipoles,0,
     .   nmultipoles,i)
         if(i.ne.0) call bummer('fmm: error, i = ',i)
         call fmmallocate(d2f,0,nmultipoles,mnmultipoles,nmultipoles,0,
     .   nmultipoles,i)
         if(i.ne.0) call bummer('fmm: error, i = ',i)
         call fmmallocate(d3f,0,nmultipoles,mnmultipoles,nmultipoles,0,
     .   nmultipoles,i)
         if(i.ne.0) call bummer('fmm: error, i = ',i)
c
         if(compute) then
            call fmmallocate(coeff1,0,nmultipoles,0,nmultipoles,i)
            if(i.ne.0) call bummer('fmm: error, i = ',i)
            call fmmallocate(coeff2,0,nmultipoles,0,nmultipoles,i)
            if(i.ne.0) call bummer('fmm: error, i = ',i)
            call fmmallocate(coeff3,mnmultipoles,nmultipoles,0,
     .      nmultipoles,0,nmultipoles,i)
            if(i.ne.0) call bummer('fmm: error, i = ',i)
            call fmmallocate(coeff4,mnmultipoles,nmultipoles,0,
     .      nmultipoles,i)
            if(i.ne.0) call bummer('fmm: error, i = ',i)
            call fmmallocate(coeff5,0,nmultipoles,0,nmultipoles,0,
     .      nmultipoles,i)
            if(i.ne.0) call bummer('fmm: error, i = ',i)
            call fmmallocate(coeff6,0,nmultipoles,0,nmultipoles,i)
            if(i.ne.0) call bummer('fmm: error, i = ',i)
c
            call coefficients(mnmultipoles,nmultipoles,coeff1,coeff2,
     .      coeff3,coeff4,coeff5,coeff6)
         endif
c
         if(depth.ge.0) then
            i = nmultipoles+1
            nsqmultipoles = iand(ishft((i*(i+1)),-1),maxint)
            i = depth+1
            i = nboxeslevel(i)+nboxesinlevel(i)
            if(i.gt.0) then
               j = nsqmultipoles*i
#ifdef FMM_TREETOGRAD
               if(homogen.eq.1) then
                  if((3*icharges).ge.(j+j)) then
                     allocomega = .false.
                     romegatree => fmmgrad(1:j)
                     iomegatree => fmmgrad((j+1):(j+j))
                  elseif((3*icharges).ge.j) then
                     call fmmallocate(omegatree,1,j,k)
                     if(k.eq.0) then
                        allocomega = .true.
                        romegatree => fmmgrad(1:j)
                        iomegatree => omegatree(1:j)
                     else
                        call bummer('fmm: error, k = ',k)
                     endif
                  else
                     call fmmallocate(omegatree,1,(j+j),k)
                     if(k.eq.0) then
                        allocomega = .true.
                        romegatree => omegatree(1:j)
                        iomegatree => omegatree((j+1):(j+j))
                     else
                        call bummer('fmm: error, k = ',k)
                     endif
                  endif
               else
                  call bummer('fmm: treetograd requested, homogen = ',
     .            homogen)
               endif
#else
               call fmmallocate(omegatree,1,(j+j),k)
               if(k.eq.0) then
                  romegatree => omegatree(1:j)
                  iomegatree => omegatree((j+1):(j+j))
               else
                  call bummer('fmm: error, k = ',k)
               endif
#endif
            else
               call bummer('fmm: error, i = ',i)
            endif
         else
            call bummer('fmm: error, depth = ',depth)
         endif
      else
         call bummer('fmm: error, nmultipoles = ',nmultipoles)
      endif
c
c      call coulomb(ncharges,q,xyzt,coul)
#ifdef FMM_INFO
#ifdef FMM_PARALLEL
      if(me.eq.0) then
#endif
c         write(6,'(1x,q50.40,1x,z32)') (sh*(coul/sf)),
c     .   (sh*(coul/sf))
#ifdef FMM_PARALLEL
      endif
#endif
#endif
c
      if(nmultipoles.le.maxnmultipoles) then
#ifdef FMM_INFO
#ifdef FMM_PARALLEL
         if(me.eq.0) then
#endif
            write(6,'('' nmultipoles = '',i13)') nmultipoles
#ifdef FMM_PARALLEL
         endif
#endif
#endif
      else
         call bummer('fmm: (nmultipoles-maxnmultipoles) = ',
     .   (nmultipoles-maxnmultipoles))
      endif
c
      n2multipoles = nmultipoles+nmultipoles
c
#ifdef FMM_INFO
#ifdef FMM_PARALLEL
      if(me.eq.0) then
#endif
c         write(6,*) sf,sh
#ifdef FMM_PARALLEL
      endif
#endif
#endif
c
      call calpow(nmultipoles,pow)
c
      call sgneg(n2multipoles,sg)
c
      if(.not.compute) call cpydtod13(mnmultipoles,nmultipoles,d2,d3,
     .d2f,d3f,wignerd,1)
c
      call pass1(ncharges,depth,ws,nbits,maxint,maxmint,bitpos,mbitpos,
     .maxnmultipoles,nmultipoles,mnmultipoles,nsqmultipoles,
     .nboxesinlevel,nboxeslevel,q,xyzt,ibox,piboxscr,ishx,ishy,mishx,
     .mishy,maskxy,buflen,rl,cmphi,smphi,fac,rfac,pow,sg,fr,coeff1,
     .coeff2,coeff3,coeff4,coeff5,coeff6,d2,d3,d2f,d3f,romegatree,
     .iomegatree,csar,car,sar,rar,isrt,kbxyzar,indar,kboxxyzar,
     .kboxindar,mi,jcar,hcar,hsar,sf,periodic,withaop,compute)
c
      if(nmultipoles.ge.0) then
        if(depth.ge.0) then
          i = depth+1
          i = nboxeslevel(i)+nboxesinlevel(i)
          if(i.gt.0) then
            jnbi = nboxesinlevel(depth+1)
#ifdef FMM_IBOXSCR
#ifdef FMM_TREETOGRAD
            if(homogen.eq.0) then
              doit = .true.
            elseif(homogen.eq.1) then
              doit = .false.
            else
              call bummer('fmm: error, homogen = ',homogen)
            endif
#else
            doit = .true.
#endif
            if(doit) then
              call fmmdeallocate(iboxscr,j)
              if(j.ne.0) call bummer('fmm: error, j = ',j)
              piboxscr => null()
#ifdef FMM_PARALLEL
              call fmmallocate(iboxscr,jcharge1,jcharge2,j)
              if(j.ne.0) call bummer('fmm: error, j = ',j)
#else
              call fmmallocate(iboxscr,1,jnbi,j)
              if(j.ne.0) call bummer('fmm: error, j = ',j)
#endif
              piboxscr => iboxscr
            endif
#endif
            j = nsqmultipoles*i
#ifdef FMM_TREETOGRAD
#ifdef FMM_EXTREMETREETOGRAD
            if(homogen.eq.1) then
              if((3*icharges).ge.(4*j)) then
                allocmu = .false.
                doallocbftreetograd = .true.
                ntreetograd = 2
                starttreetograd = 3*icharges-2*j+1
                itreetograd = 2*nsqmultipoles
                rmutree => fmmgrad((3*icharges-2*j+1):(3*icharges-j))
                imutree => fmmgrad((3*icharges-j+1):(3*icharges))
              elseif((3*icharges).ge.(3*j)) then
                call fmmallocate(mutree,1,j,k)
                if(k.eq.0) then
                  allocmu = .true.
                  doallocbftreetograd = .true.
                  ntreetograd = 1
                  starttreetograd = 3*icharges-j+1
                  itreetograd = nsqmultipoles
                  rmutree => fmmgrad((3*icharges-j+1):(3*icharges))
                  imutree => mutree(1:j)
                else
                  call bummer('fmm: error, k = ',k)
                endif
              else
                call fmmallocate(mutree,1,(j+j),k)
                if(k.eq.0) then
                  allocmu = .true.
                  ntreetograd = 0
                  rmutree => mutree(1:j)
                  imutree => mutree((j+1):(j+j))
                else
                  call bummer('fmm: error, k = ',k)
                endif
              endif
#ifdef FMM_INFO
#ifdef FMM_PARALLEL
              if(me.eq.0) then
#endif
                 write(6,*) ' ntreetograd = ',ntreetograd
#ifdef FMM_PARALLEL
              endif
#endif
#endif
            else
              call bummer('fmm: treetograd requested, homogen=',homogen)
            endif
#else
            call fmmallocate(mutree,1,(j+j),k)
            if(k.eq.0) then
              rmutree => mutree(1:j)
              imutree => mutree((j+1):(j+j))
            else
               call bummer('fmm: error, k = ',k)
            endif
#endif
#else
            call fmmallocate(mutree,1,(j+j),k)
            if(k.eq.0) then
              rmutree => mutree(1:j)
              imutree => mutree((j+1):(j+j))
            else
               call bummer('fmm: error, k = ',k)
            endif
#endif
            if(withtaylor) then
             call itaylor(gtaylor,bitpos,nbits,maxint,igtaylor,mgtaylor)
             call jptaylor(i,gtaylor,nbits,igtaylor,mgtaylor,ntaylor,j)
             if(ntaylor.gt.0) then
               if(ntaylor.le.i) then
                 if(j.ge.0) then
                   if(j.lt.nbits) then
                     j = nbits*(ntaylor-1)+(j+1)
                     if(j.eq.i) then
                       call fmmallocate(taylor,1,ntaylor,j)
                       if(j.ne.0) call bummer('fmm: error, j = ',j)
                     else
                       call bummer('fmm: (j-i) = ',(j-i))
                     endif
                   else
                     call bummer('fmm: (j-nbits) = ',(j-nbits))
                   endif
                 else
                   call bummer('fmm: error, j = ',j)
                 endif
               else
                 call bummer('fmm: (ntaylor-i) = ',(ntaylor-i))
               endif
             else
               call bummer('fmm: error, ntaylor = ',ntaylor)
             endif
            endif
          else
            call bummer('fmm: error, i = ',i)
          endif
        else
          call bummer('fmm: error, depth = ',depth)
        endif
      else
        call bummer('fmm: error, nmultipoles = ',nmultipoles)
      endif
c
      call fmmallocate(irar,1,maxwsd3,1,(depth+1),i)
      if(i.ne.0) call bummer('fmm: error, i = ',i)
c
      call fmmimn(ws,maxdepth,periodic,bitpos,ilevelmn)
c
      if((depth+1).ge.ilevelmn) then
         call fmmi(maxwsd3,ilevelmn,(depth+1),ws,bitpos,nrar,irar,
     .   periodic,.false.)
c
         i = maxdepth-ilevelmn+2
c
         if(i.gt.0) then
            if(nrar.le.maxnrar) then
               nrar = i*nrar
            else
               call bummer('fmm: (nrar-maxnrar) = ',(nrar-maxnrar))
            endif
         else
            call bummer('fmm: error, i = ',i)
         endif
      else
         nrar = 0
      endif
c
      nbfg = bfglen/3
      if(nbfg.le.0) call bummer('fmm: error, nbfg = ',nbfg)
c
      i = depth+1
      ntree = nboxeslevel(i)+nboxesinlevel(i)
      if(ntree.le.0) call bummer('fmm: error, ntree = ',ntree)
c
      call pass2(ncharges,jnbi,depth,ws,nbits,maxint,bitpos,mbitpos,
     .nmultipoles,mnmultipoles,n2multipoles,nsqmultipoles,nboxesinlevel,
     .nboxeslevel,q,xyzt,ibox,piboxscr,iboxsrt,nbfg,bfg,rl,cmphi,smphi,
     .cmphipi,smphipi,fac,rfac,pow,sg,fr,coeff1,coeff2,coeff3,coeff4,
     .coeff5,coeff6,ntree,romegatree,iomegatree,rmutree,imutree,
     .withtaylor,gtaylor,igtaylor,mgtaylor,ntaylor,taylor,buflen,flvlar,
     .csar,car,sar,rar,isrt,kbxyzar,indar,kboxxyzar,kboxindar,kbar,
     .maxwsd,mmaxwsd,ncsar,icsar,jcsar,ncar,icar,isar,maxwsd3,nrar,irar,
     .mmaxdepth,powsq,int3x,int3y,int3z,int3p,int3q,hugep,hugef,
     .periodic,ishx,maskx,ishy,masky,mishx,mishy,maskxy,maxmint,cachopt,
     .g2db,unrolled2,pgd,nbofmb,withbop,wignerd,compute)
c
      call fmmdeallocate(irar,i)
      if(i.ne.0) call bummer('fmm: error, i = ',i)
c
      if(jmp) call caljmpi(ncharges,iboxsrt,sjmp)
c
#ifndef FMM_NOPOT
      if(homogen.eq.0) then
         call fmmenpot(nbits,bitpos,bfg,bfglen,q,pfmmpot,enearfieldpot)
      elseif(homogen.ne.1) then
         call bummer('fmm: error, homogen = ',homogen)
      endif
#endif
c
      i = max((2*nmultipoles),(nmultipoles*nmultipoles))
c
      call fmmallocate(dbl,1,i,j)
c
      if(j.eq.0) then
         call caldbl(i,dbl)
      else
         call bummer('fmm: error, j = ',j)
      endif
c
      if(.not.compute) then
         if(ncsar.gt.0) then
            i = ncsar+2
            call cpydtod13(mnmultipoles,nmultipoles,d2,d3,d2f,d3f,
     .      wignerd,i)
         else
            call bummer('fmm: error, ncsar = ',ncsar)
         endif
      endif
c
      call pass3(depth,ws,nbits,nboxesinlevel(depth+1),bitpos,mbitpos,
     .maxnmultipoles,nmultipoles,mnmultipoles,nsqmultipoles,
     .nboxesinlevel,nboxeslevel,q,xyzt,ibox,bfg,bfglen,mask,buflen,
     .cmphipi,smphipi,rl,cmphi,smphi,fac,rfac,pow,sg,fr,coeff1,coeff2,
     .coeff3,coeff4,coeff5,coeff6,d2,d3,d2f,d3f,rmutree,imutree,
     .withtaylor,gtaylor,igtaylor,mgtaylor,taylor,csar,car,sar,rar,isrt,
     .kbxyzar,indar,kboxxyzar,kboxindar,mi,jcar,hcar,hsar,mmaxdepth,
     .mishx,mishy,maskxy,ipo,ishx,ishy,maxint,maxmint,powsq,periodic,
     .romegatree,iomegatree,efarfield,efarfieldpot,e1per,pfmmpot,
     .fmmgrad,dbl,sh3,withcop,compute)
c
      call fmmdeallocate(dbl,i)
      if(i.ne.0) call bummer('fmm: error, i = ',i)
c
#ifdef FMM_TREETOGRAD
      if(homogen.eq.1) then
#ifdef FMM_IBOXSCR
         call iboxscinfo5(ishx,ishy,mmaxdepth,depth,bitpos,powsq,xyzt,
     .   ibox,piboxscr)
#else
         call iboxinfo5(mmaxdepth,depth,bitpos,powsq,xyzt,ibox)
#endif
#ifdef FMM_PARALLEL
#ifdef FMM_COMPRESSION
         startbox = iand(ibox(icharge1),ibm)-1
         endbox = iand(ibox(icharge2),ibm)-1
#else
         startbox = ibox(icharge1)-1
         endbox = ibox(icharge2)-1
#endif
#else
#ifdef FMM_COMPRESSION
         startbox = iand((iand(ibox(1),ibm)-1),-8)
         endbox = ior((iand(ibox(ncharges),ibm)-1),7)
#else
         startbox = iand((ibox(1)-1),-8)
         endbox = ior((ibox(ncharges)-1),7)
#endif
#endif
         call skipeevector(icharges,ibox(icharge1))
#ifdef FMM_PARALLEL
         i = iand(startbox,-8)
         j = ior(endbox,7)
#else
         i = startbox
         j = endbox
#endif
         k = icharge1
         l = icharge2
         call pass5(ncharges,depth,ws,nbits,ishx,ishy,maxint,mishx,
     .   mishy,maskxy,bitpos,mbitpos,q,xyzt,ibox,piboxscr,iboxsrt,bfg,
     .   bfglen,enearfield,enfinbox,enfbibj,fmmgrad,pfmmpot,gb,gbsh,
     .   int3x,int3y,int3z,int3p,int3q,pagejump,pageshift,pageshiftg,
     .   pagemask,pageaddr,indsize,pagepossize,startbox,endbox,i,j,ibox,
     .   k,l,piboxscr,pages,pgd,.false.,periodic,indskpjump,nbofmb,sf,
     .   sh,linearpotential,ilinearpotential,lineardistance)
      else
         call bummer('fmm: error, homogen = ',homogen)
      endif
#endif
c
#ifdef FMM_IBOXSCR
      call fmmdeallocate(iboxscr,i)
      if(i.ne.0) call bummer('fmm: error, i = ',i)
      piboxscr => null()
#endif
      if(int3xyzd.ge.0) then
         call fmmdeallocate(int3x,i)
         if(i.ne.0) call bummer('fmm: error, i = ',i)
         call fmmdeallocate(int3y,i)
         if(i.ne.0) call bummer('fmm: error, i = ',i)
         call fmmdeallocate(int3z,i)
         if(i.ne.0) call bummer('fmm: error, i = ',i)
         call fmmdeallocate(int3p,i)
         if(i.ne.0) call bummer('fmm: error, i = ',i)
      endif
      call fmmdeallocate(int3q,i)
      if(i.ne.0) call bummer('fmm: error, i = ',i)
      call fmmdeallocate(nbofmb,i)
      if(i.ne.0) call bummer('fmm: error, i = ',i)
c
      call fmmdeallocate(nboxesinlevel,i)
      if(i.ne.0) call bummer('fmm: error, i = ',i)
      call fmmdeallocate(nboxeslevel,i)
      if(i.ne.0) call bummer('fmm: error, i = ',i)
      call fmmdeallocate(d2,i)
      if(i.ne.0) call bummer('fmm: error, i = ',i)
      call fmmdeallocate(d3,i)
      if(i.ne.0) call bummer('fmm: error, i = ',i)
      call fmmdeallocate(d2f,i)
      if(i.ne.0) call bummer('fmm: error, i = ',i)
      call fmmdeallocate(d3f,i)
      if(i.ne.0) call bummer('fmm: error, i = ',i)
      if(compute) then
         call fmmdeallocate(coeff1,i)
         if(i.ne.0) call bummer('fmm: error, i = ',i)
         call fmmdeallocate(coeff2,i)
         if(i.ne.0) call bummer('fmm: error, i = ',i)
         call fmmdeallocate(coeff3,i)
         if(i.ne.0) call bummer('fmm: error, i = ',i)
         call fmmdeallocate(coeff4,i)
         if(i.ne.0) call bummer('fmm: error, i = ',i)
         call fmmdeallocate(coeff5,i)
         if(i.ne.0) call bummer('fmm: error, i = ',i)
         call fmmdeallocate(coeff6,i)
         if(i.ne.0) call bummer('fmm: error, i = ',i)
      endif
#ifdef FMM_TREETOGRAD
      if(allocomega) then
        call fmmdeallocate(omegatree,i)
        if(i.ne.0) call bummer('fmm: error, i = ',i)
      endif
#else
      call fmmdeallocate(omegatree,i)
      if(i.ne.0) call bummer('fmm: error, i = ',i)
#endif
#ifdef FMM_TREETOGRAD
#ifdef FMM_EXTREMETREETOGRAD
      if(allocmu) then
        call fmmdeallocate(mutree,i)
        if(i.ne.0) call bummer('fmm: error, i = ',i)
      endif
#else
      call fmmdeallocate(mutree,i)
      if(i.ne.0) call bummer('fmm: error, i = ',i)
#endif
#else
      call fmmdeallocate(mutree,i)
      if(i.ne.0) call bummer('fmm: error, i = ',i)
#endif
      if(withtaylor) then
        call fmmdeallocate(taylor,i)
        if(i.ne.0) call bummer('fmm: error, i = ',i)
      endif
c
#ifdef FMM_COMPRESSION
      if(twoint) then
        call getibsrt2(nbytes,ibox,iboxsrt)
      else
        call getibsrt1(ibox)
      endif
c
      call fmmdeallocate(int1m,i)
      if(i.ne.0) call bummer('fmm: error, i = ',i)
c
      call fmmdeallocate(mint1,i)
      if(i.ne.0) call bummer('fmm: error, i = ',i)
c
      call fmmdeallocate(int1scr,i)
      if(i.ne.0) call bummer('fmm: error, i = ',i)
c
#ifdef FMM_EXTREMECOMPRESSION
      call fmmdeallocate(pinb,i)
      if(i.ne.0) call bummer('fmm: error, i = ',i)
#endif
c
      if(twoint) then
        if(nint1.ge.0) then
          call fmmdeallocate(iboxsrt,i)
          if(i.ne.0) call bummer('fmm: error, i = ',i)
        else
          call bummer('fmm: error, nint1 = ',nint1)
        endif
      elseif(nint1.ne.-1) then
        call bummer('fmm: error, nint1 = ',nint1)
      endif
c
      piboxsrt => ibox
#else
      call fmmdeallocate(ibox,i)
      if(i.ne.0) call bummer('fmm: error, i = ',i)
c
      piboxsrt => iboxsrt
#endif
c
#ifdef FMM_COMPRESSION
#ifdef FMM_EXTREMECOMPRESSION
#ifdef FMM_EXTREMEEXTREMECOMPRESSION
      call setcoordinates(xyzt)
#endif
#endif
#endif
c
      if(periodic.eq.3) then
#ifdef FMM_CORRECTION_OF_FORCES
       call calcxcycz(nbits,bitpos,bfg,bfglen,q,fmmgrad,cx,cy,cz)
#else
       cx = zero
       cy = zero
       cz = zero
#endif
      endif
c
      if(periodic.gt.0) then
       if(sh.eq.one) then
        if((dipolecorrection.ge.-1).and.(dipolecorrection.le.1)) then
         if(((dipolecorrection.eq.0).and.(negpos.eq.0)).or.
     .    (dipolecorrection.eq.1)) then
          if(periodic.eq.3) then
           call corrp(periodic,cx,cy,cz,q,xyzt,qqq,corrsx,corrsy,corrsz,
     .     corrs,efarfield,efarfieldpot,e1per,pfmmpot,fmmgrad)
          endif
         endif
        else
         call bummer('fmm: error, dipolecorrection = ',dipolecorrection)
        endif
       else
        call bummer('fmm: error, ncharges = ',ncharges)
       endif
      endif
c
      if(periodic.eq.3) then
       if(.not.plummer_potential) then
        if(.not.linearpotential) then
         call calvirialtensor3p(efarfield,enearfield,gp,virialtensor)
        endif
       endif
      elseif((.not.plummer_potential).and.(.not.linearpotential)) then
       call setbfgn(9,bfglen,bfgnlen)
c
       if(bfglen.eq.bfgnlen) then
        call calvirialtensor(nbits,bitpos,bfg,bfglen,q,xyzt,fmmgrad,cx,
     .  cy,cz,periodic,efarfield,enearfield,gp,virialtensor)
       elseif(bfglen.lt.bfgnlen) then
        call fmmallocate(bfgn,1,bfgnlen,i)
        if(i.ne.0) call bummer('fmm: error, i = ',i)
c
        call calvirialtensor(nbits,bitpos,bfgn,bfgnlen,q,xyzt,fmmgrad,
     .  cx,cy,cz,periodic,efarfield,enearfield,gp,virialtensor)
c
        call fmmdeallocate(bfgn,i)
        if(i.ne.0) call bummer('fmm: error, i = ',i)
       else
        call bummer('fmm: (bfglen-bfgnlen) = ',(bfglen-bfgnlen))
       endif
      else
#ifdef FMM_CORRECTION_OF_FORCES
       call calcxcycz(nbits,bitpos,bfg,bfglen,q,fmmgrad,cx,cy,cz)
#else
       cx = zero
       cy = zero
       cz = zero
#endif
      endif
c
      call sortback(ncharges,copyxyz,xyzt,piboxsrt,q,fmmgrad,pfmmpot)
c
#ifdef FMM_RESTORE_COORDINATES
      if(.not.copyxyz) call restorecoordinates(periodic,periodica,nbits,
     .mbitpos,xyzt)
#endif
c
#ifdef FMM_COMPRESSION
      call fmmdeallocate(ibox,i)
      if(i.eq.0) then
        piboxsrt => null()
      else
        call bummer('fmm: error, i = ',i)
      endif
#else
      call fmmdeallocate(iboxsrt,i)
      if(i.eq.0) then
        piboxsrt => null()
      else
        call bummer('fmm: error, i = ',i)
      endif
#endif
c
      if(copyxyz) then
        call fmmdeallocate(xyz,i)
        if(i.ne.0) call bummer('fmm: error, i = ',i)
      endif
c
      xyzt => null()
c
      if(periodic.gt.0) then
        if(nj.gt.0) then
          call fmmdeallocate(bitposm1,i)
          if(i.ne.0) call bummer('fmm: error, i = ',i)
          call fmmdeallocate(mskj,i)
          if(i.ne.0) call bummer('fmm: error, i = ',i)
          call fmmdeallocate(folder,i)
          if(i.ne.0) call bummer('fmm: error, i = ',i)
        elseif(nj.lt.0) then
          call bummer('fmm: error, nj = ',nj)
        endif
      elseif(periodic.ne.0) then
        call bummer('fmm: error, periodic = ',periodic)
      endif
c
      call scalefmmgr(ncharges,periodic,periodica,dipolecorrection,
     .negpos,cx,cy,cz,q,gsq,fmmgrad,gp,pfmmpot,virialtensor)
c
      if(.not.plummer_potential) then
        if(.not.linearpotential) call prtvirialtensor(virialtensor)
      endif
c
#ifndef FMM_NOPOT
      call fmmenpot(nbits,bitpos,bfg,bfglen,q,pfmmpot,energypot)
#endif
c
#ifdef FMM_DEBUG
      call fmmsgrad(nbits,bitpos,bfg,bfglen,fmmgrad,sgradx,sgrady,
     .sgradz)
c
      call fmmsagrad(nbits,bitpos,bfg,bfglen,fmmgrad,sagradx,sagrady,
     .sagradz)
#endif
c
      call fmmdeallocate(flvlar,i)
      if(i.ne.0) call bummer('fmm: error, i = ',i)
      call fmmdeallocate(powsq,i)
      if(i.ne.0) call bummer('fmm: error, i = ',i)
c
      energy = (sh*((half*efarfield)/sf))+(sh*(enearfield/sf))
c
#ifdef FMM_INFO
#ifdef FMM_PARALLEL
      if(me.eq.0) then
#endif
         write(6,*) ' efarfield     = ',(sh*((half*efarfield)/sf))
         write(6,*) ' efarfieldpot  = ',(sh*((half*efarfieldpot)/sf))
         write(6,*) ' enearfield    = ',(sh*(enearfield/sf))
#ifndef FMM_NOPOT
         if(homogen.eq.0) then
            write(6,*) ' enearfieldpot = ',(sh*(enearfieldpot/sf))
         elseif(homogen.ne.1) then
            call bummer('fmm: error, homogen = ',homogen)
         endif
#endif
         write(6,*) ' energy        = ',energy
#ifndef FMM_NOPOT
         write(6,*) ' energypot     = ',energypot
#endif
#ifdef FMM_PARALLEL
      endif
#endif
#endif
c
#ifndef FMM_UNIFORMGRID
      call calder(ierr,dem,energy,energypot,der)
#endif
c
      if(periodic.gt.0) then
#ifdef FMM_INFO
#ifdef FMM_PARALLEL
         if(me.eq.0) then
#endif
            write(6,*) ' energy_lattice     = ',(sh*(e1per/sf))
            write(6,*) ' energy_first_layer = ',(energy-(sh*(e1per/sf)))
#ifdef FMM_PARALLEL
         endif
#endif
#endif
      endif
#ifdef FMM_INFO
#ifdef FMM_PARALLEL
      if(me.eq.0) then
#endif
         write(6,*) ' efarfield          = ',(sh*((efarfield/two)/sf))
         write(6,*) ' enearfield         = ',(sh*(enearfield/sf))
         write(6,*) ' enfinbox           = ',(sh*(enfinbox/sf))
         write(6,*) ' enfbibj            = ',(sh*(enfbibj/sf))
         write(6,*) ' energy             = ',energy
#ifdef FMM_PARALLEL
      endif
#endif
#endif
c
      call fmmdeallocate(bitpos,i)
      if(i.ne.0) call bummer('fmm: error, i = ',i)
      call fmmdeallocate(mbitpos,i)
      if(i.ne.0) call bummer('fmm: error, i = ',i)
c
      if(precomputeallds) then
         if(.not.associatedwignerd) then
            call fmmdeallocatept(wignerd%wignerd,i)
            if(i.ne.0) call bummer('fmm: error, i = ',i)
         endif
      endif
c
      if(serroranalysis.lt.nerroranalysis) then
        call bummer('fmm: (serroranalysis-nerroranalysis) = ',
     .  (serroranalysis-nerroranalysis))
      elseif(snewerroranalysis.le.0) then
        call bummer('fmm: error, snewerroranalysis =',snewerroranalysis)
      else
        if(snewerroranalysis.le.nerroranalysis) then
          doit = .true.
        elseif(iplummerpotential.ne.iplummerpotentialsv) then
          doit = .true.
        elseif(aoplummer.ne.aoplummersv) then
          doit = .true.
        elseif(ilinearpotential.ne.ilinearpotentialsv) then
          doit = .true.
        elseif(ilinearpotential.eq.0) then
          doit = .false.
        elseif(ilinearpotential.eq.1) then
          if(linearodistance(3).ne.linearodistancesv(3)) then
            doit = .true.
          else
            doit = .false.
          endif
        elseif(ilinearpotential.eq.2) then
          if(linearodistance(2).ne.linearodistancesv(2)) then
            doit = .true.
          else
            doit = .false.
          endif
        elseif(ilinearpotential.eq.3) then
          if(linearodistance(2).ne.linearodistancesv(2)) then
            doit = .true.
          elseif(linearodistance(3).ne.linearodistancesv(3)) then
            doit = .true.
          else
            doit = .false.
          endif
        elseif(ilinearpotential.eq.4) then
          if(linearodistance(1).ne.linearodistancesv(1)) then
            doit = .true.
          else
            doit = .false.
          endif
        elseif(ilinearpotential.eq.5) then
          if(linearodistance(1).ne.linearodistancesv(1)) then
            doit = .true.
          elseif(linearodistance(3).ne.linearodistancesv(3)) then
            doit = .true.
          else
            doit = .false.
          endif
        elseif(ilinearpotential.eq.6) then
          if(linearodistance(1).ne.linearodistancesv(1)) then
            doit = .true.
          elseif(linearodistance(2).ne.linearodistancesv(2)) then
            doit = .true.
          else
            doit = .false.
          endif
        elseif(ilinearpotential.eq.7) then
          if(linearodistance(1).ne.linearodistancesv(1)) then
            doit = .true.
          elseif(linearodistance(2).ne.linearodistancesv(2)) then
            doit = .true.
          elseif(linearodistance(3).ne.linearodistancesv(3)) then
            doit = .true.
          else
            doit = .false.
          endif
        else
          call bummer('fmm: error, ilinearpotential =',ilinearpotential)
        endif
        if(doit) then
          ilinearpotentialsv = ilinearpotential
          linearodistancesv(1) = linearodistance(1)
          linearodistancesv(2) = linearodistance(2)
          linearodistancesv(3) = linearodistance(3)
          iplummerpotentialsv = iplummerpotential
          aoplummersv = aoplummer
c          call fmmdeallocate(hugep,i)
c          if(i.ne.0) call bummer('fmm: error, i = ',i)
c          call fmmdeallocate(hugef,i)
c          if(i.ne.0) call bummer('fmm: error, i = ',i)
          i = ltob*(maxdepth+2)+rtob*(maxdepth+1)
          if(nallocr.eq.i) then
            nallocr = 0
          else
            call bummer('fmm: (nallocr-i) = ',(nallocr-i))
          endif
          nallocst = nallocst-i
          serroranalysis = snewerroranalysis
          nerroranalysis = serroranalysis
        else
          serroranalysis = snewerroranalysis
        endif
      endif
c
#ifdef FMM_INFO
#ifdef FMM_PARALLEL
      if(me.eq.0) then
#endif
         write(6,'(''   gradx = '',d26.18)') fmmgrad(1)
         write(6,'(''   grady = '',d26.18)') fmmgrad(2)
         write(6,'(''   gradz = '',d26.18)') fmmgrad(3)
#ifndef FMM_NOPOT
         write(6,'(''     pot = '',d26.18)') pfmmpot(1)
#endif
#ifdef FMM_DEBUG
         write(6,'(''  sgradx = '',d26.18)') sgradx
         write(6,'(''  sgrady = '',d26.18)') sgrady
         write(6,'(''  sgradz = '',d26.18)') sgradz
         write(6,'('' sagradx = '',d26.18)') sagradx
         write(6,'('' sagrady = '',d26.18)') sagrady
         write(6,'('' sagradz = '',d26.18)') sagradz
#endif
         write(6,'(''  energy = '',d26.18)') energy
#ifdef FMM_PARALLEL
      endif
#endif
#endif
c
#ifdef FMM_PARALLEL
      call fmmdeallocate(gbinfo,i)
      if(i.ne.0) call bummer('fmm: error, i = ',i)
#endif
      pfmmpot => null()
#ifdef FMM_PARALLEL
#ifdef FMM_LOADSORT
      if(doload) iboxload => null()
#endif
#endif
      call edmdfmmalloc(nalloc,nallocst,'cfmm_tune')
      call prtmdfmmalloc(nalloc,maxnalloc,'  end of cfmm_tune')
c
#ifdef FMM_PARALLEL
#if defined(FMM_CPUTIME) || defined(FMM_WALLTIME)
      call c_cputime_walltime(.false.)
#endif
#endif
c
#ifdef FMM_STATISTICS
      call fmmstatistics(ncharges,.false.)
#endif

c-ik     
      FMM_internal_params%serroranalysis = serroranalysis
      FMM_internal_params%nerroranalysis = nerroranalysis
      FMM_internal_params%pgd = pgd
      FMM_internal_params%ws = ws
      FMM_internal_params%depth = depth
      FMM_internal_params%nmultipoles = nmultipoles
      FMM_internal_params%ncsar = ncsar
      FMM_internal_params%parabola = parabola
      FMM_internal_params%ilinearpotentialsv = ilinearpotentialsv
      FMM_internal_params%iplummerpotentialsv = iplummerpotentialsv
      FMM_internal_params%firsterroranalysis = firsterroranalysis
      FMM_internal_params%hugep = hugep
      FMM_internal_params%fracdepth = fracdepth
      FMM_internal_params%shmonopole = shmonopole
      FMM_internal_params%linearodistancesv = linearodistancesv
      FMM_internal_params%aoplummersv = aoplummersv
      FMM_internal_params%hugef = hugef
c-ik
      return
      end subroutine cfmm_tune
