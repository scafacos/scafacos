#include "fmm.h"
#include "fmmkinds.h"

      subroutine fmm_tune(ncharges,localcharges,q,xyzin,ierr,de,der,
     .energy,
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
      integer(kind=fmm_integer) ncharges,localcharges
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
c-ik
      real(kind=fmm_real), allocatable:: qtune(:),xyztune(:),
     .fmmpottune(:),fmmgradtune(:,:)
c
#ifdef FMM_PARALLEL
      integer(kind=fmm_integer), allocatable:: scr(:)
#endif
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
      icharge1 = 1
      icharge2 = localcharges
      icharges = localcharges

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

      if(icharge1.le.0) call bummer('fmm_tune: error, icharge1 = ',
     .icharge1)
      if(icharge2.lt.icharge1)
     .call bummer('fmm_tune: (icharge2-icharge1) =',(icharge2-icharge1))
      i = icharge2-icharge1+1
      if(i.ne.icharges) call bummer('fmm_tune: (i-icharges) = ',
     .(i-icharges))


      allocate(qtune(1:ncharges),stat = j)
      if(j.ne.0) call bummer('fmm_tune: error, j = ',j)
c
      allocate(xyztune(3*ncharges),stat = j)
      if(j.ne.0) call bummer('fmm_tune: error, j = ',j)
c
      qtune = 0.e0_fmm_real
      xyztune = 0.e0_fmm_real
c
#ifdef FMM_PARALLEL
      allocate(scr(0:(nnodes-1)),stat = i)
      if(i.ne.0) call bummer('fmm_tune: error, i = ',i)
c
      if(icharges.gt.0) then
         scr(me) = icharges
      else
         call bummer('fmm_tune: error, icharges = ',icharges)
      endif
c
      call mp_allgather(scr,1,MP_ALLNODES)
c
      i = 1
c
      if(me.gt.0) then
         do 1 j = 0,(me-1)
            i = i+scr(j)
 1       continue
      endif
c
      qtune(i:(i+icharges-1)) = q(icharge1:icharge2)
      call fmm_copy_tune((3*icharges),xyzin(1,icharge1),
     .xyztune((3*(i-1)+1)))
c
      call mp_allreduce(qtune,ncharges,MP_SUM,MP_ALLNODES)
      call mp_allreduce(xyztune,(3*ncharges),MP_SUM,MP_ALLNODES)
c
      deallocate(scr,stat = i)
      if(i.ne.0) call bummer('fmm_tune: error, i = ',i)
#else
       qtune(icharge1:icharge2) = q(icharge1:icharge2)
      call fmm_copy_tune((3*ncharges),xyzin(1,1),
     .xyztune(1))
#endif

      allocate(fmmpottune(1:ncharges),stat = i)
      if(i.ne.0) call bummer('fmm_tune: error, i = ',i)
c
      allocate(fmmgradtune(1:3,1:ncharges),stat = i)
      if(i.ne.0) call bummer('fmm_tune: error, i = ',i)
c
#ifdef FMM_PARALLEL
      MP_ALLNODES = MPI_COMM_SELF
c
      call mp_rank(MP_ALLNODES,me)
      call mp_nnodes(MP_ALLNODES,nnodes)
c
#ifdef FMM_INFO
      if(me.eq.0) write(6,*) ' nnodes: ',nnodes
#endif
c
      if(nnodes.le.0) call bummer('fmm_tune: error, nnodes = ',nnodes)
c
      if(me.lt.0) call bummer('fmm_tune: error, me = ',me)
c
      if(me.ge.nnodes) call bummer('fmm_tune: (me-nnodes)=',(me-nnodes))
#endif

      call cfmm_tune(ncharges,qtune,xyztune,ierr,
     .de,der,energy,fmmpottune,
     .fmmgradtune,periodic,periodica,periodlength,dipolecorrection,
     .ilinearpotential,linearodistance,iplummerpotential,aoplummer,
     .snewerroranalysis,homogen,maxdepth,unrolled2,balance_load,
     .FMM_internal_params)

c      call cfmm_tune(ncharges,qtune,xyztune,ierr,
c     .de,der,energy,fmmpottune,
c     .fmmgradtune,periodic,periodica,periodlength,dipolecorrection,
c     .ilinearpotential,lineardistance,iplummerpotential,aoplummer,
c     .snewerroranalysis,homogen,FMM_internal_params)

#ifdef FMM_PARALLEL
      MP_ALLNODES = MPI_COMM_WORLD
c
      call mp_rank(MP_ALLNODES,me)
      call mp_nnodes(MP_ALLNODES,nnodes)
c
#ifdef FMM_INFO
      if(me.eq.0) write(6,*) ' nnodes: ',nnodes
#endif
c
      if(nnodes.le.0) call bummer('fmm_tune: error, nnodes = ',nnodes)
c
      if(me.lt.0) call bummer('fmm_tune: error, me = ',me)
c
      if(me.ge.nnodes) call bummer('fmm_tune: (me-nnodes)=',(me-nnodes))
c
#endif

      deallocate(qtune,stat = i)
      if(i.ne.0) call bummer('fmm_tune: error, i = ',i)
      deallocate(xyztune,stat = i)
      if(i.ne.0) call bummer('fmm_tune: error, i = ',i)
      deallocate(fmmpottune,stat = i)
      if(i.ne.0) call bummer('fmm_tune: error, i = ',i)
      deallocate(fmmgradtune,stat = i)
      if(i.ne.0) call bummer('fmm_tune: error, i = ',i)

      return
      end subroutine fmm_tune

      subroutine fmm_copy_tune(n,xyz,xyztune)
c
      use fmmkinds
c
      implicit none
c
      real(kind=fmm_real) xyz(*),xyztune(*)
c
      integer(kind=fmm_integer) n
c
      if(n.gt.0) then
         xyztune(1:n) = xyz(1:n)
      else
         call bummer('fmm_copy_tune: error, n = ',n)
      endif
c
      return
      end subroutine fmm_copy_tune
