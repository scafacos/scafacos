#include "fmm.h"
module fmm_cbindings

  use fmm_fcs_binding
  use fmmkinds
  use iso_c_binding

  contains

    ! init subroutine for the C interface
    subroutine fmm_cinit(cptr) bind(c)
      use mp_wrapper, only : mp_init
      implicit none
      type(c_ptr), value :: cptr
      type(FMM_internal_params_t), pointer :: FMM_internal_params

      call c_f_pointer(cptr,FMM_internal_params)

      FMM_internal_params%firsterroranalysis = .true.

      FMM_internal_params%depth = 0

      FMM_internal_params%nerroranalysis = huge(1_8)
      FMM_internal_params%wignerd%wignerd => NULL()

      FMM_internal_params%presorted = 0
      FMM_internal_params%resort = 0
      FMM_internal_params%resort_ptr = c_null_ptr

      call mp_init()
    end subroutine fmm_cinit

    ! init loadvector subroutine for the C interface
    subroutine fmm_cinitload(cptr,loadptr,iboxloadlength) bind(c)
      implicit none
      integer(kind=fmm_integer), value :: iboxloadlength
      type(c_ptr), value :: cptr,loadptr
      type(FMM_internal_params_t), pointer :: FMM_internal_params

      call c_f_pointer(cptr,FMM_internal_params)

      FMM_internal_params%iboxloadcptr = loadptr
      FMM_internal_params%iboxloadlength = iboxloadlength

    end subroutine fmm_cinitload

    ! set loadvector subroutine for the C interface
    subroutine fmm_csetload(cptr,val) bind(c)
      implicit none
      integer(kind=fmm_integer) :: iboxloadlength
      real(kind=fmm_real), dimension(:), pointer :: iboxload
      real(kind=fmm_real), value :: val
      type(c_ptr), value :: cptr
      type(FMM_internal_params_t), pointer :: FMM_internal_params

      call c_f_pointer(cptr,FMM_internal_params)

      iboxloadlength = FMM_internal_params%iboxloadlength 
      call c_f_pointer(FMM_internal_params%iboxloadcptr,iboxload,[iboxloadlength])

      iboxload = val

    end subroutine fmm_csetload

    ! set presorted mode
    subroutine fmm_csetpresorted(cptr,presorted) bind(c)
      implicit none
      type(c_ptr), value :: cptr
      integer(kind=c_long_long), value :: presorted
      type(FMM_internal_params_t), pointer :: FMM_internal_params

      call c_f_pointer(cptr,FMM_internal_params)

      FMM_internal_params%presorted = presorted

    end subroutine fmm_csetpresorted

    ! init resort support
    subroutine fmm_cinitresort(cptr,resort_ptr) bind(c)
      implicit none
      type(c_ptr), value :: cptr,resort_ptr
      type(FMM_internal_params_t), pointer :: FMM_internal_params

      call c_f_pointer(cptr,FMM_internal_params)

      FMM_internal_params%resort_ptr = resort_ptr

    end subroutine fmm_cinitresort

    ! set resort support
    subroutine fmm_csetresort(cptr,resort) bind(c)
      implicit none
      type(c_ptr), value :: cptr
      integer(kind=c_long_long), value :: resort
      type(FMM_internal_params_t), pointer :: FMM_internal_params

      call c_f_pointer(cptr,FMM_internal_params)

      FMM_internal_params%resort = resort

    end subroutine fmm_csetresort


    ! tune subroutine for the C interface
    subroutine fmm_ctune(local_particles,&
                         positions,&
                         charges,&
                         total_particles,&
                         absrel,&
                         deltaE,&
                         dip_corr,&
                         periodicity,&
                         period_length,&
                         maxdepth,&
                         unrolled_limit,&
                         balance_load,&           
                         params,&
                         wignersize,&
                         result) bind(c)
    implicit none

        integer(kind = c_long_long), value :: local_particles
        integer(kind = c_long_long), value :: total_particles
        integer(kind = c_long_long), value :: absrel
        integer(kind = c_long_long), value :: dip_corr
        integer(kind = c_long_long), dimension(3) :: periodicity
        integer(kind = c_long_long), value :: maxdepth
        integer(kind = c_long_long), value :: unrolled_limit
        integer(kind = c_long_long), value :: balance_load

      	real(kind = c_double), value :: deltaE
		real(kind = c_double), value :: period_length
      	real(kind = c_double), dimension(3*local_particles) :: positions
      	real(kind = c_double), dimension(local_particles) :: charges

      	type(c_ptr), value :: params
      	type(FMM_internal_params_t), pointer :: FMM_internal_params

		integer(kind = 8) :: periodic, periodica,ws,nmultipoles
      	integer(kind = c_long_long) :: wignersize
		integer(kind = c_long_long)	::	result
        integer(kind = fmm_integer) :: i
        integer(kind = fmm_integer) :: poles, ncsar
        integer(kind = 8) :: fmm_temp_0 = 0, fmm_temp_1 = 1
      	real(kind = fmm_real) :: der = 0.0d0
      	real(kind = fmm_real) :: lineardistance = 0.0d0
      	real(kind = fmm_real) :: aplummer = 0.0d0
      	real(kind = fmm_real) :: energy = 0.0d0

      	call c_f_pointer(params,FMM_internal_params)

		periodic = COUNT(periodicity .eq. 1)
        select case(periodic)
          case (0)
			periodica = 0
		  case (1)
			if (periodicity(3) == 1) then
			  periodica = 0
		    else if (periodicity(2) == 1) then
			  periodica = 2
			else
			  periodica = 1
			end if
		  case (2)
			if (periodicity(1) == 1 .and. periodicity(2) == 1) then
			  periodica = 0
			else if (periodicity(1) == 1 .and. periodicity(3) == 1) then
			  periodica = 2
			else
			  periodica = 3
			end if
		  case (3)
			periodica = 0
		  case default
			result = 0
		  return
	  	end select
        
		call fmm_tune(total_particles,&
		              local_particles,&
		              charges,&
		              positions,&
		              absrel,&
		              deltaE,&
		              der,&
		              energy,&
	        		  periodic,&
	        		  periodica,&
	        		  period_length,&
	        		  dip_corr,&
	        		  fmm_temp_0,&
	        		  lineardistance,&
	        		  fmm_temp_0,&
           	          aplummer,&
           	          huge(fmm_temp_1),&
           	          fmm_temp_0,&
                      maxdepth,&
                      unrolled_limit,&
                      balance_load,&
           	          FMM_internal_params)

        ncsar = FMM_internal_params%ncsar
        poles = FMM_internal_params%nmultipoles

        wignersize = (poles+1)**2*(2*poles+1)*4*(ncsar+2)*fmm_real

      	result = 1
      end subroutine fmm_ctune


      subroutine fmm_ctunehomogen(params,wignersize,&
      				  result) bind(c)
      implicit none
      type(c_ptr), value :: params
      type(FMM_internal_params_t), pointer :: FMM_internal_params

      integer(kind = 8) :: ncsar,poles
      integer(kind = c_long_long) :: wignersize
      integer(kind = c_long_long) :: result
    
      call c_f_pointer(params,FMM_internal_params) 

      FMM_internal_params%ws = 1 
      FMM_internal_params%ncsar = 25 
 
      poles = FMM_MAXNMULTIPOLES
      ncsar = 25
      wignersize = (poles+1)**2*(2*poles+1)*4*(ncsar+2)*fmm_real
      result = 1
      end subroutine fmm_ctunehomogen
      
      ! compute wigner matrices subroutine for the C interface
      subroutine fmm_ccomputewigner(wignerptr,handle,dotune) bind(c)
      use mp_wrapper, only : diffcpointers
      implicit none
      type(c_ptr), value :: wignerptr,handle
      type(FMM_internal_params_t), pointer :: FMM_internal_params
      integer(kind=c_long_long), value :: dotune
      real(kind=fmm_real), dimension(:,:,:,:,:),pointer :: tmpwignerptr
      integer(kind=fmm_integer) :: i
      integer(kind=fmm_integer) :: ws,wsd,maxncsar,ncsar,poles
      integer(kind=fmm_integer), allocatable, dimension(:,:) :: fmmcos,icsar

      call c_f_pointer(handle,FMM_internal_params)

      if (dotune.eq.1) then
        poles = FMM_internal_params%nmultipoles
      else
        poles = FMM_MAXNMULTIPOLES
      endif

      ncsar = FMM_internal_params%ncsar
      ws = FMM_internal_params%ws

      call diffcpointers(wignerptr,c_null_ptr,i)
      if(i.ne.0) then
        call c_f_pointer(wignerptr,tmpwignerptr,&
                       [poles+1,2*poles+1,poles+1,4,ncsar+2])
      else
        return
      endif
      call remap_wignerptr(FMM_internal_params,tmpwignerptr,poles,ncsar)

      wsd = 2*ws+1
      maxncsar = ncsar

      allocate(icsar(0:wsd,0:2*wsd*wsd),stat=i)
      if(i.ne.0) call bummer('fmm_ccomputewigner: error, i = ',i)
      allocate(fmmcos(2,ncsar),stat=i)
      if(i.ne.0) call bummer('fmm_ccomputewigner: error, i = ',i)
      call fmmg(wsd,maxncsar,ws,ncsar,icsar,fmmcos,.false.)
      call calallds(ws,poles,ncsar,wsd,icsar,fmmcos,FMM_internal_params%wignerd)

      deallocate(fmmcos,stat=i)
      if(i.ne.0) call bummer('fmm_ccomputewigner: error, i = ',i)
      deallocate(icsar, stat=i)
      if(i.ne.0) call bummer('fmm_ccomputewigner: error, i = ',i)

      end subroutine fmm_ccomputewigner

      subroutine remap_wignerptr(FMM_internal_params,oldptr,poles,ncsar)
      implicit none
      integer(kind=fmm_integer) :: poles,ncsar
      real(kind=fmm_real), dimension(0:poles,-poles:poles,0:poles,1:4,1:ncsar+2), target :: oldptr
      type(FMM_internal_params_t) :: FMM_internal_params

        FMM_internal_params%wignerd%wignerd => oldptr
      end subroutine remap_wignerptr


      subroutine fmm_crun(local_particles,positions,charges,potentials,field,virial,total_particles,absrel,deltaE, &
      						dip_corr, periodicity, period_length, dotune,maxdepth,unrolled_limit,balance_load, &
                                                params, result) bind(c)

		implicit none

      	integer(kind = c_long_long), value	::	local_particles
      	integer(kind = c_long_long), value	::	total_particles
      	integer(kind = c_long_long), value	::	absrel
      	integer(kind = c_long_long), value	::	dip_corr
      	integer(kind = c_long_long), value	::	dotune
        integer(kind = c_long_long), value      ::      maxdepth
        integer(kind = c_long_long), value      ::      unrolled_limit 
        integer(kind = c_long_long), value      ::      balance_load 
      	integer(kind = c_long_long), dimension(3)	:: periodicity
      	real(kind = c_double), value			::	deltaE
		real(kind = c_double), value			::	period_length
      	real(kind = c_double), dimension(3*local_particles) :: positions
      	real(kind = c_double), dimension(local_particles)	:: charges
      	real(kind = c_double), dimension(3*local_particles) :: field
      	real(kind = c_double), dimension(local_particles)	:: potentials
      	real(kind = c_double), dimension(9)	:: virial
      	type(c_ptr), value						::	params

      	type(FMM_internal_params_t), pointer :: FMM_internal_params
      	integer(kind = 8)				::	fmm_temp_0 = 0, fmm_temp_1 = 1
	integer(kind = fmm_integer) :: tunehomogen,homogen
      	real(kind = fmm_real)		::  der = 0.0d0, lineardistance = 0.0d0, aplummer = 0.0d0, energy

		integer(kind = 8)				::	periodic, periodica

		integer(kind = c_long_long)	::	result

        	call c_f_pointer(params,FMM_internal_params)

		periodic = COUNT(periodicity .eq. 1)
		select case(periodic)
			case (0)
				periodica = 0
			case (1)
				if (periodicity(3) == 1) then
					periodica = 0
				else if (periodicity(2) == 1) then
					periodica = 2
				else
					periodica = 1
				end if
			case (2)
				if (periodicity(1) == 1 .and. periodicity(2) == 1) then
					periodica = 0
				else if (periodicity(1) == 1 .and. periodicity(3) == 1) then
					periodica = 2
				else
					periodica = 3
				end if
			case (3)
				periodica = 0
			case default
				result = 0
				return
	  	end select

	if (dotune.eq.1) then
          ! full analysis
          homogen = 0
	else
          ! homogen analysis
          homogen = 1
	endif	
        tunehomogen = huge(1_fmm_integer)

        call cfmm(total_particles,local_particles,charges,positions,absrel,deltaE,der,energy,potentials,&
              	   field,virial,periodic,periodica,period_length,dip_corr,fmm_temp_0,lineardistance,fmm_temp_0,&
                   aplummer,tunehomogen,homogen,maxdepth,unrolled_limit,balance_load,FMM_internal_params)
      	result = 1
      end subroutine fmm_crun

      subroutine fmm_cfinalize(params,dotune) bind(c)
        use mp_wrapper, only : mp_finalize
        use mwigner
	implicit none
        integer(kind=c_long_long), value :: dotune
        integer(kind=8) :: i
        integer(kind=8) :: poles
      	type(c_ptr), value :: params

      	type(FMM_internal_params_t), pointer :: FMM_internal_params

      	call c_f_pointer(params,FMM_internal_params)

        if (dotune.eq.1) then
          poles = FMM_internal_params%nmultipoles
        else
          poles = FMM_MAXNMULTIPOLES
        endif

      	FMM_internal_params%wignerd%wignerd => NULL()
        call mp_finalize()
	
      end subroutine fmm_cfinalize

      end module fmm_cbindings
