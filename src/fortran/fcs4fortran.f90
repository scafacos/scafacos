!
! Copyright (C) 2011, 2012, 2013 Rene Halver, Michael Hofmann
!
! This file is part of ScaFaCoS.
!
! ScaFaCoS is free software: you can redistribute it and/or modify
! it under the terms of the GNU Lesser Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! ScaFaCoS is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser Public License for more details.
!
! You should have received a copy of the GNU Lesser Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>. 
!



#ifdef HAVE_FCONFIG_H
#include <fconfig.h>
#endif

#include <fcs_fconfig.h>

#include "fcs4fortran_definitions.h"


module fcs_module

  use iso_c_binding

  implicit none

  ! boolean data type
  
  integer, parameter :: fcs_boolean_kind = fcs_integer_kind
  integer(kind = fcs_boolean_kind), parameter :: FCS_TRUE = 1
  integer(kind = fcs_boolean_kind), parameter :: FCS_FALSE = 0

  ! ScaFaCoS return values

  integer(kind = fcs_integer_kind_isoc), parameter ::  FCS_SUCCESS                   = FCS4FORTRAN_SUCCESS
  integer(kind = fcs_integer_kind_isoc), parameter ::  FCS_ERROR_NULL_ARGUMENT       = FCS4FORTRAN_ERROR_NULL_ARGUMENT
  integer(kind = fcs_integer_kind_isoc), parameter ::  FCS_ERROR_ALLOC_FAILED        = FCS4FORTRAN_ERROR_ALLOC_FAILED
  integer(kind = fcs_integer_kind_isoc), parameter ::  FCS_ERROR_WRONG_ARGUMENT      = FCS4FORTRAN_ERROR_WRONG_ARGUMENT
  integer(kind = fcs_integer_kind_isoc), parameter ::  FCS_ERROR_MISSING_ELEMENT     = FCS4FORTRAN_ERROR_MISSING_ELEMENT
  integer(kind = fcs_integer_kind_isoc), parameter ::  FCS_ERROR_LOGICAL_ERROR       = FCS4FORTRAN_ERROR_LOGICAL_ERROR
  integer(kind = fcs_integer_kind_isoc), parameter ::  FCS_ERROR_INCOMPATIBLE_METHOD = FCS4FORTRAN_ERROR_INCOMPATIBLE_METHOD
  integer(kind = fcs_integer_kind_isoc), parameter ::  FCS_ERROR_NOT_IMPLEMENTED     = FCS4FORTRAN_ERROR_NOT_IMPLEMENTED
  integer(kind = fcs_integer_kind_isoc), parameter ::  FCS_ERROR_FORTRAN_CALL_ERROR  = FCS4FORTRAN_ERROR_FORTRAN_CALL
  integer(kind = fcs_integer_kind_isoc), parameter ::  FCS_ERROR_RESULT_CREATE       = FCS4FORTRAN_ERROR_RESULT_CREATE

  ! definitions of method flags

  integer(kind = fcs_integer_kind_isoc), parameter ::  FCS_METHOD_NONE   = FCS4FORTRAN_METHOD_NONE
  integer(kind = fcs_integer_kind_isoc), parameter ::  FCS_METHOD_DIRECT = FCS4FORTRAN_METHOD_DIRECT
  integer(kind = fcs_integer_kind_isoc), parameter ::  FCS_METHOD_EWALD  = FCS4FORTRAN_METHOD_EWALD
  integer(kind = fcs_integer_kind_isoc), parameter ::  FCS_METHOD_FMM    = FCS4FORTRAN_METHOD_FMM
  integer(kind = fcs_integer_kind_isoc), parameter ::  FCS_METHOD_MEMD   = FCS4FORTRAN_METHOD_MEMD
  integer(kind = fcs_integer_kind_isoc), parameter ::  FCS_METHOD_MMM1D  = FCS4FORTRAN_METHOD_MMM1D
  integer(kind = fcs_integer_kind_isoc), parameter ::  FCS_METHOD_MMM2D  = FCS4FORTRAN_METHOD_MMM2D
  integer(kind = fcs_integer_kind_isoc), parameter ::  FCS_METHOD_P2NFFT = FCS4FORTRAN_METHOD_P2NFFT
  integer(kind = fcs_integer_kind_isoc), parameter ::  FCS_METHOD_P3M    = FCS4FORTRAN_METHOD_P3M
  integer(kind = fcs_integer_kind_isoc), parameter ::  FCS_METHOD_PEPC   = FCS4FORTRAN_METHOD_PEPC
  integer(kind = fcs_integer_kind_isoc), parameter ::  FCS_METHOD_PP3MG  = FCS4FORTRAN_METHOD_PP3MG
  integer(kind = fcs_integer_kind_isoc), parameter ::  FCS_METHOD_VMG    = FCS4FORTRAN_METHOD_VMG
  integer(kind = fcs_integer_kind_isoc), parameter ::  FCS_METHOD_WOLF   = FCS4FORTRAN_METHOD_WOLF

#ifdef FCS_ENABLE_FMM
  ! fmm specific parameter definition

  integer(kind = fcs_integer_kind_isoc), parameter :: FCS_FMM_COULOMB = 64
  integer(kind = fcs_integer_kind_isoc), parameter :: FCS_FMM_CUSP = 65
  integer(kind = fcs_integer_kind_isoc), parameter :: FCS_FMM_NO_DIPOLE_CORRECTION = -1
  integer(kind = fcs_integer_kind_isoc), parameter :: FCS_FMM_STANDARD_DIPOLE_CORRECTION = 0
  integer(kind = fcs_integer_kind_isoc), parameter :: FCS_FMM_ACTIVE_DIPOLE_CORRECTION = 1
  integer(kind = fcs_integer_kind_isoc), parameter :: FCS_FMM_STANDARD_ERROR = 0
  integer(kind = fcs_integer_kind_isoc), parameter :: FCS_FMM_CUSTOM_ABSOLUTE = 1
  integer(kind = fcs_integer_kind_isoc), parameter :: FCS_FMM_CUSTOM_RELATIVE = 2
#endif

  ! length of function names and description messages of the return state
  integer, parameter                               :: MAX_FUNCTION_LENGTH = 64
  integer, parameter                               :: MAX_MESSAGE_LENGTH = 512

  ! definitions of tolerance type flags
  integer(kind = fcs_integer_kind_isoc), parameter  :: FCS_TOLERANCE_TYPE_UNDEFINED     = 0
  integer(kind = fcs_integer_kind_isoc), parameter  :: FCS_TOLERANCE_TYPE_ENERGY        = 1
  integer(kind = fcs_integer_kind_isoc), parameter  :: FCS_TOLERANCE_TYPE_ENERGY_REL    = 2
  integer(kind = fcs_integer_kind_isoc), parameter  :: FCS_TOLERANCE_TYPE_POTENTIAL     = 3
  integer(kind = fcs_integer_kind_isoc), parameter  :: FCS_TOLERANCE_TYPE_POTENTIAL_REL = 4
  integer(kind = fcs_integer_kind_isoc), parameter  :: FCS_TOLERANCE_TYPE_FIELD         = 5
  integer(kind = fcs_integer_kind_isoc), parameter  :: FCS_TOLERANCE_TYPE_FIELD_REL     = 6

  ! interface containing the calls to the wrapper functions in C

  interface

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                               return value handling
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine fcs_result_destroy(res) BIND(C, name="fcs_result_destroy")
          use iso_c_binding
          type(c_ptr), value                                    ::  res
      end subroutine

      function fcs_result_get_return_code_f(res) BIND(C, name="fcs_result_get_return_code")
          use iso_c_binding
          implicit none
          type(c_ptr), value                                    ::  res
          integer(kind = fcs_integer_kind_isoc)                 ::  fcs_result_get_return_code_f
      end function
      
      function fcs_result_get_message_f(res) BIND(C, name="fcs_result_get_message")
          use iso_c_binding
          import MAX_MESSAGE_LENGTH
          implicit none
          type(c_ptr), value                                    ::  res
          type(c_ptr)                                           ::  fcs_result_get_message_f
      end function

      function fcs_result_get_function_f(res) BIND(C, name="fcs_result_get_function")
          use iso_c_binding
          import MAX_FUNCTION_LENGTH
          implicit none
          type(c_ptr), value                                    ::  res
          type(c_ptr)                                           ::  fcs_result_get_function_f
      end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                 basic ScaFaCoS functions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      function fcs_init(handle,method_name,communicator) BIND(C,name="fcs_init_f")
          use iso_c_binding
          implicit none
          type(c_ptr)                                       :: handle
          character(kind = c_char)                          :: method_name(*)
          integer, value                                    :: communicator
          type(c_ptr)                                       :: fcs_init
      end function

      function fcs_tune(handle,n_locparts,positions,charges) BIND(C,name="fcs_tune")
          use iso_c_binding
          implicit none
          type(c_ptr), value                                :: handle
          integer(kind = fcs_integer_kind_isoc),value       :: n_locparts
          real(kind = fcs_real_kind_isoc)                   :: positions(3*n_locparts)
          real(kind = fcs_real_kind_isoc)                   :: charges(n_locparts)
          type(c_ptr)                                       :: fcs_tune
      end function


      function fcs_run(handle,n_locparts,positions,charges,fields,&
                           potentials) BIND(C,name="fcs_run")
          use iso_c_binding
          implicit none
          type(c_ptr),value                                 :: handle
          integer(kind = fcs_integer_kind_isoc),value       :: n_locparts
          real(kind = fcs_real_kind_isoc)                   :: positions(3*n_locparts)
          real(kind = fcs_real_kind_isoc)                   :: charges(n_locparts)
          real(kind = fcs_real_kind_isoc)                   :: fields(3*n_locparts)
          real(kind = fcs_real_kind_isoc)                   :: potentials(n_locparts)
          type(c_ptr)                                       :: fcs_run
      end function
      
      function fcs_destroy(handle) BIND(C,name="fcs_destroy")
          use iso_c_binding
          implicit none
          type(c_ptr), value                        :: handle
          type(c_ptr)                               :: fcs_destroy
      end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                               general parameter handling
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      function fcs_get_method(handle) BIND(C,name="fcs_get_method")
          use iso_c_binding
          implicit none
          type(c_ptr), value                        :: handle
          integer(kind = fcs_integer_kind_isoc)     :: fcs_get_method
      end function

! Missing: fcs_get_method_name
      function fcs_get_communicator(handle) BIND(C,name="fcs_get_communicator")
          use iso_c_binding
          implicit none
          type(c_ptr), value                        ::  handle
          integer                                   ::  fcs_get_communicator
      end function

      function fcs_set_common_f(handle, near_field_flag, box_a, box_b, box_c, &
           box_origin, periodicity, total_parts) BIND(C,name="fcs_set_common")
          use iso_c_binding
          implicit none
          type(c_ptr),value                                 ::  handle
          integer(kind = fcs_integer_kind_isoc), value      ::  near_field_flag
          real(kind = fcs_real_kind_isoc)                   ::  box_a(3)
          real(kind = fcs_real_kind_isoc)                   ::  box_b(3)
          real(kind = fcs_real_kind_isoc)                   ::  box_c(3)
          real(kind = fcs_real_kind_isoc)                   ::  box_origin(3)
          integer(kind = fcs_integer_kind_isoc)             ::  periodicity(3)
          integer(kind = fcs_integer_kind_isoc), value      ::  total_parts
          type(c_ptr)                                       ::  fcs_set_common_f
      end function

      function fcs_set_dimensions(handle, dim_) &
            BIND(C,name="fcs_set_dimensions")
          use iso_c_binding
          implicit none
          type(c_ptr),value                                 ::  handle
          integer(kind = fcs_integer_kind_isoc), value      ::  dim_
          type(c_ptr)                                       ::  fcs_set_dimensions
      end function
      
      function fcs_get_dimensions(handle) &
                BIND(C,name="fcs_get_dimensions")
          use iso_c_binding
          implicit none
          type(c_ptr),value                                 ::  handle
          integer(kind = fcs_integer_kind_isoc)             ::  fcs_get_dimensions
      end function

      function fcs_set_near_field_flag(handle, near_field_flag) &
           BIND(C,name="fcs_set_near_field_flag")
          use iso_c_binding
          implicit none
          type(c_ptr),value                                 ::  handle
          integer(kind = fcs_integer_kind_isoc),value       ::  near_field_flag
          type(c_ptr)                                       ::  fcs_set_near_field_flag
      end function

      function fcs_get_near_field_flag(handle) &
              BIND(C,name="fcs_get_near_field_flag")
          use iso_c_binding
          implicit none
          type(c_ptr), value                                ::  handle
          integer(kind = fcs_integer_kind_isoc)             ::  fcs_get_near_field_flag
      end function

      function fcs_set_box_a(handle, box_a) BIND(C,name="fcs_set_box_a")
          use iso_c_binding
          implicit none
          type(c_ptr),value                                 ::  handle
          real(kind = fcs_real_kind_isoc)                   ::  box_a(3)
          type(c_ptr)                                       ::  fcs_set_box_a
      end function

      ! requires helper function to transform C pointer to Fortran pointer to
      ! real(3)
!      function fcs_get_box_a(handle) BIND(C,name="fcs_get_box_a")
!          use iso_c_binding
!          implicit none
!          type(c_ptr),value                                 ::  handle
!          real(kind = fcs_real_kind_isoc), dimension(3)     ::  fcs_get_box_a
!      end function

      function fcs_set_box_b(handle, box_b) BIND(C,name="fcs_set_box_b")
          use iso_c_binding
          implicit none
          type(c_ptr),value                                 ::  handle
          real(kind = fcs_real_kind_isoc)                   ::  box_b(3)
          type(c_ptr)                                       ::  fcs_set_box_b
      end function

! Missing: fcs_get_box_b
      
      function fcs_set_box_c(handle, box_c) BIND(C,name="fcs_set_box_c")
          use iso_c_binding
          implicit none
          type(c_ptr),value                                 ::  handle
          real(kind = fcs_real_kind_isoc)                   ::  box_c(3)
          type(c_ptr)                                       ::  fcs_set_box_c
      end function

! Missing: fcs_get_box_c

      function fcs_set_box_origin(handle, box_origin) BIND(C,name="fcs_set_box_origin")
          use iso_c_binding
          implicit none
          type(c_ptr),value                                 ::  handle
          real(kind = fcs_real_kind_isoc)                   ::  box_origin(3)
          type(c_ptr)                                       ::  fcs_set_box_origin
      end function

! Missing: fcs_get_box_origin

      function fcs_set_periodicity_f(handle, periodicity) &
                 BIND(C,name="fcs_set_periodicity")
          use iso_c_binding
          implicit none
          type(c_ptr),value                                 ::  handle
          integer(kind = fcs_integer_kind_isoc)             ::  periodicity(3)
          type(c_ptr)                                       ::  fcs_set_periodicity_f
      end function

! Missing: fcs_get_periodicity_f

      function fcs_set_total_particles(handle, total_particles) &
                 BIND(C,name="fcs_set_total_particles")
          use iso_c_binding
          implicit none
          type(c_ptr),value                                 ::  handle
          integer(kind = fcs_integer_kind_isoc),value       ::  total_particles
          type(c_ptr)                                       ::  fcs_set_total_particles
      end function

      function fcs_get_total_particles(handle) &
                BIND(C,name="fcs_get_total_particles")
          use iso_c_binding
          implicit none
          type(c_ptr),value                                 ::  handle
          integer(kind = fcs_integer_kind_isoc)             :: fcs_get_total_particles
      end function

      function fcs_set_max_local_particles(handle, max_local_particles) &
                 BIND(C,name="fcs_set_max_local_particles")
          use iso_c_binding
          implicit none
          type(c_ptr),value                                 ::  handle
          integer(kind = fcs_integer_kind_isoc),value       ::  max_local_particles
          type(c_ptr)                                       ::  fcs_set_max_local_particles
      end function

      function fcs_get_max_local_particles(handle) &
                BIND(C,name="fcs_get_max_local_particles")
          use iso_c_binding
          implicit none
          type(c_ptr),value                                 ::  handle
          integer(kind = fcs_integer_kind_isoc)             :: fcs_get_max_local_particles
      end function

      function fcs_set_tolerance(handle, tolerance_type, tolerance) &
                BIND(C,name="fcs_set_tolerance")
          use iso_c_binding
          implicit none
          type(c_ptr),value                                 ::  handle
          integer(kind = fcs_integer_kind_isoc),value       ::  tolerance_type
          real(kind = fcs_real_kind_isoc), value            ::  tolerance
          type(c_ptr)                                       ::  fcs_set_tolerance
      end function


      function fcs_get_tolerance(handle,tolerance_type,tolerance) &
                BIND(C,name="fcs_get_tolerance")
          use iso_c_binding
          implicit none
          type(c_ptr), value                                ::  handle
          integer(kind = fcs_integer_kind_isoc)             ::  tolerance_type
          real(kind = fcs_real_kind_isoc), value            ::  tolerance
          type(c_ptr)                                       ::  fcs_get_tolerance
      end function

      function fcs_set_r_cut(handle, r_cut) BIND(C,name="fcs_set_r_cut")
            use iso_c_binding
            implicit none
            type(c_ptr), value                              ::  handle
            real(kind = fcs_real_kind_isoc), value          ::  r_cut
            type(c_ptr)                                     ::  fcs_set_r_cut
      end function

      function fcs_unset_r_cut(handle) BIND(C,name="fcs_unset_r_cut")
            use iso_c_binding
            implicit none
            type(c_ptr), value                              ::  handle
            type(c_ptr)                                     ::  fcs_unset_r_cut
      end function

      function fcs_get_r_cut(handle, r_cut) BIND(C,name="fcs_get_r_cut")
            use iso_c_binding
            implicit none
            type(c_ptr), value                              ::  handle
            real(kind = fcs_real_kind_isoc)                 ::  r_cut
            type(c_ptr)                                     ::  fcs_get_r_cut
      end function

      function fcs_set_parameters(handle,parameters,continue_on_errors) BIND(C,name="fcs_set_parameters")
          use iso_c_binding
          implicit none
          type(c_ptr), value                                :: handle
          character(kind = c_char)                          :: parameters(*)
          type(c_ptr)                                       :: fcs_set_parameters
          integer(kind = fcs_boolean_kind_isoc)             :: continue_on_errors
      end function

      subroutine fcs_print_parameters(handle) BIND(C,name="fcs_print_parameters")
          use iso_c_binding
          implicit none
          type(c_ptr), value                        :: handle
      end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                    misc functions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      function fcs_compute_dipole_correction(handle, local_particles, positions, &
                                             charges, epsilon, field_correction, & 
                                             energy_correction) &
            BIND(C,name="fcs_compute_dipole_correction")
            use iso_c_binding
            implicit none
            type(c_ptr), value                                                      ::  handle
            integer(kind = fcs_integer_kind_isoc), value                            ::  local_particles
            real(kind = fcs_real_kind_isoc), dimension(3*local_particles)           ::  positions
            real(kind = fcs_real_kind_isoc), dimension(local_particles)             ::  charges
            real(kind = fcs_real_kind_isoc), value                                  ::  epsilon
            real(kind = fcs_real_kind_isoc), dimension(3)                           ::  field_correction
            real(kind = fcs_real_kind_isoc)                                         ::  energy_correction
            type(c_ptr)                                                             ::  fcs_compute_dipole_correction
      end function

      function fcs_get_near_field_delegation_f(handle, has_near) BIND(C,name="fcs_get_near_field_delegation")
          use iso_c_binding
          implicit none
          type(c_ptr), value                                  ::  handle
          integer(kind = fcs_integer_kind_isoc)               ::  has_near
          type(c_ptr)                                         ::  fcs_get_near_field_delegation_f
      end function

      function fcs_compute_near(handle, dist, pot, field) BIND(C,name="fcs_compute_near")
          use iso_c_binding
          implicit none
          type(c_ptr), value                                  ::  handle
          real(kind = fcs_real_kind_isoc), value              ::  dist
          real(kind = fcs_real_kind_isoc)                     ::  pot
          real(kind = fcs_real_kind_isoc)                     ::  field
          type(c_ptr)                                         ::  fcs_compute_near
      end function

      function fcs_compute_near_potential(handle, dist, pot) BIND(C,name="fcs_compute_near_potential")
          use iso_c_binding
          implicit none
          type(c_ptr), value                                  ::  handle
          real(kind = fcs_real_kind_isoc), value              ::  dist
          real(kind = fcs_real_kind_isoc)                     ::  pot
          type(c_ptr)                                         ::  fcs_compute_near_potential
      end function

      function fcs_compute_near_field(handle, dist, field) BIND(C,name="fcs_compute_near_field")
          use iso_c_binding
          implicit none
          type(c_ptr), value                                  ::  handle
          real(kind = fcs_real_kind_isoc), value              ::  dist
          real(kind = fcs_real_kind_isoc)                     ::  field
          type(c_ptr)                                         ::  fcs_compute_near_field
      end function

      function fcs_set_compute_virial_f(handle, flag) BIND(C,name="fcs_set_compute_virial")
          use iso_c_binding
          implicit none
          type(c_ptr), value                                :: handle
          integer(kind = fcs_integer_kind_isoc),value       :: flag
          type(c_ptr)                                       :: fcs_set_compute_virial_f
      end function

! Missing: fcs_get_compute_virial
      
      function fcs_get_virial(handle, virial) BIND(C,name="fcs_get_virial")
          use iso_c_binding
          implicit none
          type(c_ptr), value                                :: handle
          real(kind = fcs_real_kind_isoc)                   :: virial(9)
          type(c_ptr)                                       :: fcs_get_virial
      end function
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                  method specific setups
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      

#ifdef FCS_ENABLE_DIRECT
      function fcs_direct_setup(handle, cutoff) &
                                   BIND(C,name="fcs_direct_setup")
          use iso_c_binding
          implicit none
          type(c_ptr),value                                 ::  handle
          real(kind = fcs_real_kind_isoc),value             ::  cutoff
          type(c_ptr)                                       ::  fcs_direct_setup
      end function
#endif
#ifdef FCS_ENABLE_EWALD
      function fcs_ewald_set_tolerance_field_abs(handle, tolerance_field_abs) &
           BIND(C,name="fcs_ewald_set_tolerance_field_abs")
          use iso_c_binding
          implicit none
          type(c_ptr),value                                 ::  handle
          real(kind = fcs_real_kind_isoc),value             ::  tolerance_field_abs
          type(c_ptr)                                       ::  fcs_ewald_set_tolerance_field_abs
      end function
#endif
#ifdef FCS_ENABLE_FMM
      function fcs_fmm_setup(handle, absrel, energy_tolerance, dipole_correction, system, maxdepth, unroll_limit, balanceload) &
                                  BIND(C,name="fcs_fmm_setup")
          use iso_c_binding
          implicit none
          type(c_ptr),value                                 ::  handle
          integer(kind = fcs_integer_kind_isoc), value      ::  absrel
          real(kind = fcs_real_kind_isoc), value            ::  energy_tolerance
          integer(kind = fcs_integer_kind_isoc), value      ::  dipole_correction
          integer(kind = c_long_long), value                ::  system
          integer(kind = c_long_long), value                ::  maxdepth
          integer(kind = c_long_long), value                ::  unroll_limit
          integer(kind = c_long_long), value                ::  balanceload
          type(c_ptr)                                       ::  fcs_fmm_setup
      end function
#endif
#ifdef FCS_ENABLE_PEPC
      function fcs_pepc_setup(handle, epsilon, theta, level) &
                                  BIND(C,name="fcs_pepc_setup")
          use iso_c_binding
          implicit none
          type(c_ptr),value                                 ::  handle
          real(kind = fcs_real_kind_isoc), value            ::  epsilon
          real(kind = fcs_real_kind_isoc), value            ::  theta
          integer(kind = fcs_integer_kind_isoc), value      ::  level
          type(c_ptr)                                       ::  fcs_pepc_setup
      end function
#endif

#ifdef FCS_ENABLE_VMG
      function fcs_vmg_setup(handle, max_level, max_iterations, smooth_steps, &
                            gamma, precision, near_field_cells) &
                                   BIND(C,name="fcs_vmg_setup")
          use iso_c_binding
          implicit none
          type(c_ptr),value                                 ::  handle
          integer(kind = fcs_integer_kind_isoc), value      ::  max_level
          integer(kind = fcs_integer_kind_isoc), value      ::  max_iterations
          integer(kind = fcs_integer_kind_isoc), value      ::  smooth_steps
          integer(kind = fcs_integer_kind_isoc), value      ::  gamma
          integer(kind = fcs_integer_kind_isoc), value      ::  near_field_cells
          real(kind = fcs_real_kind_isoc), value            ::  precision
          
          type(c_ptr)                                       ::  fcs_vmg_setup
      end function
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                           method-specific getters and setters
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef FCS_ENABLE_DIRECT

      function fcs_direct_set_cutoff(handle, cutoff) &
                                   BIND(C,name="fcs_direct_set_cutoff")
          use iso_c_binding
          implicit none
          type(c_ptr), value                                ::  handle
          real(kind = fcs_real_kind_isoc), value            ::  cutoff
          type(c_ptr)                                       ::  fcs_direct_set_cutoff
      end function
   
      function fcs_direct_get_cutoff(handle, cutoff) &
                                   BIND(C,name="fcs_direct_get_cutoff")
          use iso_c_binding
          implicit none
          type(c_ptr), value                                ::  handle
          real(kind = fcs_real_kind_isoc)                   ::  cutoff
          type(c_ptr)                                       ::  fcs_direct_get_cutoff
      end function
#endif
#ifdef FCS_ENABLE_FMM
      function fcs_fmm_set_absrel(handle, absrel) &
                                   BIND(C,name="fcs_fmm_set_absrel")
          use iso_c_binding
          implicit none
          type(c_ptr), value                                ::  handle
          integer(kind = fcs_integer_kind_isoc), value      ::  absrel
          type(c_ptr)                                       ::  fcs_fmm_set_absrel
      end function
      
      function fcs_fmm_get_absrel(handle, absrel) &
                                   BIND(C,name="fcs_fmm_get_absrel")
          use iso_c_binding
          implicit none
          type(c_ptr), value                                ::  handle
          integer(kind = fcs_integer_kind_isoc)             ::  absrel
          type(c_ptr)                                       ::  fcs_fmm_get_absrel
      end function

      function fcs_fmm_set_dipole_correction(handle, dipole_correction) &
                                   BIND(C,name="fcs_fmm_set_dipole_correction")
          use iso_c_binding
          implicit none
          type(c_ptr), value                                ::  handle
          integer(kind = fcs_integer_kind_isoc), value      ::  dipole_correction
          type(c_ptr)                                       ::  fcs_fmm_set_dipole_correction
      end function
      
      function fcs_fmm_get_dipole_correction(handle, dipole_correction) &
                                   BIND(C,name="fcs_fmm_get_dipole_correction")
          use iso_c_binding
          implicit none
          type(c_ptr), value                                ::  handle
          integer(kind = fcs_integer_kind_isoc)             ::  dipole_correction
          type(c_ptr)                                       ::  fcs_fmm_get_dipole_correction
      end function

      function fcs_fmm_set_potential(handle, potential) &
                                   BIND(C,name="fcs_fmm_set_potential")
          use iso_c_binding
          implicit none
          type(c_ptr), value                                ::  handle
          integer(kind = fcs_integer_kind_isoc), value      ::  potential
          type(c_ptr)                                       ::  fcs_fmm_set_potential
      end function
      
      function fcs_fmm_get_potential(handle, potential) &
                                   BIND(C,name="fcs_fmm_get_potential")
          use iso_c_binding
          implicit none
          type(c_ptr), value                                ::  handle
          integer(kind = fcs_integer_kind_isoc)             ::  potential
          type(c_ptr)                                       ::  fcs_fmm_get_potential
      end function

      function fcs_fmm_set_cusp_radius(handle, cusp_radius) &
                                   BIND(C,name="fcs_fmm_set_cusp_radius")
          use iso_c_binding
          implicit none
          type(c_ptr), value                                ::  handle
          real(kind = fcs_real_kind_isoc), value            ::  cusp_radius
          type(c_ptr)                                       ::  fcs_fmm_set_cusp_radius
      end function
      
      function fcs_fmm_get_cusp_radius(handle, cusp_radius) &
                                   BIND(C,name="fcs_fmm_get_cusp_radius")
          use iso_c_binding
          implicit none
          type(c_ptr), value                                ::  handle
          real(kind = fcs_real_kind_isoc)                   ::  cusp_radius
          type(c_ptr)                                       ::  fcs_fmm_get_cusp_radius
      end function

      function fcs_fmm_set_tolerance_energy(handle, tolerance) &
                                   BIND(C,name="fcs_fmm_set_tolerance_energy")
          use iso_c_binding
          implicit none
          type(c_ptr), value                                ::  handle
          real(kind = fcs_real_kind_isoc), value            ::  tolerance
          type(c_ptr)                                       ::  fcs_fmm_set_tolerance_energy
      end function
      
      function fcs_fmm_get_tolerance_energy(handle, tolerance) &
                                   BIND(C,name="fcs_fmm_get_tolerance_energy")
          use iso_c_binding
          implicit none
          type(c_ptr), value                                ::  handle
          real(kind = fcs_real_kind_isoc)                   ::  tolerance
          type(c_ptr)                                       ::  fcs_fmm_get_tolerance_energy
      end function

      function fcs_fmm_set_maxdepth(handle, maxdepth) &
                                   BIND(C,name="fcs_fmm_set_maxdepth")
          use iso_c_binding
          implicit none
          type(c_ptr), value                                ::  handle
          integer(kind = c_long_long), value                ::  maxdepth
          type(c_ptr)                                       ::  fcs_fmm_set_maxdepth
      end function
      
      function fcs_fmm_get_maxdepth(handle, maxdepth) &
                                   BIND(C,name="fcs_fmm_get_maxdepth")
          use iso_c_binding
          implicit none
          type(c_ptr), value                                ::  handle
          integer(kind = c_long_long), value                ::  maxdepth
          type(c_ptr)                                       ::  fcs_fmm_get_maxdepth
      end function

      function fcs_fmm_set_unroll_limit(handle, unroll_limit) &
                                   BIND(C,name="fcs_fmm_set_unroll_limit")
          use iso_c_binding
          implicit none
          type(c_ptr), value                                ::  handle
          integer(kind = c_long_long), value                ::  unroll_limit
          type(c_ptr)                                       ::  fcs_fmm_set_unroll_limit
      end function
      
      function fcs_fmm_get_unroll_limit(handle, unroll_limit) &
                                   BIND(C,name="fcs_fmm_get_unroll_limit")
          use iso_c_binding
          implicit none
          type(c_ptr), value                                ::  handle
          integer(kind = c_long_long), value                ::  unroll_limit
          type(c_ptr)                                       ::  fcs_fmm_get_unroll_limit
      end function

      function fcs_fmm_set_balanceload(handle, balanceload) &
                                   BIND(C,name="fcs_fmm_set_balanceload")
          use iso_c_binding
          implicit none
          type(c_ptr), value                                ::  handle
          integer(kind = c_long_long), value                ::  balanceload
          type(c_ptr)                                       ::  fcs_fmm_set_balanceload
      end function
      
      function fcs_fmm_get_balanceload(handle, balanceload) &
                                   BIND(C,name="fcs_fmm_get_balanceload")
          use iso_c_binding
          implicit none
          type(c_ptr), value                                ::  handle
          integer(kind = c_long_long), value                ::  balanceload
          type(c_ptr)                                       ::  fcs_fmm_get_balanceload
      end function

      function fcs_fmm_set_internal_tuning(handle, tuning) &
                                  BIND(C,name="fcs_fmm_set_internal_tuning")
          use iso_c_binding
          implicit none
          type(c_ptr), value                                ::  handle
          integer(kind = c_long_long), value                ::  tuning
          type(c_ptr)                                       ::  fcs_fmm_set_internal_tuning
      end function

      function fcs_fmm_get_internal_tuning(handle, tuning) &
                                   BIND(C,name="fcs_fmm_get_internal_tuning")
          use iso_c_binding
          implicit none
          type(c_ptr), value                                ::  handle
          integer(kind = c_long_long), value                ::  tuning
          type(c_ptr)                                       ::  fcs_fmm_get_internal_tuning
      end function
          
#endif
#ifdef FCS_ENABLE_MEMD
      function fcs_memd_set_periodicity(handle, periodicity) &
                                   BIND(C,name="fcs_memd_set_periodicity")
          use iso_c_binding
          implicit none
          type(c_ptr), value                                ::  handle
          integer(kind = fcs_integer_kind_isoc), value      ::  periodicity
          type(c_ptr)                                       ::  fcs_memd_set_periodicity
      end function
      
      function fcs_memd_get_periodicity(handle, periodicity) &
                                   BIND(C,name="fcs_memd_get_periodicity")
          use iso_c_binding
          implicit none
          type(c_ptr), value                                ::  handle
          integer(kind = fcs_integer_kind_isoc)             ::  periodicity
          type(c_ptr)                                       ::  fcs_memd_get_periodicity
      end function
#endif
#ifdef FCS_ENABLE_MMM1D
      function fcs_mmm1d_set_far_switch_radius(handle, radius) &
                                   BIND(C,name="fcs_mmm1d_set_far_switch_radius")
          use iso_c_binding
          implicit none
          type(c_ptr), value                                ::  handle
          real(kind = fcs_real_kind_isoc), value            ::  radius
          type(c_ptr)                                       ::  fcs_mmm1d_set_far_switch_radius
      end function
      
      function fcs_mmm1d_get_far_switch_radius(handle, radius) &
                                   BIND(C,name="fcs_mmm1d_get_far_switch_radius")
          use iso_c_binding
          implicit none
          type(c_ptr), value                                ::  handle
          real(kind = fcs_real_kind_isoc)                   ::  radius
          type(c_ptr)                                       ::  fcs_mmm1d_get_far_switch_radius
      end function

      function fcs_mmm1d_set_maxPWerror(handle, error) &
                                   BIND(C,name="fcs_mmm1d_set_maxPWerror")
          use iso_c_binding
          implicit none
          type(c_ptr), value                                ::  handle
          real(kind = fcs_real_kind_isoc), value            ::  error
          type(c_ptr)                                       ::  fcs_mmm1d_set_maxPWerror
      end function
      
      function fcs_mmm1d_get_maxPWerror(handle, error) &
                                   BIND(C,name="fcs_mmm1d_get_maxPWerror")
          use iso_c_binding
          implicit none
          type(c_ptr), value                                ::  handle
          real(kind = fcs_real_kind_isoc)                   ::  error
          type(c_ptr)                                       ::  fcs_mmm1d_get_maxPWerror
      end function

      function fcs_mmm1d_set_coulomb_prefactor(handle, prefac) &
                                   BIND(C,name="fcs_mmm1d_set_coulomb_prefactor")
          use iso_c_binding
          implicit none
          type(c_ptr), value                                ::  handle
          real(kind = fcs_real_kind_isoc), value            ::  prefac
          type(c_ptr)                                       ::  fcs_mmm1d_set_coulomb_prefactor
      end function
      
      function fcs_mmm1d_get_coulomb_prefactor(handle, prefac) &
                                   BIND(C,name="fcs_mmm1d_get_coulomb_prefactor")
          use iso_c_binding
          implicit none
          type(c_ptr), value                                ::  handle
          real(kind = fcs_real_kind_isoc)                   ::  prefac
          type(c_ptr)                                       ::  fcs_mmm1d_get_coulomb_prefactor
      end function

      function fcs_mmm1d_set_bessel_cutoff(handle, cutoff) &
                                   BIND(C,name="fcs_mmm1d_set_bessel_cutoff")
          use iso_c_binding
          implicit none
          type(c_ptr), value                                ::  handle
          integer(kind = fcs_integer_kind_isoc), value      ::  cutoff
          type(c_ptr)                                       ::  fcs_mmm1d_set_bessel_cutoff
      end function
      
      function fcs_mmm1d_get_bessel_cutoff(handle, cutoff) &
                                   BIND(C,name="fcs_mmm1d_get_bessel_cutoff")
          use iso_c_binding
          implicit none
          type(c_ptr), value                                ::  handle
          integer(kind = fcs_integer_kind_isoc)             ::  cutoff
          type(c_ptr)                                       ::  fcs_mmm1d_get_bessel_cutoff
      end function
#endif
#ifdef FCS_ENABLE_P2NFFT
      function fcs_p2nfft_set_required_accuracy(handle, required_accuracy) &
           BIND(C,name="fcs_p2nfft_set_required_accuracy")
          use iso_c_binding
          implicit none
          type(c_ptr),value                                 ::  handle
          real(kind = fcs_real_kind_isoc),value             ::  required_accuracy
          type(c_ptr)                                       ::  fcs_p2nfft_set_required_accuracy
      end function
#endif

#ifdef FCS_ENABLE_P3M
      function fcs_p3m_set_tolerance_field_abs(handle, tolerance_field_abs) &
           BIND(C,name="fcs_p3m_set_tolerance_field_abs")
          use iso_c_binding
          implicit none
          type(c_ptr),value                                 ::  handle
          real(kind = fcs_real_kind_isoc),value             ::  tolerance_field_abs
          type(c_ptr)                                       ::  fcs_p3m_set_tolerance_field_abs
      end function
#endif

#ifdef FCS_ENABLE_VMG
      function fcs_vmg_set_gamma(handle, gamma) &
           BIND(C,name="fcs_vmg_set_gamma")
          use iso_c_binding
          implicit none
          type(c_ptr),value                                 ::  handle
          integer(kind = fcs_integer_kind_isoc),value       ::  gamma
          type(c_ptr)                                       ::  fcs_vmg_set_gamma
      end function

      function fcs_vmg_get_gamma(handle,gamma) &
           BIND(C,name="fcs_vmg_get_gamma")
          use iso_c_binding
          implicit none
          type(c_ptr), value                                ::  handle
          integer(kind = fcs_integer_kind_isoc)             ::  gamma
          type(c_ptr)                                       ::  fcs_vmg_get_gamma
      end function

      function fcs_vmg_set_max_iterations(handle, max_iterations) &
           BIND(C,name="fcs_vmg_set_max_iterations")
          use iso_c_binding
          implicit none
          type(c_ptr),value                                 ::  handle
          integer(kind = fcs_integer_kind_isoc),value       ::  max_iterations
          type(c_ptr)                                       ::  fcs_vmg_set_max_iterations
      end function

      function fcs_vmg_get_max_iterations(handle,max_iterations) &
           BIND(C,name="fcs_vmg_get_max_iterations")
          use iso_c_binding
          implicit none
          type(c_ptr), value                                ::  handle
          integer(kind = fcs_integer_kind_isoc)             ::  max_iterations
          type(c_ptr)                                       ::  fcs_vmg_get_max_iterations
      end function

      function fcs_vmg_set_max_level(handle, max_level) &
           BIND(C,name="fcs_vmg_set_max_level")
          use iso_c_binding
          implicit none
          type(c_ptr),value                                 ::  handle
          integer(kind = fcs_integer_kind_isoc),value       ::  max_level
          type(c_ptr)                                       ::  fcs_vmg_set_max_level
      end function

      function fcs_vmg_get_max_level(handle,max_level) &
           BIND(C,name="fcs_vmg_get_max_level")
          use iso_c_binding
          implicit none
          type(c_ptr), value                                ::  handle
          integer(kind = fcs_integer_kind_isoc)             ::  max_level
          type(c_ptr)                                       ::  fcs_vmg_get_max_level
      end function

      function fcs_vmg_set_near_field_cells(handle, near_field_cells) &
           BIND(C,name="fcs_vmg_set_near_field_cells")
          use iso_c_binding
          implicit none
          type(c_ptr),value                                 ::  handle
          integer(kind = fcs_integer_kind_isoc),value       ::  near_field_cells
          type(c_ptr)                                       ::  fcs_vmg_set_near_field_cells
      end function

      function fcs_vmg_get_near_field_cells(handle,near_field_cells) &
           BIND(C,name="fcs_vmg_get_near_field_cells")
          use iso_c_binding
          implicit none
          type(c_ptr), value                                ::  handle
          integer(kind = fcs_integer_kind_isoc)             ::  near_field_cells
          type(c_ptr)                                       ::  fcs_vmg_get_near_field_cells
      end function

      function fcs_vmg_set_precision(handle, prec) &
           BIND(C,name="fcs_vmg_set_precision")
          use iso_c_binding
          implicit none
          type(c_ptr),value                                 ::  handle
          real(kind = fcs_real_kind_isoc),value             ::  prec
          type(c_ptr)                                       ::  fcs_vmg_set_precision
      end function

      function fcs_vmg_get_precision(handle,prec) &
           BIND(C,name="fcs_vmg_get_precision")
          use iso_c_binding
          implicit none
          type(c_ptr), value                                ::  handle
          real(kind = fcs_real_kind_isoc)                   ::  prec
          type(c_ptr)                                       ::  fcs_vmg_get_precision
      end function

      function fcs_vmg_set_smoothing_steps(handle, smoothing_steps) &
           BIND(C,name="fcs_vmg_set_smoothing_steps")
          use iso_c_binding
          implicit none
          type(c_ptr),value                                 ::  handle
          integer(kind = fcs_integer_kind_isoc),value       ::  smoothing_steps
          type(c_ptr)                                       ::  fcs_vmg_set_smoothing_steps
      end function

      function fcs_vmg_get_smoothing_steps(handle,smoothing_steps) &
           BIND(C,name="fcs_vmg_get_smoothing_steps")
          use iso_c_binding
          implicit none
          type(c_ptr), value                                ::  handle
          integer(kind = fcs_integer_kind_isoc)             ::  smoothing_steps
          type(c_ptr)                                       ::  fcs_vmg_get_smoothing_steps
      end function

#endif
   end interface



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                             definitions for FORTRAN interface
!                       (the functions / subroutines where data must
!                        be converted in FORTRAN, before being passed
!                                              to C)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  contains

  function fcs_result_get_return_code(res)
    use iso_c_binding
    implicit none
    type(c_ptr), target                               ::  res
    integer(kind = fcs_integer_kind_isoc)             ::  fcs_result_get_return_code
    
    if (C_ASSOCIATED(res)) then
      fcs_result_get_return_code = fcs_result_get_return_code_f(res)
    else
      fcs_result_get_return_code = FCS_SUCCESS
    end if
    
  end function

  function fcs_result_get_message(res)
    use iso_c_binding
    implicit none
    type(c_ptr), target                                                         ::  res
    character(kind = c_char, len = MAX_MESSAGE_LENGTH)                              ::  fcs_result_get_message
    character(kind = c_char, len = MAX_MESSAGE_LENGTH), dimension(:),pointer        ::  message
    type(c_ptr)                                                                 ::  c_str
    
    if (C_ASSOCIATED(res)) then
      c_str = fcs_result_get_message_f(res)
      call c_f_pointer(cptr = c_str, fptr = message, shape = [1])
      fcs_result_get_message = message(1)
      if ( fcs_get_position_char(fcs_result_get_message,C_NULL_CHAR) == 1) then
        fcs_result_get_message = "no specific error message availiable"
      else
        fcs_result_get_message = fcs_result_get_message(1:fcs_get_position_char(fcs_result_get_message,C_NULL_CHAR)-1)
      end if
    else
      fcs_result_get_message = "call successful"
    end if
  end function

  function fcs_result_get_function(res)
    use iso_c_binding
    implicit none
    type(c_ptr), target                                                         ::  res
    character(kind = c_char, len = MAX_MESSAGE_LENGTH)                              ::  fcs_result_get_function
    character(kind = c_char, len = MAX_MESSAGE_LENGTH), dimension(:), pointer       ::  message
    type(c_ptr)                                                                 ::  c_str
    
    if (C_ASSOCIATED(res)) then
      c_str = fcs_result_get_function_f(res)
      call c_f_pointer(cptr = c_str, fptr = message, shape = [1])
      fcs_result_get_function = message(1)
      if ( fcs_get_position_char(fcs_result_get_function,C_NULL_CHAR) == 1) then
        fcs_result_get_function = "no specific error source availiable"
      else
        fcs_result_get_function = fcs_result_get_function(1:fcs_get_position_char(fcs_result_get_function,C_NULL_CHAR)-1)
      end if
    else
      fcs_result_get_function = ""
    end if
     
  end function

  function fcs_set_common(handle, near_field_flag, box_a, box_b, box_c, &
                          box_origin, periodicity, total_parts)
    use iso_c_binding
    implicit none
    type(c_ptr)                                       ::  handle
    logical                                           ::  near_field_flag
    real(kind = fcs_real_kind_isoc)                   ::  box_a(3)
    real(kind = fcs_real_kind_isoc)                   ::  box_b(3)
    real(kind = fcs_real_kind_isoc)                   ::  box_c(3)
    real(kind = fcs_real_kind_isoc)                   ::  box_origin(3)
    logical                                           ::  periodicity(3)
    integer(kind = fcs_integer_kind_isoc)             ::  total_parts
    integer(kind = fcs_integer_kind_isoc)             ::  p_c(3)
    integer(kind = fcs_integer_kind_isoc)             ::  srf_c
    type(c_ptr)                                       ::  fcs_set_common
    
    where (periodicity)
      p_c = 1
    elsewhere
      p_c = 0
    end where
    
    if (near_field_flag) then
      srf_c = 1
    else
      srf_c = 0
    end if
    
    fcs_set_common = fcs_set_common_f(handle, srf_c, box_a, box_b, box_c, box_origin, p_c, total_parts)
  end function

  function fcs_set_periodicity(handle, periodicity)
    use iso_c_binding
    implicit none
    type(c_ptr)                                       ::  handle
    integer(kind = fcs_integer_kind_isoc)             ::  p_c(3)
    logical                                           ::  periodicity(3)
    type(c_ptr)                                       ::  fcs_set_periodicity
    
    where (periodicity)
      p_c = 1
    elsewhere
      p_c = 0
    end where
    
    fcs_set_periodicity = fcs_set_periodicity_f(handle, p_c)
  end function


  function fcs_set_compute_virial(handle, flag) 
    use iso_c_binding
    implicit none
    type(c_ptr)                                       :: handle
    logical                                           :: flag
    type(c_ptr)                                       :: fcs_set_compute_virial
    integer(kind = fcs_integer_kind_isoc)             :: c_flag
    
    if (flag) then
      c_flag = 1
    else
      c_flag = 0
    end if
    
    fcs_set_compute_virial = fcs_set_compute_virial_f(handle, c_flag)
   end function

   function fcs_get_near_field_delegation(handle, has_near)
    use iso_c_binding
    implicit none
    type(c_ptr), value                                :: handle
    logical                                           :: has_near
    type(c_ptr)                                       :: fcs_get_near_field_delegation
    integer(kind = fcs_integer_kind_isoc)             :: c_has_near

    fcs_get_near_field_delegation = fcs_get_near_field_delegation_f(handle,c_has_near)
    if (c_has_near == 0) then
      has_near = .false.
    else
      has_near = .true.
    end if
     
   end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                   helper function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function fcs_get_position_char(str,c) result(idx)
    use iso_c_binding
    implicit none
    character(kind = c_char, len = *)   ::  str
    character(kind = c_char)            ::  c
    integer                             ::  idx
    integer                             ::  i
    
    do i = 1,len(str)
      if (str(i:i) == c) then
        idx = i
        exit
      end if
    end do
  end function


end module fcs_module
