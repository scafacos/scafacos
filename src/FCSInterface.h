/*
  Copyright (C) 2011, 2012, 2013, 2014 Rene Halver, Michael Hofmann
  Copyright (C) 2016 Michael Hofmann

  This file is part of ScaFaCoS.

  ScaFaCoS is free software: you can redistribute it and/or modify
  it under the terms of the GNU Lesser Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  ScaFaCoS is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU Lesser Public License for more details.

  You should have received a copy of the GNU Lesser Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>. 
*/


/**
 * @file FCSInterface.h
 * @brief supplying the C-interface routine for the ScaFaCoS library
 * @author Rene Halver, Olaf Lenz
 */


#ifndef FCS_INTERFACE_INCLUDED
#define FCS_INTERFACE_INCLUDED

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>

#include "FCSInterface_p.h"
#include "FCSResult.h"

#ifdef FCS_ENABLE_DIRECT
#include "fcs_direct.h"
#endif
#ifdef FCS_ENABLE_EWALD
#include "fcs_ewald.h"
#endif
#ifdef FCS_ENABLE_FMM
#include "fcs_fmm.h"
#endif
#ifdef FCS_ENABLE_MEMD
#include "fcs_memd.h"
#endif
#ifdef FCS_ENABLE_MMM1D
#include "fcs_mmm1d.h"
#endif
#ifdef FCS_ENABLE_MMM2D
#include "fcs_mmm2d.h"
#endif
#ifdef FCS_ENABLE_P2NFFT
#include "fcs_p2nfft.h"
#endif
#ifdef FCS_ENABLE_P3M
#include "fcs_p3m.h"
#endif
#ifdef FCS_ENABLE_PEPC
#include "fcs_pepc.h"
#endif
#ifdef FCS_ENABLE_PP3MG
#include "fcs_pp3mg.h"
#endif
#ifdef FCS_ENABLE_VMG
#include "fcs_vmg.h"
#endif
#ifdef FCS_ENABLE_WOLF
#include "fcs_wolf.h"
#endif


#define FCS_MAX_METHOD_NAME_LENGTH  32


/* fallback definition, see "6.47 Function Names as Strings" in gcc-4.9 doc */
#if __STDC_VERSION__ < 199901L
# if __GNUC__ >= 2
#  define __func__ __FUNCTION__
# else
#  define __func__ "<unknown>"
# endif
#endif

#define CHECK_HANDLE_RETURN_RESULT(_h_, _f_) do { \
  if ((_h_) == FCS_NULL) \
    return fcs_result_create(FCS_ERROR_NULL_ARGUMENT, _f_, "null handle supplied"); \
  } while (0)

#define CHECK_HANDLE_RETURN_VAL(_h_, _f_, _v_) do { \
  if ((_h_) == FCS_NULL) { \
    fprintf(stderr, "%s: null handle supplied, returning " #_v_, _f_); \
    return (_v_); \
  } } while (0)

#define CHECK_METHOD_RETURN_RESULT(_h_, _f_, _m_, _n_) do { \
  if ((_h_)->method != (_m_)) \
    return fcs_result_create(FCS_ERROR_INCOMPATIBLE_METHOD, (_f_), "handle does not represent method \"" _n_ "\""); \
} while (0)

#define CHECK_METHOD_RETURN_VAL(_h_, _f_, _m_, _n_, _v_) do { \
  if ((_h_)->method != (_m_)) { \
    fprintf(stderr, "%s: handle does not represent method '" _n_ "', returning "  #_v_, _f_); \
    return (_v_); \
  } } while (0)

#define CHECK_RESULT_RETURN(_r_) do { \
    if ((_r_) != FCS_RESULT_SUCCESS) return (_r_); \
  } while (0)

#ifdef FCS_ENABLE_DEBUG
# define FCS_DEBUG_MOP(_mop_)  do { _mop_; } while (0)
#else
# define FCS_DEBUG_MOP(_mop_)  do { } while (0)
#endif

#define FCS_DEBUG_FUNC_INTRO(_f_)       FCS_DEBUG_MOP(printf("%s\n", _f_))
#define FCS_DEBUG_FUNC_OUTRO(_f_, _r_)  FCS_DEBUG_MOP(printf("%s: return: %s\n", _f_, ((_r_) == FCS_RESULT_SUCCESS)?"success":"failed"))


/*
 * @brief data structure that is used for storing the parameters of an FCS solver
 */
typedef struct _FCS_t
{
  /* numerical identifier of the solver method */
  fcs_int method;
  /* name of the solver method */
  char method_name[FCS_MAX_METHOD_NAME_LENGTH];

  /* MPI communicator to be used for the parallel execution */
  MPI_Comm communicator;

  /* dimensions of the system */
  fcs_int dimensions;

  /* the three base vectors of the system box */
  fcs_float box_a[3], box_b[3], box_c[3];
  /* origin vector of the system box */
  fcs_float box_origin[3];

  /* periodicity of the system in each dimension (value 0: open, value 1: periodic) */
  fcs_int periodicity[3];

  /* total number of particles in the system */
  fcs_int total_particles;
  /* maximum number of particles that can be stored in the local arrays provided to fcs_run */
  fcs_int max_local_particles;

  /* whether near-field computations should be performed
     0 = done by calling routine
     1 = done by library routine*/
  fcs_int near_field_flag;

  /* structures containing the method-specific parameters */
#ifdef FCS_ENABLE_DIRECT
  fcs_direct_parameters direct_param;
#endif
#ifdef FCS_ENABLE_FMM
  fcs_fmm_parameters fmm_param;
#endif
#ifdef FCS_ENABLE_MEMD
  fcs_memd_parameters memd_param;
#endif
#ifdef FCS_ENABLE_MMM1D
  fcs_mmm1d_parameters mmm1d_param;
#endif
#ifdef FCS_ENABLE_P2NFFT
  fcs_p2nfft_parameters p2nfft_param;
#endif
#ifdef FCS_ENABLE_PEPC
  fcs_pepc_parameters pepc_param;
#endif
#ifdef FCS_ENABLE_PP3MG
  fcs_pp3mg_parameters pp3mg_param;
#endif
#ifdef FCS_ENABLE_VMG
  fcs_vmg_parameters vmg_param;
#endif
#ifdef FCS_ENABLE_WOLF
  fcs_wolf_parameters wolf_param;
#endif

  /* current instance of the method */
  void *method_context;

  fcs_int values_changed;

  fcs_int shift_positions;

  /* functions and parameters set by the solvers */
  FCSResult (*destroy)(FCS handle);

  FCSResult (*set_tolerance)(FCS handle, fcs_int tolerance_type, fcs_float tolerance);
  FCSResult (*get_tolerance)(FCS handle, fcs_int *tolerance_type, fcs_float *tolerance);

  FCSResult (*set_r_cut)(FCS handle, fcs_float r_cut);
  FCSResult (*unset_r_cut)(FCS handle);
  FCSResult (*get_r_cut)(FCS handle, fcs_float *r_cut);

  FCSResult (*set_parameter)(FCS handle, fcs_bool continue_on_errors, char **current, char **next, fcs_int *matched);
  FCSResult (*print_parameters)(FCS handle);

  FCSResult (*tune)(FCS handle, fcs_int local_particles, fcs_float *positions, fcs_float *charges);
  FCSResult (*run)(FCS handle, fcs_int local_particles, fcs_float *positions, fcs_float *charges, fcs_float *field, fcs_float *potentials);

  FCSResult (*set_compute_virial)(FCS handle, fcs_int compute_virial);
  FCSResult (*get_compute_virial)(FCS handle, fcs_int *compute_virial);
  FCSResult (*get_virial)(FCS handle, fcs_float *virial);

  FCSResult (*set_max_particle_move)(FCS handle, fcs_float max_particle_move);
  FCSResult (*set_resort)(FCS handle, fcs_int resort);
  FCSResult (*get_resort)(FCS handle, fcs_int *resort);
  FCSResult (*get_resort_availability)(FCS handle, fcs_int *availability);
  FCSResult (*get_resort_particles)(FCS handle, fcs_int *resort_particles);
  FCSResult (*resort_ints)(FCS handle, fcs_int *src, fcs_int *dst, fcs_int n, MPI_Comm comm);
  FCSResult (*resort_floats)(FCS handle, fcs_float *src, fcs_float *dst, fcs_int n, MPI_Comm comm);
  FCSResult (*resort_bytes)(FCS handle, void *src, void *dst, fcs_int n, MPI_Comm comm);

} FCS_t;


/**
 * @brief function to set the method context information
 * @param handle FCS-object representing an FCS solver
 * @param pointer to the method context informations
 * @return FCSResult-object containing the return state
 */
FCSResult fcs_set_method_context(FCS handle, void *method_context);

/**
 * @brief function to return the method context information
 * @param handle FCS handle representing an FCS solver object
 * @return pointer to the method context informations
 */
void *fcs_get_method_context(FCS handle);

/**
 * @brief function to set whether parameter values of the FCS solver have changed
 * @param handle FCS-object representing an FCS solver
 * @param values_changed whether parameter values of the FCS solver have changed
 * @return FCSResult-object containing the return state
 */
FCSResult fcs_set_values_changed(FCS handle, fcs_int values_changed);

/**
 * @brief function to return whether parameter values of the FCS solver have changed
 * @param handle FCS-object representing an FCS solver
 * @return whether parameter values of the FCS solver have changed
 */
fcs_int fcs_get_values_changed(FCS handle);

/**
 * @brief Fortran wrapper function ot initialize an FCS solver method
 * @param handle FCS-object representing an FCS solver
 * @param method string for selecting the solver method
 * @param communicator MPI communicator to be used for parallel execution
 * @return FCSResult-object containing the return state
 */
FCSResult fcs_init_f(FCS *handle, const char *method_name, MPI_Fint communicator);


/**
 * tools for parameter parsing of fcs_set_parameters
 */
#if defined(FCS_ENABLE_DEBUG)
# define FCS_PARSE_PRINT_PARAM_BEGIN(_f_)     printf("%s: calling " #_f_ "(handle", __func__)
# define FCS_PARSE_PRINT_PARAM_VAL(_f_, _v_)  printf( _f_, _v_)
# define FCS_PARSE_PRINT_PARAM_STR(_str_)     printf("%s", (_str_))
# define FCS_PARSE_PRINT_PARAM_END(_r_)       printf(") -> %p\n", _r_)
#else
# define FCS_PARSE_PRINT_PARAM_BEGIN(_f_)     do {} while (0)
# define FCS_PARSE_PRINT_PARAM_VAL(_f_, _v_)  do {} while (0)
# define FCS_PARSE_PRINT_PARAM_STR(_str_)     do {} while (0)
# define FCS_PARSE_PRINT_PARAM_END(_r_)       do {} while (0)
#endif

typedef long long fcs_long_long_t;
typedef char *fcs_p_char_t;

#define FCS_PARSE_MAKE_TYPE_FUNC(_type_, _atox_, _format_) \
  static inline _type_ *parse_##_type_(char **s, _type_ *v) { \
    if (v == NULL) return NULL; \
    *v = (_type_) _atox_(*s); \
    *s = strchr(*s, ','); \
    if (*s) { **s = 0; *s += 1; } \
    FCS_PARSE_PRINT_PARAM_VAL(_format_, *v); \
    return v; \
  } \
  static inline _type_ *const_##_type_(_type_ c, _type_ *v) { \
    if (v == NULL) return NULL; \
    *v = c; \
    FCS_PARSE_PRINT_PARAM_VAL(_format_, *v); \
    return v; \
  }

static inline fcs_bool atob(const char *nptr)
{
  const char false_str[] = "false";
  if ((strlen(nptr) == 1 && strncmp(nptr, "0", 1) == 0) || (strlen(nptr) == strlen(false_str) && strncasecmp(nptr, false_str, strlen(false_str)))) return FCS_FALSE;
  return FCS_TRUE;
}

FCS_PARSE_MAKE_TYPE_FUNC(fcs_int, atoll, "%" FCS_LMOD_INT "d")
FCS_PARSE_MAKE_TYPE_FUNC(fcs_float, atof, "%" FCS_LMOD_FLOAT "f")
FCS_PARSE_MAKE_TYPE_FUNC(fcs_bool, atob, "%" FCS_LMOD_INT "d")
FCS_PARSE_MAKE_TYPE_FUNC(fcs_long_long_t, atoll, "%lld")
FCS_PARSE_MAKE_TYPE_FUNC(fcs_p_char_t, , "%s")

#define FCS_PARSE_SEQ_MAX  3

#define FCS_PARSE_PARAM_SELECTED(_str_, _param_) \
  (strcmp(param, #_param_) == 0 || strcmp(param, _str_) == 0)

#define FCS_PARSE_IF_PARAM_INTRO(_str_, _param_) \
  if (FCS_PARSE_PARAM_SELECTED(_str_, _param_)) { \
    FCSResult _r; \
    struct { \
      void *t; \
      fcs_int v_fcs_int[FCS_PARSE_SEQ_MAX]; \
      fcs_float v_fcs_float[FCS_PARSE_SEQ_MAX]; \
      fcs_bool v_fcs_bool[FCS_PARSE_SEQ_MAX]; \
      fcs_long_long_t v_fcs_long_long_t[FCS_PARSE_SEQ_MAX]; \
      fcs_p_char_t v_fcs_p_char_t[FCS_PARSE_SEQ_MAX]; \
    } _t; \
    char *_n=NULL; \
    FCS_PARSE_PRINT_PARAM_BEGIN(_param_);

#define FCS_PARSE_IF_PARAM_EXTRO() \
    FCS_PARSE_PRINT_PARAM_END(_r); \
    if (_r != FCS_RESULT_SUCCESS && FCS_IS_FALSE(continue_on_errors)) return _r; \
    goto next_param; \
  }

#define FCS_PARSE_VAL(_type_) \
  _t.v_##_type_; \
  if (cur) { \
    FCS_PARSE_PRINT_PARAM_STR(", "); \
    parse_##_type_(&cur, &_t.v_##_type_[0]); \
    _n = cur; cur = NULL; \
  } _type_

#define FCS_PARSE_SEQ(_type_, _n_) \
  &_t.v_##_type_; \
  if (cur) { \
    FCS_PARSE_PRINT_PARAM_STR(", ["); \
    for (int _i = 0; _i < (_n_) && _i < FCS_PARSE_SEQ_MAX; ++_i) { \
      FCS_PARSE_PRINT_PARAM_STR((_i == 0)?"":", "); \
      parse_##_type_(&cur, &_t.v_##_type_[_i]); \
    } \
    _n = cur; cur = NULL; \
    FCS_PARSE_PRINT_PARAM_STR("]"); \
  } _type_ *

#define FCS_PARSE_CONST(_type_, _c_) \
  _t.v_##_type_; \
  if (cur) { \
    FCS_PARSE_PRINT_PARAM_STR(", "); \
    const_##_type_(_c_, &_t.v_##_type_[0]); \
    _n = cur; cur = NULL; \
  } _type_

#define FCS_PARSE_IF_PARAM_THEN_FUNC0_GOTO_NEXT(_str_, _param_) \
  FCS_PARSE_IF_PARAM_INTRO(_str_, _param_) \
    _r = fcs_##_param_(handle); \
  FCS_PARSE_IF_PARAM_EXTRO()

#define FCS_PARSE_IF_PARAM_THEN_FUNC1_GOTO_NEXT(_str_, _param_, _p0_) \
  FCS_PARSE_IF_PARAM_INTRO(_str_, _param_) \
    _t.t = _p0_ _v0 = *_p0_ _vv0 = _v0; cur = _n; \
    _r = fcs_##_param_(handle, _vv0); \
  FCS_PARSE_IF_PARAM_EXTRO()

#define FCS_PARSE_IF_PARAM_THEN_FUNC2_GOTO_NEXT(_str_, _param_, _p0_, _p1_) \
  FCS_PARSE_IF_PARAM_INTRO(_str_, _param_) \
    _t.t = _p0_ _v0 = *_p0_ _vv0 = _v0; cur = _n; \
    _t.t = _p1_ _v1 = *_p1_ _vv1 = _v1; cur = _n; \
    _r = fcs_##_param_(handle, _vv0, _vv1); \
  FCS_PARSE_IF_PARAM_EXTRO()

#define FCS_PARSE_IF_PARAM_THEN_FUNC3_GOTO_NEXT(_str_, _param_, _p0_, _p1_, _p2_) \
  FCS_PARSE_IF_PARAM_INTRO(_str_, _param_) \
    _t.t = _p0_ _v0 = *_p0_ _vv0 = _v0; cur = _n; \
    _t.t = _p1_ _v1 = *_p1_ _vv1 = _v1; cur = _n; \
    _t.t = _p2_ _v2 = *_p2_ _vv2 = _v2; cur = _n; \
    _r = fcs_##_param_(handle, _vv0, _vv1, _vv2); \
  FCS_PARSE_IF_PARAM_EXTRO()

#define FCS_PARSE_DUMMY() \
  if (FCS_PARSE_PARAM_SELECTED("EVEN_IF_IT_MATCHES_IT_DOES_NOTHING", EVEN_IF_IT_MATCHES_IT_DOES_NOTHING)) { \
    parse_fcs_int(NULL, NULL); \
    const_fcs_int(0, NULL); \
    parse_fcs_float(NULL, NULL); \
    const_fcs_float(0, NULL); \
    parse_fcs_bool(NULL, NULL); \
    const_fcs_bool(0, NULL); \
    parse_fcs_long_long_t(NULL, NULL); \
    const_fcs_long_long_t(0, NULL); \
    parse_fcs_p_char_t(NULL, NULL); \
    const_fcs_p_char_t(NULL, NULL); \
    goto next_param; \
  }


#endif
