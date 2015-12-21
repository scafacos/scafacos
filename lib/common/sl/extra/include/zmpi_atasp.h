/*
 *  Copyright (C) 2011, 2012, 2013, 2014, 2015 Michael Hofmann
 *  
 *  This file is part of ScaFaCoS.
 *  
 *  ScaFaCoS is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *  
 *  ScaFaCoS is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser Public License for more details.
 *  
 *  You should have received a copy of the GNU Lesser Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *  

 *  
 *  SL - Sorting Library, michael <dot> hofmann <at> informatik <dot> tu-chemnitz <dot> de
 */


#ifndef __ZMPI_ATASP_H__
#define __ZMPI_ATASP_H__





#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#ifdef HAVE_ZMPI_LOCAL_H
# include "zmpi_local.h"
#else
# error zmpi-local is required for local copying of mpi-typed data
#endif


typedef int zmpi_spec_int_t;

typedef int zmpi_spec_proc_t;

#define zmpi_SPEC_PROC_NONE  MPI_PROC_NULL

typedef void *zmpi_spec_tproc_data_t;

typedef void *zmpi_spec_elem_buf_t;

typedef struct
{
  void *buf;
  int count, max_count;
  MPI_Datatype mpi_type;

#ifdef HAVE_ZMPI_LOCAL_H
  zmpil_t zmpil_type;
#endif

} zmpi_spec_elem_t;

#if MPI_VERSION >= 3
typedef MPI_Count zmpi_spec_elem_index_t;
#else
typedef long zmpi_spec_elem_index_t;
#endif

#define zmpi_spec_elem_unset(_e_)          Z_MOP((_e_)->buf = NULL; (_e_)->count = (_e_)->max_count = 0; (_e_)->mpi_type = MPI_DATATYPE_NULL;)

#define zmpi_spec_elem_set_n(_e_, _n_)     (_e_)->count = (_n_)
#define zmpi_spec_elem_get_n(_e_)          (_e_)->count
#define zmpi_spec_elem_set_nmax(_e_, _n_)  (_e_)->max_count = (_n_)
#define zmpi_spec_elem_get_nmax(_e_)       (_e_)->max_count

#define zmpi_spec_elem_set_buf(_e_, _b_)   (_e_)->buf = (_b_)
#define zmpi_spec_elem_get_buf(_e_)        (_e_)->buf

#define zmpi_spec_elem_x_get_buf(_e_, _x_)  (_e_)->buf

#ifdef HAVE_ZMPI_LOCAL_H
# define zmpi_spec_elem_copy_at(_se_, _sat_, _de_, _dat_)  zmpil_memcpy_at((_de_)->buf, (_dat_), (_se_)->buf, (_sat_), 1, &(_se_)->zmpil_type)
#endif

#define zmpi_spec_elem_exchange_at(_s0_, _s0at_, _s1_, _s1at_, _t_) do { \
    zmpi_spec_elem_copy_at((_s0_), (_s0at_), (_t_), 0); \
    zmpi_spec_elem_copy_at((_s1_), (_s1at_), (_s0_), (_s0at_)); \
    zmpi_spec_elem_copy_at((_t_), 0, (_s1_), (_s1at_)); \
  } while (0)

#ifdef HAVE_ZMPI_LOCAL_H

# define ZMPI_FIXSIZE_DECLARE(_fx_, _fxp_)
# define ZMPI_FIXSIZE_CREATE(_fx_, _fxp_)
# define ZMPI_FIXSIZE_COPY_AT(_se_, _sat_, _de_, _dat_, _fx_, _fxp_)             memcpy(((char *) (_de_)->buf) + ((_dat_) * (_fxp_)), ((char *) (_se_)->buf) + ((_sat_) * (_fxp_)), (_fxp_))
# define ZMPI_FIXSIZE_EXCHANGE_AT(_s0_, _s0at_, _s1_, _s1at_, _t_, _fx_, _fxp_)  do { \
    ZMPI_FIXSIZE_COPY_AT((_s0_), (_s0at_), (_t_), 0, _fx_, _fxp_); \
    ZMPI_FIXSIZE_COPY_AT((_s1_), (_s1at_), (_s0_), (_s0at_), _fx_, _fxp_); \
    ZMPI_FIXSIZE_COPY_AT((_t_), 0, (_s1_), (_s1at_), _fx_, _fxp_); \
  } while (0)
# define ZMPI_FIXSIZE_DESTROY(_fx_, _fxp_)

# define ZMPI_FIXTYPE_DECLARE(_fx_, _fxp_)
# define ZMPI_FIXTYPE_CREATE(_fx_, _fxp_)
# define ZMPI_FIXTYPE_COPY_AT(_se_, _sat_, _de_, _dat_, _fx_, _fxp_)             ((_fxp_ *) (_de_)->buf)[_dat_] = ((_fxp_ *) (_se_)->buf)[_sat_]
# define ZMPI_FIXTYPE_EXCHANGE_AT(_s0_, _s0at_, _s1_, _s1at_, _t_, _fx_, _fxp_)  do { \
    ZMPI_FIXTYPE_COPY_AT((_s0_), (_s0at_), (_t_), 0, _fx_, _fxp_); \
    ZMPI_FIXTYPE_COPY_AT((_s1_), (_s1at_), (_s0_), (_s0at_), _fx_, _fxp_); \
    ZMPI_FIXTYPE_COPY_AT((_t_), 0, (_s1_), (_s1at_), _fx_, _fxp_); \
  } while (0)
# define ZMPI_FIXTYPE_DESTROY(_fx_, _fxp_)

#endif





/* tproc count */

/* sp_macro zmpi_SPEC_DECLARE_TPROC_COUNT_DB */
#define zmpi_SPEC_DECLARE_TPROC_COUNT_DB \
  struct { zmpi_spec_elem_index_t n, i; zmpi_spec_proc_t p; } spec0cd;

/* sp_macro zmpi_SPEC_DO_TPROC_COUNT_DB */
#define zmpi_SPEC_DO_TPROC_COUNT_DB(_tp_, _tpd_, _b_, _cs_)  do { \
  spec0cd.n = zmpi_spec_elem_get_n(_b_); \
  for (spec0cd.i = 0; spec0cd.i < spec0cd.n; ++spec0cd.i) { \
    spec0cd.p = (_tp_)(zmpi_spec_elem_get_buf(_b_), spec0cd.i, _tpd_); \
    if (spec0cd.p == zmpi_SPEC_PROC_NONE) continue; \
    ++(_cs_)[spec0cd.p]; \
  } } while (0)

/* sp_macro zmpi_SPEC_FUNC_TPROC_COUNT_DB */
#define zmpi_SPEC_FUNC_TPROC_COUNT_DB(_name_, _tp_, _s_...) \
_s_ void _name_##_tproc_count_db(zmpi_spec_tproc_data_t tproc_data, zmpi_spec_elem_t *s, int *counts) \
{ \
  zmpi_SPEC_DECLARE_TPROC_COUNT_DB \
  zmpi_SPEC_DO_TPROC_COUNT_DB(_tp_, tproc_data, s, counts); \
}

/* sp_macro zmpi_SPEC_DECLARE_TPROC_COUNT_IP */
#define zmpi_SPEC_DECLARE_TPROC_COUNT_IP \
  struct { zmpi_spec_elem_index_t n, t, i; zmpi_spec_proc_t p; } spec0ci;

/* sp_macro zmpi_SPEC_DO_TPROC_COUNT_IP */
#define zmpi_SPEC_DO_TPROC_COUNT_IP(_tp_, _tpd_, _b_, _cs_)  do { \
  spec0ci.n = zmpi_spec_elem_get_n(_b_); \
  spec0ci.t = 0; \
  for (spec0ci.i = 0; spec0ci.i < spec0ci.n; ++spec0ci.i) { \
    spec0ci.p = (_tp_)(zmpi_spec_elem_get_buf(_b_), spec0ci.i, _tpd_); \
    if (spec0ci.p == zmpi_SPEC_PROC_NONE) continue; \
    ++(_cs_)[spec0ci.p]; \
    if (spec0ci.t < spec0ci.i) zmpi_spec_elem_copy_at((_b_), spec0ci.i, (_b_), spec0ci.t); \
    ++spec0ci.t; \
  } \
  zmpi_spec_elem_set_n(_b_, spec0ci.t); \
} while (0)

/* sp_macro zmpi_SPEC_FUNC_TPROC_COUNT_IP */
#define zmpi_SPEC_FUNC_TPROC_COUNT_IP(_name_, _tp_, _s_...) \
_s_ void _name_##_tproc_count_ip(zmpi_spec_tproc_data_t tproc_data, zmpi_spec_elem_t *s, int *counts) \
{ \
  zmpi_SPEC_DECLARE_TPROC_COUNT_IP \
  zmpi_SPEC_DO_TPROC_COUNT_IP(_tp_, tproc_data, s, counts); \
}


/* tproc_mod count */

/* sp_macro zmpi_SPEC_DECLARE_TPROC_MOD_COUNT_DB */
#define zmpi_SPEC_DECLARE_TPROC_MOD_COUNT_DB \
  struct { zmpi_spec_elem_index_t n, i; zmpi_spec_proc_t p; } spec1cd;

/* sp_macro zmpi_SPEC_DO_TPROC_MOD_COUNT_DB */
#define zmpi_SPEC_DO_TPROC_MOD_COUNT_DB(_tp_, _tpd_, _b_, _cs_)  do { \
  spec1cd.n = zmpi_spec_elem_get_n(_b_); \
  for (spec1cd.i = 0; spec1cd.i < spec1cd.n; ++spec1cd.i) { \
    spec1cd.p = (_tp_)(zmpi_spec_elem_get_buf(_b_), spec1cd.i, _tpd_, NULL); \
    if (spec1cd.p == zmpi_SPEC_PROC_NONE) continue; \
    ++(_cs_)[spec1cd.p]; \
  } } while (0)

/* sp_macro zmpi_SPEC_FUNC_TPROC_MOD_COUNT_DB */
#define zmpi_SPEC_FUNC_TPROC_MOD_COUNT_DB(_name_, _tp_, _s_...) \
_s_ void _name_##_tproc_mod_count_db(zmpi_spec_tproc_data_t tproc_data, zmpi_spec_elem_t *s, int *counts) \
{ \
  zmpi_SPEC_DECLARE_TPROC_MOD_COUNT_DB \
  zmpi_SPEC_DO_TPROC_MOD_COUNT_DB(_tp_, tproc_data, s, counts); \
}

/* sp_macro zmpi_SPEC_DECLARE_TPROC_MOD_COUNT_IP */
#define zmpi_SPEC_DECLARE_TPROC_MOD_COUNT_IP \
  struct { zmpi_spec_elem_index_t n, t, i; zmpi_spec_proc_t p; } spec1ci;

/* sp_macro zmpi_SPEC_DO_TPROC_MOD_COUNT_IP */
#define zmpi_SPEC_DO_TPROC_MOD_COUNT_IP(_tp_, _tpd_, _b_, _cs_)  do { \
  spec1ci.n = zmpi_spec_elem_get_n(_b_); \
  spec1ci.t = 0; \
  for (spec1ci.i = 0; spec1ci.i < spec1ci.n; ++spec1ci.i) { \
    spec1ci.p = (_tp_)(zmpi_spec_elem_get_buf(_b_), spec1ci.i, _tpd_, NULL); \
    if (spec1ci.p == zmpi_SPEC_PROC_NONE) continue; \
    ++(_cs_)[spec1ci.p]; \
    if (spec1ci.t < spec1ci.i) zmpi_spec_elem_copy_at((_b_), spec1ci.i, (_b_), spec1ci.t); \
    ++spec1ci.t; \
  } \
  zmpi_spec_elem_set_n(_b_, spec1ci.t); \
} while (0)

/* sp_macro zmpi_SPEC_FUNC_TPROC_MOD_COUNT_IP */
#define zmpi_SPEC_FUNC_TPROC_MOD_COUNT_IP(_name_, _tp_, _s_...) \
_s_ void _name_##_tproc_mod_count_ip(zmpi_spec_tproc_data_t tproc_data, zmpi_spec_elem_t *s, int *counts) \
{ \
  zmpi_SPEC_DECLARE_TPROC_MOD_COUNT_IP \
  zmpi_SPEC_DO_TPROC_MOD_COUNT_IP(_tp_, tproc_data, s, counts); \
}


/* tprocs count */

/* sp_macro zmpi_SPEC_DECLARE_TPROCS_COUNT_DB */
#define zmpi_SPEC_DECLARE_TPROCS_COUNT_DB \
  struct { zmpi_spec_elem_index_t n, i; zmpi_spec_int_t j, m; } spec2cd;

/* sp_macro zmpi_SPEC_DO_TPROCS_COUNT_DB */
#define zmpi_SPEC_DO_TPROCS_COUNT_DB(_tp_, _tpd_, _b_, _cs_, _ps_)  do { \
  spec2cd.n = zmpi_spec_elem_get_n(_b_); \
  for (spec2cd.i = 0; spec2cd.i < spec2cd.n; ++spec2cd.i) { \
    (_tp_)(zmpi_spec_elem_get_buf(_b_), spec2cd.i, (_tpd_), &spec2cd.m, (_ps_)); \
    for (spec2cd.j = 0; spec2cd.j < spec2cd.m; ++spec2cd.j) ++(_cs_)[(_ps_)[spec2cd.j]]; \
  } } while (0)

/* sp_macro zmpi_SPEC_FUNC_TPROCS_COUNT_DB */
#define zmpi_SPEC_FUNC_TPROCS_COUNT_DB(_name_, _tp_, _s_...) \
_s_ void _name_##_tprocs_count_db(zmpi_spec_tproc_data_t tproc_data, zmpi_spec_elem_t *s, int *counts, zmpi_spec_proc_t *procs) \
{ \
  zmpi_SPEC_DECLARE_TPROCS_COUNT_DB \
  zmpi_SPEC_DO_TPROCS_COUNT_DB(_tp_, tproc_data, s, counts, procs); \
}

/* sp_macro zmpi_SPEC_DECLARE_TPROCS_COUNT_IP */
#define zmpi_SPEC_DECLARE_TPROCS_COUNT_IP \
  struct { zmpi_spec_elem_index_t n, t, i; zmpi_spec_int_t j, m; } spec2ci;

/* sp_macro zmpi_SPEC_DO_TPROCS_COUNT_IP */
#define zmpi_SPEC_DO_TPROCS_COUNT_IP(_tp_, _tpd_, _b_, _cs_, _ps_)  do { \
  spec2ci.n = zmpi_spec_elem_get_n(_b_); \
  spec2ci.t = 0; \
  for (spec2ci.i = 0; spec2ci.i < spec2ci.n; ++spec2ci.i) { \
    (_tp_)(zmpi_spec_elem_get_buf(_b_), spec2ci.i, (_tpd_), &spec2ci.m, (_ps_)); \
    if (spec2ci.m <= 0) continue; \
    for (spec2ci.j = 0; spec2ci.j < spec2ci.m; ++spec2ci.j) ++(_cs_)[(_ps_)[spec2ci.j]]; \
    if (spec2ci.t < spec2ci.i) zmpi_spec_elem_copy_at((_b_), spec2ci.i, (_b_), spec2ci.t); \
    ++spec2ci.t; \
  } \
  zmpi_spec_elem_set_n(_b_, spec2ci.t); \
} while (0)

/* sp_macro zmpi_SPEC_FUNC_TPROCS_COUNT_IP */
#define zmpi_SPEC_FUNC_TPROCS_COUNT_IP(_name_, _tp_, _s_...) \
_s_ void _name_##_tprocs_count_ip(zmpi_spec_tproc_data_t tproc_data, zmpi_spec_elem_t *s, int *counts, zmpi_spec_proc_t *procs) \
{ \
  zmpi_SPEC_DECLARE_TPROCS_COUNT_IP \
  zmpi_SPEC_DO_TPROCS_COUNT_IP(_tp_, tproc_data, s, counts, procs); \
}


/* tprocs_mod count */

/* sp_macro zmpi_SPEC_DECLARE_TPROCS_MOD_COUNT_DB */
#define zmpi_SPEC_DECLARE_TPROCS_MOD_COUNT_DB \
  struct { zmpi_spec_elem_index_t n, i; zmpi_spec_int_t j, m; } spec3cd;

/* sp_macro zmpi_SPEC_DO_TPROCS_MOD_COUNT_DB */
#define zmpi_SPEC_DO_TPROCS_MOD_COUNT_DB(_tp_, _tpd_, _b_, _cs_, _ps_)  do { \
  spec3cd.n = zmpi_spec_elem_get_n(_b_); \
  for (spec3cd.i = 0; spec3cd.i < spec3cd.n; ++spec3cd.i) { \
    (_tp_)(zmpi_spec_elem_get_buf(_b_), spec3cd.i, (_tpd_), &spec3cd.m, (_ps_), NULL); \
    for (spec3cd.j = 0; spec3cd.j < spec3cd.m; ++spec3cd.j) ++(_cs_)[(_ps_)[spec3cd.j]]; \
  } } while (0)

/* sp_macro zmpi_SPEC_FUNC_TPROCS_MOD_COUNT_DB */
#define zmpi_SPEC_FUNC_TPROCS_MOD_COUNT_DB(_name_, _tp_, _s_...) \
_s_ void _name_##_tprocs_mod_count_db(zmpi_spec_tproc_data_t tproc_data, zmpi_spec_elem_t *s, int *counts, zmpi_spec_proc_t *procs) \
{ \
  zmpi_SPEC_DECLARE_TPROCS_MOD_COUNT_DB \
  zmpi_SPEC_DO_TPROCS_MOD_COUNT_DB(_tp_, tproc_data, s, counts, procs); \
}

/* sp_macro zmpi_SPEC_DECLARE_TPROCS_MOD_COUNT_IP */
#define zmpi_SPEC_DECLARE_TPROCS_MOD_COUNT_IP \
  struct { zmpi_spec_elem_index_t n, t, i; zmpi_spec_int_t j, m; } spec3ci;

/* sp_macro zmpi_SPEC_DO_TPROCS_MOD_COUNT_IP */
#define zmpi_SPEC_DO_TPROCS_MOD_COUNT_IP(_tp_, _tpd_, _b_, _cs_, _ps_)  do { \
  spec3ci.n = zmpi_spec_elem_get_n(_b_); \
  spec3ci.t = 0; \
  for (spec3ci.i = 0; spec3ci.i < spec3ci.n; ++spec3ci.i) { \
    (_tp_)(zmpi_spec_elem_get_buf(_b_), spec3ci.i, (_tpd_), &spec3ci.m, (_ps_), NULL); \
    if (spec3ci.m <= 0) continue; \
    for (spec3ci.j = 0; spec3ci.j < spec3ci.m; ++spec3ci.j) ++(_cs_)[(_ps_)[spec3ci.j]]; \
    if (spec3ci.t < spec3ci.i) zmpi_spec_elem_copy_at((_b_), spec3ci.i, (_b_), spec3ci.t); \
    ++spec3ci.t; \
  } \
  zmpi_spec_elem_set_n(_b_, spec3ci.t); \
} while (0)

/* sp_macro zmpi_SPEC_FUNC_TPROCS_MOD_COUNT_IP */
#define zmpi_SPEC_FUNC_TPROCS_MOD_COUNT_IP(_name_, _tp_, _s_...) \
_s_ void _name_##_tprocs_mod_count_ip(zmpi_spec_tproc_data_t tproc_data, zmpi_spec_elem_t *s, int *counts, zmpi_spec_proc_t *procs) \
{ \
  zmpi_SPEC_DECLARE_TPROCS_MOD_COUNT_IP \
  zmpi_SPEC_DO_TPROCS_MOD_COUNT_IP(_tp_, tproc_data, s, counts, procs); \
}


/* un-fixed macros, sp_macro zmpi_spec_fixed_default_declare zmpi_spec_fixed_default_create zmpi_spec_fixed_default_copy_at zmpi_spec_fixed_default_exchange_at zmpi_spec_fixed_default_destroy */
#define zmpi_spec_fixed_default_declare(_fx_, _fxp_)
#define zmpi_spec_fixed_default_create(_fx_, _fxp_)
#define zmpi_spec_fixed_default_copy_at(_se_, _sat_, _de_, _dat_, _fx_, _fxp_)             zmpi_spec_elem_copy_at(_se_, _sat_, _de_, _dat_)
#define zmpi_spec_fixed_default_exchange_at(_s0_, _s0at_, _s1_, _s1at_, _t_, _fx_, _fxp_)  zmpi_spec_elem_exchange_at(_s0_, _s0at_, _s1_, _s1at_, _t_)
#define zmpi_spec_fixed_default_destroy(_fx_, _fxp_)


/* tproc rearrange */

/* sp_macro zmpi_SPEC_DECLARE_FIXED_TPROC_REARRANGE_DB */
#define zmpi_SPEC_DECLARE_FIXED_TPROC_REARRANGE_DB(_fxdcl_, _fxp_) \
  struct { zmpi_spec_elem_index_t n, i; zmpi_spec_proc_t p; _fxdcl_(fx, _fxp_) } spec0d;

/* sp_macro zmpi_SPEC_DO_FIXED_TPROC_REARRANGE_DB */
#define zmpi_SPEC_DO_FIXED_TPROC_REARRANGE_DB(_fxp_, _fxc_, _fxca_, _fxxa_, _fxd_, _tp_, _tpd_, _sb_, _db_, _ds_)  do { \
  _fxc_(spec0d.fx, _fxp_); \
  spec0d.n = zmpi_spec_elem_get_n(_sb_); \
  for (spec0d.i = 0; spec0d.i < spec0d.n; ++spec0d.i) { \
    spec0d.p = (_tp_)(zmpi_spec_elem_get_buf(_sb_), spec0d.i, _tpd_); \
    if (spec0d.p == zmpi_SPEC_PROC_NONE) continue; \
    _fxca_((_sb_), spec0d.i, (_db_), (_ds_)[spec0d.p], spec0d.fx, _fxp_); \
    ++(_ds_)[spec0d.p]; \
  } \
  _fxd_(spec0d.fx, _fxp_); \
  } while (0)

/* sp_macro zmpi_SPEC_FUNC_FIXED_TPROC_REARRANGE_DB */
#define zmpi_SPEC_FUNC_FIXED_TPROC_REARRANGE_DB(_name_, _fxdcl_, _fxp_, _fxc_, _fxca_, _fxxa_, _fxd_, _tp_, _s_...) \
_s_ void _name_##_tproc_rearrange_db(zmpi_spec_tproc_data_t tproc_data, zmpi_spec_elem_t *s, zmpi_spec_elem_t *d, int *displs) \
{ \
  zmpi_SPEC_DECLARE_FIXED_TPROC_REARRANGE_DB(_fxdcl_, _fxp_) \
  zmpi_SPEC_DO_FIXED_TPROC_REARRANGE_DB(_fxp_, _fxc_, _fxca_, _fxxa_, _fxd_, _tp_, tproc_data, s, d, displs); \
}

/* sp_macro zmpi_SPEC_DECLARE_TPROC_REARRANGE_DB */
#define zmpi_SPEC_DECLARE_TPROC_REARRANGE_DB \
  zmpi_SPEC_DECLARE_FIXED_TPROC_REARRANGE_DB(zmpi_spec_fixed_default_declare, NOPARAM)

/* sp_macro zmpi_SPEC_DO_TPROC_REARRANGE_DB */
#define zmpi_SPEC_DO_TPROC_REARRANGE_DB(_tp_, _tpd_, _sb_, _db_, _ds_) \
  zmpi_SPEC_DO_FIXED_TPROC_REARRANGE_DB(NOPARAM, zmpi_spec_fixed_default_create, zmpi_spec_fixed_default_copy_at, zmpi_spec_fixed_default_exchange_at, zmpi_spec_fixed_default_destroy, _tp_, _tpd_, _sb_, _db_, _ds_)

/* sp_macro zmpi_SPEC_FUNC_TPROC_REARRANGE_DB */
#define zmpi_SPEC_FUNC_TPROC_REARRANGE_DB(_name_, _tp_, _s_...) \
  zmpi_SPEC_FUNC_FIXED_TPROC_REARRANGE_DB(_name_, zmpi_spec_fixed_default_declare, NOPARAM, zmpi_spec_fixed_default_create, zmpi_spec_fixed_default_copy_at, zmpi_spec_fixed_default_exchange_at, zmpi_spec_fixed_default_destroy, _tp_, _s_)

/* sp_macro zmpi_SPEC_DECLARE_TPROC_REARRANGE_IP */
#define zmpi_SPEC_DECLARE_TPROC_REARRANGE_IP \
  struct { zmpi_spec_elem_index_t e, i, j; zmpi_spec_proc_t p, np; } spec0i;

/* sp_macro zmpi_SPEC_DO_TPROC_REARRANGE_IP */
#define zmpi_SPEC_DO_TPROC_REARRANGE_IP(_tp_, _tpd_, _b_, _xb_, _ds_, _cs_, _n_)  do { \
  for (spec0i.e = 0, spec0i.i = 0; spec0i.i < (_n_); ++spec0i.i) { \
    spec0i.e += (_cs_)[spec0i.i]; \
    spec0i.j = (_ds_)[spec0i.i]; \
    while (spec0i.j < spec0i.e) { \
      spec0i.p = (_tp_)(zmpi_spec_elem_get_buf(_b_), spec0i.j, _tpd_); \
      while (spec0i.p != spec0i.i) { \
        spec0i.np = (_tp_)(zmpi_spec_elem_get_buf(_b_), (_ds_)[spec0i.p], _tpd_); \
        if (spec0i.np != spec0i.p) zmpi_spec_elem_exchange_at((_b_), (_ds_)[spec0i.p], (_b_), spec0i.j, (_xb_)); \
        ++(_ds_)[spec0i.p]; \
        spec0i.p = spec0i.np; \
      } \
      ++spec0i.j; \
    } \
  } } while (0)

/* sp_macro zmpi_SPEC_FUNC_TPROC_REARRANGE_IP */
#define zmpi_SPEC_FUNC_TPROC_REARRANGE_IP(_name_, _tp_, _s_) \
_s_ void _name_##_tproc_rearrange_ip(zmpi_spec_tproc_data_t tproc_data, zmpi_spec_elem_t *s, zmpi_spec_elem_t *x, int *displs, int *counts, zmpi_spec_int_t n) \
{ \
  zmpi_SPEC_DECLARE_TPROC_REARRANGE_IP \
  zmpi_SPEC_DO_TPROC_REARRANGE_IP(_tp_, tproc_data, s, x, displs, counts, n); \
}


/* tproc_mod rearrange */

/* sp_macro zmpi_SPEC_DECLARE_FIXED_TPROC_MOD_REARRANGE_DB */
#define zmpi_SPEC_DECLARE_FIXED_TPROC_MOD_REARRANGE_DB(_fxdcl_, _fxp_) \
  struct { zmpi_spec_elem_index_t n, i; zmpi_spec_proc_t p; _fxdcl_(fx, _fxp_) } spec1d;

/* sp_macro zmpi_SPEC_DO_FIXED_TPROC_MOD_REARRANGE_DB */
#define zmpi_SPEC_DO_FIXED_TPROC_MOD_REARRANGE_DB(_fxp_, _fxc_, _fxca_, _fxxa_, _fxd_, _tp_, _tpd_, _sb_, _db_, _ds_, _ib_)  do { \
  spec1d.n = zmpi_spec_elem_get_n(_sb_); \
  _fxc_(spec0d.fx, _fxp_); \
  for (spec1d.i = 0; spec1d.i < spec1d.n; ++spec1d.i) { \
    spec1d.p = (_tp_)(zmpi_spec_elem_get_buf(_sb_), spec1d.i, _tpd_, zmpi_spec_elem_get_buf(_ib_)); \
    if (spec1d.p == zmpi_SPEC_PROC_NONE) continue; \
    _fxca_((_ib_), 0, (_db_), (_ds_)[spec1d.p], spec1d.fx, _fxp_); \
    ++(_ds_)[spec1d.p]; \
  } \
  _fxd_(spec0d.fx, _fxp_); \
  } while (0)

/* sp_macro zmpi_SPEC_FUNC_FIXED_TPROC_MOD_REARRANGE_DB */
#define zmpi_SPEC_FUNC_FIXED_TPROC_MOD_REARRANGE_DB(_name_, _fxdcl_, _fxp_, _fxc_, _fxca_, _fxxa_, _fxd_, _tp_, _s_...) \
_s_ void _name_##_tproc_mod_rearrange_db(zmpi_spec_tproc_data_t tproc_data, zmpi_spec_elem_t *s, zmpi_spec_elem_t *d, int *displs, zmpi_spec_elem_t *mod) \
{ \
  zmpi_SPEC_DECLARE_FIXED_TPROC_MOD_REARRANGE_DB(_fxdcl_, _fxp_) \
  zmpi_SPEC_DO_FIXED_TPROC_MOD_REARRANGE_DB(_fxp_, _fxc_, _fxca_, _fxxa_, _fxd_, _tp_, tproc_data, s, d, displs, mod); \
}

/* sp_macro zmpi_SPEC_DECLARE_TPROC_MOD_REARRANGE_DB */
#define zmpi_SPEC_DECLARE_TPROC_MOD_REARRANGE_DB \
  zmpi_SPEC_DECLARE_FIXED_TPROC_MOD_REARRANGE_DB(zmpi_spec_fixed_default_declare, NOPARAM)

/* sp_macro zmpi_SPEC_DO_TPROC_MOD_REARRANGE_DB */
#define zmpi_SPEC_DO_TPROC_MOD_REARRANGE_DB(_tp_, _tpd_, _sb_, _db_, _ds_, _ib_) \
  zmpi_SPEC_DO_FIXED_TPROC_MOD_REARRANGE_DB(NOPARAM, zmpi_spec_fixed_default_create, zmpi_spec_fixed_default_copy_at, zmpi_spec_fixed_default_exchange_at, zmpi_spec_fixed_default_destroy, _tp_, _tpd_, _sb_, _db_, _ds_, _ib_)

/* sp_macro zmpi_SPEC_FUNC_TPROC_MOD_REARRANGE_DB */
#define zmpi_SPEC_FUNC_TPROC_MOD_REARRANGE_DB(_name_, _tp_, _s_...) \
  zmpi_SPEC_FUNC_FIXED_TPROC_MOD_REARRANGE_DB(_name_, zmpi_spec_fixed_default_declare, NOPARAM, zmpi_spec_fixed_default_create, zmpi_spec_fixed_default_copy_at, zmpi_spec_fixed_default_exchange_at, zmpi_spec_fixed_default_destroy, _tp_, _s_)

/* sp_macro zmpi_SPEC_DECLARE_TPROC_MOD_REARRANGE_IP */
#define zmpi_SPEC_DECLARE_TPROC_MOD_REARRANGE_IP \
  struct { zmpi_spec_elem_index_t e, i, j; zmpi_spec_proc_t p, np; } spec1i;

/* sp_macro zmpi_SPEC_DO_TPROC_MOD_REARRANGE_IP */
#define zmpi_SPEC_DO_TPROC_MOD_REARRANGE_IP(_tp_, _tpd_, _b_, _xb_, _ds_, _cs_, _n_, _ib_)  do { \
  for (spec1i.e = 0, spec1i.i = 0; spec1i.i < (_n_); ++spec1i.i) { \
    spec1i.e += (_cs_)[spec1i.i]; \
    spec1i.j = (_ds_)[spec1i.i]; \
    while (spec1i.j < spec1i.e) { \
      spec1i.p = (_tp_)(zmpi_spec_elem_get_buf(_b_), spec1i.j, _tpd_, zmpi_spec_elem_get_buf(_ib_)); \
      zmpi_spec_elem_copy_at((_ib_), 0, (_b_), spec1i.j); \
      while (spec1i.p != spec1i.i) { \
        spec1i.np = (_tp_)(zmpi_spec_elem_get_buf(_b_), (_ds_)[spec1i.p], _tpd_, zmpi_spec_elem_get_buf(_ib_)); \
        if (spec1i.np != spec1i.p) { \
          zmpi_spec_elem_copy_at((_b_), spec1i.j, (_b_), (_ds_)[spec1i.p]); \
          zmpi_spec_elem_copy_at((_ib_), 0, (_b_), spec1i.j); \
        } else zmpi_spec_elem_copy_at((_ib_), 0, (_b_), (_ds_)[spec1i.p]); \
        ++(_ds_)[spec1i.p]; \
        spec1i.p = spec1i.np; \
      } \
      ++spec1i.j; \
    } \
  } } while (0)

/* sp_macro zmpi_SPEC_FUNC_TPROC_MOD_REARRANGE_IP */
#define zmpi_SPEC_FUNC_TPROC_MOD_REARRANGE_IP(_name_, _tp_, _s_) \
_s_ void _name_##_tproc_mod_rearrange_ip(zmpi_spec_tproc_data_t tproc_data, zmpi_spec_elem_t *s, zmpi_spec_elem_t *x, int *displs, int *counts, zmpi_spec_int_t n, zmpi_spec_elem_t *mod) \
{ \
  zmpi_SPEC_DECLARE_TPROC_MOD_REARRANGE_IP \
  zmpi_SPEC_DO_TPROC_MOD_REARRANGE_IP(_tp_, tproc_data, s, x, displs, counts, n, mod); \
}


/* tprocs rearrange */

/* sp_macro zmpi_SPEC_DECLARE_FIXED_TPROCS_REARRANGE_DB */
#define zmpi_SPEC_DECLARE_FIXED_TPROCS_REARRANGE_DB(_fxdcl_, _fxp_) \
  struct { zmpi_spec_elem_index_t n, i; zmpi_spec_int_t j, m; _fxdcl_(fx, _fxp_) } spec2d;

/* sp_macro zmpi_SPEC_DO_FIXED_TPROCS_REARRANGE_DB */
#define zmpi_SPEC_DO_FIXED_TPROCS_REARRANGE_DB(_fxp_, _fxc_, _fxca_, _fxxa_, _fxd_, _tp_, _tpd_, _sb_, _db_, _ds_, _ps_)  do { \
  _fxc_(spec2d.fx, _fxp_); \
  spec2d.n = zmpi_spec_elem_get_n(_sb_); \
  for (spec2d.i = 0; spec2d.i < spec2d.n; ++spec2d.i) { \
    (_tp_)(zmpi_spec_elem_get_buf(_sb_), spec2d.i, (_tpd_), &spec2d.m, (_ps_)); \
    for (spec2d.j = 0; spec2d.j < spec2d.m; ++spec2d.j) { \
      _fxca_((_sb_), spec2d.i, (_db_), (_ds_)[(_ps_)[spec2d.j]], spec2d.fx, _fxp_); \
      ++(_ds_)[(_ps_)[spec2d.j]]; \
    } \
  } \
  _fxd_(spec2d.fx, _fxp_); \
  } while (0)

/* sp_macro zmpi_SPEC_FUNC_FIXED_TPROCS_REARRANGE_DB */
#define zmpi_SPEC_FUNC_FIXED_TPROCS_REARRANGE_DB(_name_, _fxdcl_, _fxp_, _fxc_, _fxca_, _fxxa_, _fxd_, _tp_, _s_...) \
_s_ void _name_##_tprocs_rearrange_db(zmpi_spec_tproc_data_t tproc_data, zmpi_spec_elem_t *s, zmpi_spec_elem_t *d, int *displs, zmpi_spec_proc_t *procs) \
{ \
  zmpi_SPEC_DECLARE_FIXED_TPROCS_REARRANGE_DB(_fxdcl_, _fxp_) \
  zmpi_SPEC_DO_FIXED_TPROCS_REARRANGE_DB(_fxp_, _fxc_, _fxca_, _fxxa_, _fxd_, _tp_, tproc_data, s, d, displs, procs); \
}

/* sp_macro zmpi_SPEC_DECLARE_TPROCS_REARRANGE_DB */
#define zmpi_SPEC_DECLARE_TPROCS_REARRANGE_DB \
  zmpi_SPEC_DECLARE_FIXED_TPROCS_REARRANGE_DB(zmpi_spec_fixed_default_declare, NOPARAM)

/* sp_macro zmpi_SPEC_DO_TPROCS_REARRANGE_DB */
#define zmpi_SPEC_DO_TPROCS_REARRANGE_DB(_tp_, _tpd_, _sb_, _db_, _ds_, _ps_) \
  zmpi_SPEC_DO_FIXED_TPROCS_REARRANGE_DB(NOPARAM, zmpi_spec_fixed_default_create, zmpi_spec_fixed_default_copy_at, zmpi_spec_fixed_default_exchange_at, zmpi_spec_fixed_default_destroy, _tp_, _tpd_, _sb_, _db_, _ds_, _ps_)

/* sp_macro zmpi_SPEC_FUNC_TPROCS_REARRANGE_DB */
#define zmpi_SPEC_FUNC_TPROCS_REARRANGE_DB(_name_, _tp_, _s_...) \
  zmpi_SPEC_FUNC_FIXED_TPROCS_REARRANGE_DB(_name_, zmpi_spec_fixed_default_declare, NOPARAM, zmpi_spec_fixed_default_create, zmpi_spec_fixed_default_copy_at, zmpi_spec_fixed_default_exchange_at, zmpi_spec_fixed_default_destroy, _tp_, _s_)

/* sp_macro zmpi_SPEC_DECLARE_TPROCS_REARRANGE_IP */
#define zmpi_SPEC_DECLARE_TPROCS_REARRANGE_IP \
  struct { zmpi_spec_elem_index_t e, j, fe, fc, le, lc; zmpi_spec_int_t i, n, f, l, o; } spec2i;

/* sp_macro zmpi_SPEC_DO_TPROCS_REARRANGE_IP */
#define zmpi_SPEC_DO_TPROCS_REARRANGE_IP(_tp_, _tpd_, _b_, _xb_, _ds_, _cs_, _n_, _ps_)  do { \
  spec2i.f = 0; spec2i.fe = (_cs_)[0]; spec2i.fc = zmpi_spec_elem_get_n(_b_); \
  while (spec2i.f + 1 < (_n_) && spec2i.fc >= spec2i.fe) { ++spec2i.f; spec2i.fe += (_cs_)[spec2i.f]; } \
  spec2i.l = 0; spec2i.le = (_cs_)[0]; spec2i.lc = zmpi_spec_elem_get_n(_b_) - 1; \
  while (spec2i.lc >= spec2i.le) { ++spec2i.l; spec2i.le += (_cs_)[spec2i.l]; } \
  for (spec2i.e = 0, spec2i.i = 0; spec2i.i < (_n_); ++spec2i.i) { \
    spec2i.e += (_cs_)[spec2i.i]; \
    spec2i.j = (_ds_)[spec2i.i]; \
    while (spec2i.j < spec2i.e) { \
      (_tp_)(zmpi_spec_elem_get_buf(_b_), spec2i.j, (_tpd_), &spec2i.n, (_ps_)); \
      spec2i.o = -1; \
      while (spec2i.n > 0) { \
        --spec2i.n; \
        if ((_ps_)[spec2i.n] == spec2i.i && spec2i.o < 0) spec2i.o = spec2i.n; \
        else if ((_ds_)[(_ps_)[spec2i.n]] < spec2i.fc) { \
          spec2i.l = spec2i.f; spec2i.le = spec2i.fe; spec2i.lc = spec2i.fc; \
          if (spec2i.fc < spec2i.fe) { \
            zmpi_spec_elem_copy_at((_b_), (_ds_)[(_ps_)[spec2i.n]], (_b_), spec2i.fc); \
            ++spec2i.fc; \
          } else zmpi_spec_elem_copy_at((_b_), (_ds_)[(_ps_)[spec2i.n]], (_xb_), 0); \
        } else if ((_ds_)[(_ps_)[spec2i.n]] == spec2i.fc) ++spec2i.fc; \
        if (spec2i.j != (_ds_)[(_ps_)[spec2i.n]]) zmpi_spec_elem_copy_at((_b_), spec2i.j, (_b_), (_ds_)[(_ps_)[spec2i.n]]); \
        ++(_ds_)[(_ps_)[spec2i.n]]; \
        while (spec2i.f + 1 < (_n_) && spec2i.fc >= spec2i.fe) { ++spec2i.f; spec2i.fe += (_cs_)[spec2i.f]; spec2i.fc = (_ds_)[spec2i.f]; } \
      } \
      if (spec2i.o < 0) { \
        if (spec2i.lc < spec2i.le) {  \
          zmpi_spec_elem_copy_at((_b_), spec2i.lc, (_b_), spec2i.j); \
          spec2i.f = spec2i.l; spec2i.fe = spec2i.le; spec2i.fc = spec2i.lc; \
          --spec2i.lc; \
          while (spec2i.l > 0 && spec2i.lc < (_ds_)[spec2i.l]) { spec2i.le -= (_cs_)[spec2i.l]; spec2i.lc = spec2i.le - 1; --spec2i.l; } \
        } else zmpi_spec_elem_copy_at((_xb_), 0, (_b_), spec2i.j); \
      } \
      spec2i.j = (_ds_)[spec2i.i]; \
    } \
  } } while (0)

/* sp_macro zmpi_SPEC_FUNC_TPROCS_REARRANGE_IP */
#define zmpi_SPEC_FUNC_TPROCS_REARRANGE_IP(_name_, _tp_, _s_...) \
_s_ void _name_##_tprocs_rearrange_ip(zmpi_spec_tproc_data_t tproc_data, zmpi_spec_elem_t *s, zmpi_spec_elem_t *d, int *displs, int *counts, zmpi_spec_int_t n, zmpi_spec_proc_t *procs) \
{ \
  zmpi_SPEC_DECLARE_TPROCS_REARRANGE_IP \
  zmpi_SPEC_DO_TPROCS_REARRANGE_IP(_tp_, tproc_data, s, d, displs, counts, n, procs); \
}


/* tprocs_mod rearrange */

/* sp_macro zmpi_SPEC_DECLARE_FIXED_TPROCS_MOD_REARRANGE_DB */
#define zmpi_SPEC_DECLARE_FIXED_TPROCS_MOD_REARRANGE_DB(_fxdcl_, _fxp_) \
  struct { zmpi_spec_elem_index_t n, i; zmpi_spec_int_t j, m; _fxdcl_(fx, _fxp_) } spec3d;

/* sp_macro zmpi_SPEC_DO_FIXED_TPROCS_MOD_REARRANGE_DB */
#define zmpi_SPEC_DO_FIXED_TPROCS_MOD_REARRANGE_DB(_fxp_, _fxc_, _fxca_, _fxxa_, _fxd_, _tp_, _tpd_, _sb_, _db_, _ds_, _ps_, _ib_)  do { \
  _fxc_(spec3d.fx, _fxp_); \
  spec3d.n = zmpi_spec_elem_get_n(_sb_); \
  for (spec3d.i = 0; spec3d.i < spec3d.n; ++spec3d.i) { \
    (_tp_)(zmpi_spec_elem_get_buf(_sb_), spec3d.i, (_tpd_), &spec3d.m, (_ps_), zmpi_spec_elem_get_buf(_ib_)); \
    for (spec3d.j = 0; spec3d.j < spec3d.m; ++spec3d.j) { \
      _fxca_((_ib_), spec3d.j, (_db_), (_ds_)[(_ps_)[spec3d.j]], spec3d.fx, _fxp_); \
      ++(_ds_)[(_ps_)[spec3d.j]]; \
    } \
  } \
  _fxd_(spec3d.fx, _fxp_); \
  } while (0)

/* sp_macro zmpi_SPEC_FUNC_FIXED_TPROCS_MOD_REARRANGE_DB */
#define zmpi_SPEC_FUNC_FIXED_TPROCS_MOD_REARRANGE_DB(_name_, _fxdcl_, _fxp_, _fxc_, _fxca_, _fxxa_, _fxd_, _tp_, _s_...) \
_s_ void _name_##_tprocs_mod_rearrange_db(zmpi_spec_tproc_data_t tproc_data, zmpi_spec_elem_t *s, zmpi_spec_elem_t *d, int *displs, zmpi_spec_proc_t *procs, zmpi_spec_elem_t *mod) \
{ \
  zmpi_SPEC_DECLARE_FIXED_TPROCS_MOD_REARRANGE_DB(_fxdcl_, _fxp_) \
  zmpi_SPEC_DO_FIXED_TPROCS_MOD_REARRANGE_DB(_fxp_, _fxc_, _fxca_, _fxxa_, _fxd_, _tp_, tproc_data, s, d, displs, procs, mod); \
}

/* sp_macro zmpi_SPEC_DECLARE_TPROCS_MOD_REARRANGE_DB */
#define zmpi_SPEC_DECLARE_TPROCS_MOD_REARRANGE_DB \
  zmpi_SPEC_DECLARE_FIXED_TPROCS_MOD_REARRANGE_DB(zmpi_spec_fixed_default_declare, NOPARAM)

/* sp_macro zmpi_SPEC_DO_TPROCS_MOD_REARRANGE_DB */
#define zmpi_SPEC_DO_TPROCS_MOD_REARRANGE_DB(_tp_, _tpd_, _sb_, _db_, _ds_, _ps_, _ib_) \
  zmpi_SPEC_DO_FIXED_TPROCS_MOD_REARRANGE_DB(NOPARAM, zmpi_spec_fixed_default_create, zmpi_spec_fixed_default_copy_at, zmpi_spec_fixed_default_exchange_at, zmpi_spec_fixed_default_destroy, _tp_, _tpd_, _sb_, _db_, _ds_, _ps_, _ib_)

/* sp_macro zmpi_SPEC_FUNC_TPROCS_MOD_REARRANGE_DB */
#define zmpi_SPEC_FUNC_TPROCS_MOD_REARRANGE_DB(_name_, _tp_, _s_...) \
  zmpi_SPEC_FUNC_FIXED_TPROCS_MOD_REARRANGE_DB(_name_, zmpi_spec_fixed_default_declare, NOPARAM, zmpi_spec_fixed_default_create, zmpi_spec_fixed_default_copy_at, zmpi_spec_fixed_default_exchange_at, zmpi_spec_fixed_default_destroy, _tp_, _s_)

/* sp_macro zmpi_SPEC_DECLARE_TPROCS_MOD_REARRANGE_IP */
#define zmpi_SPEC_DECLARE_TPROCS_MOD_REARRANGE_IP \
  struct { zmpi_spec_elem_index_t e, j, fe, fc, le, lc; zmpi_spec_int_t i, n, f, l, o; } spec3i;

/* sp_macro zmpi_SPEC_DO_TPROCS_MOD_REARRANGE_IP */
#define zmpi_SPEC_DO_TPROCS_MOD_REARRANGE_IP(_tp_, _tpd_, _b_, _xb_, _ds_, _cs_, _n_, _ps_, _ib_)  do { \
  spec3i.f = 0; spec3i.fe = (_cs_)[0]; spec3i.fc = zmpi_spec_elem_get_n(_b_); \
  while (spec3i.f + 1 < (_n_) && spec3i.fc >= spec3i.fe) { ++spec3i.f; spec3i.fe += (_cs_)[spec3i.f]; } \
  spec3i.l = 0; spec3i.le = (_cs_)[0]; spec3i.lc = zmpi_spec_elem_get_n(_b_) - 1; \
  while (spec3i.lc >= spec3i.le) { ++spec3i.l; spec3i.le += (_cs_)[spec3i.l]; } \
  for (spec3i.e = 0, spec3i.i = 0; spec3i.i < (_n_); ++spec3i.i) { \
    spec3i.e += (_cs_)[spec3i.i]; \
    spec3i.j = (_ds_)[spec3i.i]; \
    while (spec3i.j < spec3i.e) { \
      (_tp_)(zmpi_spec_elem_get_buf(_b_), spec3i.j, (_tpd_), &spec3i.n, (_ps_), zmpi_spec_elem_get_buf(_ib_)); \
      spec3i.o = -1; \
      while (spec3i.n > 0) { \
        --spec3i.n; \
        if ((_ps_)[spec3i.n] == spec3i.i && spec3i.o < 0) spec3i.o = spec3i.n; \
        else if ((_ds_)[(_ps_)[spec3i.n]] < spec3i.fc) { \
          spec3i.l = spec3i.f; spec3i.le = spec3i.fe; spec3i.lc = spec3i.fc; \
          if (spec3i.fc < spec3i.fe) { \
            zmpi_spec_elem_copy_at((_b_), (_ds_)[(_ps_)[spec3i.n]], (_b_), spec3i.fc); \
            ++spec3i.fc; \
          } else zmpi_spec_elem_copy_at((_b_), (_ds_)[(_ps_)[spec3i.n]], (_xb_), 0); \
        } else if ((_ds_)[(_ps_)[spec3i.n]] == spec3i.fc) ++spec3i.fc; \
        zmpi_spec_elem_copy_at((_ib_), spec3i.n, (_b_), (_ds_)[(_ps_)[spec3i.n]]); \
        ++(_ds_)[(_ps_)[spec3i.n]]; \
        while (spec3i.f + 1 < (_n_) && spec3i.fc >= spec3i.fe) { ++spec3i.f; spec3i.fe += (_cs_)[spec3i.f]; spec3i.fc = (_ds_)[spec3i.f]; } \
      } \
      if (spec3i.o < 0) { \
        if (spec3i.lc < spec3i.le) {  \
          zmpi_spec_elem_copy_at((_b_), spec3i.lc, (_b_), spec3i.j); \
          spec3i.f = spec3i.l; spec3i.fe = spec3i.le; spec3i.fc = spec3i.lc; \
          --spec3i.lc; \
          while (spec3i.l > 0 && spec3i.lc < (_ds_)[spec3i.l]) { spec3i.le -= (_cs_)[spec3i.l]; spec3i.lc = spec3i.le - 1; --spec3i.l; } \
        } else zmpi_spec_elem_copy_at((_xb_), 0, (_b_), spec3i.j); \
      } \
      spec3i.j = (_ds_)[spec3i.i]; \
    } \
  } } while (0)

/* sp_macro zmpi_SPEC_FUNC_TPROCS_MOD_REARRANGE_IP */
#define zmpi_SPEC_FUNC_TPROCS_MOD_REARRANGE_IP(_name_, _tp_, _s_...) \
_s_ void _name_##_tprocs_mod_rearrange_ip(zmpi_spec_tproc_data_t tproc_data, zmpi_spec_elem_t *s, zmpi_spec_elem_t *x, int *displs, int *counts, zmpi_spec_int_t n, zmpi_spec_proc_t *procs, zmpi_spec_elem_t *mod) \
{ \
  zmpi_SPEC_DECLARE_TPROCS_MOD_REARRANGE_IP \
  zmpi_SPEC_DO_TPROCS_MOD_REARRANGE_IP(_tp_, tproc_data, s, x, displs, counts, n, procs, mod); \
}


/* tproc indices */

/* sp_macro zmpi_SPEC_DECLARE_TPROC_INDICES_DB */
#define zmpi_SPEC_DECLARE_TPROC_INDICES_DB \
  struct { zmpi_spec_elem_index_t i; zmpi_spec_proc_t p; } spec0xd;

/* sp_macro zmpi_SPEC_DO_TPROC_INDICES_DB */
#define zmpi_SPEC_DO_TPROC_INDICES_DB(_tp_, _tpd_, _b_, _ix_, _id_)  do { \
  for (spec0xd.i = 0; spec0xd.i < zmpi_spec_elem_get_n(_b_); ++spec0xd.i) { \
    spec0xd.p = (_tp_)(zmpi_spec_elem_get_buf(_b_), spec0xd.i, (_tpd_)); \
    if (spec0xd.p == zmpi_SPEC_PROC_NONE) continue; \
    (_ix_)[(_id_)[spec0xd.p]] = spec0xd.i; \
    ++(_id_)[spec0xd.p]; \
  } } while (0)

/* sp_macro zmpi_SPEC_FUNC_TPROC_INDICES_DB */
#define zmpi_SPEC_FUNC_TPROC_INDICES_DB(_name_, _tp_, _s_...) \
_s_ void _name_##_tproc_indices_db(zmpi_spec_tproc_data_t tproc_data, zmpi_spec_elem_t *s, int *indices, int *idispls) \
{ \
  zmpi_SPEC_DECLARE_TPROC_INDICES_DB \
  zmpi_SPEC_DO_TPROC_INDICES_DB(_tp_, tproc_data, s, indices, idispls); \
}


/* tproc_mod indices */

/* sp_macro zmpi_SPEC_DECLARE_TPROC_MOD_INDICES_DB */
#define zmpi_SPEC_DECLARE_TPROC_MOD_INDICES_DB \
  struct { zmpi_spec_elem_index_t i, k; zmpi_spec_proc_t p; } spec1xd;

/* sp_macro zmpi_SPEC_DO_TPROC_MOD_INDICES_DB */
#define zmpi_SPEC_DO_TPROC_MOD_INDICES_DB(_tp_, _tpd_, _b_, _ix_, _id_, _ib_, _d_)  do { \
  spec1xd.k = 0; \
  for (spec1xd.i = 0; spec1xd.i < zmpi_spec_elem_get_n(_b_); ++spec1xd.i) { \
    spec1xd.p = (_tp_)(zmpi_spec_elem_get_buf(_b_), spec1xd.i, (_tpd_), zmpi_spec_elem_get_buf(_ib_)); \
    if (spec1xd.p == zmpi_SPEC_PROC_NONE) continue; \
    zmpi_spec_elem_copy_at((_ib_), 0, (_d_), spec1xd.k); \
    (_ix_)[(_id_)[spec1xd.p]] = spec1xd.k; \
    ++spec1xd.k; \
    ++(_id_)[spec1xd.p]; \
  } } while (0)

/* sp_macro zmpi_SPEC_FUNC_TPROC_MOD_INDICES_DB */
#define zmpi_SPEC_FUNC_TPROC_MOD_INDICES_DB(_name_, _tp_, _s_...) \
_s_ void _name_##_tproc_mod_indices_db(zmpi_spec_tproc_data_t tproc_data, zmpi_spec_elem_t *s, int *indices, int *idispls, zmpi_spec_elem_t *mod, zmpi_spec_elem_t *d) \
{ \
  zmpi_SPEC_DECLARE_TPROC_MOD_INDICES_DB \
  zmpi_SPEC_DO_TPROC_MOD_INDICES_DB(_tp_, tproc_data, s, indices, idispls, mod, d); \
}


/* tprocs indices */

/* sp_macro zmpi_SPEC_DECLARE_TPROCS_INDICES_DB */
#define zmpi_SPEC_DECLARE_TPROCS_INDICES_DB \
  struct { zmpi_spec_elem_index_t i; zmpi_spec_int_t j, n; } spec2xd;

/* sp_macro zmpi_SPEC_DO_TPROCS_INDICES_DB */
#define zmpi_SPEC_DO_TPROCS_INDICES_DB(_tp_, _tpd_, _b_, _ix_, _id_, _ps_)  do { \
  for (spec2xd.i = 0; spec2xd.i < zmpi_spec_elem_get_n(_b_); ++spec2xd.i) { \
    (_tp_)(zmpi_spec_elem_get_buf(_b_), spec2xd.i, (_tpd_), &spec2xd.n, (_ps_)); \
    for (spec2xd.j = 0; spec2xd.j < spec2xd.n; ++spec2xd.j) { \
      (_ix_)[(_id_)[(_ps_)[spec2xd.j]]] = spec2xd.i; \
      ++(_id_)[(_ps_)[spec2xd.j]]; \
    } \
  } } while (0)

/* sp_macro zmpi_SPEC_FUNC_TPROCS_INDICES_DB */
#define zmpi_SPEC_FUNC_TPROCS_INDICES_DB(_name_, _tp_, _s_...) \
_s_ void _name_##_tprocs_indices_db(zmpi_spec_tproc_data_t tproc_data, zmpi_spec_elem_t *s, int *indices, int *idispls, zmpi_spec_proc_t *procs) \
{ \
  zmpi_SPEC_DECLARE_TPROCS_INDICES_DB \
  zmpi_SPEC_DO_TPROCS_INDICES_DB(_tp_, tproc_data, s, indices, idispls, procs); \
}


/* tprocs_mod indices */

/* sp_macro zmpi_SPEC_DECLARE_TPROCS_MOD_INDICES_DB */
#define zmpi_SPEC_DECLARE_TPROCS_MOD_INDICES_DB \
  struct { zmpi_spec_elem_index_t i, k; zmpi_spec_int_t j, n; } spec3xd;

/* sp_macro zmpi_SPEC_DO_TPROCS_MOD_INDICES_DB */
#define zmpi_SPEC_DO_TPROCS_MOD_INDICES_DB(_tp_, _tpd_, _b_, _ix_, _id_, _ps_, _ib_, _d_)  do { \
  spec3xd.k = 0; \
  for (spec3xd.i = 0; spec3xd.i < zmpi_spec_elem_get_n(_b_); ++spec3xd.i) { \
    (_tp_)(zmpi_spec_elem_get_buf(_b_), spec3xd.i, (_tpd_), &spec3xd.n, (_ps_), zmpi_spec_elem_get_buf(_ib_)); \
    for (spec3xd.j = 0; spec3xd.j < spec3xd.n; ++spec3xd.j) { \
      zmpi_spec_elem_copy_at((_ib_), spec3xd.j, (_d_), spec3xd.k); \
      (_ix_)[(_id_)[(_ps_)[spec3xd.j]]] = spec3xd.k; \
      ++spec3xd.k; \
      ++(_id_)[(_ps_)[spec3xd.j]]; \
    } \
  } } while (0)

/* sp_macro zmpi_SPEC_FUNC_TPROCS_MOD_INDICES_DB */
#define zmpi_SPEC_FUNC_TPROCS_MOD_INDICES_DB(_name_, _tp_, _s_...) \
_s_ void _name_##_tprocs_mod_indices_db(zmpi_spec_tproc_data_t tproc_data, zmpi_spec_elem_t *s, int *indices, int *idispls, zmpi_spec_proc_t *procs, zmpi_spec_elem_t *mod, zmpi_spec_elem_t *d) \
{ \
  zmpi_SPEC_DECLARE_TPROCS_MOD_INDICES_DB \
  zmpi_SPEC_DO_TPROCS_MOD_INDICES_DB(_tp_, tproc_data, s, indices, idispls, procs, mod, d); \
}


/* tproc sendrecv */

/* sp_macro zmpi_SPEC_DECLARE_FIXED_TPROC_SENDRECV_DB */
#define zmpi_SPEC_DECLARE_FIXED_TPROC_SENDRECV_DB(_fxdcl_, _fxp_)
/*#define zmpi_SPEC_DECLARE_FIXED_TPROC_SENDRECV_DB(_fxdcl_, _fxp_) \
  struct { _fxdcl_(fx, _fxp_) } spec0srd;*/

/* sp_macro zmpi_SPEC_DO_FIXED_TPROC_SENDRECV_DB */
#define zmpi_SPEC_DO_FIXED_TPROC_SENDRECV_DB(_fxp_, _fxc_, _fxca_, _fxxa_, _fxd_, _tp_, _tpd_, _sb_, _rb_, _sc_, _sd_, _rd_, _ab_, _ad_, _as_, _aq_, _aqn_, _aqs_, _r_, _p_)  do { \
  _fxc_(spec0srd.fx, _fxp_); \
  while (*(_sd_) < (_sc_)) { \
    if ((_p_) == zmpi_SPEC_PROC_NONE) (_p_) = (_tp_)(zmpi_spec_elem_get_buf(_sb_), *(_sd_), (_tpd_)); \
    if ((_p_) != zmpi_SPEC_PROC_NONE) { \
      if ((_p_) == (_r_)) { \
        _fxca_((_sb_), *(_sd_), (_rb_), *(_rd_), spec0srd.fx, _fxp_); \
        ++(*(_rd_)); \
      } else { \
        if ((_ad_)[_p_] >= ((_p_) + 1) * (_as_)) break; \
        _fxca_((_sb_), *(_sd_), (_ab_), (_ad_)[_p_], spec0srd.fx, _fxp_); \
        ++(_ad_)[_p_]; \
        if ((_ad_)[_p_] >= ((_p_) + 1) * (_as_)) { \
          (_aq_)[*(_aqn_)] = (_p_); ++(*(_aqn_)); *(_aqn_) %= (_aqs_); \
        } \
      } \
    } \
    (_p_) = zmpi_SPEC_PROC_NONE; \
    ++(*(_sd_)); \
 } \
 _fxd_(spec0srd.fx, _fxp_); \
 } while (0)

/* sp_macro zmpi_SPEC_FUNC_FIXED_TPROC_SENDRECV_DB */
#define zmpi_SPEC_FUNC_FIXED_TPROC_SENDRECV_DB(_name_, _fxdcl_, _fxp_, _fxc_, _fxca_, _fxxa_, _fxd_, _tp_, _s_...) \
_s_ zmpi_spec_proc_t _name_##_tproc_sendrecv_db(zmpi_spec_tproc_data_t tproc_data, zmpi_spec_elem_t *sb, zmpi_spec_elem_t *rb, zmpi_spec_int_t scount, zmpi_spec_int_t *sdispl, zmpi_spec_int_t *rdispl, zmpi_spec_elem_t *ax, zmpi_spec_int_t *aux_displs, zmpi_spec_int_t aux_size_max, zmpi_spec_int_t *aux_queue, zmpi_spec_int_t *aux_queue_next, zmpi_spec_int_t aux_queue_size, zmpi_spec_proc_t rank, zmpi_spec_proc_t p) \
{ \
  zmpi_SPEC_DECLARE_FIXED_TPROC_SENDRECV_DB(_fxdcl_, _fxp_) \
  zmpi_SPEC_DO_FIXED_TPROC_SENDRECV_DB(_fxp_, _fxc_, _fxca_, _fxxa_, _fxd_, _tp_, tproc_data, sb, rb, scount, sdispl, rdispl, ax, aux_displs, aux_size_max, aux_queue, aux_queue_next, aux_queue_size, rank, p); \
  return p; \
}

/* sp_macro zmpi_SPEC_DECLARE_TPROC_SENDRECV_DB */
#define zmpi_SPEC_DECLARE_TPROC_SENDRECV_DB \
  zmpi_SPEC_DECLARE_FIXED_TPROC_SENDRECV_DB(zmpi_spec_fixed_default_declare, NOPARAM)

/* sp_macro zmpi_SPEC_DO_TPROC_SENDRECV_DB */
#define zmpi_SPEC_DO_TPROC_SENDRECV_DB(_tp_, _tpd_, _sb_, _rb_, _sc_, _sd_, _rd_, _ab_, _ad_, _as_, _aq_, _aqn_, _aqs_, _r_, _p_) \
  zmpi_SPEC_DO_FIXED_TPROC_SENDRECV_DB(NOPARAM, zmpi_spec_fixed_default_create, zmpi_spec_fixed_default_copy_at, zmpi_spec_fixed_default_exchange_at, zmpi_spec_fixed_default_destroy, _tp_, _tpd_, _sb_, _rb_, _sc_, _sd_, _rd_, _ab_, _ad_, _as_, _aq_, _aqn_, _aqs_, _r_, _p_)

/* sp_macro zmpi_SPEC_FUNC_TPROC_SENDRECV_DB */
#define zmpi_SPEC_FUNC_TPROC_SENDRECV_DB(_name_, _tp_, _s_...) \
  zmpi_SPEC_FUNC_FIXED_TPROC_SENDRECV_DB(_name_, zmpi_spec_fixed_default_declare, NOPARAM, zmpi_spec_fixed_default_create, zmpi_spec_fixed_default_copy_at, zmpi_spec_fixed_default_exchange_at, zmpi_spec_fixed_default_destroy, _tp_, _s_)


/* sp_macro zmpi_SPEC_DEFINE_TPROC */
#define zmpi_SPEC_DEFINE_TPROC(_name_, _tp_, _s_...) \
  zmpi_SPEC_FUNC_TPROC_COUNT_DB(_name_, _tp_, _s_) \
  zmpi_SPEC_FUNC_TPROC_COUNT_IP(_name_, _tp_, _s_) \
  zmpi_SPEC_FUNC_TPROC_REARRANGE_DB(_name_, _tp_, _s_) \
  zmpi_SPEC_FUNC_TPROC_REARRANGE_IP(_name_, _tp_, _s_) \
  zmpi_SPEC_FUNC_TPROC_INDICES_DB(_name_, _tp_, _s_) \
  zmpi_SPEC_FUNC_TPROC_SENDRECV_DB(_name_, _tp_, _s_)

/* sp_macro zmpi_SPEC_DEFINE_TPROC_MOD */
#define zmpi_SPEC_DEFINE_TPROC_MOD(_name_, _tp_, _s_...) \
  zmpi_SPEC_FUNC_TPROC_MOD_COUNT_DB(_name_, _tp_, _s_) \
  zmpi_SPEC_FUNC_TPROC_MOD_COUNT_IP(_name_, _tp_, _s_) \
  zmpi_SPEC_FUNC_TPROC_MOD_REARRANGE_DB(_name_, _tp_, _s_) \
  zmpi_SPEC_FUNC_TPROC_MOD_REARRANGE_IP(_name_, _tp_, _s_) \
  zmpi_SPEC_FUNC_TPROC_MOD_INDICES_DB(_name_, _tp_, _s_)

/* sp_macro zmpi_SPEC_DEFINE_TPROCS */
#define zmpi_SPEC_DEFINE_TPROCS(_name_, _tp_, _s_...) \
  zmpi_SPEC_FUNC_TPROCS_COUNT_DB(_name_, _tp_, _s_) \
  zmpi_SPEC_FUNC_TPROCS_COUNT_IP(_name_, _tp_, _s_) \
  zmpi_SPEC_FUNC_TPROCS_REARRANGE_DB(_name_, _tp_, _s_) \
  zmpi_SPEC_FUNC_TPROCS_REARRANGE_IP(_name_, _tp_, _s_) \
  zmpi_SPEC_FUNC_TPROCS_INDICES_DB(_name_, _tp_, _s_)

/* sp_macro zmpi_SPEC_DEFINE_TPROCS_MOD */
#define zmpi_SPEC_DEFINE_TPROCS_MOD(_name_, _tp_, _s_...) \
  zmpi_SPEC_FUNC_TPROCS_MOD_COUNT_DB(_name_, _tp_, _s_) \
  zmpi_SPEC_FUNC_TPROCS_MOD_COUNT_IP(_name_, _tp_, _s_) \
  zmpi_SPEC_FUNC_TPROCS_MOD_REARRANGE_DB(_name_, _tp_, _s_) \
  zmpi_SPEC_FUNC_TPROCS_MOD_REARRANGE_IP(_name_, _tp_, _s_) \
  zmpi_SPEC_FUNC_TPROCS_MOD_INDICES_DB(_name_, _tp_, _s_)

/* sp_macro zmpi_SPEC_DEFINE_FIXED_TPROC */
#define zmpi_SPEC_DEFINE_FIXED_TPROC(_name_, _fxdcl_, _fxp_, _fxc_, _fxca_, _fxxa_, _fxd_, _tp_, _s_...) \
  zmpi_SPEC_FUNC_TPROC_COUNT_DB(_name_, _tp_, _s_) \
  zmpi_SPEC_FUNC_TPROC_COUNT_IP(_name_, _tp_, _s_) \
  zmpi_SPEC_FUNC_FIXED_TPROC_REARRANGE_DB(_name_, _fxdcl_, _fxp_, _fxc_, _fxca_, _fxxa_, _fxd_, _tp_, _s_) \
  zmpi_SPEC_FUNC_TPROC_REARRANGE_IP(_name_, _tp_, _s_) \
  zmpi_SPEC_FUNC_TPROC_INDICES_DB(_name_, _tp_, _s_) \
  zmpi_SPEC_FUNC_FIXED_TPROC_SENDRECV_DB(_name_, _fxdcl_, _fxp_, _fxc_, _fxca_, _fxxa_, _fxd_, _tp_, _s_)

/* sp_macro zmpi_SPEC_DEFINE_FIXED_TPROC_MOD */
#define zmpi_SPEC_DEFINE_FIXED_TPROC_MOD(_name_, _fxdcl_, _fxp_, _fxc_, _fxca_, _fxxa_, _fxd_, _tp_, _s_...) \
  zmpi_SPEC_FUNC_TPROC_MOD_COUNT_DB(_name_, _tp_, _s_) \
  zmpi_SPEC_FUNC_TPROC_MOD_COUNT_IP(_name_, _tp_, _s_) \
  zmpi_SPEC_FUNC_FIXED_TPROC_MOD_REARRANGE_DB(_name_, _fxdcl_, _fxp_, _fxc_, _fxca_, _fxxa_, _fxd_, _tp_, _s_) \
  zmpi_SPEC_FUNC_TPROC_MOD_REARRANGE_IP(_name_, _tp_, _s_) \
  zmpi_SPEC_FUNC_TPROC_MOD_INDICES_DB(_name_, _tp_, _s_)

/* sp_macro zmpi_SPEC_DEFINE_FIXED_TPROCS */
#define zmpi_SPEC_DEFINE_FIXED_TPROCS(_name_, _fxdcl_, _fxp_, _fxc_, _fxca_, _fxxa_, _fxd_, _tp_, _s_...) \
  zmpi_SPEC_FUNC_TPROCS_COUNT_DB(_name_, _tp_, _s_) \
  zmpi_SPEC_FUNC_TPROCS_COUNT_IP(_name_, _tp_, _s_) \
  zmpi_SPEC_FUNC_FIXED_TPROCS_REARRANGE_DB(_name_, _fxdcl_, _fxp_, _fxc_, _fxca_, _fxxa_, _fxd_, _tp_, _s_) \
  zmpi_SPEC_FUNC_TPROCS_REARRANGE_IP(_name_, _tp_, _s_) \
  zmpi_SPEC_FUNC_TPROCS_INDICES_DB(_name_, _tp_, _s_)

/* sp_macro zmpi_SPEC_DEFINE_FIXED_TPROCS_MOD */
#define zmpi_SPEC_DEFINE_FIXED_TPROCS_MOD(_name_, _fxdcl_, _fxp_, _fxc_, _fxca_, _fxxa_, _fxd_, _tp_, _s_...) \
  zmpi_SPEC_FUNC_TPROCS_MOD_COUNT_DB(_name_, _tp_, _s_) \
  zmpi_SPEC_FUNC_TPROCS_MOD_COUNT_IP(_name_, _tp_, _s_) \
  zmpi_SPEC_FUNC_FIXED_TPROCS_MOD_REARRANGE_DB(_name_, _fxdcl_, _fxp_, _fxc_, _fxca_, _fxxa_, _fxd_, _tp_, _s_) \
  zmpi_SPEC_FUNC_TPROCS_MOD_REARRANGE_IP(_name_, _tp_, _s_) \
  zmpi_SPEC_FUNC_TPROCS_MOD_INDICES_DB(_name_, _tp_, _s_)

/* sp_type zmpi_spec_tproc_f zmpi_spec_tproc_count_f zmpi_spec_tproc_rearrange_db_f zmpi_spec_tproc_rearrange_ip_f zmpi_spec_tproc_indices_db_f zmpi_spec_tproc_sendrecv_db_f */
typedef zmpi_spec_proc_t zmpi_spec_tproc_f(zmpi_spec_elem_buf_t b, zmpi_spec_elem_index_t x, zmpi_spec_tproc_data_t tproc_data);
typedef void zmpi_spec_tproc_count_f(zmpi_spec_tproc_data_t tproc_data, zmpi_spec_elem_t *s, int *counts);
typedef void zmpi_spec_tproc_rearrange_db_f(zmpi_spec_tproc_data_t tproc_data, zmpi_spec_elem_t *s, zmpi_spec_elem_t *d, int *displs);
typedef void zmpi_spec_tproc_rearrange_ip_f(zmpi_spec_tproc_data_t tproc_data, zmpi_spec_elem_t *s, zmpi_spec_elem_t *x, int *displs, int *counts, zmpi_spec_int_t n);
typedef void zmpi_spec_tproc_indices_db_f(zmpi_spec_tproc_data_t tproc_data, zmpi_spec_elem_t *s, int *indices, int *idispls);
typedef zmpi_spec_proc_t zmpi_spec_tproc_sendrecv_db_f(zmpi_spec_tproc_data_t tproc_data, zmpi_spec_elem_t *sb, zmpi_spec_elem_t *rb, zmpi_spec_int_t scount, zmpi_spec_int_t *sdispl, zmpi_spec_int_t *rdispl, zmpi_spec_elem_t *ax, zmpi_spec_int_t *aux_displs, zmpi_spec_int_t aux_size_max, zmpi_spec_int_t *aux_queue, zmpi_spec_int_t *aux_queue_next, zmpi_spec_int_t aux_queue_size, zmpi_spec_proc_t rank, zmpi_spec_proc_t p);

/* sp_type zmpi_spec_tproc_mod_f zmpi_spec_tproc_mod_count_f zmpi_spec_tproc_mod_rearrange_db_f zmpi_spec_tproc_mod_rearrange_ip_f zmpi_spec_tproc_mod_indices_db_f */
typedef zmpi_spec_proc_t zmpi_spec_tproc_mod_f(zmpi_spec_elem_buf_t b, zmpi_spec_elem_index_t x, zmpi_spec_tproc_data_t tproc_data, zmpi_spec_elem_buf_t mod);
typedef void zmpi_spec_tproc_mod_count_f(zmpi_spec_tproc_data_t tproc_data, zmpi_spec_elem_t *s, int *counts);
typedef void zmpi_spec_tproc_mod_rearrange_db_f(zmpi_spec_tproc_data_t tproc_data, zmpi_spec_elem_t *s, zmpi_spec_elem_t *d, int *displs, zmpi_spec_elem_t *mod);
typedef void zmpi_spec_tproc_mod_rearrange_ip_f(zmpi_spec_tproc_data_t tproc_data, zmpi_spec_elem_t *s, zmpi_spec_elem_t *x, int *displs, int *counts, zmpi_spec_int_t n, zmpi_spec_elem_t *mod);
typedef void zmpi_spec_tproc_mod_indices_db_f(zmpi_spec_tproc_data_t tproc_data, zmpi_spec_elem_t *s, int *indices, int *idispls, zmpi_spec_elem_t *mod, zmpi_spec_elem_t *d);

/* sp_type zmpi_spec_tprocs_f zmpi_spec_tprocs_count_f zmpi_spec_tprocs_rearrange_db_f zmpi_spec_tprocs_rearrange_ip_f zmpi_spec_tprocs_indices_db_f */
typedef void zmpi_spec_tprocs_f(zmpi_spec_elem_buf_t b, zmpi_spec_elem_index_t x, zmpi_spec_tproc_data_t tproc_data, zmpi_spec_int_t *nprocs, zmpi_spec_proc_t *procs);
typedef void zmpi_spec_tprocs_count_f(zmpi_spec_tproc_data_t tproc_data, zmpi_spec_elem_t *s, int *counts, zmpi_spec_proc_t *procs);
typedef void zmpi_spec_tprocs_rearrange_db_f(zmpi_spec_tproc_data_t tproc_data, zmpi_spec_elem_t *s, zmpi_spec_elem_t *d, int *displs, zmpi_spec_proc_t *procs);
typedef void zmpi_spec_tprocs_rearrange_ip_f(zmpi_spec_tproc_data_t tproc_data, zmpi_spec_elem_t *s, zmpi_spec_elem_t *x, int *displs, int *counts, zmpi_spec_int_t n, zmpi_spec_proc_t *procs);
typedef void zmpi_spec_tprocs_indices_db_f(zmpi_spec_tproc_data_t tproc_data, zmpi_spec_elem_t *s, int *indices, int *idispls, zmpi_spec_proc_t *procs);

/* sp_type zmpi_spec_tprocs_mod_f zmpi_spec_tprocs_mod_count_f zmpi_spec_tprocs_mod_rearrange_db_f zmpi_spec_tprocs_mod_rearrange_ip_f zmpi_spec_tprocs_mod_indices_db_f */
typedef void zmpi_spec_tprocs_mod_f(zmpi_spec_elem_buf_t b, zmpi_spec_elem_index_t x, zmpi_spec_tproc_data_t tproc_data, zmpi_spec_int_t *nprocs, zmpi_spec_proc_t *procs, zmpi_spec_elem_buf_t mod);
typedef void zmpi_spec_tprocs_mod_count_f(zmpi_spec_tproc_data_t tproc_data, zmpi_spec_elem_t *s, int *counts, zmpi_spec_proc_t *procs);
typedef void zmpi_spec_tprocs_mod_rearrange_db_f(zmpi_spec_tproc_data_t tproc_data, zmpi_spec_elem_t *s, zmpi_spec_elem_t *d, int *displs, zmpi_spec_proc_t *procs, zmpi_spec_elem_t *mod);
typedef void zmpi_spec_tprocs_mod_rearrange_ip_f(zmpi_spec_tproc_data_t tproc_data, zmpi_spec_elem_t *s, zmpi_spec_elem_t *x, int *displs, int *counts, zmpi_spec_int_t n, zmpi_spec_proc_t *procs, zmpi_spec_elem_t *mod);
typedef void zmpi_spec_tprocs_mod_indices_db_f(zmpi_spec_tproc_data_t tproc_data, zmpi_spec_elem_t *s, int *indices, int *idispls, zmpi_spec_proc_t *procs, zmpi_spec_elem_t *mod, zmpi_spec_elem_t *d);

/* sp_type zmpi_spec_tproc_reset_f */
typedef void zmpi_spec_tproc_reset_f(zmpi_spec_tproc_data_t tproc_data);


/* sp_type zmpi__spec_tproc_ext_t zmpi_spec_tproc_ext_t */
typedef struct zmpi__spec_tproc_ext_t
{
  zmpi_spec_tproc_count_f *count_db, *count_ip;
  zmpi_spec_tproc_rearrange_db_f *rearrange_db;
  zmpi_spec_tproc_rearrange_ip_f *rearrange_ip;
  zmpi_spec_tproc_indices_db_f *indices_db;
  zmpi_spec_tproc_sendrecv_db_f *sendrecv_db;

} zmpi_spec_tproc_ext_t;

/* sp_type zmpi__spec_tproc_mod_ext_tproc_t zmpi_spec_tproc_mod_ext_t */
typedef struct zmpi__spec_tproc_mod_ext_tproc_t
{
  zmpi_spec_tproc_mod_count_f *count_db, *count_ip;
  zmpi_spec_tproc_mod_rearrange_db_f *rearrange_db;
  zmpi_spec_tproc_mod_rearrange_ip_f *rearrange_ip;
  zmpi_spec_tproc_mod_indices_db_f *indices_db;

} zmpi_spec_tproc_mod_ext_t;

/* sp_type zmpi__spec_tprocs_ext_t zmpi_spec_tprocs_ext_t */
typedef struct zmpi__spec_tprocs_ext_t
{
  zmpi_spec_tprocs_count_f *count_db, *count_ip;
  zmpi_spec_tprocs_rearrange_db_f *rearrange_db;
  zmpi_spec_tprocs_rearrange_ip_f *rearrange_ip;
  zmpi_spec_tprocs_indices_db_f *indices_db;

} zmpi_spec_tprocs_ext_t;

/* sp_type zmpi__spec_tprocs_mod_ext_t zmpi_spec_tprocs_mod_ext_t */
typedef struct zmpi__spec_tprocs_mod_ext_t
{
  zmpi_spec_tprocs_mod_count_f *count_db, *count_ip;
  zmpi_spec_tprocs_mod_rearrange_db_f *rearrange_db;
  zmpi_spec_tprocs_mod_rearrange_ip_f *rearrange_ip;
  zmpi_spec_tprocs_mod_indices_db_f *indices_db;

} zmpi_spec_tprocs_mod_ext_t;

/* sp_macro zmpi_SPEC_EXT_PARAM_TPROC zmpi_SPEC_EXT_PARAM_TPROC_NULL zmpi_SPEC_EXT_PARAM_TPROC_MOD zmpi_SPEC_EXT_PARAM_TPROC_MOD_NULL zmpi_SPEC_EXT_PARAM_TPROCS zmpi_SPEC_EXT_PARAM_TPROCS_NULL zmpi_SPEC_EXT_PARAM_TPROCS_MOD zmpi_SPEC_EXT_PARAM_TPROCS_MOD_NULL */
#define zmpi_SPEC_EXT_PARAM_TPROC(_name_)       { _name_##_tproc_count_db, _name_##_tproc_count_ip, _name_##_tproc_rearrange_db, _name_##_tproc_rearrange_ip, _name_##_tproc_indices_db, _name_##_tproc_sendrecv_db }
#define zmpi_SPEC_EXT_PARAM_TPROC_NULL          { NULL, NULL, NULL, NULL, NULL, NULL }
#define zmpi_SPEC_EXT_PARAM_TPROC_MOD(_name_)   { _name_##_tproc_mod_count_db, _name_##_tproc_mod_count_ip, _name_##_tproc_mod_rearrange_db, _name_##_tproc_mod_rearrange_ip, _name_##_tproc_mod_indices_db }
#define zmpi_SPEC_EXT_PARAM_TPROC_MOD_NULL      { NULL, NULL, NULL, NULL, NULL }
#define zmpi_SPEC_EXT_PARAM_TPROCS(_name_)      { _name_##_tprocs_count_db, _name_##_tprocs_count_ip, _name_##_tprocs_rearrange_db, _name_##_tprocs_rearrange_ip, _name_##_tprocs_indices_db }
#define zmpi_SPEC_EXT_PARAM_TPROCS_NULL         { NULL, NULL, NULL, NULL, NULL }
#define zmpi_SPEC_EXT_PARAM_TPROCS_MOD(_name_)  { _name_##_tprocs_mod_count_db, _name_##_tprocs_mod_count_ip, _name_##_tprocs_mod_rearrange_db, _name_##_tprocs_mod_rearrange_ip, _name_##_tprocs_mod_indices_db }
#define zmpi_SPEC_EXT_PARAM_TPROCS_MOD_NULL     { NULL, NULL, NULL, NULL, NULL }


/* enable tloc features */
#ifdef zmpi_SPEC_TLOC

/* sp_macro zmpi_SPEC_TLOC zmpi_SPEC_LOC_NONE */


/* tloc rearrange */

/* sp_macro zmpi_SPEC_DECLARE_TLOC_REARRANGE_DB */
#define zmpi_SPEC_DECLARE_TLOC_REARRANGE_DB \
  struct { zmpi_spec_int_t i, p; } spec0d;

/* sp_macro zmpi_SPEC_DO_TLOC_REARRANGE_DB */
#define zmpi_SPEC_DO_TLOC_REARRANGE_DB(_tl_, _tld_, _sb_, _db_)  do { \
  for (spec0d.i = 0; spec0d.i < zmpi_spec_elem_get_n(_sb_); ++spec0d.i) { \
    spec0d.p = (_tl_)(zmpi_spec_elem_get_buf(_sb_), spec0d.i, _tld_); \
    if (spec0d.p == zmpi_SPEC_LOC_NONE) continue; \
    zmpi_spec_elem_copy_at((_sb_), spec0d.i, (_db_), spec0d.p); \
  } } while (0)

/* sp_macro zmpi_SPEC_FUNC_TLOC_REARRANGE_DB */
#define zmpi_SPEC_FUNC_TLOC_REARRANGE_DB(_name_, _tl_, _s_...) \
_s_ void _name_##_tloc_rearrange_db(zmpi_spec_tloc_data_t tloc_data, zmpi_spec_elem_t *s, zmpi_spec_elem_t *d) \
{ \
  zmpi_SPEC_DECLARE_TLOC_REARRANGE_DB \
  zmpi_SPEC_DO_TLOC_REARRANGE_DB(_tl_, tloc_data, s, d); \
}

/* sp_macro zmpi_SPEC_DECLARE_TLOC_REARRANGE_IP */
#define zmpi_SPEC_DECLARE_TLOC_REARRANGE_IP \
  struct { zmpi_spec_int_t i, p, np; } spec0i;

/* sp_macro zmpi_SPEC_DO_TLOC_REARRANGE_IP */
#define zmpi_SPEC_DO_TLOC_REARRANGE_IP(_tl_, _tld_, _b_, _xb_)  do { \
  for (spec0i.i = 0; spec0i.i < zmpi_spec_elem_get_n(_b_); ++spec0i.i) { \
    spec0i.p = (_tl_)(zmpi_spec_elem_get_buf(_b_), spec0i.i, _tld_); \
    if (spec0i.p == zmpi_SPEC_LOC_NONE) continue; \
    while (spec0i.i != spec0i.p) { \
      spec0i.np = (_tl_)(zmpi_spec_elem_get_buf(_b_), spec0i.p, _tld_); \
      if (spec0i.np == zmpi_SPEC_LOC_NONE) { zmpi_spec_elem_copy_at((_b_), spec0i.i, (_b_), spec0i.p); break; } \
      zmpi_spec_elem_exchange_at((_b_), spec0i.i, (_b_), spec0i.p, (_xb_)); \
      spec0i.p = spec0i.np; \
    } \
  } } while (0)

/* sp_macro zmpi_SPEC_FUNC_TLOC_REARRANGE_IP */
#define zmpi_SPEC_FUNC_TLOC_REARRANGE_IP(_name_, _tl_, _s_) \
_s_ void _name_##_tloc_rearrange_ip(zmpi_spec_elem_t *s, zmpi_spec_elem_t *x, zmpi_spec_tloc_data_t tloc_data) \
{ \
  zmpi_SPEC_DECLARE_TLOC_REARRANGE_IP \
  zmpi_SPEC_DO_TLOC_REARRANGE_IP(_tl_, tloc_data, s, x); \
}


/* tloc_mod_mod rearrange */

/* sp_macro zmpi_SPEC_DECLARE_TLOC_MOD_REARRANGE_DB */
#define zmpi_SPEC_DECLARE_TLOC_MOD_REARRANGE_DB \
  struct { zmpi_spec_int_t i, p; } spec1d;

/* sp_macro zmpi_SPEC_DO_TLOC_MOD_REARRANGE_DB */
#define zmpi_SPEC_DO_TLOC_MOD_REARRANGE_DB(_tl_, _tld_, _sb_, _db_, _ib_)  do { \
  if (_ib_) { \
    for (spec1d.i = 0; spec1d.i < zmpi_spec_elem_get_n(_sb_); ++spec1d.i) { \
      spec1d.p = (_tl_)(zmpi_spec_elem_get_buf(_sb_), spec1d.i, _tld_, zmpi_spec_elem_get_buf(_ib_)); \
      if (spec1d.p == zmpi_SPEC_LOC_NONE) continue; \
      zmpi_spec_elem_copy_at((_ib_), 0, (_db_), spec1d.p); \
    } \
  } else { \
    for (spec1d.i = 0; spec1d.i < zmpi_spec_elem_get_n(_sb_); ++spec1d.i) { \
      spec1d.p = (_tl_)(zmpi_spec_elem_get_buf(_sb_), spec1d.i, _tld_, NULL); \
      if (spec1d.p == zmpi_SPEC_LOC_NONE) continue; \
      zmpi_spec_elem_copy_at((_sb_), spec1d.i, (_db_), spec1d.p); \
    } \
  } } while (0) 

/* sp_macro zmpi_SPEC_FUNC_TLOC_MOD_REARRANGE_DB */
#define zmpi_SPEC_FUNC_TLOC_MOD_REARRANGE_DB(_name_, _tl_, _s_...) \
_s_ void _name_##_tloc_mod_rearrange_db(zmpi_spec_elem_t *s, zmpi_spec_elem_t *d, zmpi_spec_tloc_data_t tloc_data, zmpi_spec_elem_t *mod) \
{ \
  zmpi_SPEC_DECLARE_TLOC_MOD_REARRANGE_DB \
  zmpi_SPEC_DO_TLOC_MOD_REARRANGE_DB(_tl_, tloc_data, s, d, mod); \
}

/* sp_macro zmpi_SPEC_DECLARE_TLOC_MOD_REARRANGE_IP */
#define zmpi_SPEC_DECLARE_TLOC_MOD_REARRANGE_IP \
  struct { zmpi_spec_int_t i, p, np; } spec1i;

/* sp_macro zmpi_SPEC_DO_TLOC_MOD_REARRANGE_IP */
#define zmpi_SPEC_DO_TLOC_MOD_REARRANGE_IP(_tl_, _tld_, _b_, _xb_, _ib_)  do { \
  if (_ib_) { \
    for (spec1i.i = 0; spec1i.i < zmpi_spec_elem_get_n(_b_); ++spec1i.i) { \
      spec1i.p = (_tl_)(zmpi_spec_elem_get_buf(_b_), spec1i.i, _tld_, zmpi_spec_elem_get_buf(_ib_)); \
      if (spec1i.p == zmpi_SPEC_LOC_NONE) continue; \
      while (spec1i.i != spec1i.p) { \
        spec1i.np = (_tl_)(zmpi_spec_elem_get_buf(_b_), spec1i.p, _tld_, zmpi_spec_elem_get_buf(_xb_)); \
        if (spec1i.np == zmpi_SPEC_LOC_NONE) break; \
        zmpi_spec_elem_copy_at((_ib_), 0, (_b_), spec1i.p); \
        zmpi_spec_elem_copy_at((_xb_), 0, (_ib_), 0); \
        spec1i.p = spec1i.np; \
      } \
      zmpi_spec_elem_copy_at((_ib_), 0, (_b_), spec1i.i); \
    } \
  } else { \
    for (spec1i.i = 0; spec1i.i < zmpi_spec_elem_get_n(_b_); ++spec1i.i) { \
      spec1i.p = (_tl_)(zmpi_spec_elem_get_buf(_b_), spec1i.i, _tld_, NULL); \
      if (spec1i.p == zmpi_SPEC_LOC_NONE) continue; \
      while (spec1i.i != spec1i.p) { \
        spec1i.np = (_tl_)(zmpi_spec_elem_get_buf(_b_), spec1i.p, _tld_, NULL); \
        if (spec1i.np == zmpi_SPEC_LOC_NONE) { zmpi_spec_elem_copy_at((_b_), spec1i.i, (_b_), spec1i.p); break; } \
        zmpi_spec_elem_exchange_at((_b_), spec1i.i, (_b_), spec1i.p, (_xb_)); \
        spec1i.p = spec1i.np; \
      } \
    } \
 } } while (0) 

/* sp_macro zmpi_SPEC_FUNC_TLOC_MOD_REARRANGE_IP */
#define zmpi_SPEC_FUNC_TLOC_MOD_REARRANGE_IP(_name_, _tl_, _s_) \
_s_ void _name_##_tloc_mod_rearrange_ip(zmpi_spec_elem_t *s, zmpi_spec_elem_t *x, zmpi_spec_tloc_data_t tloc_data, zmpi_spec_elem_t *mod) \
{ \
  zmpi_SPEC_DECLARE_TLOC_MOD_REARRANGE_IP \
  zmpi_SPEC_DO_TLOC_MOD_REARRANGE_IP(_tl_, tloc_data, s, x, mod); \
}

/* sp_macro zmpi_SPEC_DEFINE_TLOC */
#define zmpi_SPEC_DEFINE_TLOC(_name_, _tl_, _s_...) \
  zmpi_SPEC_FUNC_TLOC_REARRANGE_DB(_name_, _tl_, _s_) \
  zmpi_SPEC_FUNC_TLOC_REARRANGE_IP(_name_, _tl_, _s_)

/* sp_macro zmpi_SPEC_DEFINE_TLOC_MOD */
#define zmpi_SPEC_DEFINE_TLOC_MOD(_name_, _tl_, _s_...) \
  zmpi_SPEC_FUNC_TLOC_MOD_REARRANGE_DB(_name_, _tl_, _s_) \
  zmpi_SPEC_FUNC_TLOC_MOD_REARRANGE_IP(_name_, _tl_, _s_)

/* sp_macro zmpi_SPEC_EXT_PARAM_TLOC zmpi_SPEC_EXT_PARAM_TLOC_NULL zmpi_SPEC_EXT_PARAM_TLOC_MOD zmpi_SPEC_EXT_PARAM_TLOC_MOD_NULL */
#define zmpi_SPEC_EXT_PARAM_TLOC(_name_)      _name_##_tloc_rearrange_db, _name_##_tloc_rearrange_ip
#define zmpi_SPEC_EXT_PARAM_TLOC_NULL         NULL, NULL
#define zmpi_SPEC_EXT_PARAM_TLOC_MOD(_name_)  _name_##_tloc_mod_rearrange_db, _name_##_tloc_mod_rearrange_ip
#define zmpi_SPEC_EXT_PARAM_TLOC_MOD_NULL     NULL, NULL


/* sp_type zmpi_spec_tloc_f zmpi_spec_tloc_rearrange_db_f zmpi_spec_tloc_rearrange_ip_f */
typedef zmpi_spec_elem_index_t zmpi_spec_tloc_f(zmpi_spec_elem_buf_t b, zmpi_spec_elem_index_t x, zmpi_spec_tloc_data_t tloc_data);
typedef void zmpi_spec_tloc_rearrange_db_f(zmpi_spec_elem_t *s, zmpi_spec_elem_t *d, zmpi_spec_tloc_data_t tloc_data);
typedef void zmpi_spec_tloc_rearrange_ip_f(zmpi_spec_elem_t *s, zmpi_spec_elem_t *x, zmpi_spec_tloc_data_t tloc_data);

/* sp_type zmpi_spec_tloc_mod_f zmpi_spec_tloc_mod_rearrange_db_f zmpi_spec_tloc_mod_rearrange_ip_f */
typedef zmpi_spec_elem_index_t zmpi_spec_tloc_mod_f(zmpi_spec_elem_buf_t b, zmpi_spec_elem_index_t x, zmpi_spec_tloc_data_t tloc_data, zmpi_spec_elem_buf_t mod);
typedef void zmpi_spec_tloc_mod_rearrange_db_f(zmpi_spec_elem_t *s, zmpi_spec_elem_t *d, zmpi_spec_tloc_data_t tloc_data, zmpi_spec_elem_t *mod);
typedef void zmpi_spec_tloc_mod_rearrange_ip_f(zmpi_spec_elem_t *s, zmpi_spec_elem_t *x, zmpi_spec_tloc_data_t tloc_data, zmpi_spec_elem_t *mod);


#endif /* zmpi_SPEC_TLOC */




#ifdef ZMPI_PREFIX
# include "zmpi_atasp_rename.h"
#endif


#if MPI_VERSION >= 3
# define IF_ELSE_MPI_VERSION_3(_if_, _else_)  _if_
#else
# define IF_ELSE_MPI_VERSION_3(_if_, _else_)  _else_
#endif

typedef struct zmpi__spec_tproc_t *ZMPI_Tproc;

#define ZMPI_TPROC_NULL  NULL

#if MPI_VERSION < 3
typedef long ZMPI_Count;
#else
typedef MPI_Count ZMPI_Count;
#endif

typedef int ZMPI_TPROC_FN(void *b, IF_ELSE_MPI_VERSION_3(MPI_Count, ZMPI_Count) x, void *tproc_data);
typedef int ZMPI_TPROC_MOD_FN(void *b, IF_ELSE_MPI_VERSION_3(MPI_Count, ZMPI_Count) x, void *tproc_data, void *mod);
typedef void ZMPI_TPROCS_FN(void *b, IF_ELSE_MPI_VERSION_3(MPI_Count, ZMPI_Count) x, void *tproc_data, int *nprocs, int *procs);
typedef void ZMPI_TPROCS_MOD_FN(void *b, IF_ELSE_MPI_VERSION_3(MPI_Count, ZMPI_Count) x, void *tproc_data, int *nprocs, int *procs, void *mod);

typedef void ZMPI_TPROC_RESET_FN(void *tproc_data);

#define ZMPI_TPROC_RESET_NULL  NULL

typedef struct _ZMPI_Tproc_exdef
{
  int type;

  zmpi_spec_tproc_ext_t tproc_ext;
  zmpi_spec_tproc_mod_ext_t tproc_mod_ext;
  zmpi_spec_tprocs_ext_t tprocs_ext;
  zmpi_spec_tprocs_mod_ext_t tprocs_mod_ext;

} const *ZMPI_Tproc_exdef;

#define ZMPI_TPROC_EXDEF_NULL  NULL

/* default */
#define ZMPI_TPROC_EXDEF_DEFINE_TPROC(_name_, _tp_, _s_...) \
  zmpi_SPEC_DEFINE_TPROC(_name_, _tp_, _s_) \
  _s_ const struct _ZMPI_Tproc_exdef _##_name_ = { 1, zmpi_SPEC_EXT_PARAM_TPROC(_name_), zmpi_SPEC_EXT_PARAM_TPROC_MOD_NULL, zmpi_SPEC_EXT_PARAM_TPROCS_NULL, zmpi_SPEC_EXT_PARAM_TPROCS_MOD_NULL }, *_name_ = &_##_name_;

#define ZMPI_TPROC_EXDEF_DEFINE_TPROC_MOD(_name_, _tp_, _s_...) \
  zmpi_SPEC_DEFINE_TPROC_MOD(_name_, _tp_, _s_) \
  _s_ const struct _ZMPI_Tproc_exdef _##_name_ = { 2, zmpi_SPEC_EXT_PARAM_TPROC_NULL, zmpi_SPEC_EXT_PARAM_TPROC_MOD(_name_), zmpi_SPEC_EXT_PARAM_TPROCS_NULL, zmpi_SPEC_EXT_PARAM_TPROCS_MOD_NULL }, *_name_ = &_##_name_;

#define ZMPI_TPROC_EXDEF_DEFINE_TPROCS(_name_, _tp_, _s_...) \
  zmpi_SPEC_DEFINE_TPROCS(_name_, _tp_, _s_) \
  _s_ const struct _ZMPI_Tproc_exdef _##_name_ = { 3, zmpi_SPEC_EXT_PARAM_TPROC_NULL, zmpi_SPEC_EXT_PARAM_TPROC_MOD_NULL, zmpi_SPEC_EXT_PARAM_TPROCS(_name_), zmpi_SPEC_EXT_PARAM_TPROCS_MOD_NULL }, *_name_ = &_##_name_;

#define ZMPI_TPROC_EXDEF_DEFINE_TPROCS_MOD(_name_, _tp_, _s_...) \
  zmpi_SPEC_DEFINE_TPROCS_MOD(_name_, _tp_, _s_) \
  _s_ const struct _ZMPI_Tproc_exdef _##_name_ = { 4, zmpi_SPEC_EXT_PARAM_TPROC_NULL, zmpi_SPEC_EXT_PARAM_TPROC_MOD_NULL, zmpi_SPEC_EXT_PARAM_TPROCS_NULL, zmpi_SPEC_EXT_PARAM_TPROCS_MOD(_name_) }, *_name_ = &_##_name_;

/* fixtype */
#define ZMPI_TPROC_EXDEF_DEFINE_FIXTYPE_TPROC(_name_, _fxt_, _tp_, _s_...) \
  zmpi_SPEC_DEFINE_FIXED_TPROC(_name_, ZMPI_FIXTYPE_DECLARE, _fxt_, ZMPI_FIXTYPE_CREATE, ZMPI_FIXTYPE_COPY_AT, ZMPI_FIXTYPE_EXCHANGE_AT, ZMPI_FIXTYPE_DESTROY, _tp_, _s_) \
  _s_ const struct _ZMPI_Tproc_exdef _##_name_ = { 1, zmpi_SPEC_EXT_PARAM_TPROC(_name_), zmpi_SPEC_EXT_PARAM_TPROC_MOD_NULL, zmpi_SPEC_EXT_PARAM_TPROCS_NULL, zmpi_SPEC_EXT_PARAM_TPROCS_MOD_NULL }, *_name_ = &_##_name_;

#define ZMPI_TPROC_EXDEF_DEFINE_FIXTYPE_TPROC_MOD(_name_, _fxt_, _tp_, _s_...) \
  zmpi_SPEC_DEFINE_FIXED_TPROC_MOD(_name_, ZMPI_FIXTYPE_DECLARE, _fxt_, ZMPI_FIXTYPE_CREATE, ZMPI_FIXTYPE_COPY_AT, ZMPI_FIXTYPE_EXCHANGE_AT, ZMPI_FIXTYPE_DESTROY, _tp_, _s_) \
  _s_ const struct _ZMPI_Tproc_exdef _##_name_ = { 2, zmpi_SPEC_EXT_PARAM_TPROC_NULL, zmpi_SPEC_EXT_PARAM_TPROC_MOD(_name_), zmpi_SPEC_EXT_PARAM_TPROCS_NULL, zmpi_SPEC_EXT_PARAM_TPROCS_MOD_NULL }, *_name_ = &_##_name_;

#define ZMPI_TPROC_EXDEF_DEFINE_FIXTYPE_TPROCS(_name_, _fxt_, _tp_, _s_...) \
  zmpi_SPEC_DEFINE_FIXED_TPROCS(_name_, ZMPI_FIXTYPE_DECLARE, _fxt_, ZMPI_FIXTYPE_CREATE, ZMPI_FIXTYPE_COPY_AT, ZMPI_FIXTYPE_EXCHANGE_AT, ZMPI_FIXTYPE_DESTROY, _tp_, _s_) \
  _s_ const struct _ZMPI_Tproc_exdef _##_name_ = { 3, zmpi_SPEC_EXT_PARAM_TPROC_NULL, zmpi_SPEC_EXT_PARAM_TPROC_MOD_NULL, zmpi_SPEC_EXT_PARAM_TPROCS(_name_), zmpi_SPEC_EXT_PARAM_TPROCS_MOD_NULL }, *_name_ = &_##_name_;

#define ZMPI_TPROC_EXDEF_DEFINE_FIXTYPE_TPROCS_MOD(_name_, _fxt_, _tp_, _s_...) \
  zmpi_SPEC_DEFINE_FIXED_TPROCS_MOD(_name_, ZMPI_FIXTYPE_DECLARE, _fxt_, ZMPI_FIXTYPE_CREATE, ZMPI_FIXTYPE_COPY_AT, ZMPI_FIXTYPE_EXCHANGE_AT, ZMPI_FIXTYPE_DESTROY, _tp_, _s_) \
  _s_ const struct _ZMPI_Tproc_exdef _##_name_ = { 4, zmpi_SPEC_EXT_PARAM_TPROC_NULL, zmpi_SPEC_EXT_PARAM_TPROC_MOD_NULL, zmpi_SPEC_EXT_PARAM_TPROCS_NULL, zmpi_SPEC_EXT_PARAM_TPROCS_MOD(_name_) }, *_name_ = &_##_name_;

/* fixsize */
#define ZMPI_TPROC_EXDEF_DEFINE_FIXSIZE_TPROC(_name_, _fxs_, _tp_, _s_...) \
  _s_ const int _name_##_params = _fxs_; \
  zmpi_SPEC_DEFINE_FIXED_TPROC(_name_, ZMPI_FIXSIZE_DECLARE, _name_##_params, ZMPI_FIXSIZE_CREATE, ZMPI_FIXSIZE_COPY_AT, ZMPI_FIXSIZE_EXCHANGE_AT, ZMPI_FIXSIZE_DESTROY, _tp_, _s_) \
  _s_ const struct _ZMPI_Tproc_exdef _##_name_ = { 1, zmpi_SPEC_EXT_PARAM_TPROC(_name_), zmpi_SPEC_EXT_PARAM_TPROC_MOD_NULL, zmpi_SPEC_EXT_PARAM_TPROCS_NULL, zmpi_SPEC_EXT_PARAM_TPROCS_MOD_NULL }, *_name_ = &_##_name_;

int ZMPI_Tproc_create_tproc(ZMPI_Tproc *tproc, ZMPI_TPROC_FN *tfn, ZMPI_TPROC_RESET_FN *rfn, const ZMPI_Tproc_exdef exdef);
int ZMPI_Tproc_create_tproc_mod(ZMPI_Tproc *tproc, ZMPI_TPROC_MOD_FN *tfn, ZMPI_TPROC_RESET_FN *rfn, const ZMPI_Tproc_exdef exdef);
int ZMPI_Tproc_create_tprocs(ZMPI_Tproc *tproc, int max_tprocs, ZMPI_TPROCS_FN *tfn, ZMPI_TPROC_RESET_FN *rfn, const ZMPI_Tproc_exdef exdef);
int ZMPI_Tproc_create_tprocs_mod(ZMPI_Tproc *tproc, int max_tprocs, ZMPI_TPROCS_MOD_FN *tfn, ZMPI_TPROC_RESET_FN *rfn, const ZMPI_Tproc_exdef exdef);
int ZMPI_Tproc_free(ZMPI_Tproc *tproc);

int ZMPI_Tproc_set_neighbors(ZMPI_Tproc tproc, int nneighbors, int *neighbors, MPI_Comm comm);
int ZMPI_Tproc_set_proclists(ZMPI_Tproc tproc, int ndstprocs, int *dstprocs, int nsrcprocs, int *srcprocs, MPI_Comm comm);

#define ZMPI_ALLTOALL_SPECIFIC_TYPE_ALLTOALLV    0
#define ZMPI_ALLTOALL_SPECIFIC_TYPE_ALLTOALLW    1
#define ZMPI_ALLTOALL_SPECIFIC_TYPE_PUT          2
#define ZMPI_ALLTOALL_SPECIFIC_TYPE_PUT_2PHASES  3
#define ZMPI_ALLTOALL_SPECIFIC_TYPE_SENDRECV     4

#define ZMPI_ALLTOALL_SPECIFIC_TYPE_DEFAULT  ZMPI_ALLTOALL_SPECIFIC_TYPE_ALLTOALLV

extern int ZMPI_Alltoall_specific_type;

#define ZMPI_NEIGHBOR_ALLTOALL_SPECIFIC_TYPE_DEFAULT  ZMPI_ALLTOALL_SPECIFIC_TYPE_ALLTOALLV

extern int ZMPI_Neighbor_alltoall_specific_type;

#if MPI_VERSION < 3
typedef int ZMPI_Status;
int ZMPI_Get_elements(const ZMPI_Status *status, MPI_Datatype datatype, int *count);
# define ZMPI_STATUS_IGNORE  NULL
#else
typedef MPI_Status ZMPI_Status;
# define ZMPI_Get_elements   MPI_Get_elements
# define ZMPI_STATUS_IGNORE  MPI_STATUS_IGNORE
#endif

typedef int ZMPI_ALLTOALL_SPECIFIC_FN(void *sbuf, int scount, MPI_Datatype stype, void *rbuf, int rcount, MPI_Datatype rtype, ZMPI_Tproc tproc, void *tproc_data, MPI_Comm comm, IF_ELSE_MPI_VERSION_3(MPI_Status, ZMPI_Status) *status);
int ZMPI_Alltoall_specific(void *sbuf, int scount, MPI_Datatype stype, void *rbuf, int rcount, MPI_Datatype rtype, ZMPI_Tproc tproc, void *tproc_data, MPI_Comm comm, IF_ELSE_MPI_VERSION_3(MPI_Status, ZMPI_Status) *status);
int ZMPI_Neighbor_alltoall_specific(void *sbuf, int scount, MPI_Datatype stype, void *rbuf, int rcount, MPI_Datatype rtype, ZMPI_Tproc tproc, void *tproc_data, MPI_Comm comm, IF_ELSE_MPI_VERSION_3(MPI_Status, ZMPI_Status) *status);

#undef IF_ELSE_MPI_VERSION_3


#endif /* __ZMPI_ATASP_H__ */
