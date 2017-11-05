/*
 *  Copyright (C) 2011, 2012, 2013, 2014, 2015, 2016 Michael Hofmann
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


#ifndef __SPEC_PUBLIC_H__
#define __SPEC_PUBLIC_H__


/* tproc count */

/* sp_macro SPEC_DECLARE_TPROC_COUNT_DB */
#define SPEC_DECLARE_TPROC_COUNT_DB \
  struct { spec_elem_index_t n, i; spec_proc_t p; } spec0cd;

/* sp_macro SPEC_DO_TPROC_COUNT_DB */
#define SPEC_DO_TPROC_COUNT_DB(_tp_, _tpd_, _b_, _cs_)  do { \
  spec0cd.n = spec_elem_get_n(_b_); \
  for (spec0cd.i = 0; spec0cd.i < spec0cd.n; ++spec0cd.i) { \
    spec0cd.p = (_tp_)(spec_elem_get_buf(_b_), spec0cd.i, _tpd_); \
    if (spec0cd.p == SPEC_PROC_NONE) continue; \
    ++(_cs_)[spec0cd.p]; \
  } } while (0)

/* sp_macro SPEC_FUNC_TPROC_COUNT_DB */
#define SPEC_FUNC_TPROC_COUNT_DB(_name_, _tp_, _s_...) \
_s_ void _name_##_tproc_count_db(spec_tproc_data_t tproc_data, spec_elem_t *s, int *counts) \
{ \
  SPEC_DECLARE_TPROC_COUNT_DB \
  SPEC_DO_TPROC_COUNT_DB(_tp_, tproc_data, s, counts); \
}

/* sp_macro SPEC_DECLARE_TPROC_COUNT_IP */
#define SPEC_DECLARE_TPROC_COUNT_IP \
  struct { spec_elem_index_t n, t, i; spec_proc_t p; } spec0ci;

/* sp_macro SPEC_DO_TPROC_COUNT_IP */
#define SPEC_DO_TPROC_COUNT_IP(_tp_, _tpd_, _b_, _cs_)  do { \
  spec0ci.n = spec_elem_get_n(_b_); \
  spec0ci.t = 0; \
  for (spec0ci.i = 0; spec0ci.i < spec0ci.n; ++spec0ci.i) { \
    spec0ci.p = (_tp_)(spec_elem_get_buf(_b_), spec0ci.i, _tpd_); \
    if (spec0ci.p == SPEC_PROC_NONE) continue; \
    ++(_cs_)[spec0ci.p]; \
    if (spec0ci.t < spec0ci.i) spec_elem_copy_at((_b_), spec0ci.i, (_b_), spec0ci.t); \
    ++spec0ci.t; \
  } \
  spec_elem_set_n(_b_, spec0ci.t); \
} while (0)

/* sp_macro SPEC_FUNC_TPROC_COUNT_IP */
#define SPEC_FUNC_TPROC_COUNT_IP(_name_, _tp_, _s_...) \
_s_ void _name_##_tproc_count_ip(spec_tproc_data_t tproc_data, spec_elem_t *s, int *counts) \
{ \
  SPEC_DECLARE_TPROC_COUNT_IP \
  SPEC_DO_TPROC_COUNT_IP(_tp_, tproc_data, s, counts); \
}


/* tproc_mod count */

/* sp_macro SPEC_DECLARE_TPROC_MOD_COUNT_DB */
#define SPEC_DECLARE_TPROC_MOD_COUNT_DB \
  struct { spec_elem_index_t n, i; spec_proc_t p; } spec1cd;

/* sp_macro SPEC_DO_TPROC_MOD_COUNT_DB */
#define SPEC_DO_TPROC_MOD_COUNT_DB(_tp_, _tpd_, _b_, _cs_)  do { \
  spec1cd.n = spec_elem_get_n(_b_); \
  for (spec1cd.i = 0; spec1cd.i < spec1cd.n; ++spec1cd.i) { \
    spec1cd.p = (_tp_)(spec_elem_get_buf(_b_), spec1cd.i, _tpd_, NULL); \
    if (spec1cd.p == SPEC_PROC_NONE) continue; \
    ++(_cs_)[spec1cd.p]; \
  } } while (0)

/* sp_macro SPEC_FUNC_TPROC_MOD_COUNT_DB */
#define SPEC_FUNC_TPROC_MOD_COUNT_DB(_name_, _tp_, _s_...) \
_s_ void _name_##_tproc_mod_count_db(spec_tproc_data_t tproc_data, spec_elem_t *s, int *counts) \
{ \
  SPEC_DECLARE_TPROC_MOD_COUNT_DB \
  SPEC_DO_TPROC_MOD_COUNT_DB(_tp_, tproc_data, s, counts); \
}

/* sp_macro SPEC_DECLARE_TPROC_MOD_COUNT_IP */
#define SPEC_DECLARE_TPROC_MOD_COUNT_IP \
  struct { spec_elem_index_t n, t, i; spec_proc_t p; } spec1ci;

/* sp_macro SPEC_DO_TPROC_MOD_COUNT_IP */
#define SPEC_DO_TPROC_MOD_COUNT_IP(_tp_, _tpd_, _b_, _cs_)  do { \
  spec1ci.n = spec_elem_get_n(_b_); \
  spec1ci.t = 0; \
  for (spec1ci.i = 0; spec1ci.i < spec1ci.n; ++spec1ci.i) { \
    spec1ci.p = (_tp_)(spec_elem_get_buf(_b_), spec1ci.i, _tpd_, NULL); \
    if (spec1ci.p == SPEC_PROC_NONE) continue; \
    ++(_cs_)[spec1ci.p]; \
    if (spec1ci.t < spec1ci.i) spec_elem_copy_at((_b_), spec1ci.i, (_b_), spec1ci.t); \
    ++spec1ci.t; \
  } \
  spec_elem_set_n(_b_, spec1ci.t); \
} while (0)

/* sp_macro SPEC_FUNC_TPROC_MOD_COUNT_IP */
#define SPEC_FUNC_TPROC_MOD_COUNT_IP(_name_, _tp_, _s_...) \
_s_ void _name_##_tproc_mod_count_ip(spec_tproc_data_t tproc_data, spec_elem_t *s, int *counts) \
{ \
  SPEC_DECLARE_TPROC_MOD_COUNT_IP \
  SPEC_DO_TPROC_MOD_COUNT_IP(_tp_, tproc_data, s, counts); \
}


/* tprocs count */

/* sp_macro SPEC_DECLARE_TPROCS_COUNT_DB */
#define SPEC_DECLARE_TPROCS_COUNT_DB \
  struct { spec_elem_index_t n, i; spec_int_t j, m; } spec2cd;

/* sp_macro SPEC_DO_TPROCS_COUNT_DB */
#define SPEC_DO_TPROCS_COUNT_DB(_tp_, _tpd_, _b_, _cs_, _ps_)  do { \
  spec2cd.n = spec_elem_get_n(_b_); \
  for (spec2cd.i = 0; spec2cd.i < spec2cd.n; ++spec2cd.i) { \
    (_tp_)(spec_elem_get_buf(_b_), spec2cd.i, (_tpd_), &spec2cd.m, (_ps_)); \
    for (spec2cd.j = 0; spec2cd.j < spec2cd.m; ++spec2cd.j) ++(_cs_)[(_ps_)[spec2cd.j]]; \
  } } while (0)

/* sp_macro SPEC_FUNC_TPROCS_COUNT_DB */
#define SPEC_FUNC_TPROCS_COUNT_DB(_name_, _tp_, _s_...) \
_s_ void _name_##_tprocs_count_db(spec_tproc_data_t tproc_data, spec_elem_t *s, int *counts, spec_proc_t *procs) \
{ \
  SPEC_DECLARE_TPROCS_COUNT_DB \
  SPEC_DO_TPROCS_COUNT_DB(_tp_, tproc_data, s, counts, procs); \
}

/* sp_macro SPEC_DECLARE_TPROCS_COUNT_IP */
#define SPEC_DECLARE_TPROCS_COUNT_IP \
  struct { spec_elem_index_t n, t, i; spec_int_t j, m; } spec2ci;

/* sp_macro SPEC_DO_TPROCS_COUNT_IP */
#define SPEC_DO_TPROCS_COUNT_IP(_tp_, _tpd_, _b_, _cs_, _ps_)  do { \
  spec2ci.n = spec_elem_get_n(_b_); \
  spec2ci.t = 0; \
  for (spec2ci.i = 0; spec2ci.i < spec2ci.n; ++spec2ci.i) { \
    (_tp_)(spec_elem_get_buf(_b_), spec2ci.i, (_tpd_), &spec2ci.m, (_ps_)); \
    if (spec2ci.m <= 0) continue; \
    for (spec2ci.j = 0; spec2ci.j < spec2ci.m; ++spec2ci.j) ++(_cs_)[(_ps_)[spec2ci.j]]; \
    if (spec2ci.t < spec2ci.i) spec_elem_copy_at((_b_), spec2ci.i, (_b_), spec2ci.t); \
    ++spec2ci.t; \
  } \
  spec_elem_set_n(_b_, spec2ci.t); \
} while (0)

/* sp_macro SPEC_FUNC_TPROCS_COUNT_IP */
#define SPEC_FUNC_TPROCS_COUNT_IP(_name_, _tp_, _s_...) \
_s_ void _name_##_tprocs_count_ip(spec_tproc_data_t tproc_data, spec_elem_t *s, int *counts, spec_proc_t *procs) \
{ \
  SPEC_DECLARE_TPROCS_COUNT_IP \
  SPEC_DO_TPROCS_COUNT_IP(_tp_, tproc_data, s, counts, procs); \
}


/* tprocs_mod count */

/* sp_macro SPEC_DECLARE_TPROCS_MOD_COUNT_DB */
#define SPEC_DECLARE_TPROCS_MOD_COUNT_DB \
  struct { spec_elem_index_t n, i; spec_int_t j, m; } spec3cd;

/* sp_macro SPEC_DO_TPROCS_MOD_COUNT_DB */
#define SPEC_DO_TPROCS_MOD_COUNT_DB(_tp_, _tpd_, _b_, _cs_, _ps_)  do { \
  spec3cd.n = spec_elem_get_n(_b_); \
  for (spec3cd.i = 0; spec3cd.i < spec3cd.n; ++spec3cd.i) { \
    (_tp_)(spec_elem_get_buf(_b_), spec3cd.i, (_tpd_), &spec3cd.m, (_ps_), NULL); \
    for (spec3cd.j = 0; spec3cd.j < spec3cd.m; ++spec3cd.j) ++(_cs_)[(_ps_)[spec3cd.j]]; \
  } } while (0)

/* sp_macro SPEC_FUNC_TPROCS_MOD_COUNT_DB */
#define SPEC_FUNC_TPROCS_MOD_COUNT_DB(_name_, _tp_, _s_...) \
_s_ void _name_##_tprocs_mod_count_db(spec_tproc_data_t tproc_data, spec_elem_t *s, int *counts, spec_proc_t *procs) \
{ \
  SPEC_DECLARE_TPROCS_MOD_COUNT_DB \
  SPEC_DO_TPROCS_MOD_COUNT_DB(_tp_, tproc_data, s, counts, procs); \
}

/* sp_macro SPEC_DECLARE_TPROCS_MOD_COUNT_IP */
#define SPEC_DECLARE_TPROCS_MOD_COUNT_IP \
  struct { spec_elem_index_t n, t, i; spec_int_t j, m; } spec3ci;

/* sp_macro SPEC_DO_TPROCS_MOD_COUNT_IP */
#define SPEC_DO_TPROCS_MOD_COUNT_IP(_tp_, _tpd_, _b_, _cs_, _ps_)  do { \
  spec3ci.n = spec_elem_get_n(_b_); \
  spec3ci.t = 0; \
  for (spec3ci.i = 0; spec3ci.i < spec3ci.n; ++spec3ci.i) { \
    (_tp_)(spec_elem_get_buf(_b_), spec3ci.i, (_tpd_), &spec3ci.m, (_ps_), NULL); \
    if (spec3ci.m <= 0) continue; \
    for (spec3ci.j = 0; spec3ci.j < spec3ci.m; ++spec3ci.j) ++(_cs_)[(_ps_)[spec3ci.j]]; \
    if (spec3ci.t < spec3ci.i) spec_elem_copy_at((_b_), spec3ci.i, (_b_), spec3ci.t); \
    ++spec3ci.t; \
  } \
  spec_elem_set_n(_b_, spec3ci.t); \
} while (0)

/* sp_macro SPEC_FUNC_TPROCS_MOD_COUNT_IP */
#define SPEC_FUNC_TPROCS_MOD_COUNT_IP(_name_, _tp_, _s_...) \
_s_ void _name_##_tprocs_mod_count_ip(spec_tproc_data_t tproc_data, spec_elem_t *s, int *counts, spec_proc_t *procs) \
{ \
  SPEC_DECLARE_TPROCS_MOD_COUNT_IP \
  SPEC_DO_TPROCS_MOD_COUNT_IP(_tp_, tproc_data, s, counts, procs); \
}


/* un-fixed macros, sp_macro spec_fixed_default_declare spec_fixed_default_create spec_fixed_default_copy_at spec_fixed_default_exchange_at spec_fixed_default_destroy */
#define spec_fixed_default_declare(_fx_, _fxp_)
#define spec_fixed_default_create(_fx_, _fxp_)
#define spec_fixed_default_copy_at(_se_, _sat_, _de_, _dat_, _fx_, _fxp_)             spec_elem_copy_at(_se_, _sat_, _de_, _dat_)
#define spec_fixed_default_exchange_at(_s0_, _s0at_, _s1_, _s1at_, _t_, _fx_, _fxp_)  spec_elem_exchange_at(_s0_, _s0at_, _s1_, _s1at_, _t_)
#define spec_fixed_default_destroy(_fx_, _fxp_)


/* tproc rearrange */

/* sp_macro SPEC_DECLARE_FIXED_TPROC_REARRANGE_DB */
#define SPEC_DECLARE_FIXED_TPROC_REARRANGE_DB(_fxdcl_, _fxp_) \
  struct { spec_elem_index_t n, i; spec_proc_t p; _fxdcl_(fx, _fxp_) } spec0d;

/* sp_macro SPEC_DO_FIXED_TPROC_REARRANGE_DB */
#define SPEC_DO_FIXED_TPROC_REARRANGE_DB(_fxp_, _fxc_, _fxca_, _fxxa_, _fxd_, _tp_, _tpd_, _sb_, _db_, _ds_)  do { \
  _fxc_(spec0d.fx, _fxp_); \
  spec0d.n = spec_elem_get_n(_sb_); \
  for (spec0d.i = 0; spec0d.i < spec0d.n; ++spec0d.i) { \
    spec0d.p = (_tp_)(spec_elem_get_buf(_sb_), spec0d.i, _tpd_); \
    if (spec0d.p == SPEC_PROC_NONE) continue; \
    _fxca_((_sb_), spec0d.i, (_db_), (_ds_)[spec0d.p], spec0d.fx, _fxp_); \
    ++(_ds_)[spec0d.p]; \
  } \
  _fxd_(spec0d.fx, _fxp_); \
  } while (0)

/* sp_macro SPEC_FUNC_FIXED_TPROC_REARRANGE_DB */
#define SPEC_FUNC_FIXED_TPROC_REARRANGE_DB(_name_, _fxdcl_, _fxp_, _fxc_, _fxca_, _fxxa_, _fxd_, _tp_, _s_...) \
_s_ void _name_##_tproc_rearrange_db(spec_tproc_data_t tproc_data, spec_elem_t *s, spec_elem_t *d, int *displs) \
{ \
  SPEC_DECLARE_FIXED_TPROC_REARRANGE_DB(_fxdcl_, _fxp_) \
  SPEC_DO_FIXED_TPROC_REARRANGE_DB(_fxp_, _fxc_, _fxca_, _fxxa_, _fxd_, _tp_, tproc_data, s, d, displs); \
}

/* sp_macro SPEC_DECLARE_TPROC_REARRANGE_DB */
#define SPEC_DECLARE_TPROC_REARRANGE_DB \
  SPEC_DECLARE_FIXED_TPROC_REARRANGE_DB(spec_fixed_default_declare, NOPARAM)

/* sp_macro SPEC_DO_TPROC_REARRANGE_DB */
#define SPEC_DO_TPROC_REARRANGE_DB(_tp_, _tpd_, _sb_, _db_, _ds_) \
  SPEC_DO_FIXED_TPROC_REARRANGE_DB(NOPARAM, spec_fixed_default_create, spec_fixed_default_copy_at, spec_fixed_default_exchange_at, spec_fixed_default_destroy, _tp_, _tpd_, _sb_, _db_, _ds_)

/* sp_macro SPEC_FUNC_TPROC_REARRANGE_DB */
#define SPEC_FUNC_TPROC_REARRANGE_DB(_name_, _tp_, _s_...) \
  SPEC_FUNC_FIXED_TPROC_REARRANGE_DB(_name_, spec_fixed_default_declare, NOPARAM, spec_fixed_default_create, spec_fixed_default_copy_at, spec_fixed_default_exchange_at, spec_fixed_default_destroy, _tp_, _s_)

/* sp_macro SPEC_DECLARE_TPROC_REARRANGE_IP */
#define SPEC_DECLARE_TPROC_REARRANGE_IP \
  struct { spec_elem_index_t e, i, j; spec_proc_t p, np; } spec0i;

/* sp_macro SPEC_DO_TPROC_REARRANGE_IP */
#define SPEC_DO_TPROC_REARRANGE_IP(_tp_, _tpd_, _b_, _xb_, _ds_, _cs_, _n_)  do { \
  for (spec0i.e = 0, spec0i.i = 0; spec0i.i < (_n_); ++spec0i.i) { \
    spec0i.e += (_cs_)[spec0i.i]; \
    spec0i.j = (_ds_)[spec0i.i]; \
    while (spec0i.j < spec0i.e) { \
      spec0i.p = (_tp_)(spec_elem_get_buf(_b_), spec0i.j, _tpd_); \
      while (spec0i.p != spec0i.i) { \
        spec0i.np = (_tp_)(spec_elem_get_buf(_b_), (_ds_)[spec0i.p], _tpd_); \
        if (spec0i.np != spec0i.p) spec_elem_exchange_at((_b_), (_ds_)[spec0i.p], (_b_), spec0i.j, (_xb_)); \
        ++(_ds_)[spec0i.p]; \
        spec0i.p = spec0i.np; \
      } \
      ++spec0i.j; \
    } \
  } } while (0)

/* sp_macro SPEC_FUNC_TPROC_REARRANGE_IP */
#define SPEC_FUNC_TPROC_REARRANGE_IP(_name_, _tp_, _s_) \
_s_ void _name_##_tproc_rearrange_ip(spec_tproc_data_t tproc_data, spec_elem_t *s, spec_elem_t *x, int *displs, int *counts, spec_int_t n) \
{ \
  SPEC_DECLARE_TPROC_REARRANGE_IP \
  SPEC_DO_TPROC_REARRANGE_IP(_tp_, tproc_data, s, x, displs, counts, n); \
}


/* tproc_mod rearrange */

/* sp_macro SPEC_DECLARE_FIXED_TPROC_MOD_REARRANGE_DB */
#define SPEC_DECLARE_FIXED_TPROC_MOD_REARRANGE_DB(_fxdcl_, _fxp_) \
  struct { spec_elem_index_t n, i; spec_proc_t p; _fxdcl_(fx, _fxp_) } spec1d;

/* sp_macro SPEC_DO_FIXED_TPROC_MOD_REARRANGE_DB */
#define SPEC_DO_FIXED_TPROC_MOD_REARRANGE_DB(_fxp_, _fxc_, _fxca_, _fxxa_, _fxd_, _tp_, _tpd_, _sb_, _db_, _ds_, _ib_)  do { \
  spec1d.n = spec_elem_get_n(_sb_); \
  _fxc_(spec0d.fx, _fxp_); \
  for (spec1d.i = 0; spec1d.i < spec1d.n; ++spec1d.i) { \
    spec1d.p = (_tp_)(spec_elem_get_buf(_sb_), spec1d.i, _tpd_, spec_elem_get_buf(_ib_)); \
    if (spec1d.p == SPEC_PROC_NONE) continue; \
    _fxca_((_ib_), 0, (_db_), (_ds_)[spec1d.p], spec1d.fx, _fxp_); \
    ++(_ds_)[spec1d.p]; \
  } \
  _fxd_(spec0d.fx, _fxp_); \
  } while (0)

/* sp_macro SPEC_FUNC_FIXED_TPROC_MOD_REARRANGE_DB */
#define SPEC_FUNC_FIXED_TPROC_MOD_REARRANGE_DB(_name_, _fxdcl_, _fxp_, _fxc_, _fxca_, _fxxa_, _fxd_, _tp_, _s_...) \
_s_ void _name_##_tproc_mod_rearrange_db(spec_tproc_data_t tproc_data, spec_elem_t *s, spec_elem_t *d, int *displs, spec_elem_t *mod) \
{ \
  SPEC_DECLARE_FIXED_TPROC_MOD_REARRANGE_DB(_fxdcl_, _fxp_) \
  SPEC_DO_FIXED_TPROC_MOD_REARRANGE_DB(_fxp_, _fxc_, _fxca_, _fxxa_, _fxd_, _tp_, tproc_data, s, d, displs, mod); \
}

/* sp_macro SPEC_DECLARE_TPROC_MOD_REARRANGE_DB */
#define SPEC_DECLARE_TPROC_MOD_REARRANGE_DB \
  SPEC_DECLARE_FIXED_TPROC_MOD_REARRANGE_DB(spec_fixed_default_declare, NOPARAM)

/* sp_macro SPEC_DO_TPROC_MOD_REARRANGE_DB */
#define SPEC_DO_TPROC_MOD_REARRANGE_DB(_tp_, _tpd_, _sb_, _db_, _ds_, _ib_) \
  SPEC_DO_FIXED_TPROC_MOD_REARRANGE_DB(NOPARAM, spec_fixed_default_create, spec_fixed_default_copy_at, spec_fixed_default_exchange_at, spec_fixed_default_destroy, _tp_, _tpd_, _sb_, _db_, _ds_, _ib_)

/* sp_macro SPEC_FUNC_TPROC_MOD_REARRANGE_DB */
#define SPEC_FUNC_TPROC_MOD_REARRANGE_DB(_name_, _tp_, _s_...) \
  SPEC_FUNC_FIXED_TPROC_MOD_REARRANGE_DB(_name_, spec_fixed_default_declare, NOPARAM, spec_fixed_default_create, spec_fixed_default_copy_at, spec_fixed_default_exchange_at, spec_fixed_default_destroy, _tp_, _s_)

/* sp_macro SPEC_DECLARE_TPROC_MOD_REARRANGE_IP */
#define SPEC_DECLARE_TPROC_MOD_REARRANGE_IP \
  struct { spec_elem_index_t e, i, j; spec_proc_t p, np; } spec1i;

/* sp_macro SPEC_DO_TPROC_MOD_REARRANGE_IP */
#define SPEC_DO_TPROC_MOD_REARRANGE_IP(_tp_, _tpd_, _b_, _xb_, _ds_, _cs_, _n_, _ib_)  do { \
  for (spec1i.e = 0, spec1i.i = 0; spec1i.i < (_n_); ++spec1i.i) { \
    spec1i.e += (_cs_)[spec1i.i]; \
    spec1i.j = (_ds_)[spec1i.i]; \
    while (spec1i.j < spec1i.e) { \
      spec1i.p = (_tp_)(spec_elem_get_buf(_b_), spec1i.j, _tpd_, spec_elem_get_buf(_ib_)); \
      spec_elem_copy_at((_ib_), 0, (_b_), spec1i.j); \
      while (spec1i.p != spec1i.i) { \
        spec1i.np = (_tp_)(spec_elem_get_buf(_b_), (_ds_)[spec1i.p], _tpd_, spec_elem_get_buf(_ib_)); \
        if (spec1i.np != spec1i.p) { \
          spec_elem_copy_at((_b_), spec1i.j, (_b_), (_ds_)[spec1i.p]); \
          spec_elem_copy_at((_ib_), 0, (_b_), spec1i.j); \
        } else spec_elem_copy_at((_ib_), 0, (_b_), (_ds_)[spec1i.p]); \
        ++(_ds_)[spec1i.p]; \
        spec1i.p = spec1i.np; \
      } \
      ++spec1i.j; \
    } \
  } } while (0)

/* sp_macro SPEC_FUNC_TPROC_MOD_REARRANGE_IP */
#define SPEC_FUNC_TPROC_MOD_REARRANGE_IP(_name_, _tp_, _s_) \
_s_ void _name_##_tproc_mod_rearrange_ip(spec_tproc_data_t tproc_data, spec_elem_t *s, spec_elem_t *x, int *displs, int *counts, spec_int_t n, spec_elem_t *mod) \
{ \
  SPEC_DECLARE_TPROC_MOD_REARRANGE_IP \
  SPEC_DO_TPROC_MOD_REARRANGE_IP(_tp_, tproc_data, s, x, displs, counts, n, mod); \
}


/* tprocs rearrange */

/* sp_macro SPEC_DECLARE_FIXED_TPROCS_REARRANGE_DB */
#define SPEC_DECLARE_FIXED_TPROCS_REARRANGE_DB(_fxdcl_, _fxp_) \
  struct { spec_elem_index_t n, i; spec_int_t j, m; _fxdcl_(fx, _fxp_) } spec2d;

/* sp_macro SPEC_DO_FIXED_TPROCS_REARRANGE_DB */
#define SPEC_DO_FIXED_TPROCS_REARRANGE_DB(_fxp_, _fxc_, _fxca_, _fxxa_, _fxd_, _tp_, _tpd_, _sb_, _db_, _ds_, _ps_)  do { \
  _fxc_(spec2d.fx, _fxp_); \
  spec2d.n = spec_elem_get_n(_sb_); \
  for (spec2d.i = 0; spec2d.i < spec2d.n; ++spec2d.i) { \
    (_tp_)(spec_elem_get_buf(_sb_), spec2d.i, (_tpd_), &spec2d.m, (_ps_)); \
    for (spec2d.j = 0; spec2d.j < spec2d.m; ++spec2d.j) { \
      _fxca_((_sb_), spec2d.i, (_db_), (_ds_)[(_ps_)[spec2d.j]], spec2d.fx, _fxp_); \
      ++(_ds_)[(_ps_)[spec2d.j]]; \
    } \
  } \
  _fxd_(spec2d.fx, _fxp_); \
  } while (0)

/* sp_macro SPEC_FUNC_FIXED_TPROCS_REARRANGE_DB */
#define SPEC_FUNC_FIXED_TPROCS_REARRANGE_DB(_name_, _fxdcl_, _fxp_, _fxc_, _fxca_, _fxxa_, _fxd_, _tp_, _s_...) \
_s_ void _name_##_tprocs_rearrange_db(spec_tproc_data_t tproc_data, spec_elem_t *s, spec_elem_t *d, int *displs, spec_proc_t *procs) \
{ \
  SPEC_DECLARE_FIXED_TPROCS_REARRANGE_DB(_fxdcl_, _fxp_) \
  SPEC_DO_FIXED_TPROCS_REARRANGE_DB(_fxp_, _fxc_, _fxca_, _fxxa_, _fxd_, _tp_, tproc_data, s, d, displs, procs); \
}

/* sp_macro SPEC_DECLARE_TPROCS_REARRANGE_DB */
#define SPEC_DECLARE_TPROCS_REARRANGE_DB \
  SPEC_DECLARE_FIXED_TPROCS_REARRANGE_DB(spec_fixed_default_declare, NOPARAM)

/* sp_macro SPEC_DO_TPROCS_REARRANGE_DB */
#define SPEC_DO_TPROCS_REARRANGE_DB(_tp_, _tpd_, _sb_, _db_, _ds_, _ps_) \
  SPEC_DO_FIXED_TPROCS_REARRANGE_DB(NOPARAM, spec_fixed_default_create, spec_fixed_default_copy_at, spec_fixed_default_exchange_at, spec_fixed_default_destroy, _tp_, _tpd_, _sb_, _db_, _ds_, _ps_)

/* sp_macro SPEC_FUNC_TPROCS_REARRANGE_DB */
#define SPEC_FUNC_TPROCS_REARRANGE_DB(_name_, _tp_, _s_...) \
  SPEC_FUNC_FIXED_TPROCS_REARRANGE_DB(_name_, spec_fixed_default_declare, NOPARAM, spec_fixed_default_create, spec_fixed_default_copy_at, spec_fixed_default_exchange_at, spec_fixed_default_destroy, _tp_, _s_)

/* sp_macro SPEC_DECLARE_TPROCS_REARRANGE_IP */
#define SPEC_DECLARE_TPROCS_REARRANGE_IP \
  struct { spec_elem_index_t e, j, fe, fc, le, lc; spec_int_t i, n, f, l, o; } spec2i;

/* sp_macro SPEC_DO_TPROCS_REARRANGE_IP */
#define SPEC_DO_TPROCS_REARRANGE_IP(_tp_, _tpd_, _b_, _xb_, _ds_, _cs_, _n_, _ps_)  do { \
  spec2i.f = 0; spec2i.fe = (_cs_)[0]; spec2i.fc = spec_elem_get_n(_b_); \
  while (spec2i.f + 1 < (_n_) && spec2i.fc >= spec2i.fe) { ++spec2i.f; spec2i.fe += (_cs_)[spec2i.f]; } \
  spec2i.l = 0; spec2i.le = (_cs_)[0]; spec2i.lc = spec_elem_get_n(_b_) - 1; \
  while (spec2i.lc >= spec2i.le) { ++spec2i.l; spec2i.le += (_cs_)[spec2i.l]; } \
  for (spec2i.e = 0, spec2i.i = 0; spec2i.i < (_n_); ++spec2i.i) { \
    spec2i.e += (_cs_)[spec2i.i]; \
    spec2i.j = (_ds_)[spec2i.i]; \
    while (spec2i.j < spec2i.e) { \
      (_tp_)(spec_elem_get_buf(_b_), spec2i.j, (_tpd_), &spec2i.n, (_ps_)); \
      spec2i.o = -1; \
      while (spec2i.n > 0) { \
        --spec2i.n; \
        if ((_ps_)[spec2i.n] == spec2i.i && spec2i.o < 0) spec2i.o = spec2i.n; \
        else if ((_ds_)[(_ps_)[spec2i.n]] < spec2i.fc) { \
          spec2i.l = spec2i.f; spec2i.le = spec2i.fe; spec2i.lc = spec2i.fc; \
          if (spec2i.fc < spec2i.fe) { \
            spec_elem_copy_at((_b_), (_ds_)[(_ps_)[spec2i.n]], (_b_), spec2i.fc); \
            ++spec2i.fc; \
          } else spec_elem_copy_at((_b_), (_ds_)[(_ps_)[spec2i.n]], (_xb_), 0); \
        } else if ((_ds_)[(_ps_)[spec2i.n]] == spec2i.fc) ++spec2i.fc; \
        if (spec2i.j != (_ds_)[(_ps_)[spec2i.n]]) spec_elem_copy_at((_b_), spec2i.j, (_b_), (_ds_)[(_ps_)[spec2i.n]]); \
        ++(_ds_)[(_ps_)[spec2i.n]]; \
        while (spec2i.f + 1 < (_n_) && spec2i.fc >= spec2i.fe) { ++spec2i.f; spec2i.fe += (_cs_)[spec2i.f]; spec2i.fc = (_ds_)[spec2i.f]; } \
      } \
      if (spec2i.o < 0) { \
        if (spec2i.lc < spec2i.le) {  \
          spec_elem_copy_at((_b_), spec2i.lc, (_b_), spec2i.j); \
          spec2i.f = spec2i.l; spec2i.fe = spec2i.le; spec2i.fc = spec2i.lc; \
          --spec2i.lc; \
          while (spec2i.l > 0 && spec2i.lc < (_ds_)[spec2i.l]) { spec2i.le -= (_cs_)[spec2i.l]; spec2i.lc = spec2i.le - 1; --spec2i.l; } \
        } else spec_elem_copy_at((_xb_), 0, (_b_), spec2i.j); \
      } \
      spec2i.j = (_ds_)[spec2i.i]; \
    } \
  } } while (0)

/* sp_macro SPEC_FUNC_TPROCS_REARRANGE_IP */
#define SPEC_FUNC_TPROCS_REARRANGE_IP(_name_, _tp_, _s_...) \
_s_ void _name_##_tprocs_rearrange_ip(spec_tproc_data_t tproc_data, spec_elem_t *s, spec_elem_t *d, int *displs, int *counts, spec_int_t n, spec_proc_t *procs) \
{ \
  SPEC_DECLARE_TPROCS_REARRANGE_IP \
  SPEC_DO_TPROCS_REARRANGE_IP(_tp_, tproc_data, s, d, displs, counts, n, procs); \
}


/* tprocs_mod rearrange */

/* sp_macro SPEC_DECLARE_FIXED_TPROCS_MOD_REARRANGE_DB */
#define SPEC_DECLARE_FIXED_TPROCS_MOD_REARRANGE_DB(_fxdcl_, _fxp_) \
  struct { spec_elem_index_t n, i; spec_int_t j, m; _fxdcl_(fx, _fxp_) } spec3d;

/* sp_macro SPEC_DO_FIXED_TPROCS_MOD_REARRANGE_DB */
#define SPEC_DO_FIXED_TPROCS_MOD_REARRANGE_DB(_fxp_, _fxc_, _fxca_, _fxxa_, _fxd_, _tp_, _tpd_, _sb_, _db_, _ds_, _ps_, _ib_)  do { \
  _fxc_(spec3d.fx, _fxp_); \
  spec3d.n = spec_elem_get_n(_sb_); \
  for (spec3d.i = 0; spec3d.i < spec3d.n; ++spec3d.i) { \
    (_tp_)(spec_elem_get_buf(_sb_), spec3d.i, (_tpd_), &spec3d.m, (_ps_), spec_elem_get_buf(_ib_)); \
    for (spec3d.j = 0; spec3d.j < spec3d.m; ++spec3d.j) { \
      _fxca_((_ib_), spec3d.j, (_db_), (_ds_)[(_ps_)[spec3d.j]], spec3d.fx, _fxp_); \
      ++(_ds_)[(_ps_)[spec3d.j]]; \
    } \
  } \
  _fxd_(spec3d.fx, _fxp_); \
  } while (0)

/* sp_macro SPEC_FUNC_FIXED_TPROCS_MOD_REARRANGE_DB */
#define SPEC_FUNC_FIXED_TPROCS_MOD_REARRANGE_DB(_name_, _fxdcl_, _fxp_, _fxc_, _fxca_, _fxxa_, _fxd_, _tp_, _s_...) \
_s_ void _name_##_tprocs_mod_rearrange_db(spec_tproc_data_t tproc_data, spec_elem_t *s, spec_elem_t *d, int *displs, spec_proc_t *procs, spec_elem_t *mod) \
{ \
  SPEC_DECLARE_FIXED_TPROCS_MOD_REARRANGE_DB(_fxdcl_, _fxp_) \
  SPEC_DO_FIXED_TPROCS_MOD_REARRANGE_DB(_fxp_, _fxc_, _fxca_, _fxxa_, _fxd_, _tp_, tproc_data, s, d, displs, procs, mod); \
}

/* sp_macro SPEC_DECLARE_TPROCS_MOD_REARRANGE_DB */
#define SPEC_DECLARE_TPROCS_MOD_REARRANGE_DB \
  SPEC_DECLARE_FIXED_TPROCS_MOD_REARRANGE_DB(spec_fixed_default_declare, NOPARAM)

/* sp_macro SPEC_DO_TPROCS_MOD_REARRANGE_DB */
#define SPEC_DO_TPROCS_MOD_REARRANGE_DB(_tp_, _tpd_, _sb_, _db_, _ds_, _ps_, _ib_) \
  SPEC_DO_FIXED_TPROCS_MOD_REARRANGE_DB(NOPARAM, spec_fixed_default_create, spec_fixed_default_copy_at, spec_fixed_default_exchange_at, spec_fixed_default_destroy, _tp_, _tpd_, _sb_, _db_, _ds_, _ps_, _ib_)

/* sp_macro SPEC_FUNC_TPROCS_MOD_REARRANGE_DB */
#define SPEC_FUNC_TPROCS_MOD_REARRANGE_DB(_name_, _tp_, _s_...) \
  SPEC_FUNC_FIXED_TPROCS_MOD_REARRANGE_DB(_name_, spec_fixed_default_declare, NOPARAM, spec_fixed_default_create, spec_fixed_default_copy_at, spec_fixed_default_exchange_at, spec_fixed_default_destroy, _tp_, _s_)

/* sp_macro SPEC_DECLARE_TPROCS_MOD_REARRANGE_IP */
#define SPEC_DECLARE_TPROCS_MOD_REARRANGE_IP \
  struct { spec_elem_index_t e, j, fe, fc, le, lc; spec_int_t i, n, f, l, o; } spec3i;

/* sp_macro SPEC_DO_TPROCS_MOD_REARRANGE_IP */
#define SPEC_DO_TPROCS_MOD_REARRANGE_IP(_tp_, _tpd_, _b_, _xb_, _ds_, _cs_, _n_, _ps_, _ib_)  do { \
  spec3i.f = 0; spec3i.fe = (_cs_)[0]; spec3i.fc = spec_elem_get_n(_b_); \
  while (spec3i.f + 1 < (_n_) && spec3i.fc >= spec3i.fe) { ++spec3i.f; spec3i.fe += (_cs_)[spec3i.f]; } \
  spec3i.l = 0; spec3i.le = (_cs_)[0]; spec3i.lc = spec_elem_get_n(_b_) - 1; \
  while (spec3i.lc >= spec3i.le) { ++spec3i.l; spec3i.le += (_cs_)[spec3i.l]; } \
  for (spec3i.e = 0, spec3i.i = 0; spec3i.i < (_n_); ++spec3i.i) { \
    spec3i.e += (_cs_)[spec3i.i]; \
    spec3i.j = (_ds_)[spec3i.i]; \
    while (spec3i.j < spec3i.e) { \
      (_tp_)(spec_elem_get_buf(_b_), spec3i.j, (_tpd_), &spec3i.n, (_ps_), spec_elem_get_buf(_ib_)); \
      spec3i.o = -1; \
      while (spec3i.n > 0) { \
        --spec3i.n; \
        if ((_ps_)[spec3i.n] == spec3i.i && spec3i.o < 0) spec3i.o = spec3i.n; \
        else if ((_ds_)[(_ps_)[spec3i.n]] < spec3i.fc) { \
          spec3i.l = spec3i.f; spec3i.le = spec3i.fe; spec3i.lc = spec3i.fc; \
          if (spec3i.fc < spec3i.fe) { \
            spec_elem_copy_at((_b_), (_ds_)[(_ps_)[spec3i.n]], (_b_), spec3i.fc); \
            ++spec3i.fc; \
          } else spec_elem_copy_at((_b_), (_ds_)[(_ps_)[spec3i.n]], (_xb_), 0); \
        } else if ((_ds_)[(_ps_)[spec3i.n]] == spec3i.fc) ++spec3i.fc; \
        spec_elem_copy_at((_ib_), spec3i.n, (_b_), (_ds_)[(_ps_)[spec3i.n]]); \
        ++(_ds_)[(_ps_)[spec3i.n]]; \
        while (spec3i.f + 1 < (_n_) && spec3i.fc >= spec3i.fe) { ++spec3i.f; spec3i.fe += (_cs_)[spec3i.f]; spec3i.fc = (_ds_)[spec3i.f]; } \
      } \
      if (spec3i.o < 0) { \
        if (spec3i.lc < spec3i.le) {  \
          spec_elem_copy_at((_b_), spec3i.lc, (_b_), spec3i.j); \
          spec3i.f = spec3i.l; spec3i.fe = spec3i.le; spec3i.fc = spec3i.lc; \
          --spec3i.lc; \
          while (spec3i.l > 0 && spec3i.lc < (_ds_)[spec3i.l]) { spec3i.le -= (_cs_)[spec3i.l]; spec3i.lc = spec3i.le - 1; --spec3i.l; } \
        } else spec_elem_copy_at((_xb_), 0, (_b_), spec3i.j); \
      } \
      spec3i.j = (_ds_)[spec3i.i]; \
    } \
  } } while (0)

/* sp_macro SPEC_FUNC_TPROCS_MOD_REARRANGE_IP */
#define SPEC_FUNC_TPROCS_MOD_REARRANGE_IP(_name_, _tp_, _s_...) \
_s_ void _name_##_tprocs_mod_rearrange_ip(spec_tproc_data_t tproc_data, spec_elem_t *s, spec_elem_t *x, int *displs, int *counts, spec_int_t n, spec_proc_t *procs, spec_elem_t *mod) \
{ \
  SPEC_DECLARE_TPROCS_MOD_REARRANGE_IP \
  SPEC_DO_TPROCS_MOD_REARRANGE_IP(_tp_, tproc_data, s, x, displs, counts, n, procs, mod); \
}


/* tproc indices */

/* sp_macro SPEC_DECLARE_TPROC_INDICES_DB */
#define SPEC_DECLARE_TPROC_INDICES_DB \
  struct { spec_elem_index_t i; spec_proc_t p; } spec0xd;

/* sp_macro SPEC_DO_TPROC_INDICES_DB */
#define SPEC_DO_TPROC_INDICES_DB(_tp_, _tpd_, _b_, _ix_, _id_)  do { \
  for (spec0xd.i = 0; spec0xd.i < spec_elem_get_n(_b_); ++spec0xd.i) { \
    spec0xd.p = (_tp_)(spec_elem_get_buf(_b_), spec0xd.i, (_tpd_)); \
    if (spec0xd.p == SPEC_PROC_NONE) continue; \
    (_ix_)[(_id_)[spec0xd.p]] = spec0xd.i; \
    ++(_id_)[spec0xd.p]; \
  } } while (0)

/* sp_macro SPEC_FUNC_TPROC_INDICES_DB */
#define SPEC_FUNC_TPROC_INDICES_DB(_name_, _tp_, _s_...) \
_s_ void _name_##_tproc_indices_db(spec_tproc_data_t tproc_data, spec_elem_t *s, int *indices, int *idispls) \
{ \
  SPEC_DECLARE_TPROC_INDICES_DB \
  SPEC_DO_TPROC_INDICES_DB(_tp_, tproc_data, s, indices, idispls); \
}


/* tproc_mod indices */

/* sp_macro SPEC_DECLARE_TPROC_MOD_INDICES_DB */
#define SPEC_DECLARE_TPROC_MOD_INDICES_DB \
  struct { spec_elem_index_t i, k; spec_proc_t p; } spec1xd;

/* sp_macro SPEC_DO_TPROC_MOD_INDICES_DB */
#define SPEC_DO_TPROC_MOD_INDICES_DB(_tp_, _tpd_, _b_, _ix_, _id_, _ib_, _d_)  do { \
  spec1xd.k = 0; \
  for (spec1xd.i = 0; spec1xd.i < spec_elem_get_n(_b_); ++spec1xd.i) { \
    spec1xd.p = (_tp_)(spec_elem_get_buf(_b_), spec1xd.i, (_tpd_), spec_elem_get_buf(_ib_)); \
    if (spec1xd.p == SPEC_PROC_NONE) continue; \
    spec_elem_copy_at((_ib_), 0, (_d_), spec1xd.k); \
    (_ix_)[(_id_)[spec1xd.p]] = spec1xd.k; \
    ++spec1xd.k; \
    ++(_id_)[spec1xd.p]; \
  } } while (0)

/* sp_macro SPEC_FUNC_TPROC_MOD_INDICES_DB */
#define SPEC_FUNC_TPROC_MOD_INDICES_DB(_name_, _tp_, _s_...) \
_s_ void _name_##_tproc_mod_indices_db(spec_tproc_data_t tproc_data, spec_elem_t *s, int *indices, int *idispls, spec_elem_t *mod, spec_elem_t *d) \
{ \
  SPEC_DECLARE_TPROC_MOD_INDICES_DB \
  SPEC_DO_TPROC_MOD_INDICES_DB(_tp_, tproc_data, s, indices, idispls, mod, d); \
}


/* tprocs indices */

/* sp_macro SPEC_DECLARE_TPROCS_INDICES_DB */
#define SPEC_DECLARE_TPROCS_INDICES_DB \
  struct { spec_elem_index_t i; spec_int_t j, n; } spec2xd;

/* sp_macro SPEC_DO_TPROCS_INDICES_DB */
#define SPEC_DO_TPROCS_INDICES_DB(_tp_, _tpd_, _b_, _ix_, _id_, _ps_)  do { \
  for (spec2xd.i = 0; spec2xd.i < spec_elem_get_n(_b_); ++spec2xd.i) { \
    (_tp_)(spec_elem_get_buf(_b_), spec2xd.i, (_tpd_), &spec2xd.n, (_ps_)); \
    for (spec2xd.j = 0; spec2xd.j < spec2xd.n; ++spec2xd.j) { \
      (_ix_)[(_id_)[(_ps_)[spec2xd.j]]] = spec2xd.i; \
      ++(_id_)[(_ps_)[spec2xd.j]]; \
    } \
  } } while (0)

/* sp_macro SPEC_FUNC_TPROCS_INDICES_DB */
#define SPEC_FUNC_TPROCS_INDICES_DB(_name_, _tp_, _s_...) \
_s_ void _name_##_tprocs_indices_db(spec_tproc_data_t tproc_data, spec_elem_t *s, int *indices, int *idispls, spec_proc_t *procs) \
{ \
  SPEC_DECLARE_TPROCS_INDICES_DB \
  SPEC_DO_TPROCS_INDICES_DB(_tp_, tproc_data, s, indices, idispls, procs); \
}


/* tprocs_mod indices */

/* sp_macro SPEC_DECLARE_TPROCS_MOD_INDICES_DB */
#define SPEC_DECLARE_TPROCS_MOD_INDICES_DB \
  struct { spec_elem_index_t i, k; spec_int_t j, n; } spec3xd;

/* sp_macro SPEC_DO_TPROCS_MOD_INDICES_DB */
#define SPEC_DO_TPROCS_MOD_INDICES_DB(_tp_, _tpd_, _b_, _ix_, _id_, _ps_, _ib_, _d_)  do { \
  spec3xd.k = 0; \
  for (spec3xd.i = 0; spec3xd.i < spec_elem_get_n(_b_); ++spec3xd.i) { \
    (_tp_)(spec_elem_get_buf(_b_), spec3xd.i, (_tpd_), &spec3xd.n, (_ps_), spec_elem_get_buf(_ib_)); \
    for (spec3xd.j = 0; spec3xd.j < spec3xd.n; ++spec3xd.j) { \
      spec_elem_copy_at((_ib_), spec3xd.j, (_d_), spec3xd.k); \
      (_ix_)[(_id_)[(_ps_)[spec3xd.j]]] = spec3xd.k; \
      ++spec3xd.k; \
      ++(_id_)[(_ps_)[spec3xd.j]]; \
    } \
  } } while (0)

/* sp_macro SPEC_FUNC_TPROCS_MOD_INDICES_DB */
#define SPEC_FUNC_TPROCS_MOD_INDICES_DB(_name_, _tp_, _s_...) \
_s_ void _name_##_tprocs_mod_indices_db(spec_tproc_data_t tproc_data, spec_elem_t *s, int *indices, int *idispls, spec_proc_t *procs, spec_elem_t *mod, spec_elem_t *d) \
{ \
  SPEC_DECLARE_TPROCS_MOD_INDICES_DB \
  SPEC_DO_TPROCS_MOD_INDICES_DB(_tp_, tproc_data, s, indices, idispls, procs, mod, d); \
}


/* tproc sendrecv */

/* sp_macro SPEC_DECLARE_FIXED_TPROC_SENDRECV_DB */
#define SPEC_DECLARE_FIXED_TPROC_SENDRECV_DB(_fxdcl_, _fxp_)
/*#define SPEC_DECLARE_FIXED_TPROC_SENDRECV_DB(_fxdcl_, _fxp_) \
  struct { _fxdcl_(fx, _fxp_) } spec0srd;*/

/* sp_macro SPEC_DO_FIXED_TPROC_SENDRECV_DB */
#define SPEC_DO_FIXED_TPROC_SENDRECV_DB(_fxp_, _fxc_, _fxca_, _fxxa_, _fxd_, _tp_, _tpd_, _sb_, _rb_, _sc_, _sd_, _rd_, _ab_, _ad_, _as_, _aq_, _aqn_, _aqs_, _r_, _p_)  do { \
  _fxc_(spec0srd.fx, _fxp_); \
  while (*(_sd_) < (_sc_)) { \
    if ((_p_) == SPEC_PROC_NONE) (_p_) = (_tp_)(spec_elem_get_buf(_sb_), *(_sd_), (_tpd_)); \
    if ((_p_) != SPEC_PROC_NONE) { \
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
    (_p_) = SPEC_PROC_NONE; \
    ++(*(_sd_)); \
 } \
 _fxd_(spec0srd.fx, _fxp_); \
 } while (0)

/* sp_macro SPEC_FUNC_FIXED_TPROC_SENDRECV_DB */
#define SPEC_FUNC_FIXED_TPROC_SENDRECV_DB(_name_, _fxdcl_, _fxp_, _fxc_, _fxca_, _fxxa_, _fxd_, _tp_, _s_...) \
_s_ spec_proc_t _name_##_tproc_sendrecv_db(spec_tproc_data_t tproc_data, spec_elem_t *sb, spec_elem_t *rb, spec_int_t scount, spec_int_t *sdispl, spec_int_t *rdispl, spec_elem_t *ax, spec_int_t *aux_displs, spec_int_t aux_size_max, spec_int_t *aux_queue, spec_int_t *aux_queue_next, spec_int_t aux_queue_size, spec_proc_t rank, spec_proc_t p) \
{ \
  SPEC_DECLARE_FIXED_TPROC_SENDRECV_DB(_fxdcl_, _fxp_) \
  SPEC_DO_FIXED_TPROC_SENDRECV_DB(_fxp_, _fxc_, _fxca_, _fxxa_, _fxd_, _tp_, tproc_data, sb, rb, scount, sdispl, rdispl, ax, aux_displs, aux_size_max, aux_queue, aux_queue_next, aux_queue_size, rank, p); \
  return p; \
}

/* sp_macro SPEC_DECLARE_TPROC_SENDRECV_DB */
#define SPEC_DECLARE_TPROC_SENDRECV_DB \
  SPEC_DECLARE_FIXED_TPROC_SENDRECV_DB(spec_fixed_default_declare, NOPARAM)

/* sp_macro SPEC_DO_TPROC_SENDRECV_DB */
#define SPEC_DO_TPROC_SENDRECV_DB(_tp_, _tpd_, _sb_, _rb_, _sc_, _sd_, _rd_, _ab_, _ad_, _as_, _aq_, _aqn_, _aqs_, _r_, _p_) \
  SPEC_DO_FIXED_TPROC_SENDRECV_DB(NOPARAM, spec_fixed_default_create, spec_fixed_default_copy_at, spec_fixed_default_exchange_at, spec_fixed_default_destroy, _tp_, _tpd_, _sb_, _rb_, _sc_, _sd_, _rd_, _ab_, _ad_, _as_, _aq_, _aqn_, _aqs_, _r_, _p_)

/* sp_macro SPEC_FUNC_TPROC_SENDRECV_DB */
#define SPEC_FUNC_TPROC_SENDRECV_DB(_name_, _tp_, _s_...) \
  SPEC_FUNC_FIXED_TPROC_SENDRECV_DB(_name_, spec_fixed_default_declare, NOPARAM, spec_fixed_default_create, spec_fixed_default_copy_at, spec_fixed_default_exchange_at, spec_fixed_default_destroy, _tp_, _s_)


/* sp_macro SPEC_DEFINE_TPROC */
#define SPEC_DEFINE_TPROC(_name_, _tp_, _s_...) \
  SPEC_FUNC_TPROC_COUNT_DB(_name_, _tp_, _s_) \
  SPEC_FUNC_TPROC_COUNT_IP(_name_, _tp_, _s_) \
  SPEC_FUNC_TPROC_REARRANGE_DB(_name_, _tp_, _s_) \
  SPEC_FUNC_TPROC_REARRANGE_IP(_name_, _tp_, _s_) \
  SPEC_FUNC_TPROC_INDICES_DB(_name_, _tp_, _s_) \
  SPEC_FUNC_TPROC_SENDRECV_DB(_name_, _tp_, _s_)

/* sp_macro SPEC_DEFINE_TPROC_MOD */
#define SPEC_DEFINE_TPROC_MOD(_name_, _tp_, _s_...) \
  SPEC_FUNC_TPROC_MOD_COUNT_DB(_name_, _tp_, _s_) \
  SPEC_FUNC_TPROC_MOD_COUNT_IP(_name_, _tp_, _s_) \
  SPEC_FUNC_TPROC_MOD_REARRANGE_DB(_name_, _tp_, _s_) \
  SPEC_FUNC_TPROC_MOD_REARRANGE_IP(_name_, _tp_, _s_) \
  SPEC_FUNC_TPROC_MOD_INDICES_DB(_name_, _tp_, _s_)

/* sp_macro SPEC_DEFINE_TPROCS */
#define SPEC_DEFINE_TPROCS(_name_, _tp_, _s_...) \
  SPEC_FUNC_TPROCS_COUNT_DB(_name_, _tp_, _s_) \
  SPEC_FUNC_TPROCS_COUNT_IP(_name_, _tp_, _s_) \
  SPEC_FUNC_TPROCS_REARRANGE_DB(_name_, _tp_, _s_) \
  SPEC_FUNC_TPROCS_REARRANGE_IP(_name_, _tp_, _s_) \
  SPEC_FUNC_TPROCS_INDICES_DB(_name_, _tp_, _s_)

/* sp_macro SPEC_DEFINE_TPROCS_MOD */
#define SPEC_DEFINE_TPROCS_MOD(_name_, _tp_, _s_...) \
  SPEC_FUNC_TPROCS_MOD_COUNT_DB(_name_, _tp_, _s_) \
  SPEC_FUNC_TPROCS_MOD_COUNT_IP(_name_, _tp_, _s_) \
  SPEC_FUNC_TPROCS_MOD_REARRANGE_DB(_name_, _tp_, _s_) \
  SPEC_FUNC_TPROCS_MOD_REARRANGE_IP(_name_, _tp_, _s_) \
  SPEC_FUNC_TPROCS_MOD_INDICES_DB(_name_, _tp_, _s_)

/* sp_macro SPEC_DEFINE_FIXED_TPROC */
#define SPEC_DEFINE_FIXED_TPROC(_name_, _fxdcl_, _fxp_, _fxc_, _fxca_, _fxxa_, _fxd_, _tp_, _s_...) \
  SPEC_FUNC_TPROC_COUNT_DB(_name_, _tp_, _s_) \
  SPEC_FUNC_TPROC_COUNT_IP(_name_, _tp_, _s_) \
  SPEC_FUNC_FIXED_TPROC_REARRANGE_DB(_name_, _fxdcl_, _fxp_, _fxc_, _fxca_, _fxxa_, _fxd_, _tp_, _s_) \
  SPEC_FUNC_TPROC_REARRANGE_IP(_name_, _tp_, _s_) \
  SPEC_FUNC_TPROC_INDICES_DB(_name_, _tp_, _s_) \
  SPEC_FUNC_FIXED_TPROC_SENDRECV_DB(_name_, _fxdcl_, _fxp_, _fxc_, _fxca_, _fxxa_, _fxd_, _tp_, _s_)

/* sp_macro SPEC_DEFINE_FIXED_TPROC_MOD */
#define SPEC_DEFINE_FIXED_TPROC_MOD(_name_, _fxdcl_, _fxp_, _fxc_, _fxca_, _fxxa_, _fxd_, _tp_, _s_...) \
  SPEC_FUNC_TPROC_MOD_COUNT_DB(_name_, _tp_, _s_) \
  SPEC_FUNC_TPROC_MOD_COUNT_IP(_name_, _tp_, _s_) \
  SPEC_FUNC_FIXED_TPROC_MOD_REARRANGE_DB(_name_, _fxdcl_, _fxp_, _fxc_, _fxca_, _fxxa_, _fxd_, _tp_, _s_) \
  SPEC_FUNC_TPROC_MOD_REARRANGE_IP(_name_, _tp_, _s_) \
  SPEC_FUNC_TPROC_MOD_INDICES_DB(_name_, _tp_, _s_)

/* sp_macro SPEC_DEFINE_FIXED_TPROCS */
#define SPEC_DEFINE_FIXED_TPROCS(_name_, _fxdcl_, _fxp_, _fxc_, _fxca_, _fxxa_, _fxd_, _tp_, _s_...) \
  SPEC_FUNC_TPROCS_COUNT_DB(_name_, _tp_, _s_) \
  SPEC_FUNC_TPROCS_COUNT_IP(_name_, _tp_, _s_) \
  SPEC_FUNC_FIXED_TPROCS_REARRANGE_DB(_name_, _fxdcl_, _fxp_, _fxc_, _fxca_, _fxxa_, _fxd_, _tp_, _s_) \
  SPEC_FUNC_TPROCS_REARRANGE_IP(_name_, _tp_, _s_) \
  SPEC_FUNC_TPROCS_INDICES_DB(_name_, _tp_, _s_)

/* sp_macro SPEC_DEFINE_FIXED_TPROCS_MOD */
#define SPEC_DEFINE_FIXED_TPROCS_MOD(_name_, _fxdcl_, _fxp_, _fxc_, _fxca_, _fxxa_, _fxd_, _tp_, _s_...) \
  SPEC_FUNC_TPROCS_MOD_COUNT_DB(_name_, _tp_, _s_) \
  SPEC_FUNC_TPROCS_MOD_COUNT_IP(_name_, _tp_, _s_) \
  SPEC_FUNC_FIXED_TPROCS_MOD_REARRANGE_DB(_name_, _fxdcl_, _fxp_, _fxc_, _fxca_, _fxxa_, _fxd_, _tp_, _s_) \
  SPEC_FUNC_TPROCS_MOD_REARRANGE_IP(_name_, _tp_, _s_) \
  SPEC_FUNC_TPROCS_MOD_INDICES_DB(_name_, _tp_, _s_)

/* sp_type spec_tproc_f spec_tproc_count_f spec_tproc_rearrange_db_f spec_tproc_rearrange_ip_f spec_tproc_indices_db_f spec_tproc_sendrecv_db_f */
typedef spec_proc_t spec_tproc_f(spec_elem_buf_t b, spec_elem_index_t x, spec_tproc_data_t tproc_data);
typedef void spec_tproc_count_f(spec_tproc_data_t tproc_data, spec_elem_t *s, int *counts);
typedef void spec_tproc_rearrange_db_f(spec_tproc_data_t tproc_data, spec_elem_t *s, spec_elem_t *d, int *displs);
typedef void spec_tproc_rearrange_ip_f(spec_tproc_data_t tproc_data, spec_elem_t *s, spec_elem_t *x, int *displs, int *counts, spec_int_t n);
typedef void spec_tproc_indices_db_f(spec_tproc_data_t tproc_data, spec_elem_t *s, int *indices, int *idispls);
typedef spec_proc_t spec_tproc_sendrecv_db_f(spec_tproc_data_t tproc_data, spec_elem_t *sb, spec_elem_t *rb, spec_int_t scount, spec_int_t *sdispl, spec_int_t *rdispl, spec_elem_t *ax, spec_int_t *aux_displs, spec_int_t aux_size_max, spec_int_t *aux_queue, spec_int_t *aux_queue_next, spec_int_t aux_queue_size, spec_proc_t rank, spec_proc_t p);

/* sp_type spec_tproc_mod_f spec_tproc_mod_count_f spec_tproc_mod_rearrange_db_f spec_tproc_mod_rearrange_ip_f spec_tproc_mod_indices_db_f */
typedef spec_proc_t spec_tproc_mod_f(spec_elem_buf_t b, spec_elem_index_t x, spec_tproc_data_t tproc_data, spec_elem_buf_t mod);
typedef void spec_tproc_mod_count_f(spec_tproc_data_t tproc_data, spec_elem_t *s, int *counts);
typedef void spec_tproc_mod_rearrange_db_f(spec_tproc_data_t tproc_data, spec_elem_t *s, spec_elem_t *d, int *displs, spec_elem_t *mod);
typedef void spec_tproc_mod_rearrange_ip_f(spec_tproc_data_t tproc_data, spec_elem_t *s, spec_elem_t *x, int *displs, int *counts, spec_int_t n, spec_elem_t *mod);
typedef void spec_tproc_mod_indices_db_f(spec_tproc_data_t tproc_data, spec_elem_t *s, int *indices, int *idispls, spec_elem_t *mod, spec_elem_t *d);

/* sp_type spec_tprocs_f spec_tprocs_count_f spec_tprocs_rearrange_db_f spec_tprocs_rearrange_ip_f spec_tprocs_indices_db_f */
typedef void spec_tprocs_f(spec_elem_buf_t b, spec_elem_index_t x, spec_tproc_data_t tproc_data, spec_int_t *nprocs, spec_proc_t *procs);
typedef void spec_tprocs_count_f(spec_tproc_data_t tproc_data, spec_elem_t *s, int *counts, spec_proc_t *procs);
typedef void spec_tprocs_rearrange_db_f(spec_tproc_data_t tproc_data, spec_elem_t *s, spec_elem_t *d, int *displs, spec_proc_t *procs);
typedef void spec_tprocs_rearrange_ip_f(spec_tproc_data_t tproc_data, spec_elem_t *s, spec_elem_t *x, int *displs, int *counts, spec_int_t n, spec_proc_t *procs);
typedef void spec_tprocs_indices_db_f(spec_tproc_data_t tproc_data, spec_elem_t *s, int *indices, int *idispls, spec_proc_t *procs);

/* sp_type spec_tprocs_mod_f spec_tprocs_mod_count_f spec_tprocs_mod_rearrange_db_f spec_tprocs_mod_rearrange_ip_f spec_tprocs_mod_indices_db_f */
typedef void spec_tprocs_mod_f(spec_elem_buf_t b, spec_elem_index_t x, spec_tproc_data_t tproc_data, spec_int_t *nprocs, spec_proc_t *procs, spec_elem_buf_t mod);
typedef void spec_tprocs_mod_count_f(spec_tproc_data_t tproc_data, spec_elem_t *s, int *counts, spec_proc_t *procs);
typedef void spec_tprocs_mod_rearrange_db_f(spec_tproc_data_t tproc_data, spec_elem_t *s, spec_elem_t *d, int *displs, spec_proc_t *procs, spec_elem_t *mod);
typedef void spec_tprocs_mod_rearrange_ip_f(spec_tproc_data_t tproc_data, spec_elem_t *s, spec_elem_t *x, int *displs, int *counts, spec_int_t n, spec_proc_t *procs, spec_elem_t *mod);
typedef void spec_tprocs_mod_indices_db_f(spec_tproc_data_t tproc_data, spec_elem_t *s, int *indices, int *idispls, spec_proc_t *procs, spec_elem_t *mod, spec_elem_t *d);

/* sp_type spec_tproc_reset_f */
typedef void spec_tproc_reset_f(spec_tproc_data_t tproc_data);


/* sp_type _spec_tproc_ext_t spec_tproc_ext_t */
typedef struct _spec_tproc_ext_t
{
  spec_tproc_count_f *count_db, *count_ip;
  spec_tproc_rearrange_db_f *rearrange_db;
  spec_tproc_rearrange_ip_f *rearrange_ip;
  spec_tproc_indices_db_f *indices_db;
  spec_tproc_sendrecv_db_f *sendrecv_db;

} spec_tproc_ext_t;

/* sp_type _spec_tproc_mod_ext_tproc_t spec_tproc_mod_ext_t */
typedef struct _spec_tproc_mod_ext_tproc_t
{
  spec_tproc_mod_count_f *count_db, *count_ip;
  spec_tproc_mod_rearrange_db_f *rearrange_db;
  spec_tproc_mod_rearrange_ip_f *rearrange_ip;
  spec_tproc_mod_indices_db_f *indices_db;

} spec_tproc_mod_ext_t;

/* sp_type _spec_tprocs_ext_t spec_tprocs_ext_t */
typedef struct _spec_tprocs_ext_t
{
  spec_tprocs_count_f *count_db, *count_ip;
  spec_tprocs_rearrange_db_f *rearrange_db;
  spec_tprocs_rearrange_ip_f *rearrange_ip;
  spec_tprocs_indices_db_f *indices_db;

} spec_tprocs_ext_t;

/* sp_type _spec_tprocs_mod_ext_t spec_tprocs_mod_ext_t */
typedef struct _spec_tprocs_mod_ext_t
{
  spec_tprocs_mod_count_f *count_db, *count_ip;
  spec_tprocs_mod_rearrange_db_f *rearrange_db;
  spec_tprocs_mod_rearrange_ip_f *rearrange_ip;
  spec_tprocs_mod_indices_db_f *indices_db;

} spec_tprocs_mod_ext_t;

/* sp_macro SPEC_EXT_PARAM_TPROC SPEC_EXT_PARAM_TPROC_NULL SPEC_EXT_PARAM_TPROC_MOD SPEC_EXT_PARAM_TPROC_MOD_NULL SPEC_EXT_PARAM_TPROCS SPEC_EXT_PARAM_TPROCS_NULL SPEC_EXT_PARAM_TPROCS_MOD SPEC_EXT_PARAM_TPROCS_MOD_NULL */
#define SPEC_EXT_PARAM_TPROC(_name_)       { _name_##_tproc_count_db, _name_##_tproc_count_ip, _name_##_tproc_rearrange_db, _name_##_tproc_rearrange_ip, _name_##_tproc_indices_db, _name_##_tproc_sendrecv_db }
#define SPEC_EXT_PARAM_TPROC_NULL          { NULL, NULL, NULL, NULL, NULL, NULL }
#define SPEC_EXT_PARAM_TPROC_MOD(_name_)   { _name_##_tproc_mod_count_db, _name_##_tproc_mod_count_ip, _name_##_tproc_mod_rearrange_db, _name_##_tproc_mod_rearrange_ip, _name_##_tproc_mod_indices_db }
#define SPEC_EXT_PARAM_TPROC_MOD_NULL      { NULL, NULL, NULL, NULL, NULL }
#define SPEC_EXT_PARAM_TPROCS(_name_)      { _name_##_tprocs_count_db, _name_##_tprocs_count_ip, _name_##_tprocs_rearrange_db, _name_##_tprocs_rearrange_ip, _name_##_tprocs_indices_db }
#define SPEC_EXT_PARAM_TPROCS_NULL         { NULL, NULL, NULL, NULL, NULL }
#define SPEC_EXT_PARAM_TPROCS_MOD(_name_)  { _name_##_tprocs_mod_count_db, _name_##_tprocs_mod_count_ip, _name_##_tprocs_mod_rearrange_db, _name_##_tprocs_mod_rearrange_ip, _name_##_tprocs_mod_indices_db }
#define SPEC_EXT_PARAM_TPROCS_MOD_NULL     { NULL, NULL, NULL, NULL, NULL }


/* enable tloc features */
#ifdef SPEC_TLOC

/* sp_macro SPEC_TLOC SPEC_LOC_NONE */


/* tloc rearrange */

/* sp_macro SPEC_DECLARE_TLOC_REARRANGE_DB */
#define SPEC_DECLARE_TLOC_REARRANGE_DB \
  struct { spec_int_t i, p; } spec0d;

/* sp_macro SPEC_DO_TLOC_REARRANGE_DB */
#define SPEC_DO_TLOC_REARRANGE_DB(_tl_, _tld_, _sb_, _db_)  do { \
  for (spec0d.i = 0; spec0d.i < spec_elem_get_n(_sb_); ++spec0d.i) { \
    spec0d.p = (_tl_)(spec_elem_get_buf(_sb_), spec0d.i, _tld_); \
    if (spec0d.p == SPEC_LOC_NONE) continue; \
    spec_elem_copy_at((_sb_), spec0d.i, (_db_), spec0d.p); \
  } } while (0)

/* sp_macro SPEC_FUNC_TLOC_REARRANGE_DB */
#define SPEC_FUNC_TLOC_REARRANGE_DB(_name_, _tl_, _s_...) \
_s_ void _name_##_tloc_rearrange_db(spec_tloc_data_t tloc_data, spec_elem_t *s, spec_elem_t *d) \
{ \
  SPEC_DECLARE_TLOC_REARRANGE_DB \
  SPEC_DO_TLOC_REARRANGE_DB(_tl_, tloc_data, s, d); \
}

/* sp_macro SPEC_DECLARE_TLOC_REARRANGE_IP */
#define SPEC_DECLARE_TLOC_REARRANGE_IP \
  struct { spec_int_t i, p, np; } spec0i;

/* sp_macro SPEC_DO_TLOC_REARRANGE_IP */
#define SPEC_DO_TLOC_REARRANGE_IP(_tl_, _tld_, _b_, _xb_)  do { \
  for (spec0i.i = 0; spec0i.i < spec_elem_get_n(_b_); ++spec0i.i) { \
    spec0i.p = (_tl_)(spec_elem_get_buf(_b_), spec0i.i, _tld_); \
    if (spec0i.p == SPEC_LOC_NONE) continue; \
    while (spec0i.i != spec0i.p) { \
      spec0i.np = (_tl_)(spec_elem_get_buf(_b_), spec0i.p, _tld_); \
      if (spec0i.np == SPEC_LOC_NONE) { spec_elem_copy_at((_b_), spec0i.i, (_b_), spec0i.p); break; } \
      spec_elem_exchange_at((_b_), spec0i.i, (_b_), spec0i.p, (_xb_)); \
      spec0i.p = spec0i.np; \
    } \
  } } while (0)

/* sp_macro SPEC_FUNC_TLOC_REARRANGE_IP */
#define SPEC_FUNC_TLOC_REARRANGE_IP(_name_, _tl_, _s_) \
_s_ void _name_##_tloc_rearrange_ip(spec_elem_t *s, spec_elem_t *x, spec_tloc_data_t tloc_data) \
{ \
  SPEC_DECLARE_TLOC_REARRANGE_IP \
  SPEC_DO_TLOC_REARRANGE_IP(_tl_, tloc_data, s, x); \
}


/* tloc_mod_mod rearrange */

/* sp_macro SPEC_DECLARE_TLOC_MOD_REARRANGE_DB */
#define SPEC_DECLARE_TLOC_MOD_REARRANGE_DB \
  struct { spec_int_t i, p; } spec1d;

/* sp_macro SPEC_DO_TLOC_MOD_REARRANGE_DB */
#define SPEC_DO_TLOC_MOD_REARRANGE_DB(_tl_, _tld_, _sb_, _db_, _ib_)  do { \
  if (_ib_) { \
    for (spec1d.i = 0; spec1d.i < spec_elem_get_n(_sb_); ++spec1d.i) { \
      spec1d.p = (_tl_)(spec_elem_get_buf(_sb_), spec1d.i, _tld_, spec_elem_get_buf(_ib_)); \
      if (spec1d.p == SPEC_LOC_NONE) continue; \
      spec_elem_copy_at((_ib_), 0, (_db_), spec1d.p); \
    } \
  } else { \
    for (spec1d.i = 0; spec1d.i < spec_elem_get_n(_sb_); ++spec1d.i) { \
      spec1d.p = (_tl_)(spec_elem_get_buf(_sb_), spec1d.i, _tld_, NULL); \
      if (spec1d.p == SPEC_LOC_NONE) continue; \
      spec_elem_copy_at((_sb_), spec1d.i, (_db_), spec1d.p); \
    } \
  } } while (0) 

/* sp_macro SPEC_FUNC_TLOC_MOD_REARRANGE_DB */
#define SPEC_FUNC_TLOC_MOD_REARRANGE_DB(_name_, _tl_, _s_...) \
_s_ void _name_##_tloc_mod_rearrange_db(spec_elem_t *s, spec_elem_t *d, spec_tloc_data_t tloc_data, spec_elem_t *mod) \
{ \
  SPEC_DECLARE_TLOC_MOD_REARRANGE_DB \
  SPEC_DO_TLOC_MOD_REARRANGE_DB(_tl_, tloc_data, s, d, mod); \
}

/* sp_macro SPEC_DECLARE_TLOC_MOD_REARRANGE_IP */
#define SPEC_DECLARE_TLOC_MOD_REARRANGE_IP \
  struct { spec_int_t i, p, np; } spec1i;

/* sp_macro SPEC_DO_TLOC_MOD_REARRANGE_IP */
#define SPEC_DO_TLOC_MOD_REARRANGE_IP(_tl_, _tld_, _b_, _xb_, _ib_)  do { \
  if (_ib_) { \
    for (spec1i.i = 0; spec1i.i < spec_elem_get_n(_b_); ++spec1i.i) { \
      spec1i.p = (_tl_)(spec_elem_get_buf(_b_), spec1i.i, _tld_, spec_elem_get_buf(_ib_)); \
      if (spec1i.p == SPEC_LOC_NONE) continue; \
      while (spec1i.i != spec1i.p) { \
        spec1i.np = (_tl_)(spec_elem_get_buf(_b_), spec1i.p, _tld_, spec_elem_get_buf(_xb_)); \
        if (spec1i.np == SPEC_LOC_NONE) break; \
        spec_elem_copy_at((_ib_), 0, (_b_), spec1i.p); \
        spec_elem_copy_at((_xb_), 0, (_ib_), 0); \
        spec1i.p = spec1i.np; \
      } \
      spec_elem_copy_at((_ib_), 0, (_b_), spec1i.i); \
    } \
  } else { \
    for (spec1i.i = 0; spec1i.i < spec_elem_get_n(_b_); ++spec1i.i) { \
      spec1i.p = (_tl_)(spec_elem_get_buf(_b_), spec1i.i, _tld_, NULL); \
      if (spec1i.p == SPEC_LOC_NONE) continue; \
      while (spec1i.i != spec1i.p) { \
        spec1i.np = (_tl_)(spec_elem_get_buf(_b_), spec1i.p, _tld_, NULL); \
        if (spec1i.np == SPEC_LOC_NONE) { spec_elem_copy_at((_b_), spec1i.i, (_b_), spec1i.p); break; } \
        spec_elem_exchange_at((_b_), spec1i.i, (_b_), spec1i.p, (_xb_)); \
        spec1i.p = spec1i.np; \
      } \
    } \
 } } while (0) 

/* sp_macro SPEC_FUNC_TLOC_MOD_REARRANGE_IP */
#define SPEC_FUNC_TLOC_MOD_REARRANGE_IP(_name_, _tl_, _s_) \
_s_ void _name_##_tloc_mod_rearrange_ip(spec_elem_t *s, spec_elem_t *x, spec_tloc_data_t tloc_data, spec_elem_t *mod) \
{ \
  SPEC_DECLARE_TLOC_MOD_REARRANGE_IP \
  SPEC_DO_TLOC_MOD_REARRANGE_IP(_tl_, tloc_data, s, x, mod); \
}

/* sp_macro SPEC_DEFINE_TLOC */
#define SPEC_DEFINE_TLOC(_name_, _tl_, _s_...) \
  SPEC_FUNC_TLOC_REARRANGE_DB(_name_, _tl_, _s_) \
  SPEC_FUNC_TLOC_REARRANGE_IP(_name_, _tl_, _s_)

/* sp_macro SPEC_DEFINE_TLOC_MOD */
#define SPEC_DEFINE_TLOC_MOD(_name_, _tl_, _s_...) \
  SPEC_FUNC_TLOC_MOD_REARRANGE_DB(_name_, _tl_, _s_) \
  SPEC_FUNC_TLOC_MOD_REARRANGE_IP(_name_, _tl_, _s_)

/* sp_macro SPEC_EXT_PARAM_TLOC SPEC_EXT_PARAM_TLOC_NULL SPEC_EXT_PARAM_TLOC_MOD SPEC_EXT_PARAM_TLOC_MOD_NULL */
#define SPEC_EXT_PARAM_TLOC(_name_)      _name_##_tloc_rearrange_db, _name_##_tloc_rearrange_ip
#define SPEC_EXT_PARAM_TLOC_NULL         NULL, NULL
#define SPEC_EXT_PARAM_TLOC_MOD(_name_)  _name_##_tloc_mod_rearrange_db, _name_##_tloc_mod_rearrange_ip
#define SPEC_EXT_PARAM_TLOC_MOD_NULL     NULL, NULL


/* sp_type spec_tloc_f spec_tloc_rearrange_db_f spec_tloc_rearrange_ip_f */
typedef spec_elem_index_t spec_tloc_f(spec_elem_buf_t b, spec_elem_index_t x, spec_tloc_data_t tloc_data);
typedef void spec_tloc_rearrange_db_f(spec_elem_t *s, spec_elem_t *d, spec_tloc_data_t tloc_data);
typedef void spec_tloc_rearrange_ip_f(spec_elem_t *s, spec_elem_t *x, spec_tloc_data_t tloc_data);

/* sp_type spec_tloc_mod_f spec_tloc_mod_rearrange_db_f spec_tloc_mod_rearrange_ip_f */
typedef spec_elem_index_t spec_tloc_mod_f(spec_elem_buf_t b, spec_elem_index_t x, spec_tloc_data_t tloc_data, spec_elem_buf_t mod);
typedef void spec_tloc_mod_rearrange_db_f(spec_elem_t *s, spec_elem_t *d, spec_tloc_data_t tloc_data, spec_elem_t *mod);
typedef void spec_tloc_mod_rearrange_ip_f(spec_elem_t *s, spec_elem_t *x, spec_tloc_data_t tloc_data, spec_elem_t *mod);


#endif /* SPEC_TLOC */


#endif /* __SPEC_PUBLIC_H__ */
