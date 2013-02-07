/*
 *  Copyright (C) 2011, 2012, 2013 Michael Hofmann
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


/* sp_macro SPEC_PROC_NULL */


/* tproc count */

/* sp_macro SPEC_DECLARE_TPROC_COUNT_DB */
#define SPEC_DECLARE_TPROC_COUNT_DB \
  struct { spint_t i, p; } spec0cd;

/* sp_macro SPEC_DO_TPROC_COUNT_DB */
#define SPEC_DO_TPROC_COUNT_DB(_tp_, _tpd_, _b_, _cs_)  do { \
  for (spec0cd.i = 0; spec0cd.i < spec_elem_get_n(_b_); ++spec0cd.i) { \
    spec0cd.p = (_tp_)(spec_elem_get_buf(_b_), spec0cd.i, _tpd_); \
    if (spec0cd.p == SPEC_PROC_NULL) continue; \
    ++(_cs_)[spec0cd.p]; \
  } } while (0)

/* sp_macro SPEC_FUNC_TPROC_COUNT_DB */
#define SPEC_FUNC_TPROC_COUNT_DB(_name_, _tp_, _s_...) \
_s_ void _name_##_tproc_count_db(spec_elem_t *s, spec_tproc_data_t tproc_data, int *counts) \
{ \
  SPEC_DECLARE_TPROC_COUNT_DB \
  SPEC_DO_TPROC_COUNT_DB(_tp_, tproc_data, s, counts); \
}

/* sp_macro SPEC_DECLARE_TPROC_COUNT_IP */
#define SPEC_DECLARE_TPROC_COUNT_IP \
  struct { spint_t i, p, t; } spec0ci;

/* sp_macro SPEC_DO_TPROC_COUNT_IP */
#define SPEC_DO_TPROC_COUNT_IP(_tp_, _tpd_, _b_, _cs_)  do { \
  spec0ci.t = 0; \
  for (spec0ci.i = 0; spec0ci.i < spec_elem_get_n(_b_); ++spec0ci.i) { \
    spec0ci.p = (_tp_)(spec_elem_get_buf(_b_), spec0ci.i, _tpd_); \
    if (spec0ci.p == SPEC_PROC_NULL) continue; \
    ++(_cs_)[spec0ci.p]; \
    if (spec0ci.t < spec0ci.i) spec_elem_copy_at((_b_), spec0ci.i, (_b_), spec0ci.t); \
    ++spec0ci.t; \
  } \
  spec_elem_set_n(_b_, spec0ci.t); \
} while (0)

/* sp_macro SPEC_FUNC_TPROC_COUNT_IP */
#define SPEC_FUNC_TPROC_COUNT_IP(_name_, _tp_, _s_...) \
_s_ void _name_##_tproc_count_ip(spec_elem_t *s, spec_tproc_data_t tproc_data, int *counts) \
{ \
  SPEC_DECLARE_TPROC_COUNT_IP \
  SPEC_DO_TPROC_COUNT_IP(_tp_, tproc_data, s, counts); \
}


/* tproc_mod count */

/* sp_macro SPEC_DECLARE_TPROC_MOD_COUNT_DB */
#define SPEC_DECLARE_TPROC_MOD_COUNT_DB \
  struct { spint_t i, p; } spec1cd;

/* sp_macro SPEC_DO_TPROC_MOD_COUNT_DB */
#define SPEC_DO_TPROC_MOD_COUNT_DB(_tp_, _tpd_, _b_, _cs_)  do { \
  for (spec1cd.i = 0; spec1cd.i < spec_elem_get_n(_b_); ++spec1cd.i) { \
    spec1cd.p = (_tp_)(spec_elem_get_buf(_b_), spec1cd.i, _tpd_, NULL); \
    if (spec1cd.p == SPEC_PROC_NULL) continue; \
    ++(_cs_)[spec1cd.p]; \
  } } while (0)

/* sp_macro SPEC_FUNC_TPROC_MOD_COUNT_DB */
#define SPEC_FUNC_TPROC_MOD_COUNT_DB(_name_, _tp_, _s_...) \
_s_ void _name_##_tproc_mod_count_db(spec_elem_t *s, spec_tproc_data_t tproc_data, int *counts) \
{ \
  SPEC_DECLARE_TPROC_MOD_COUNT_DB \
  SPEC_DO_TPROC_MOD_COUNT_DB(_tp_, tproc_data, s, counts); \
}

/* sp_macro SPEC_DECLARE_TPROC_MOD_COUNT_IP */
#define SPEC_DECLARE_TPROC_MOD_COUNT_IP \
  struct { spint_t i, p, t; } spec1ci;

/* sp_macro SPEC_DO_TPROC_MOD_COUNT_IP */
#define SPEC_DO_TPROC_MOD_COUNT_IP(_tp_, _tpd_, _b_, _cs_)  do { \
  spec1ci.t = 0; \
  for (spec1ci.i = 0; spec1ci.i < spec_elem_get_n(_b_); ++spec1ci.i) { \
    spec1ci.p = (_tp_)(spec_elem_get_buf(_b_), spec1ci.i, _tpd_, NULL); \
    if (spec1ci.p == SPEC_PROC_NULL) continue; \
    ++(_cs_)[spec1ci.p]; \
    if (spec1ci.t < spec1ci.i) spec_elem_copy_at((_b_), spec1ci.i, (_b_), spec1ci.t); \
    ++spec1ci.t; \
  } \
  spec_elem_set_n(_b_, spec1ci.t); \
} while (0)

/* sp_macro SPEC_FUNC_TPROC_MOD_COUNT_IP */
#define SPEC_FUNC_TPROC_MOD_COUNT_IP(_name_, _tp_, _s_...) \
_s_ void _name_##_tproc_mod_count_ip(spec_elem_t *s, spec_tproc_data_t tproc_data, int *counts) \
{ \
  SPEC_DECLARE_TPROC_MOD_COUNT_IP \
  SPEC_DO_TPROC_MOD_COUNT_IP(_tp_, tproc_data, s, counts); \
}


/* tprocs count */

/* sp_macro SPEC_DECLARE_TPROCS_COUNT_DB */
#define SPEC_DECLARE_TPROCS_COUNT_DB \
  struct { spint_t i, j, n; } spec2cd;

/* sp_macro SPEC_DO_TPROCS_COUNT_DB */
#define SPEC_DO_TPROCS_COUNT_DB(_tp_, _tpd_, _b_, _cs_, _ps_)  do { \
  for (spec2cd.i = 0; spec2cd.i < spec_elem_get_n(_b_); ++spec2cd.i) { \
    spec2cd.n = (_tp_)(spec_elem_get_buf(_b_), spec2cd.i, (_tpd_), (_ps_)); \
    for (spec2cd.j = 0; spec2cd.j < spec2cd.n; ++spec2cd.j) ++(_cs_)[(_ps_)[spec2cd.j]]; \
  } } while (0)

/* sp_macro SPEC_FUNC_TPROCS_COUNT_DB */
#define SPEC_FUNC_TPROCS_COUNT_DB(_name_, _tp_, _s_...) \
_s_ void _name_##_tprocs_count_db(spec_elem_t *s, spec_tproc_data_t tproc_data, int *counts, int *procs) \
{ \
  SPEC_DECLARE_TPROCS_COUNT_DB \
  SPEC_DO_TPROCS_COUNT_DB(_tp_, tproc_data, s, counts, procs); \
}

/* sp_macro SPEC_DECLARE_TPROCS_COUNT_IP */
#define SPEC_DECLARE_TPROCS_COUNT_IP \
  struct { spint_t i, j, n, t; } spec2ci;

/* sp_macro SPEC_DO_TPROCS_COUNT_IP */
#define SPEC_DO_TPROCS_COUNT_IP(_tp_, _tpd_, _b_, _cs_, _ps_)  do { \
  spec2ci.t = 0; \
  for (spec2ci.i = 0; spec2ci.i < spec_elem_get_n(_b_); ++spec2ci.i) { \
    spec2ci.n = (_tp_)(spec_elem_get_buf(_b_), spec2ci.i, (_tpd_), (_ps_)); \
    if (spec2ci.n <= 0) continue; \
    for (spec2ci.j = 0; spec2ci.j < spec2ci.n; ++spec2ci.j) ++(_cs_)[(_ps_)[spec2ci.j]]; \
    if (spec2ci.t < spec2ci.i) spec_elem_copy_at((_b_), spec2ci.i, (_b_), spec2ci.t); \
    ++spec2ci.t; \
  } \
  spec_elem_set_n(_b_, spec2ci.t); \
} while (0)

/* sp_macro SPEC_FUNC_TPROCS_COUNT_IP */
#define SPEC_FUNC_TPROCS_COUNT_IP(_name_, _tp_, _s_...) \
_s_ void _name_##_tprocs_count_ip(spec_elem_t *s, spec_tproc_data_t tproc_data, int *counts, int *procs) \
{ \
  SPEC_DECLARE_TPROCS_COUNT_IP \
  SPEC_DO_TPROCS_COUNT_IP(_tp_, tproc_data, s, counts, procs); \
}


/* tprocs_mod count */

/* sp_macro SPEC_DECLARE_TPROCS_MOD_COUNT_DB */
#define SPEC_DECLARE_TPROCS_MOD_COUNT_DB \
  struct { spint_t i, j, n; } spec3cd;

/* sp_macro SPEC_DO_TPROCS_MOD_COUNT_DB */
#define SPEC_DO_TPROCS_MOD_COUNT_DB(_tp_, _tpd_, _b_, _cs_, _ps_)  do { \
  for (spec3cd.i = 0; spec3cd.i < spec_elem_get_n(_b_); ++spec3cd.i) \
  { \
    spec3cd.n = (_tp_)(spec_elem_get_buf(_b_), spec3cd.i, (_tpd_), (_ps_), NULL); \
    for (spec3cd.j = 0; spec3cd.j < spec3cd.n; ++spec3cd.j) ++(_cs_)[(_ps_)[spec3cd.j]]; \
  } } while (0)

/* sp_macro SPEC_FUNC_TPROCS_MOD_COUNT_DB */
#define SPEC_FUNC_TPROCS_MOD_COUNT_DB(_name_, _tp_, _s_...) \
_s_ void _name_##_tprocs_mod_count_db(spec_elem_t *s, spec_tproc_data_t tproc_data, int *counts, int *procs) \
{ \
  SPEC_DECLARE_TPROCS_MOD_COUNT_DB \
  SPEC_DO_TPROCS_MOD_COUNT_DB(_tp_, tproc_data, s, counts, procs); \
}

/* sp_macro SPEC_DECLARE_TPROCS_MOD_COUNT_IP */
#define SPEC_DECLARE_TPROCS_MOD_COUNT_IP \
  struct { spint_t i, j, n, t; } spec3ci;

/* sp_macro SPEC_DO_TPROCS_MOD_COUNT_IP */
#define SPEC_DO_TPROCS_MOD_COUNT_IP(_tp_, _tpd_, _b_, _cs_, _ps_)  do { \
  spec3ci.t = 0; \
  for (spec3ci.i = 0; spec3ci.i < spec_elem_get_n(_b_); ++spec3ci.i) { \
    spec3ci.n = (_tp_)(spec_elem_get_buf(_b_), spec3ci.i, (_tpd_), (_ps_), NULL); \
    if (spec3ci.n <= 0) continue; \
    for (spec3ci.j = 0; spec3ci.j < spec3ci.n; ++spec3ci.j) ++(_cs_)[(_ps_)[spec3ci.j]]; \
    if (spec3ci.t < spec3ci.i) spec_elem_copy_at((_b_), spec3ci.i, (_b_), spec3ci.t); \
    ++spec3ci.t; \
  } \
  spec_elem_set_n(_b_, spec3ci.t); \
} while (0)

/* sp_macro SPEC_FUNC_TPROCS_MOD_COUNT_IP */
#define SPEC_FUNC_TPROCS_MOD_COUNT_IP(_name_, _tp_, _s_...) \
_s_ void _name_##_tprocs_mod_count_ip(spec_elem_t *s, spec_tproc_data_t tproc_data, int *counts, int *procs) \
{ \
  SPEC_DECLARE_TPROCS_MOD_COUNT_IP \
  SPEC_DO_TPROCS_MOD_COUNT_IP(_tp_, tproc_data, s, counts, procs); \
}


/* tproc rearrange */

/* sp_macro SPEC_DECLARE_TPROC_REARRANGE_DB */
#define SPEC_DECLARE_TPROC_REARRANGE_DB \
  struct { spint_t i, p; } spec0d;

/* sp_macro SPEC_DO_TPROC_REARRANGE_DB */
#define SPEC_DO_TPROC_REARRANGE_DB(_tp_, _tpd_, _sb_, _db_, _ds_)  do { \
  for (spec0d.i = 0; spec0d.i < spec_elem_get_n(_sb_); ++spec0d.i) { \
    spec0d.p = (_tp_)(spec_elem_get_buf(_sb_), spec0d.i, _tpd_); \
    if (spec0d.p == SPEC_PROC_NULL) continue; \
    spec_elem_copy_at((_sb_), spec0d.i, (_db_), (_ds_)[spec0d.p]); \
    ++(_ds_)[spec0d.p]; \
  } } while (0)

/* sp_macro SPEC_FUNC_TPROC_REARRANGE_DB */
#define SPEC_FUNC_TPROC_REARRANGE_DB(_name_, _tp_, _s_...) \
_s_ void _name_##_tproc_rearrange_db(spec_elem_t *s, spec_elem_t *d, spec_tproc_data_t tproc_data, int *displs) \
{ \
  SPEC_DECLARE_TPROC_REARRANGE_DB \
  SPEC_DO_TPROC_REARRANGE_DB(_tp_, tproc_data, s, d, displs); \
}

/* sp_macro SPEC_DECLARE_TPROC_REARRANGE_IP */
#define SPEC_DECLARE_TPROC_REARRANGE_IP \
  struct { spint_t e, i, j, p, np; } spec0i;

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
_s_ void _name_##_tproc_rearrange_ip(spec_elem_t *s, spec_elem_t *x, spec_tproc_data_t tproc_data, int *displs, int *counts, int n) \
{ \
  SPEC_DECLARE_TPROC_REARRANGE_IP \
  SPEC_DO_TPROC_REARRANGE_IP(_tp_, tproc_data, s, x, displs, counts, n); \
}


/* tproc_mod rearrange */

/* sp_macro SPEC_DECLARE_TPROC_MOD_REARRANGE_DB */
#define SPEC_DECLARE_TPROC_MOD_REARRANGE_DB \
  struct { spint_t i, p; } spec1d;

/* sp_macro SPEC_DO_TPROC_MOD_REARRANGE_DB */
#define SPEC_DO_TPROC_MOD_REARRANGE_DB(_tp_, _tpd_, _sb_, _db_, _ds_, _ib_)  do { \
  if (_ib_) { \
    for (spec1d.i = 0; spec1d.i < spec_elem_get_n(_sb_); ++spec1d.i) { \
      spec1d.p = (_tp_)(spec_elem_get_buf(_sb_), spec1d.i, _tpd_, spec_elem_get_buf(_ib_)); \
      if (spec1d.p == SPEC_PROC_NULL) continue; \
      spec_elem_copy_at((_ib_), 0, (_db_), (_ds_)[spec1d.p]); \
      ++(_ds_)[spec1d.p]; \
    } \
  } else { \
    for (spec1d.i = 0; spec1d.i < spec_elem_get_n(_sb_); ++spec1d.i) { \
      spec1d.p = (_tp_)(spec_elem_get_buf(_sb_), spec1d.i, _tpd_, NULL); \
      if (spec1d.p == SPEC_PROC_NULL) continue; \
      spec_elem_copy_at((_sb_), spec1d.i, (_db_), (_ds_)[spec1d.p]); \
      ++(_ds_)[spec1d.p]; \
    } \
  } } while (0)

/* sp_macro SPEC_FUNC_TPROC_MOD_REARRANGE_DB */
#define SPEC_FUNC_TPROC_MOD_REARRANGE_DB(_name_, _tp_, _s_...) \
_s_ void _name_##_tproc_mod_rearrange_db(spec_elem_t *s, spec_elem_t *d, spec_tproc_data_t tproc_data, int *displs, spec_elem_t *mod) \
{ \
  SPEC_DECLARE_TPROC_MOD_REARRANGE_DB \
  SPEC_DO_TPROC_MOD_REARRANGE_DB(_tp_, tproc_data, s, d, displs, mod); \
}

/* sp_macro SPEC_DECLARE_TPROC_MOD_REARRANGE_IP */
#define SPEC_DECLARE_TPROC_MOD_REARRANGE_IP \
  struct { spint_t e, i, j, p, np; } spec1i;

/* sp_macro SPEC_DO_TPROC_MOD_REARRANGE_IP */
#define SPEC_DO_TPROC_MOD_REARRANGE_IP(_tp_, _tpd_, _b_, _xb_, _ds_, _cs_, _n_, _ib_)  do { \
  if (_ib_) { \
    for (spec1i.e = 0, spec1i.i = 0; spec1i.i < (_n_); ++spec1i.i) { \
      spec1i.e += (_cs_)[spec1i.i]; \
      spec1i.j = (_ds_)[spec1i.i]; \
      while (spec1i.j < spec1i.e) { \
        spec1i.p = (_tp_)(spec_elem_get_buf(_b_), spec1i.j, _tpd_, spec_elem_get_buf(_ib_)); \
        spec_elem_copy_at((_ib_), 0, (_b_), spec1i.j); \
        while (spec1i.p != spec1i.i) { \
          spec1i.np = (_tp_)(spec_elem_get_buf(_b_), (_ds_)[spec1i.p], _tpd_, (_ib_)); \
          if (spec1i.np != spec1i.p) { \
            spec_elem_copy_at((_b_), spec1i.j, (_b_), (_ds_)[spec1i.p]); \
            spec_elem_copy_at((_ib_), 0, (_b_), spec1i.j); \
          } else spec_elem_copy_at((_ib_), 0, (_b_), (_ds_)[spec1i.p]); \
          ++(_ds_)[spec1i.p]; \
          spec1i.p = spec1i.np; \
        } \
        ++spec1i.j; \
      } \
    } \
  } else { \
    for (spec1i.e = 0, spec1i.i = 0; spec1i.i < (_n_); ++spec1i.i) { \
      spec1i.e += (_cs_)[spec1i.i]; \
      spec1i.j = (_ds_)[spec1i.i]; \
      while (spec1i.j < spec1i.e) { \
        spec1i.p = (_tp_)(spec_elem_get_buf(_b_), spec1i.j, _tpd_, NULL); \
        while (spec1i.p != spec1i.i) { \
          spec1i.np = (_tp_)(spec_elem_get_buf(_b_), (_ds_)[spec1i.p], _tpd_, NULL); \
          if (spec1i.np != spec1i.p) spec_elem_exchange_at((_b_), (_ds_)[spec1i.p], (_b_), spec1i.j, (_xb_)); \
          ++(_ds_)[spec1i.p]; \
          spec1i.p = spec1i.np; \
        } \
        ++spec1i.j; \
      } \
    } \
  } } while (0)

/* sp_macro SPEC_FUNC_TPROC_MOD_REARRANGE_IP */
#define SPEC_FUNC_TPROC_MOD_REARRANGE_IP(_name_, _tp_, _s_) \
_s_ void _name_##_tproc_mod_rearrange_ip(spec_elem_t *s, spec_elem_t *x, spec_tproc_data_t tproc_data, int *displs, int *counts, int n, spec_elem_t *mod) \
{ \
  SPEC_DECLARE_TPROC_MOD_REARRANGE_IP \
  SPEC_DO_TPROC_MOD_REARRANGE_IP(_tp_, tproc_data, s, x, displs, counts, n, mod); \
}


/* tprocs rearrange */

/* sp_macro SPEC_DECLARE_TPROCS_REARRANGE_DB */
#define SPEC_DECLARE_TPROCS_REARRANGE_DB \
  struct { spint_t i, j, n; } spec2d;

/* sp_macro SPEC_DO_TPROCS_REARRANGE_DB */
#define SPEC_DO_TPROCS_REARRANGE_DB(_tp_, _tpd_, _sb_, _db_, _ds_, _ps_)  do { \
  for (spec2d.i = 0; spec2d.i < spec_elem_get_n(_sb_); ++spec2d.i) { \
    spec2d.n = (_tp_)(spec_elem_get_buf(_sb_), spec2d.i, (_tpd_), (_ps_)); \
    for (spec2d.j = 0; spec2d.j < spec2d.n; ++spec2d.j) { \
      spec_elem_copy_at((_sb_), spec2d.i, (_db_), (_ds_)[(_ps_)[spec2d.j]]); \
      ++(_ds_)[(_ps_)[spec2d.j]]; \
    } \
  } } while (0)

/* sp_macro SPEC_FUNC_TPROCS_REARRANGE_DB */
#define SPEC_FUNC_TPROCS_REARRANGE_DB(_name_, _tp_, _s_...) \
_s_ void _name_##_tprocs_rearrange_db(spec_elem_t *s, spec_elem_t *d, spec_tproc_data_t tproc_data, int *displs, int *procs) \
{ \
  SPEC_DECLARE_TPROCS_REARRANGE_DB \
  SPEC_DO_TPROCS_REARRANGE_DB(_tp_, tproc_data, s, d, displs, procs); \
}

/* sp_macro SPEC_DECLARE_TPROCS_REARRANGE_IP */
#define SPEC_DECLARE_TPROCS_REARRANGE_IP \
  struct { spint_t e, i, j, n, f, fe, fc, l, le, lc, o; } spec2i;

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
      spec2i.n = (_tp_)(spec_elem_get_buf(_b_), spec2i.j, (_tpd_), (_ps_)); \
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
_s_ void _name_##_tprocs_rearrange_ip(spec_elem_t *s, spec_elem_t *d, spec_tproc_data_t tproc_data, int *displs, int *counts, int n, int *procs) \
{ \
  SPEC_DECLARE_TPROCS_REARRANGE_IP \
  SPEC_DO_TPROCS_REARRANGE_IP(_tp_, tproc_data, s, d, displs, counts, n, procs); \
}


/* tprocs_mod rearrange */

/* sp_macro SPEC_DECLARE_TPROCS_MOD_REARRANGE_DB */
#define SPEC_DECLARE_TPROCS_MOD_REARRANGE_DB \
  struct { spint_t i, j, n; } spec3d;

/* sp_macro SPEC_DO_TPROCS_MOD_REARRANGE_DB */
#define SPEC_DO_TPROCS_MOD_REARRANGE_DB(_tp_, _tpd_, _sb_, _db_, _ds_, _ps_, _ib_)  do { \
  if (_ib_) { \
    for (spec3d.i = 0; spec3d.i < spec_elem_get_n(_sb_); ++spec3d.i) { \
      spec3d.n = (_tp_)(spec_elem_get_buf(_sb_), spec3d.i, (_tpd_), (_ps_), spec_elem_get_buf(_ib_)); \
      for (spec3d.j = 0; spec3d.j < spec3d.n; ++spec3d.j) { \
        spec_elem_copy_at((_ib_), spec3d.j, (_db_), (_ds_)[(_ps_)[spec3d.j]]); \
        ++(_ds_)[(_ps_)[spec3d.j]]; \
      } \
    } \
  } else { \
    for (spec3d.i = 0; spec3d.i < spec_elem_get_n(_sb_); ++spec3d.i) { \
      spec3d.n = (_tp_)(spec_elem_get_buf(_sb_), spec3d.i, (_tpd_), (_ps_), NULL); \
      for (spec3d.j = 0; spec3d.j < spec3d.n; ++spec3d.j) { \
        spec_elem_copy_at((_sb_), spec3d.i, (_db_), (_ds_)[(_ps_)[spec3d.j]]); \
        ++(_ds_)[(_ps_)[spec3d.j]]; \
      } \
    } \
  } } while (0)

/* sp_macro SPEC_FUNC_TPROCS_MOD_REARRANGE_DB */
#define SPEC_FUNC_TPROCS_MOD_REARRANGE_DB(_name_, _tp_, _s_...) \
_s_ void _name_##_tprocs_mod_rearrange_db(spec_elem_t *s, spec_elem_t *d, spec_tproc_data_t tproc_data, int *displs, int *procs, spec_elem_t *mod) \
{ \
  SPEC_DECLARE_TPROCS_MOD_REARRANGE_DB \
  SPEC_DO_TPROCS_MOD_REARRANGE_DB(_tp_, tproc_data, s, d, displs, procs, mod); \
}

/* sp_macro SPEC_DECLARE_TPROCS_MOD_REARRANGE_IP */
#define SPEC_DECLARE_TPROCS_MOD_REARRANGE_IP \
  struct { spint_t e, i, j, n, o, f, fe, fc, l, le, lc; } spec3i;

/* sp_macro SPEC_DO_TPROCS_MOD_REARRANGE_IP */
#define SPEC_DO_TPROCS_MOD_REARRANGE_IP(_tp_, _tpd_, _b_, _xb_, _ds_, _cs_, _n_, _ps_, _ib_)  do { \
  if (_ib_) { \
    spec3i.f = 0; spec3i.fe = (_cs_)[0]; spec3i.fc = spec_elem_get_n(_b_); \
    while (spec3i.f + 1 < (_n_) && spec3i.fc >= spec3i.fe) { ++spec3i.f; spec3i.fe += (_cs_)[spec3i.f]; } \
    spec3i.l = 0; spec3i.le = (_cs_)[0]; spec3i.lc = spec_elem_get_n(_b_) - 1; \
    while (spec3i.lc >= spec3i.le) { ++spec3i.l; spec3i.le += (_cs_)[spec3i.l]; } \
    for (spec3i.e = 0, spec3i.i = 0; spec3i.i < (_n_); ++spec3i.i) { \
      spec3i.e += (_cs_)[spec3i.i]; \
      spec3i.j = (_ds_)[spec3i.i]; \
      while (spec3i.j < spec3i.e) { \
        spec3i.n = (_tp_)(spec_elem_get_buf(_b_), spec3i.j, (_tpd_), (_ps_), spec_elem_get_buf(_ib_)); \
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
    } \
  } else { \
    spec3i.f = 0; spec3i.fe = (_cs_)[0]; spec3i.fc = spec_elem_get_n(_b_); \
    while (spec3i.f + 1 < (_n_) && spec3i.fc >= spec3i.fe) { ++spec3i.f; spec3i.fe += (_cs_)[spec3i.f]; } \
    spec3i.l = 0; spec3i.le = (_cs_)[0]; spec3i.lc = spec_elem_get_n(_b_) - 1; \
    while (spec3i.lc >= spec3i.le) { ++spec3i.l; spec3i.le += (_cs_)[spec3i.l]; } \
    for (spec3i.e = 0, spec3i.i = 0; spec3i.i < (_n_); ++spec3i.i) { \
      spec3i.e += (_cs_)[spec3i.i]; \
      spec3i.j = (_ds_)[spec3i.i]; \
      while (spec3i.j < spec3i.e) { \
        spec3i.n = (_tp_)(spec_elem_get_buf(_b_), spec3i.j, (_tpd_), (_ps_), NULL); \
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
          if (spec3i.j != (_ds_)[(_ps_)[spec3i.n]]) spec_elem_copy_at((_b_), spec3i.j, (_b_), (_ds_)[(_ps_)[spec3i.n]]); \
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
    } \
  } } while (0)

/* sp_macro SPEC_FUNC_TPROCS_MOD_REARRANGE_IP */
#define SPEC_FUNC_TPROCS_MOD_REARRANGE_IP(_name_, _tp_, _s_...) \
_s_ void _name_##_tprocs_mod_rearrange_ip(spec_elem_t *s, spec_elem_t *x, spec_tproc_data_t tproc_data, int *displs, int *counts, int n, int *procs, spec_elem_t *mod) \
{ \
  SPEC_DECLARE_TPROCS_MOD_REARRANGE_IP \
  SPEC_DO_TPROCS_MOD_REARRANGE_IP(_tp_, tproc_data, s, x, displs, counts, n, procs, mod); \
}

/* sp_macro SPEC_DEFINE_TPROC */
#define SPEC_DEFINE_TPROC(_name_, _tp_, _s_...) \
  SPEC_FUNC_TPROC_COUNT_DB(_name_, _tp_, _s_) \
  SPEC_FUNC_TPROC_COUNT_IP(_name_, _tp_, _s_) \
  SPEC_FUNC_TPROC_REARRANGE_DB(_name_, _tp_, _s_) \
  SPEC_FUNC_TPROC_REARRANGE_IP(_name_, _tp_, _s_)

/* sp_macro SPEC_DEFINE_TPROC_MOD */
#define SPEC_DEFINE_TPROC_MOD(_name_, _tp_, _s_...) \
  SPEC_FUNC_TPROC_MOD_COUNT_DB(_name_, _tp_, _s_) \
  SPEC_FUNC_TPROC_MOD_COUNT_IP(_name_, _tp_, _s_) \
  SPEC_FUNC_TPROC_MOD_REARRANGE_DB(_name_, _tp_, _s_) \
  SPEC_FUNC_TPROC_MOD_REARRANGE_IP(_name_, _tp_, _s_)

/* sp_macro SPEC_DEFINE_TPROCS */
#define SPEC_DEFINE_TPROCS(_name_, _tp_, _s_...) \
  SPEC_FUNC_TPROCS_COUNT_DB(_name_, _tp_, _s_) \
  SPEC_FUNC_TPROCS_COUNT_IP(_name_, _tp_, _s_) \
  SPEC_FUNC_TPROCS_REARRANGE_DB(_name_, _tp_, _s_) \
  SPEC_FUNC_TPROCS_REARRANGE_IP(_name_, _tp_, _s_)

/* sp_macro SPEC_DEFINE_TPROCS_MOD */
#define SPEC_DEFINE_TPROCS_MOD(_name_, _tp_, _s_...) \
  SPEC_FUNC_TPROCS_MOD_COUNT_DB(_name_, _tp_, _s_) \
  SPEC_FUNC_TPROCS_MOD_COUNT_IP(_name_, _tp_, _s_) \
  SPEC_FUNC_TPROCS_MOD_REARRANGE_DB(_name_, _tp_, _s_) \
  SPEC_FUNC_TPROCS_MOD_REARRANGE_IP(_name_, _tp_, _s_)

/* sp_macro SPEC_EXT_PARAM_TPROC SPEC_EXT_PARAM_TPROC_NULL SPEC_EXT_PARAM_TPROC_MOD SPEC_EXT_PARAM_TPROC_MOD_NULL SPEC_EXT_PARAM_TPROCS SPEC_EXT_PARAM_TPROCS_NULL SPEC_EXT_PARAM_TPROCS_MOD SPEC_EXT_PARAM_TPROCS_MOD_NULL */
#define SPEC_EXT_PARAM_TPROC(_name_)       _name_##_tproc_count_db, _name_##_tproc_count_ip, _name_##_tproc_rearrange_db, _name_##_tproc_rearrange_ip
#define SPEC_EXT_PARAM_TPROC_NULL          NULL, NULL, NULL, NULL
#define SPEC_EXT_PARAM_TPROC_MOD(_name_)   _name_##_tproc_mod_count_db, _name_##_tproc_mod_count_ip, _name_##_tproc_mod_rearrange_db, _name_##_tproc_mod_rearrange_ip
#define SPEC_EXT_PARAM_TPROC_MOD_NULL      NULL, NULL, NULL, NULL
#define SPEC_EXT_PARAM_TPROCS(_name_)      _name_##_tprocs_count_db, _name_##_tprocs_count_ip, _name_##_tprocs_rearrange_db, _name_##_tprocs_rearrange_ip
#define SPEC_EXT_PARAM_TPROCS_NULL         NULL, NULL, NULL, NULL
#define SPEC_EXT_PARAM_TPROCS_MOD(_name_)  _name_##_tprocs_mod_count_db, _name_##_tprocs_mod_count_ip, _name_##_tprocs_mod_rearrange_db, _name_##_tprocs_mod_rearrange_ip
#define SPEC_EXT_PARAM_TPROCS_MOD_NULL     NULL, NULL, NULL, NULL


/* sp_type spec_tproc_f spec_tproc_count_f spec_tproc_rearrange_db_f spec_tproc_rearrange_ip_f */
typedef int spec_tproc_f(spec_elem_buf_t b, spec_elem_index_t x, spec_tproc_data_t tproc_data);
typedef void spec_tproc_count_f(spec_elem_t *s, spec_tproc_data_t tproc_data, int *counts);
typedef void spec_tproc_rearrange_db_f(spec_elem_t *s, spec_elem_t *d, spec_tproc_data_t tproc_data, int *displs);
typedef void spec_tproc_rearrange_ip_f(spec_elem_t *s, spec_elem_t *x, spec_tproc_data_t tproc_data, int *displs, int *counts, int n);

/* sp_type spec_tproc_mod_f spec_tproc_mod_count_f spec_tproc_mod_rearrange_db_f spec_tproc_mod_rearrange_ip_f */
typedef int spec_tproc_mod_f(spec_elem_buf_t b, spec_elem_index_t x, spec_tproc_data_t tproc_data, spec_elem_buf_t mod);
typedef void spec_tproc_mod_count_f(spec_elem_t *s, spec_tproc_data_t tproc_data, int *counts);
typedef void spec_tproc_mod_rearrange_db_f(spec_elem_t *s, spec_elem_t *d, spec_tproc_data_t tproc_data, int *displs, spec_elem_t *mod);
typedef void spec_tproc_mod_rearrange_ip_f(spec_elem_t *s, spec_elem_t *x, spec_tproc_data_t tproc_data, int *displs, int *counts, int n, spec_elem_t *mod);

/* sp_type spec_tprocs_f spec_tprocs_count_f spec_tprocs_rearrange_db_f spec_tprocs_rearrange_ip_f */
typedef int spec_tprocs_f(spec_elem_buf_t b, spec_elem_index_t x, spec_tproc_data_t tproc_data, int *procs);
typedef void spec_tprocs_count_f(spec_elem_t *s, spec_tproc_data_t tproc_data, int *procs, int *counts);
typedef void spec_tprocs_rearrange_db_f(spec_elem_t *s, spec_elem_t *d, spec_tproc_data_t tproc_data, int *displs, int *procs);
typedef void spec_tprocs_rearrange_ip_f(spec_elem_t *s, spec_elem_t *x, spec_tproc_data_t tproc_data, int *displs, int *counts, int n, int *procs);

/* sp_type spec_tprocs_mod_f spec_tprocs_mod_count_f spec_tprocs_mod_rearrange_db_f spec_tprocs_mod_rearrange_ip_f */
typedef int spec_tprocs_mod_f(spec_elem_buf_t b, spec_elem_index_t x, spec_tproc_data_t tproc_data, int *procs, spec_elem_buf_t mod);
typedef void spec_tprocs_mod_count_f(spec_elem_t *s, spec_tproc_data_t tproc_data, int *procs, int *counts);
typedef void spec_tprocs_mod_rearrange_db_f(spec_elem_t *s, spec_elem_t *d, spec_tproc_data_t tproc_data, int *displs, int *procs, spec_elem_t *mod);
typedef void spec_tprocs_mod_rearrange_ip_f(spec_elem_t *s, spec_elem_t *x, spec_tproc_data_t tproc_data, int *displs, int *counts, int n, int *procs, spec_elem_t *mod);

/* sp_type spec_tproc_reset_f */
typedef void spec_tproc_reset_f(spec_tproc_data_t tproc_data);


#endif /* __SPEC_PUBLIC_H__ */
