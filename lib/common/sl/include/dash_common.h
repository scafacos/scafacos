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


#ifndef __DASH_COMMON_H__
#define __DASH_COMMON_H__


#ifdef DASH_TIMING
# define DS_TIMING_CMD(_cmd_)  Z_MOP(_cmd_)
#else
# define DS_TIMING_CMD(_cmd_)  Z_NOP()
#endif


#define DASH_RADIXSORT_SORT_DECLARE(_t_) \
  struct { \
    dsint_t o, n, n0, rh, l, h, sp; \
    struct { dsint_t o, n, rh; } s[sizeof(_t_) * 8]; \
    _t_ bm; \
  } dsi;

#define DASH_RADIXSORT_SORT(_n_, _rl_, _rh_, _get_, _xchg_) \
  dsi.n = _n_; dsi.rh = _rh_; dsi.sp = 0; dsi.o = 0; \
  while (1) { \
    while (dsi.n > 1) { \
      dsi.bm = 1 << dsi.rh; \
      dsi.l = dsi.o; dsi.h = dsi.o + dsi.n - 1; \
      while (1) { \
        while (dsi.l < dsi.h) if (_get_(dsi.l) & dsi.bm) break; else ++dsi.l; \
        while (dsi.l < dsi.h) if (_get_(dsi.h) & dsi.bm) --dsi.h; else break; \
        if (dsi.l >= dsi.h) break; \
        _xchg_(dsi.l, dsi.h); \
        ++dsi.l; --dsi.h; \
      } \
      dsi.n0 = dsi.l + ((_get_(dsi.l) & dsi.bm) == 0) - dsi.o; \
      if (dsi.rh <= _rl_) break; \
      --dsi.rh; \
      dsi.s[dsi.sp].o = dsi.o; dsi.s[dsi.sp].n = dsi.n0; dsi.s[dsi.sp].rh = dsi.rh; \
      ++dsi.sp; \
      dsi.o += dsi.n0; dsi.n -= dsi.n0; \
    } \
    if (dsi.sp <= 0) break; \
    --dsi.sp; \
    dsi.o = dsi.s[dsi.sp].o; dsi.n = dsi.s[dsi.sp].n; dsi.rh = dsi.s[dsi.sp].rh; \
  }


void ds_sort_dsints(dsint_t *ints, dsint_t n, dsint_t x);


#define DASH_BINARY_SEARCH_DECLARE() \
  struct { \
    dsint_t l, h, m; \
  } dbs;

#define DASH_BINARY_SEARCH_LT(_l_, _h_, _k_, _get_, _r_) do { \
  dbs.l = _l_; dbs.h = _h_; \
  if (dbs.l <= dbs.h) { \
    while (dbs.l < dbs.h) { \
      dbs.m = (dbs.l + dbs.h) / 2; \
      if (_get_(dbs.m) < (_k_)) dbs.l = dbs.m + 1; else dbs.h = dbs.m; \
    } \
    (_r_) = dbs.l - (!(_get_(dbs.l) < (_k_))); \
  } else (_r_) = dbs.l - 1; \
} while (0)

#define DASH_BINARY_SEARCH_LE(_l_, _h_, _k_, _get_, _r_) do { \
  dbs.l = _l_; dbs.h = _h_; \
  if (dbs.l <= dbs.h) { \
    while (dbs.l < dbs.h) { \
      dbs.m = (dbs.l + dbs.h) / 2; \
      if (_get_(dbs.m) <= (_k_)) dbs.l = dbs.m + 1; else dbs.h = dbs.m; \
    } \
    (_r_) = dbs.l - (!(_get_(dbs.l) <= (_k_))); \
  } else (_r_) = dbs.l - 1; \
} while (0)

#define DASH_BINARY_SEARCH_GT(_l_, _h_, _k_, _get_, _r_) do { \
  dbs.l = _l_; dbs.h = _h_; \
  if (dbs.l <= dbs.h) { \
    while (dbs.l < dbs.h) { \
      dbs.m = (dbs.l + dbs.h) / 2; \
      if (_get_(dbs.m) > (_k_)) dbs.h = dbs.m; else dbs.l = dbs.m + 1;\
    } \
    (_r_) = dbs.l + (!(_get_(dbs.l) > (_k_))); \
  } else (_r_) = dbs.l + 1; \
} while (0)

#define DASH_BINARY_SEARCH_GE(_l_, _h_, _k_, _get_, _r_) do { \
  dbs.l = _l_; dbs.h = _h_; \
  if (dbs.l <= dbs.h) { \
    while (dbs.l < dbs.h) { \
      dbs.m = (dbs.l + dbs.h) / 2; \
      if (_get_(dbs.m) >= (_k_)) dbs.h = dbs.m; else dbs.l = dbs.m + 1;\
    } \
    (_r_) = dbs.l + (!(_get_(dbs.l) >= (_k_))); \
  } else (_r_) = dbs.l + 1; \
} while (0)


#endif /* __DASH_COMMON_H__ */
