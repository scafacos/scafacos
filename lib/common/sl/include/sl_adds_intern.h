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


#ifndef __SL_ADDS_INTERN_H__
#define __SL_ADDS_INTERN_H__


/* slint_t math functions */
#define slint_ceil(_x_)   ceill(_x_)
#define slint_fabs(_x_)   fabsl(_x_)
#define slint_floor(_x_)  floorl(_x_)
#define slint_sqrt(_x_)   sqrtl(_x_)


inline static void elem_npack_at(elements_t *_s_, slint_t _sat_, packed_elements_t *_d_, slint_t _dat_, slint_t _n_)
{
  slint_t i;

  for (i = 0; i < _n_; ++i) elem_pack_at(_s_, _sat_ + i, _d_, _dat_ + i);
}


inline static void pelem_nunpack_at(packed_elements_t *_s_, slint_t _sat_, elements_t *_d_, slint_t _dat_, slint_t _n_)
{
  slint_t i;

  for (i = 0; i < _n_; ++i) pelem_unpack_at(_s_, _sat_ + i, _d_, _dat_ + i);
}


typedef void *tune_data_t;

typedef void (*tune_f)(void *tune_data_t);

#define TUNE_RESULT_STREAM    stdout
#define TUNE_PROGRESS_STREAM  stderr

/* src/tune/tune_common.c */
double intersect_poly2(slint_t n, double *x, double *y0, double *y1, double xmin, double xmax);

/* src/tune/tune_auto.c */
void tune_auto_begin(const char *name);
void tune_auto_end(const char *name);
void tune_auto_output_result(const char *key, const char *val);
void tune_auto_output_result_slint(const char *key, slint_t val);

/* src/tune/tune_sort_radix.c */
void tune_sort_radix_threshold(tune_data_t tune_data);
void tune_sort_radix_width(tune_data_t tune_data);


#endif /* __SL_ADDS_INTERN_H__ */
