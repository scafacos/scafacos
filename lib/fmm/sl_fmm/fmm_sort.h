/*
 *  Copyright (C) 2011, 2012 Michael Hofmann
 *  
 *  This file is part of ScaFaCoS/FMM.
 *  
 *  ScaFaCoS/FMM is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *  
 *  ScaFaCoS/FMM is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser Public License for more details.
 *  
 *  You should have received a copy of the GNU Lesser Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *  
 */


#ifndef __FMM_SORT_H__
#define __FMM_SORT_H__


#ifndef NOT_sl_front_xqsa0
# include "sl_front_xqsa0.h"
# define front_xqsa0(_x_)  _x_
#else
# define front_xqsa0(_x_)
#endif

#ifndef NOT_sl_front_xqsaI
# include "sl_front_xqsaI.h"
# define front_xqsaI(_x_)  _x_
#else
# define front_xqsaI(_x_)
#endif

#ifndef NOT_sl_front_xqsaX
# include "sl_front_xqsaX.h"
# define front_xqsaX(_x_)  _x_
#else
# define front_xqsaX(_x_)
#endif

#ifndef NOT_sl_front_xq_a0
# include "sl_front_xq_a0.h"
# define front_xq_a0(_x_)  _x_
#else
# define front_xq_a0(_x_)
#endif

#ifndef NOT_sl_front_xq_aI
# include "sl_front_xq_aI.h"
# define front_xq_aI(_x_)  _x_
#else
# define front_xq_aI(_x_)
#endif

#ifndef NOT_sl_front_xq_aX
# include "sl_front_xq_aX.h"
# define front_xq_aX(_x_)  _x_
#else
# define front_xq_aX(_x_)
#endif

#ifndef NOT_sl_back_qxpg
# include "sl_back_qxpg.h"
# define back_qxpg(_x_)  _x_
#else
# define back_qxpg(_x_)
#endif

#ifndef NOT_sl_back_qx_g
# include "sl_back_qx_g.h"
# define back_qx_g(_x_)  _x_
#else
# define back_qx_g(_x_)
#endif

#ifndef NOT_sl_back_q_pg
# include "sl_back_q_pg.h"
# define back_q_pg(_x_)  _x_
#else
# define back_q_pg(_x_)
#endif

#ifndef NOT_sl_back_q__g
# include "sl_back_q__g.h"
# define back_q__g(_x_)  _x_
#else
# define back_q__g(_x_)
#endif


#if !defined(NOT_sl_front_xqsa0)
# define front_(_x_)      front_xqsa0_##_x_
# define front_slint_fmt  front_xqsa0_slint_fmt
#elif !defined(NOT_sl_front_xqsaI)
# define front_(_x_)      front_xqsaI_##_x_
# define front_slint_fmt  front_xqsaI_slint_fmt
#elif !defined(NOT_sl_front_xqsaX)
# define front_(_x_)      front_xqsaX_##_x_
# define front_slint_fmt  front_xqsaX_slint_fmt
#elif !defined(NOT_sl_front_xq_a0)
# define front_(_x_)      front_xq_a0_##_x_
# define front_slint_fmt  front_xq_a0_slint_fmt
#elif !defined(NOT_sl_front_xq_aI)
# define front_(_x_)      front_xq_aI_##_x_
# define front_slint_fmt  front_xq_aI_slint_fmt
#elif !defined(NOT_sl_front_xq_aX)
# define front_(_x_)      front_xq_aX_##_x_
# define front_slint_fmt  front_xq_aX_slint_fmt
#else
# define NO_SL_FRONT
#endif

#if !defined(NOT_sl_back_qxpg)
# define back_(_x_)       back_qxpg_##_x_
# define back_slint_fmt   back_qxpg_slint_fmt
#elif !defined(NOT_sl_back_qx_g)
# define back_(_x_)       back_qx_g_##_x_
# define back_slint_fmt   back_qx_g_slint_fmt
#elif !defined(NOT_sl_back_q_pg)
# define back_(_x_)       back_q_pg_##_x_
# define back_slint_fmt   back_q_pg_slint_fmt
#elif !defined(NOT_sl_back_q__g)
# define back_(_x_)       back_q__g_##_x_
# define back_slint_fmt   back_q__g_slint_fmt
#else
# define NO_SL_BACK
#endif


#ifndef PINT_T
# define PINT_T
typedef PARAM_INTEGER_C pint_t;
#endif


#ifndef NO_SL_FRONT
void fmm_sort_front_mem(void *mem0, void *mem1, pint_t *mem_sizes, pint_t *depth, pint_t *subx, pint_t *n, front_(slkey_t) *ibox, front_(sldata0_t) *xyz, front_(sldata1_t) *q, pint_t *addr_desc, void *addr, front_(sldata2_t) *scr, pint_t *type);
void fmm_sort_front(                                               pint_t *depth, pint_t *subx, pint_t *n, front_(slkey_t) *ibox, front_(sldata0_t) *xyz, front_(sldata1_t) *q, pint_t *addr_desc, void *addr, front_(sldata2_t) *scr, pint_t *type);

void fmm_sort_front_3bit_mem(void *mem0, void *mem1, pint_t *mem_sizes, pint_t *subx, pint_t *n, front_(slkey_t) *ibox, front_(sldata0_t) *xyz, front_(sldata1_t) *q, pint_t *addr_desc, void *addr, front_(sldata2_t) *scr, pint_t *type);
void fmm_sort_front_3bit(                                               pint_t *subx, pint_t *n, front_(slkey_t) *ibox, front_(sldata0_t) *xyz, front_(sldata1_t) *q, pint_t *addr_desc, void *addr, front_(sldata2_t) *scr, pint_t *type);
#endif

#ifndef NO_SL_BACK
void fmm_sort_back_mem(void *mem0, void *mem1, pint_t *mem_sizes, pint_t *n, back_(slkey_t) *addr, back_(sldata0_t) *q, back_(sldata1_t) *xyz, back_(sldata2_t) *pot, back_(sldata3_t) *grad, pint_t *type);
void fmm_sort_back(                                               pint_t *n, back_(slkey_t) *addr, back_(sldata0_t) *q, back_(sldata1_t) *xyz, back_(sldata2_t) *pot, back_(sldata3_t) *grad, pint_t *type);
#endif


/* deprecated */
#ifdef DEPRECATED

void fmm_sort_blue_mem(void *mem, pint_t *mem_size, pint_t *depth, pint_t *subx, pint_t *n, front_xq_aI_slkey_t *ibox, front_xq_aI_sldata0_t *xyz, front_xq_aX_sldata1_t *q, front_xq_aX_sldata3_t *addr);
void fmm_sort_blue(pint_t *depth, pint_t *subx, pint_t *n, front_xq_aI_slkey_t *ibox, front_xq_aI_sldata0_t *xyz, front_xq_aX_sldata1_t *q, front_xq_aX_sldata3_t *addr);

void fmm_sort_blue_3bit_mem(void *mem, pint_t *mem_size, pint_t *subx, pint_t *n, front_xq_aI_slkey_t *ibox, front_xq_aI_sldata0_t *xyz, front_xq_aX_sldata1_t *q, front_xq_aX_sldata3_t *addr);
void fmm_sort_blue_3bit(pint_t *subx, pint_t *n, front_xq_aI_slkey_t *ibox, front_xq_aI_sldata0_t *xyz, front_xq_aX_sldata1_t *q, front_xq_aX_sldata3_t *addr);

void fmm_sort_green_mem(void *mem, pint_t *mem_size, pint_t *depth, pint_t *subx, pint_t *n, front_xqsaX_slkey_t *ibox, front_xqsaX_sldata0_t *xyz, front_xqsaX_sldata1_t *q, front_xqsaX_sldata3_t *addr, front_xqsaX_sldata2_t *scr);
void fmm_sort_green(pint_t *depth, pint_t *subx, pint_t *n, front_xqsaX_slkey_t *ibox, front_xqsaX_sldata0_t *xyz, front_xqsaX_sldata1_t *q, front_xqsaX_sldata3_t *addr, front_xqsaX_sldata2_t *scr);

void fmm_sort_green_3bit_mem(void *mem, pint_t *mem_size, pint_t *subx, pint_t *n, front_xqsaX_slkey_t *ibox, front_xqsaX_sldata0_t *xyz, front_xqsaX_sldata1_t *q, front_xqsaX_sldata3_t *addr, front_xqsaX_sldata2_t *scr);
void fmm_sort_green_3bit(pint_t *subx, pint_t *n, front_xqsaX_slkey_t *ibox, front_xqsaX_sldata0_t *xyz, front_xqsaX_sldata1_t *q, front_xqsaX_sldata3_t *addr, front_xqsaX_sldata2_t *scr);

void fmm_sort_red_mem(void *mem, pint_t *mem_size, pint_t *n, back_qxpg_slkey_t *addr, back_qxpg_sldata3_t *grad, back_qxpg_sldata2_t *pot, back_qxpg_sldata0_t *q);
void fmm_sort_red(pint_t *n, back_qxpg_slkey_t *addr, back_qxpg_sldata3_t *grad, back_qxpg_sldata2_t *pot, back_qxpg_sldata0_t *q);

void fmm_sort_red_without_pot_mem(void *mem, pint_t *mem_size, pint_t *n, back_qxpg_slkey_t *addr, back_qxpg_sldata3_t *grad, back_qxpg_sldata0_t *q);
void fmm_sort_red_without_pot(pint_t *n, back_qxpg_slkey_t *addr, back_qxpg_sldata3_t *grad, back_qxpg_sldata0_t *q);

#endif


#endif /* __FMM_SORT_H__ */
