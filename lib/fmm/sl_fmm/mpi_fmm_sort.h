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


#ifndef __MPI_FMM_SORT_H__
#define __MPI_FMM_SORT_H__


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


#ifndef NOT_sl_front_xqsaIl
# include "sl_front_xqsaIl.h"
# define front_xqsaIl(_x_)  _x_
#else
# define front_xqsaIl(_x_)
#endif

#ifndef NOT_sl_front_xq_aIl
# include "sl_front_xq_aIl.h"
# define front_xq_aIl(_x_)  _x_
#else
# define front_xq_aIl(_x_)
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

#ifndef NOT_sl_back_qxpgl
# include "sl_back_qxpgl.h"
# define back_qxpgl(_x_)  _x_
#else
# define back_qxpgl(_x_)
#endif

#ifndef NOT_sl_back_qx_gl
# include "sl_back_qx_gl.h"
# define back_qx_gl(_x_)  _x_
#else
# define back_qx_gl(_x_)
#endif

#ifndef NOT_sl_back_q_pgl
# include "sl_back_q_pgl.h"
# define back_q_pgl(_x_)  _x_
#else
# define back_q_pgl(_x_)
#endif

#ifndef NOT_sl_back_q__gl
# include "sl_back_q__gl.h"
# define back_q__gl(_x_)  _x_
#else
# define back_q__gl(_x_)
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

#if defined(front_)
typedef front_(slint_t) front_slint_t;
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
#elif !defined(NOT_sl_back_qxpgl)
# define back_(_x_)       back_qxpgl_##_x_
# define back_slint_fmt   back_qxpgl_slint_fmt
#elif !defined(NOT_sl_back_qx_gl)
# define back_(_x_)       back_qx_gl_##_x_
# define back_slint_fmt   back_qx_gl_slint_fmt
#elif !defined(NOT_sl_back_q_pgl)
# define back_(_x_)       back_q_pgl_##_x_
# define back_slint_fmt   back_q_pgl_slint_fmt
#elif !defined(NOT_sl_back_q__gl)
# define back_(_x_)       back_q__gl_##_x_
# define back_slint_fmt   back_q__gl_slint_fmt
#else
# define NO_SL_BACK
#endif

#if defined(back_)
typedef back_(slint_t) back_slint_t;
#endif

#define WITH_SORT_FRONT_LOAD
#define WITH_FCOMM

typedef SL_INTEGER_C slint_t;
#define slint_fmt  SL_INTEGER_FMT

#ifndef PINT_T
# define PINT_T
typedef PARAM_INTEGER_C pint_t;
#endif


#ifndef NO_SL_FRONT

void mpi_fmm_sort_front_mem(
 void *mem0, void *mem1, pint_t *mem_sizes,
 pint_t *depth,
 pint_t *subx, pint_t *n, front_(slkey_t) *ibox, front_(sldata0_t) *xyz, front_(sldata1_t) *q, pint_t *addr_desc, void *addr, front_(sldata2_t) *scr, pint_t *type
#ifdef WITH_SORT_FRONT_LOAD
 , front_xqsaIl_sldata4_t *load, REAL_C *imba, pint_t *nmin, pint_t *nmax
#endif
#ifdef WITH_FCOMM
 , FINT8_TYPE_C *fcomm
#endif
);
void mpi_fmm_sort_front(
 pint_t *depth,
 pint_t *subx, pint_t *n, front_(slkey_t) *ibox, front_(sldata0_t) *xyz, front_(sldata1_t) *q, pint_t *addr_desc, void *addr, front_(sldata2_t) *scr, pint_t *type
#ifdef WITH_SORT_FRONT_LOAD
 , front_xqsaIl_sldata4_t *load, REAL_C *imba, pint_t *nmin, pint_t *nmax
#endif
#ifdef WITH_FCOMM
 , FINT8_TYPE_C *fcomm
#endif
);

void mpi_fmm_sort_front_3bit_mem(
 void *mem0, void *mem1, pint_t *mem_sizes,
 pint_t *subx, pint_t *n, front_(slkey_t) *ibox, front_(sldata0_t) *xyz, front_(sldata1_t) *q, pint_t *addr_desc, void *addr, front_(sldata2_t) *scr, pint_t *type
#ifdef WITH_SORT_FRONT_LOAD
 , front_xqsaIl_sldata4_t *load, REAL_C *imba, pint_t *nmin, pint_t *nmax
#endif
#ifdef WITH_FCOMM
 , FINT8_TYPE_C *fcomm
#endif
);
void mpi_fmm_sort_front_3bit(
 pint_t *subx, pint_t *n, front_(slkey_t) *ibox, front_(sldata0_t) *xyz, front_(sldata1_t) *q, pint_t *addr_desc, void *addr, front_(sldata2_t) *scr, pint_t *type
#ifdef WITH_SORT_FRONT_LOAD
 , front_xqsaIl_sldata4_t *load, REAL_C *imba, pint_t *nmin, pint_t *nmax
#endif
#ifdef WITH_FCOMM
 , FINT8_TYPE_C *fcomm
#endif
);

void mpi_fmm_sort_front_rebalance_mem(
 void *mem0, void *mem1, pint_t *mem_sizes,
 pint_t *depth,
 pint_t *subx, pint_t *n, front_(slkey_t) *ibox, front_(sldata0_t) *xyz, front_(sldata1_t) *q, pint_t *addr_desc, void *addr, front_(sldata2_t) *scr, pint_t *type
#ifdef WITH_SORT_FRONT_LOAD
 , front_xqsaIl_sldata4_t *load, REAL_C *imba, pint_t *nmin, pint_t *nmax
#endif
#ifdef WITH_FCOMM
 , FINT8_TYPE_C *fcomm
#endif
);
void mpi_fmm_sort_front_rebalance(
 pint_t *depth,
 pint_t *subx, pint_t *n, front_(slkey_t) *ibox, front_(sldata0_t) *xyz, front_(sldata1_t) *q, pint_t *addr_desc, void *addr, front_(sldata2_t) *scr, pint_t *type
#ifdef WITH_SORT_FRONT_LOAD
 , front_xqsaIl_sldata4_t *load, REAL_C *imba, pint_t *nmin, pint_t *nmax
#endif
#ifdef WITH_FCOMM
 , FINT8_TYPE_C *fcomm
#endif
);

#endif /* NO_SL_FRONT */

#ifndef NO_SL_BACK

void mpi_fmm_sort_back_mem(
 void *mem0, void *mem1, pint_t *mem_sizes,
 pint_t *ntotal, pint_t *nin, pint_t *nout, back_(slkey_t) *addr, back_(sldata0_t) *q, back_(sldata1_t) *xyz, back_(sldata2_t) *pot, back_(sldata3_t) *grad, back_(sldata4_t) *load, pint_t *type
#ifdef WITH_FCOMM
 , FINT8_TYPE_C *fcomm
#endif
);
void mpi_fmm_sort_back(
 pint_t *ntotal, pint_t *nin, pint_t *nout, back_(slkey_t) *addr, back_(sldata0_t) *q, back_(sldata1_t) *xyz, back_(sldata2_t) *pot, back_(sldata3_t) *grad, back_(sldata4_t) *load, pint_t *type
#ifdef WITH_FCOMM
 , FINT8_TYPE_C *fcomm
#endif
);

#endif /* NO_SL_BACK */


#endif /* __MPI_FMM_SORT_H__ */
