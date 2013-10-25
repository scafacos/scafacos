/*
 *  Copyright (C) 2011, 2012, 2013 Michael Hofmann
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


#ifndef __COMMON_H__
#define __COMMON_H__


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

#ifndef NOT_sl_front_xqsaIl
# include "sl_front_xqsaIl.h"
# define front_xqsaIl(_x_)  _x_
#else
# define front_xqsaIl(_x_)
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
#elif !defined(NOT_sl_front_xqsaIl)
# define front_(_x_)      front_xqsaIl_##_x_
# define front_slint_fmt  front_xqsaIl_slint_fmt
#elif !defined(NOT_sl_front_xq_a0)
# define front_(_x_)      front_xq_a0_##_x_
# define front_slint_fmt  front_xq_a0_slint_fmt
#elif !defined(NOT_sl_front_xq_aI)
# define front_(_x_)      front_xq_aI_##_x_
# define front_slint_fmt  front_xq_aI_slint_fmt
#elif !defined(NOT_sl_front_xq_aX)
# define front_(_x_)      front_xq_aX_##_x_
# define front_slint_fmt  front_xq_aX_slint_fmt
#elif !defined(NOT_sl_front_xq_aIl)
# define front_(_x_)      front_xq_aIl_##_x_
# define front_slint_fmt  front_xq_aIl_slint_fmt
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


#define Z_MOP(_mop_)  do { _mop_ } while (0)
#define Z_NOP()       Z_MOP()

#define z_max(_a_, _b_)  (((_a_)>(_b_))?(_a_):(_b_))
#define z_min(_a_, _b_)  (((_a_)<(_b_))?(_a_):(_b_))


#ifdef DO_DEBUG
# define DEBUG_CMD(_cmd_)  Z_MOP(_cmd_)
#else
# define DEBUG_CMD(_cmd_)  Z_NOP()
#endif

#ifdef DO_INFO
# define INFO_CMD(_cmd_)  Z_MOP(_cmd_)
# error
#else
# define INFO_CMD(_cmd_)  Z_NOP()
#endif

#ifdef DO_TIMING
# define TIMING_DECL(_decl_)       _decl_
# define TIMING_CMD(_cmd_)         Z_MOP(_cmd_)
#else
# define TIMING_DECL(_decl_)
# define TIMING_CMD(_cmd_)         Z_NOP()
#endif
#ifdef DO_TIMING_SYNC
# define TIMING_SYNC(_c_)          TIMING_CMD(MPI_Barrier(_c_);)
#else
# define TIMING_SYNC(_c_)          Z_NOP()
#endif
#define TIMING_START(_t_)          TIMING_CMD(((_t_) = MPI_Wtime());)
#define TIMING_STOP(_t_)           TIMING_CMD(((_t_) = MPI_Wtime() - (_t_));)
#define TIMING_STOP_ADD(_t_, _r_)  TIMING_CMD(((_r_) += MPI_Wtime() - (_t_));)


typedef SL_INTEGER_C slint_t;
#define slint_fmt  SL_INTEGER_FMT

#ifndef PINT_T
# define PINT_T
typedef PARAM_INTEGER_C pint_t;
#endif


#endif /* __COMMON_H__ */
