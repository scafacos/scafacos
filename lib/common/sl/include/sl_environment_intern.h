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


#ifndef __SL_ENVIRONMENT_INTERN_H__
#define __SL_ENVIRONMENT_INTERN_H__


#ifndef ENV_TRACE_IF
# ifdef GLOBAL_TRACE_IF
#  define ENV_TRACE_IF  GLOBAL_TRACE_IF
# else
/*#  define ENV_TRACE_IF  (sl_default_context.mpi.rank == -1)*/
#  define ENV_TRACE_IF  (SL_PROC_RANK == -1)
# endif
#endif


#define z_alloc_hook(_n_, _s_, _p_, _fi_, _li_, _fu_)  Z_MOP( \
  Z_TRACE_IF_(_fi_, _li_, _fu_, ENV_TRACE_IF, "z_alloc: %" slint_fmt " * %" slint_fmt " -> %p", (_n_), (_s_), (_p_)); \
  rti_minc_alloc(); \
  rti_malloc((_n_) * (_s_));)

#define z_free_hook(_p_)  Z_MOP( \
  Z_TRACE_IF(ENV_TRACE_IF, "z_free: %p", (_p_)); \
  rti_minc_free();)

#define z_alloca_hook(_n_, _s_, _p_, _fi_, _li_, _fu_) \
  Z_TRACE_IF_(_fi_, _li_, _fu_, ENV_TRACE_IF, "z_alloca: %" slint_fmt " * %" slint_fmt " -> %p", (_n_), (_s_), (_p_))

#define z_freea_hook(_p_) \
  Z_TRACE_IF(ENV_TRACE_IF, "z_freea: %p", (_p_))


#endif /* __SL_ENVIRONMENT_INTERN_H__ */
