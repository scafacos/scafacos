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


#ifndef __DASH_CONF_H__
#define __DASH_CONF_H__


#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#ifdef SLDEBUG
# ifndef SLDEBUG_OUTPUT
#  define SLDEBUG_OUTPUT  SLDEBUG
# endif
# ifndef SLDEBUG_ALLOC
#  define SLDEBUG_ALLOC   SLDEBUG
# endif
#endif

#ifdef SLDEBUG_OUTPUT_NOT
# undef SLDEBUG_OUTPUT
#endif

#ifdef SLDEBUG_ALLOC_NOT
# undef SLDEBUG_ALLOC
#endif


#include "sl_common.h"

/*#include "sl_rename.h"

#include "sl_config.h"
#include "sl_config_intern.h"

#define SL_PROC_RANK  SL_DEFCON(mpi.rank)*/

#define z_mpi_rank  SL_PROC_RANK

#define DS_RENAME

#ifdef SL_PREFIX
# define DS_PREFIX  SL_PREFIX
#endif

typedef sl_int_type_c dsint_t;
#define dsint_fmt  sl_int_type_fmt
#define MPI_DSINT  sl_int_type_mpi

typedef int dspint_t;
#define dspint_fmt  "d"
#define MPI_DSPINT  MPI_INT


#define DASH_MAX_NBUFFERS  2

#define DASH_SYMMETRIC

#define DASH_SCHED_A2AV_AUX_STATIC
#define DASH_SCHED_A2AV_AUX_HEAP

#define DASH_SCHED_A2AV_OVERLAP


#endif /* __DASH_CONF_H__ */
