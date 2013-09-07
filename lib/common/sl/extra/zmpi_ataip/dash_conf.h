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


#define DS_RENAME

#ifdef ZMPI_PREFIX
# define DS_PREFIX  DS_CONCAT(ZMPI_PREFIX, zmpi_)
#else
# define DS_PREFIX  zmpi_
#endif


typedef long dsint_t;
#define dsint_fmt   "ld"
#define MPI_DSINT   MPI_LONG

typedef int dspint_t;
#define dspint_fmt  "d"
#define MPI_DSPINT  MPI_INT


#define DASH_MAX_NBUFFERS  2

#define DASH_SYMMETRIC

#define DASH_SCHED_A2AV_AUX_STATIC
#define DASH_SCHED_A2AV_AUX_HEAP

#define DASH_SCHED_A2AV_OVERLAP


#endif /* __DASH_CONF_H__ */
