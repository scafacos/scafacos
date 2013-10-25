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


#ifndef __CONFIG_TYPES_H__
#define __CONFIG_TYPES_H__


/* integer size (default: INTEGER64) */
/*#define INTEGER32
#define INTEGER64*/

/* integer size of parameters (default: equal to integer size) */
/*#define PARAM_INTEGER32
#define PARAM_INTEGER64*/

/* sorting library internal integer size (default: integer size of parameters) */
/*#define SL_INTEGER32
#define SL_INTEGER64*/

/* real size (default: REAL64) */
/*#define REAL32
#define REAL64
#define REAL128*/


/* make defaults */
#if !defined(INTEGER32) && !defined(INTEGER64)
# define INTEGER64
#endif

#if !defined(PARAM_INTEGER32) && !defined(PARAM_INTEGER64)
# ifdef INTEGER32
#  define PARAM_INTEGER32
# else
#  define PARAM_INTEGER64
# endif
#endif

#if !defined(SL_INTEGER32) && !defined(SL_INTEGER64)
# ifdef PARAM_INTEGER32
#  define SL_INTEGER32
# else
#  define SL_INTEGER64
# endif
#endif

#if !defined(REAL32) && !defined(REAL64) && !defined(REAL128)
# define REAL64
#endif


/* make types */
#ifdef INTEGER32
# define INTEGER_C    FINT4_TYPE_C
# define INTEGER_MPI  FINT4_TYPE_MPI
# define INTEGER_FMT  FINT4_TYPE_FMT
#endif
#ifdef INTEGER64
# define INTEGER_C    FINT8_TYPE_C
# define INTEGER_MPI  FINT8_TYPE_MPI
# define INTEGER_FMT  FINT8_TYPE_FMT
#endif

#ifdef PARAM_INTEGER32
# define PARAM_INTEGER_C    FINT4_TYPE_C
# define PARAM_INTEGER_MPI  FINT4_TYPE_MPI
# define PARAM_INTEGER_FMT  FINT4_TYPE_FMT
#endif
#ifdef PARAM_INTEGER64
# define PARAM_INTEGER_C    FINT8_TYPE_C
# define PARAM_INTEGER_MPI  FINT8_TYPE_MPI
# define PARAM_INTEGER_FMT  FINT8_TYPE_FMT
#endif

#ifdef SL_INTEGER32
# define SL_INTEGER_C    FINT4_TYPE_C
# define SL_INTEGER_MPI  FINT4_TYPE_MPI
# define SL_INTEGER_FMT  FINT4_TYPE_FMT
#endif
#ifdef SL_INTEGER64
# define SL_INTEGER_C    FINT8_TYPE_C
# define SL_INTEGER_MPI  FINT8_TYPE_MPI
# define SL_INTEGER_FMT  FINT8_TYPE_FMT
#endif

#ifdef REAL32
# define REAL_C    FREAL4_TYPE_C
# define REAL_MPI  FREAL4_TYPE_MPI
# define REAL_FMT  FREAL4_TYPE_FMT
#endif
#ifdef REAL64
# define REAL_C    FREAL8_TYPE_C
# define REAL_MPI  FREAL8_TYPE_MPI
# define REAL_FMT  FREAL8_TYPE_FMT
#endif
#ifdef REAL128
# define REAL_C    FREAL16_TYPE_C
# define REAL_MPI  FREAL16_TYPE_MPI
# define REAL_FMT  FREAL16_TYPE_FMT
#endif


#endif /* __CONFIG_TYPES_H__ */
