/*
 *  Copyright (C) 2011, 2012, 2013, 2014, 2015, 2016 Michael Hofmann
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


#ifndef __Z_CONFIG_H__
#define __Z_CONFIG_H__


#ifndef HAVE_CONFIG_H

/* C standard library */
#define HAVE_ASSERT_H  1
#define HAVE_CTYPE_H  1
#define HAVE_ERRNO_H  1
#define HAVE_FLOAT_H  1
#define HAVE_LIMITS_H  1
#define HAVE_LOCALE_H  1
#define HAVE_MATH_H  1
#define HAVE_SETJMP_H  1
#define HAVE_SIGNAL_H  1
#define HAVE_STDARG_H  1
#define HAVE_STDDEF_H  1
#define HAVE_STDIO_H  1
#define HAVE_STDLIB_H  1
#define HAVE_STRING_H  1
#define HAVE_TIME_H  1

/* AMD1 extensions */
#if __STDC_VERSION__ >= 199409L
# define HAVE_ISO646_H  1
# define HAVE_WCHAR_H  1
# define HAVE_WCTYPE_H  1
#endif

/* C99 extensions */
#if __STDC_VERSION__ >= 199901L
# define HAVE_COMPLEX_H  1
# define HAVE_FENV_H  1
# define HAVE_INTTYPES_H  1
# define HAVE_STDBOOL_H  1
# define HAVE_STDINT_H  1
# define HAVE_TGMATH_H  1
#endif

/* C11 extensions */
#if __STDC_VERSION__ >= 201112L
# define HAVE_STDALIGN_H  1
# define HAVE_STDATOMIC_H  1
# define HAVE_STDNORETURN_H  1
# define HAVE_THREADS_H  1
# define HAVE_UCHAR_H  1
#endif


/* features of the GNU C Compiler */
#ifdef __GNUC__
# ifdef __STDC_VERSION__
#  define HAVE_ROUND  1
# else
#  define HAVE_RANDOM  1
#  define HAVE_SRANDOM  1
# endif
#endif


#ifdef __bgp__
# define HAVE_SPI_KERNEL_INTERFACE_H  1
# define HAVE_COMMON_BGP_PERSONALITY_H  1
# define HAVE_COMMON_BGP_PERSONALITY_INLINES_H  1
# define HAVE__BGP_PERSONALITY_T  1
#endif


#ifdef __bgq__
# define HAVE_MPIX_H  1
# define HAVE_MPIX_HARDWARE_T  1
#endif


#endif /* HAVE_CONFIG_H */


#if !defined(HAVE_MPI_IN_PLACE) && !defined(IGNORE_MPI_IN_PLACE)
# if defined(MPI_VERSION) && (MPI_VERSION >= 2)
#  define HAVE_MPI_IN_PLACE  1
# endif
#endif


#endif
