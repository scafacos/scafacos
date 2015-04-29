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


#ifndef __Z_CONFIG_H__
#define __Z_CONFIG_H__


#ifndef HAVE_CONFIG_H


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
