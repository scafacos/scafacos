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


#ifndef __CONFIG_FMM_SORT_H__
#define __CONFIG_FMM_SORT_H__


#include "fortran2c_types.h"

/* FCS mode? */
#if defined(fcs_int) && defined(fcs_float)

/* fixed 8-byte integer required in fmmkinds.h */
#define INTEGER_C          FINT8_TYPE_C
#define INTEGER_MPI        FINT8_TYPE_MPI
#define INTEGER_FMT        FINT8_TYPE_FMT

#define PARAM_INTEGER_C    FINT8_TYPE_C
#define PARAM_INTEGER_MPI  FINT8_TYPE_MPI
#define PARAM_INTEGER_FMT  FINT8_TYPE_FMT

/* use FCS int as sl-internal integer */
#define SL_INTEGER_C       fcs_int
#define SL_INTEGER_MPI     FCS_MPI_INT
#define SL_INTEGER_FMT     FCS_LMOD_INT "d"

/* use FCS float as real data type */
#define REAL_C             fcs_float
#define REAL_MPI           FCS_MPI_FLOAT
#define REAL_FMT           FCS_LMOD_FLOAT "f"

#define SL_FMM_PREFIX  fcs_

#else /* FCS mode? */

# include "config_types.h"

#endif /* FCS mode? */


#define SL_FMM_CONFIG_CONCAT(_a_, _b_)   SL_FMM_CONFIG_CONCAT_(_a_, _b_)
#define SL_FMM_CONFIG_CONCAT_(_a_, _b_)  _a_##_b_

#ifdef SL_FMM_PREFIX
# define SL_FMM_CONFIG_VAR(_v_)    SL_FMM_CONFIG_CONCAT(SL_FMM_PREFIX, _v_)
#else
# define SL_FMM_CONFIG_VAR(_v_)    _v_
#endif


extern INTEGER_C SL_FMM_CONFIG_VAR(fmm_front_key_mask);
extern int SL_FMM_CONFIG_VAR(fmm_front_aX);


#if defined(FCS_ENABLE_DEBUG_SL_FMM)
# define DO_DEBUG
#endif
#define DEBUG_PRINT_PREFIX  "SL_FMM_DEBUG: "

#if defined(FCS_ENABLE_INFO_SL_FMM)
# define DO_INFO
# error
#endif
#define INFO_PRINT_PREFIX  "SL_FMM_INFO: "

#if defined(FCS_ENABLE_TIMING_SL_FMM)
# define DO_TIMING
#endif
#define TIMING_PRINT_PREFIX  "SL_FMM_TIMING: "


#define WITH_SORT_FRONT_LOAD
#define WITH_FCOMM


#endif /* __CONFIG_FMM_SORT_H__ */
