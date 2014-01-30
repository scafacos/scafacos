/*
  Copyright (C) 2014 Olaf Lenz
  
  This file is part of ScaFaCoS.
  
  ScaFaCoS is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  ScaFaCoS is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>. 
*/
#ifndef _P3M_CONFIG_HPP
#define _P3M_CONFIG_HPP

#include <config.h>

#ifdef FCS_ENABLE_DEBUG
#define P3M_ENABLE_DEBUG 1
#endif

#ifdef FCS_ENABLE_INFO
#define P3M_ENABLE_INFO 1
#endif

#define P3M_LMOD_FLOAT FCS_LMOD_FLOAT
#define P3M_LMOD_INT FCS_LMOD_INT
#define FFLOAT "%" FCS_LMOD_FLOAT "f"
#define FFLOATE "%" FCS_LMOD_FLOAT "e"
#define F3FLOAT "(" FFLOAT ", " FFLOAT ", " FFLOAT ")"
#define FINT "%" FCS_LMOD_INT "d"
#define F3INT "(" FINT ", " FINT ", " FINT ")"

#define P3M_MPI_FLOAT FCS_MPI_FLOAT
#define P3M_MPI_INT FCS_MPI_INT

#define p3m_float fcs_float
#define p3m_int fcs_int

/** maximal precision */
#ifdef FCS_FLOAT_IS_DOUBLE
#define ROUND_ERROR_PREC 1.0e-14
#else
#define ROUND_ERROR_PREC 1.0e-6
#endif

#endif
