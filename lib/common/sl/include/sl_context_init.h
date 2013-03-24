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

/* src/base/base.c */
  {
 -2,
  },
#ifdef SL_USE_RTI
 { },
#endif
  {
 sort_radix_ip_threshold,
 sort_radix_db_threshold,
 sort_radix_db_threshold,
  },
  {
 sort_radix_iter_threshold,
  },
/* src/base_mpi/base_mpi.c */
#ifdef SL_USE_MPI
  {
 MPI_DATATYPE_NULL,
 MPI_DATATYPE_NULL,
 MPI_DATATYPE_NULL,
 MPI_DATATYPE_NULL,
 MPI_DATATYPE_NULL,
 MPI_DATATYPE_NULL,
{
#define xelem_call_data      MPI_DATATYPE_NULL,
#define xelem_call_data_not  MPI_DATATYPE_NULL,
#include "sl_xelem_call.h"
MPI_DATATYPE_NULL
},
 -2,
  },
#endif
#ifdef SL_USE_MPI
  {
 NULL,
 0,
 -1,
  },
#endif
#ifdef SL_USE_MPI
  {
 { 0 },
 8,
 0,
 0,
 0,
  },
#endif
#ifdef SL_USE_MPI
  {
 0,
 0,
 0,
  },
#endif
#ifdef SL_USE_MPI
  {
#ifdef MSEG_ROOT
 -1,
#endif
#ifdef MSEG_BORDER_UPDATE_REDUCTION
 0.0,
 0.0,
#endif
#ifdef MSEG_FORWARD_ONLY
 0,
#endif
#ifdef MSEG_INFO
 0,
 NULL,
 0,
 0,
#endif
 -1,
 SL_MSEG_FM_EXACT,
  },
#endif
#ifdef SL_USE_MPI
  {
#ifdef MSS_ROOT
 -1,
#endif
  },
#endif
#ifdef SL_USE_MPI
  {
 { 0 },
 0,
  },
#endif
#ifdef SL_USE_MPI
  {
 { 0 },
 0,
 NULL,
  },
#endif
#ifdef SL_USE_MPI
  {
 { 0 },
 { 0 },
 { 0 },
 0,
 0,
 0,
 0,
 0,
  },
#endif
