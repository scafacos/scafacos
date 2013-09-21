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


#ifndef __MPI_FMM_RESORT_H__
#define __MPI_FMM_RESORT_H__


#include "../../common/resort/resort.h"


typedef struct _fcs_fmm_resort_t
{
  MPI_Comm comm;

  fcs_resort_t resort;

} *fcs_fmm_resort_t;

#define FCS_FMM_RESORT_NULL  NULL


void fcs_fmm_resort_create(fcs_fmm_resort_t *fr, fcs_int nparticles, MPI_Comm comm);
void fcs_fmm_resort_destroy(fcs_fmm_resort_t *fr);


#endif /* __MPI_FMM_RESORT_H__ */
