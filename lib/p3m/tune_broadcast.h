/*
  Copyright (C) 2013 Olaf Lenz
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
#ifndef _P3M_TUNE_BROADCAST_H 
#define _P3M_TUNE_BROADCAST_H 
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "types.h"
#include "FCSResult.h"

/* Events during tuning */
#define FAILED -1
#define FINISHED 0
#define COMPUTE_ERROR_ESTIMATE 1
#define TEST_RUN 2

void
ifcs_p3m_tune_broadcast_command
(ifcs_p3m_data_struct *d, fcs_int command);

FCSResult
ifcs_p3m_tune_broadcast_slave
(ifcs_p3m_data_struct *d, fcs_int num_particles, fcs_int max_particles,
 fcs_float *positions, fcs_float *charges);

#endif
