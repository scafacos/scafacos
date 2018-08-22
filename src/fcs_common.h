/*
  Copyright (C) 2011, 2012, 2013, 2014 Michael Hofmann

  This file is part of ScaFaCoS.

  ScaFaCoS is free software: you can redistribute it and/or modify
  it under the terms of the GNU Lesser Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  ScaFaCoS is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU Lesser Public License for more details.

  You should have received a copy of the GNU Lesser Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>. 
*/


#ifndef __FCS_COMMON_H__
#define __FCS_COMMON_H__


#include <mpi.h>

#include "fcs_result.h"
#include "fcs_interface.h"


#ifdef __cplusplus
extern "C" {
#endif


FCSResult fcs_common_set_parameter(FCS handle, fcs_bool continue_on_errors, char **current, char **next, fcs_int *matched);
FCSResult fcs_common_print_parameters(FCS handle);


#ifdef __cplusplus
}
#endif


#endif /* __FCS_COMMON_H__ */
