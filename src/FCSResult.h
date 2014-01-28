/*
  Copyright (C) 2011, 2012, 2013 Lidia Westphal, Rene Halver, Michael Hofmann

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


/**
 * @file FCSResult.h
 * @brief internal interface definitions for the FCSResult-object that is used
 * for handling the return state of the ScaFaCoS library functions
 * @author Lidia Westphal, Rene Halver, Michael Hofmann
 */


#ifndef FCS_RESULT_INCLUDED
#define FCS_RESULT_INCLUDED

#include "FCSResult_p.h"


/**
 * @brief function to create an FCSResult-object for storing the return state
 * @param code return code to associate with the return state
 * @param function name of a function to associate with the return state
 * @param message description message to associate with the return state
 * @return FCSResult-object containing the return state
 */
FCSResult fcs_result_create(fcs_int code, const char *function, const char *message, ...);


#endif
