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
 * @file FCSResult_p.h
 * @brief public interface definitions for the FCSResult-object that is used
 * for handling the return state of the ScaFaCoS library functions
 * @author Lidia Westphal, Rene Halver, Michael Hofmann
 */


#ifndef _FCSRESULT_P_H
#define _FCSRESULT_P_H

#include "FCSDefinitions.h"


#ifdef __cplusplus
extern "C" {
#endif


/**
 * @brief FCSResult-object that is used for handling the return state of the
 * ScaFaCoS library functions
 */
typedef struct FCSResult_t *FCSResult;

#define FCS_RESULT_SUCCESS  NULL

#define FCS_RESULT_MAX_FUNCTION_LENGTH   64
#define FCS_RESULT_MAX_MESSAGE_LENGTH   512

/**
 * @brief function to destroy an FCSResult-object
 * @param result FCSResult-object containing the return state
 */
void fcs_result_destroy(FCSResult result);

/**
 * @brief function to return the return code associated with an return state
 * @param result FCSResult-object containing the return state
 * @return return code associated with the return state
 */
fcs_int fcs_result_get_return_code(FCSResult result);

/**
 * @brief function to return the function name associated with an return state
 * @param result FCSResult-object containing the return state
 * @return function name associated with the return state
 */
const char *fcs_result_get_function(FCSResult result);

/**
 * @brief function to return the description message associated with an return state
 * @param result FCSResult-object containing the return state
 * @return description message associated with the return state
 */
const char *fcs_result_get_message(FCSResult result);

/**
 * @brief function to print the return state to stdout
 * @param result FCSResult-object containing the return state
 */
void fcs_result_print_result(FCSResult result);


#ifdef FCS_ENABLE_DEPRECATED
#define fcsResult_destroy          fcs_result_destroy
#define fcsResult_getReturnCode    fcs_result_get_return_code
#define fcsResult_getErrorMessage  fcs_result_get_message
#define fcsResult_getErrorSource   fcs_result_get_function
#define fcsResult_printResult      fcs_result_print_result
#endif


#ifdef __cplusplus
}
#endif


#endif
