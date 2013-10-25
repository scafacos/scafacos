/*
  Copyright (C) 2011-2012 Rene Halver

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



#ifndef FCS_RESULT_P_INCLUDED
#define FCS_RESULT_P_INCLUDED

#include "FCSDefinitions.h"


/**
 * @file FCSResult.h
 * @brief FCSResult-object is used for the error handling in the
 * ScaFaCos-Library.
 * @author Lidia Westphal, Rene Halver
 */

/**
 * @brief a handle for the FCSResult-object which should be used instead of the
 * data structure itself
 */
typedef struct FCSResult_t *FCSResult;

#define FCS_RESULT_SUCCESS  NULL


#ifdef __cplusplus
extern "C" {
#endif

/**
* @brief initializes an FCSResult object
* @param code is the return code
* @param where is the name of the function the error object was made in
* @param message is the error message
* @return pointer to the FCSResult-object
*/
FCSResult fcsResult_create(fcs_int code, const char *origin, const char *message, ...);

/**
* @brief destroys an FCSResult object
* @param err is the pointer to the FCSResult object to destroy
* @return 0 if successful, otherwise -1
*/
fcs_int fcsResult_destroy(FCSResult err);

/**
 * @brief returns error code
 * @param err is the pointer to the FCSResult object
 * @return fcs_int as an return code
 */
fcs_int fcsResult_getReturnCode(FCSResult err);

/**
 * @brief returns an error message
 * @param err is the pointer to the FCSResult object
 * @return a string containing the corresonding error message
 */
const char *fcsResult_getErrorMessage(FCSResult err);

/**
 * @brief returns a source name of the error (function name)
 * @param err is the pointer to the FCSResult object
 * @return a string containing the name of the function where
 * the error occured
 */
const char *fcsResult_getErrorSource(FCSResult err);

/**
 * @brief prints return code, message and source to stdout
 * @param err is the pofcs_inter to the FCSResult object
 */
void fcsResult_printResult(FCSResult err);

#ifdef __cplusplus
}
#endif

#endif
