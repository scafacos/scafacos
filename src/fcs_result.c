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



#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdarg.h>

#include "fcs_result.h"


/*
 * @brief data structure that is used for storing the return state of the
 * ScaFaCoS library functions
 */
typedef struct FCSResult_t {
  fcs_int error_code;
  char function[FCS_RESULT_MAX_FUNCTION_LENGTH];
  char message[FCS_RESULT_MAX_MESSAGE_LENGTH];
} FCSResult_t;


const static FCSResult_t fcs_result_create_error = { FCS_ERROR_RESULT_CREATE, "fcs_result_create", "failed to create an FCSResult object" };


/**
 * create an FCSResult-object for storing the return state
 */
FCSResult fcs_result_create(fcs_int code, const char *function, const char *message, ...)
{
  FCSResult result;


  result = malloc(sizeof(FCSResult_t));

  if (result == NULL) return (FCSResult) &fcs_result_create_error;

  result->error_code = code;

  if (function)
    strncpy(result->function, function, FCS_RESULT_MAX_FUNCTION_LENGTH);
  else
    result->function[0] = '\0';

  if (message)
  {
    va_list argp;
    va_start(argp, message);
    vsnprintf(result->message, FCS_RESULT_MAX_MESSAGE_LENGTH, message, argp);
    va_end(argp);

  } else result->message[0] = '\0';

  return result;
}


/**
 * destroy an FCSResult-object
 */
void fcs_result_destroy(FCSResult result)
{
  if (result == FCS_RESULT_SUCCESS) return;

  if (result != &fcs_result_create_error) free(result);
}


/**
 * return the return code associated with an return state
 */
fcs_int fcs_result_get_return_code(FCSResult result)
{
  if (result == FCS_RESULT_SUCCESS) return FCS_SUCCESS;

  return result->error_code;
}


/**
 * return the function name associated with an return state
 */
const char *fcs_result_get_function(FCSResult result)
{
  if (result == FCS_RESULT_SUCCESS) return NULL;

  return result->function;
}


/**
 * return the description message associated with an return state
 */
const char *fcs_result_get_message(FCSResult result)
{
  if (result == FCS_RESULT_SUCCESS) return NULL;

  return result->message;
}


/**
 * print the return state to stdout
 */
void fcs_result_print_result(FCSResult result)
{
  switch (fcs_result_get_return_code(result))
  {
    case FCS_SUCCESS:
      printf("return code: FCS_SUCCESS\n");
      break;
    case FCS_ERROR_ALLOC_FAILED:
      printf("return code: FCS_ALLOC_FAILED\n");
      break;
    case FCS_ERROR_NULL_ARGUMENT:
      printf("return code: FCS_ERROR_NULL_ARGUMENT\n");
      break;
    case FCS_ERROR_WRONG_ARGUMENT:
      printf("return code: FCS_ERROR_WRONG_ARGUMENT\n");
      break;
    case FCS_ERROR_MISSING_ELEMENT:
      printf("return code: FCS_ERROR_MISSING_ELEMENT\n");
      break;
    case FCS_ERROR_INCOMPATIBLE_METHOD:
      printf("return code: FCS_ERROR_INCOMPATIBLE_METHOD\n");
      break;
    case FCS_ERROR_RESULT_CREATE:
      printf("return code: FCS_ERROR_RESULT_CREATE\n");
      break;
    default:
      printf("return code: %" FCS_LMOD_INT "d\n", fcs_result_get_return_code(result));
      break;
  }

  if (result == FCS_RESULT_SUCCESS || fcs_result_get_function(result)[0] == '\0')
    printf("function: not available\n");
  else
    printf("function: %s\n", fcs_result_get_function(result));

  if (result == FCS_RESULT_SUCCESS || fcs_result_get_message(result)[0] == '\0')
    printf("message: not available\n");
  else
    printf("message: %s\n", fcs_result_get_message(result));
}
