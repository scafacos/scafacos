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

#include "FCSResult.h"


/*
 * creates an FCSResult object
 */
FCSResult fcsResult_create(fcs_int code, const char *origin, const char *message, ...) {

#define FCS_RESULT_MAX_MESSAGE_LENGTH  512
  char printed_message[FCS_RESULT_MAX_MESSAGE_LENGTH];
  va_list argp;
  fcs_int origin_size, message_size;

  va_start(argp, message);
  vsnprintf(printed_message, FCS_RESULT_MAX_MESSAGE_LENGTH, message, argp);
  va_end(argp);

  origin_size = strlen(origin) + 1;
  message_size = strlen(printed_message) + 1;

	FCSResult err = (FCSResult) malloc(sizeof(FCSResult_t) + origin_size + message_size);
	if(err == NULL)
		return err;

	err->code = code;
	err->origin = ((char *) err) + sizeof(FCSResult_t);
	err->message = err->origin + origin_size;

  strcpy(err->origin, origin);
  strcpy(err->message, printed_message);

	return err;
}


/*
 * destroys an FCSResult object, should be used instead of free
 */
fcs_int fcsResult_destroy(FCSResult err) {

	if(err == NULL) {
		return 0;
	}

	free(err);

	return 0;
}


fcs_int fcsResult_getReturnCode(FCSResult err) {
	if (err == NULL) return -1;
	return err->code;
}

const char *fcsResult_getErrorMessage(FCSResult err) {
	if(err == NULL)
		return NULL;
	return err->message;
}

const char *fcsResult_getErrorSource(FCSResult err) {
	if(err == NULL)
		return NULL;
	return err->origin;
}

void fcsResult_printResult(FCSResult err)
{
	if(err == NULL)
	{
		printf("OK\n");
		return;
	}

	switch (fcsResult_getReturnCode(err))
	{
		case FCS_SUCCESS:
			printf("Return code: FCS_SUCCESS\n");
			break;
		case FCS_ALLOC_FAILED:
			printf("Return code: FCS_ALLOC_FAILED\n");
			break;
		case FCS_NULL_ARGUMENT:
			printf("Return code: FCS_NULL_ARGUMENT\n");
			break;
		case FCS_WRONG_ARGUMENT:
			printf("Return code: FCS_WRONG_ARGUMENT\n");
			break;
		case FCS_MISSING_ELEMENT:
			printf("Return code: FCS_MISSING_ELEMENT\n");
			break;
		case FCS_LOGICAL_ERROR:
			printf("Return code: FCS_LOGICAL_ERROR\n");
			break;
		case FCS_INCOMPATIBLE_METHOD:
			printf("Return code: FCS_INCOMPATIBLE_METHOD\n");
			break;
		case FCS_MPI_ERROR:
			printf("Return code: FCS_MPI_ERROR\n");
			break;
		default:
			printf("Return code: %" FCS_LMOD_INT "d\n", fcsResult_getReturnCode(err));
			break;
	}

	if (fcsResult_getErrorSource(err) == NULL)
		printf("Result message: no error source\n");
	else
		printf("Result source: %s\n", fcsResult_getErrorSource(err));

	if (fcsResult_getErrorMessage(err) == NULL)
		printf("Result message: no error message\n");
	else
		printf("Result message: %s\n", fcsResult_getErrorMessage(err));
}
