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



#ifndef FCS_RESULT_INCLUDED
#define FCS_RESULT_INCLUDED

#include "FCSResult_p.h"

/**
 * @file FCSResult.h
 * @brief FCSResult-object is used for the error handling in the
 * ScaFaCos-Library.
 * @author Lidia Westphal, Rene Halver
 */

/*
 * FCSResult structure


 */
typedef struct FCSResult_t {
	fcs_int code;
	char *origin;
	char *message;
} FCSResult_t;


#endif
