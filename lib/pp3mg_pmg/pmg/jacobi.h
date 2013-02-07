/*
 *  jacobi.h
 *
 *  Copyright 2006 Matthias Bolten. All rights reserved.
 *
 */

#ifndef _JACOBI__H_
#define _JACOBI__H_

#include "mg.h"

double jacobi( double*** v, double*** f, double*** r, mg_data* data, int level, 
	       int maxiter );


#endif /* ifndef _JACOBI__H_ */
