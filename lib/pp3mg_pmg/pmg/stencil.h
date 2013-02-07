/*
 *  stencil.h
 *
 *  Copyright 2006 Matthias Bolten. All rights reserved.
 *
 */

#ifndef _STENCIL__H_
#define _STENCIL__H_

#include "mg.h"

void find_min_max( double a, double b, double c, double d, double e, double f, double g, double* min, double* max );
double calculate_omega( double a, double b, double c, double d, double e, double f, double g, double min, double max );
void set_stencil( mg_data* data, int level, int dim );

#endif /* ifndef _STENCIL__H_ */
