/*
 *  galerkin.h
 *
 *  Copyright 2006 Matthias Bolten. All rights reserved.
 *
 */

#ifndef _GALERKIN__H_
#define _GALERKIN__H_

#include "mg.h"

double** calculate_R( double* r, int dim );
void part( double*** A, double** R, int k, int l, int dim, int flag, double * );
double*** galerkin( mg_data* data, int level, int dim );

#endif /* ifndef _GALERKIN__H_ */
