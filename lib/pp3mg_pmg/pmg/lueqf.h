/*
 *  lueqf.h
 *
 *  Copyright 2006 Matthias Bolten. All rights reserved.
 *
 */

#ifndef _LUEQF__H_
#define _LUEQF__H_

#include "mg.h"

double lueqf_res( double*** v, double*** f, double*** r, mg_data* data, int level );

// double lueqf_res(double ***v, double ***f, double ***r, int m, int n, int o,
// 		 MPI_Comm cart_comm);
void lueqf_invd( double ***r, mg_data* data, int level );

#endif /* ifndef _LUEQF__H_ */
