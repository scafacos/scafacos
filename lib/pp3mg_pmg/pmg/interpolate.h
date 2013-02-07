/*
 *  interpolate.h
 *
 *  Copyright 2006 Matthias Bolten. All rights reserved.
 *
 */

#ifndef _INTERPOLATE__H_
#define _INTERPOLATE__H_

#include "mg.h"

void interpolate_prepare(double ***fine, double ***coarse, mg_data* data, int level );

void interpolate_finish(double ***fine, mg_data* data, int level );

#endif /* ifndef _INTERPOLATE__H_ */
