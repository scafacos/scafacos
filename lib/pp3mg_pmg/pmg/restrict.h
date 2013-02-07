/*
 *  restrict.h
 *
 *  Copyright 2006 Matthias Bolten. All rights reserved.
 *
 */

#ifndef _RESTRICT__H_
#define _RESTRICT__H_

#include "mg.h"

void restrict_inj( double ***fine, double ***coarse, mg_data* data, int level );

void restrict_fw(double ***fine, double ***coarse, mg_data* data, int level );

#endif /* ifndef _RESTRICT__H_ */
