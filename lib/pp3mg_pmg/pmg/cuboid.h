/*
 *  cuboid.h
 *
 *  Copyright 2006 Matthias Bolten. All rights reserved.
 *
 */

#ifndef _CUBOID__H_
#define _CUBOID__H_

double*** cuboid_alloc(int m, int n, int o);
double*** cube_alloc(int n);

void cuboid_free(double ***cuboid, int m, int n, int o);
void cube_free(double ***cuboid, int n);

#endif /* ifndef _CUBOID__H_ */
