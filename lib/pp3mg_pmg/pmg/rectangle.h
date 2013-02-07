/*
 *  rectangle.h
 *
 *  Copyright 2006 Matthias Bolten. All rights reserved.
 *
 */

#ifndef _RECTANGLE__H_
#define _RECTANGLE__H_

double** rectangle_alloc(int m, int n);
double** square_alloc(int n);

void rectangle_free(double **rectangle, int m, int n);
void square_free(double **square, int n);

#endif /* ifndef _RECTANGLE__H_ */
