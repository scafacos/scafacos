/*
 *  rectangle.c
 *
 *  Copyright 2006 Matthias Bolten. All rights reserved.
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <stdio.h>
#include <stdlib.h>
#include "rectangle.h"

double** rectangle_alloc(int m, int n)
{
	double *field;
	double **rectangle;
	int i;
	
	field = (double*) malloc(sizeof(double)*m*n);
	
	rectangle = (double**) malloc(sizeof(double**)*m);
	for (i=0;i<m;i++) {
	  rectangle[i] = &field[i*n];
	}
	
	return rectangle;
}

double** square_alloc(int n)
{
	return rectangle_alloc(n,n);
}

void rectangle_free(double **rectangle, int m, int n)
{
	free(rectangle[0]);
	free(rectangle);

	return;
}

void square_free(double **square,int n)
{
	rectangle_free(square,n,n);
	
	return;
}
