/*
 *  cuboid.c
 *
 *  Copyright 2006 Matthias Bolten. All rights reserved.
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <stdio.h>
#include <stdlib.h>
#include "cuboid.h"

double*** cuboid_alloc(int m, int n, int o)
{
	double *field;
	double ***cuboid;
	int i,j;
	
	field = (double*) malloc(sizeof(double)*m*n*o);
	if (field==NULL) 
	  {
	    printf("Malloc failed!\n");
	    exit(1);
	  }
	cuboid = (double***) malloc(sizeof(double**)*m);
	if (cuboid==NULL) 
	  {
	    printf("Malloc failed!\n");
	    exit(1);
	  }
	for (i=0;i<m;i++) {
		cuboid[i] = (double**) malloc(sizeof(double*)*n);
		if (cuboid[i]==NULL) 
		  {
		    printf("Malloc failed!\n");
		    exit(1);
		  }
		for (j=0;j<n;j++) {
			cuboid[i][j] = &field[i*n*o+j*o];
		}
	}
	
	return cuboid;
}

double*** cube_alloc(int n)
{
	return cuboid_alloc(n,n,n);
}

void cuboid_free(double ***cuboid, int m, int n, int o)
{
	int i;
	
	free(cuboid[0][0]);
	for (i=0;i<m;i++) {
		free(cuboid[i]);
	}
	free(cuboid);

	return;
}

void cube_free(double ***cube,int n)
{
	cuboid_free(cube,n,n,n);
	
	return;
}
