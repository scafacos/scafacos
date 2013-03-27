/*
 * ghosted_grid.h
 *
 * This file contains definitions of functions found in the file
 * "ghosted_grid.c".
 *
 * Based on Fortran code by Matthias Bolten.
 *
 * Authors: Matthias Bolten, Stephanie Friedhoff
 * Created: 2009/06/01
 *
 * Copyright 2009, 2010, 2011, 2012 Matthias Bolten, Stephanie Friedhoff
 * All rights reserved.
 *
 */

#ifndef _GHOSTED_GRID__H_
#define _GHOSTED_GRID__H_

void pp3mg_update_ghosts( double*** u, int m, int n, int o, int ghosts, MPI_Comm mpi_comm_cart );

#endif  /* ifndef _GHOSTED_GRID__H_ */
