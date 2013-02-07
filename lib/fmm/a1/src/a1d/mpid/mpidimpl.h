/* -*- Mode: C; c-basic-offset:4 ; -*- */
/*
 *  Copyright (C) 2010 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#include "a1.h"
#include "a1u.h"
#include "a1d.h"
#include "mpi.h"
#include <assert.h>

/*************************************************
 *                 Constants                     *
 ************************************************/



/*************************************************
 *                  Macros                       *
 *************************************************/



/*************************************************
 *             Data Structures                   *
 *************************************************/

typedef struct
{
   size_t my_rank;
   size_t num_ranks;
} A1D_Process_info_t;


/*************************************************
 *             Global variables                  *
 ************************************************/



/************************************************* 
 *             Function Prototypes               *
 ************************************************/

