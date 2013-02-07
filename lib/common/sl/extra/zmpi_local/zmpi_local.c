/*
 *  Copyright (C) 2011, 2012, 2013 Michael Hofmann
 *  
 *  This file is part of ScaFaCoS.
 *  
 *  ScaFaCoS is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *  
 *  ScaFaCoS is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser Public License for more details.
 *  
 *  You should have received a copy of the GNU Lesser Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *  

 *  
 *  SL - Sorting Library, michael <dot> hofmann <at> informatik <dot> tu-chemnitz <dot> de
 */


#include "z_pack.h"

#include "zmpi_local.h"
#include "zmpil_simple.h"


int zmpil_simple_create(zmpil_simple_t *mpil, MPI_Datatype type) /* zmpi_func zmpil_simple_create */
{
/*  mpil->type = type;*/

  MPI_Type_get_true_extent(type, &mpil->true_lb, &mpil->true_extent);

  Z_ASSERT(mpil->true_lb == 0);

  return 0;
}


void zmpil_simple_destroy(zmpil_simple_t *mpil) /* zmpi_func zmpil_simple_destroy */
{
}
