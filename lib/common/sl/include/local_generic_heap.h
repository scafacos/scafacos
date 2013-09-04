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


#ifndef __LOCAL_GENERIC_HEAP_H__
#define __LOCAL_GENERIC_HEAP_H__


#include "local_generic_heap_conf.h"


#ifdef LGH_RENAME
# include "local_generic_heap_rename.h"
#endif


typedef struct _lgh_segment_t
{
  struct _lgh_segment_t *next;

  lghint_t offset, size;

} lgh_segment_t;


typedef struct _lgh_t
{
  lghint_t nallocs;

  lgh_segment_t *segments;

} lgh_t;


void lgh_create(lgh_t *lgh, lghint_t size);
void lgh_destroy(lgh_t *lgh);

lgh_segment_t *lgh_alloc(lgh_t *lgh, lghint_t size);
lgh_segment_t *lgh_alloc_minmax(lgh_t *lgh, lghint_t min, lghint_t max);
void lgh_free(lgh_t *lgh, lgh_segment_t *seg);


#endif /* __LOCAL_GENERIC_HEAP_H__ */
