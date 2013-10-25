/*
 *  Copyright (C) 2011, 2012, 2013 Michael Hofmann
 *  
 *  This file is part of ScaFaCoS/FMM.
 *  
 *  ScaFaCoS/FMM is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *  
 *  ScaFaCoS/FMM is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser Public License for more details.
 *  
 *  You should have received a copy of the GNU Lesser Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *  
 */


#ifndef __FMM_SORT_H__
#define __FMM_SORT_H__


#include "config_fmm_sort.h"
#include "rename_fmm_sort.h"

#include "common.h"


#ifndef NO_SL_FRONT
void fmm_sort_front_mem(void *mem0, void *mem1, pint_t *mem_sizes, pint_t *depth, pint_t *subx, pint_t *n, front_(slkey_t) *ibox, front_(sldata0_t) *xyz, front_(sldata1_t) *q, pint_t *addr_desc, void *addr, front_(sldata2_t) *scr, pint_t *type);
void fmm_sort_front(                                               pint_t *depth, pint_t *subx, pint_t *n, front_(slkey_t) *ibox, front_(sldata0_t) *xyz, front_(sldata1_t) *q, pint_t *addr_desc, void *addr, front_(sldata2_t) *scr, pint_t *type);

void fmm_sort_front_3bit_mem(void *mem0, void *mem1, pint_t *mem_sizes, pint_t *subx, pint_t *n, front_(slkey_t) *ibox, front_(sldata0_t) *xyz, front_(sldata1_t) *q, pint_t *addr_desc, void *addr, front_(sldata2_t) *scr, pint_t *type);
void fmm_sort_front_3bit(                                               pint_t *subx, pint_t *n, front_(slkey_t) *ibox, front_(sldata0_t) *xyz, front_(sldata1_t) *q, pint_t *addr_desc, void *addr, front_(sldata2_t) *scr, pint_t *type);
#endif

#ifndef NO_SL_BACK
void fmm_sort_back_mem(void *mem0, void *mem1, pint_t *mem_sizes, pint_t *n, back_(slkey_t) *addr, back_(sldata0_t) *q, back_(sldata1_t) *xyz, back_(sldata2_t) *pot, back_(sldata3_t) *grad, pint_t *type);
void fmm_sort_back(                                               pint_t *n, back_(slkey_t) *addr, back_(sldata0_t) *q, back_(sldata1_t) *xyz, back_(sldata2_t) *pot, back_(sldata3_t) *grad, pint_t *type);
#endif


#endif /* __FMM_SORT_H__ */
