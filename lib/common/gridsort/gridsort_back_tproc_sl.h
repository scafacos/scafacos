/*
  Copyright (C) 2011, 2012, 2013 Michael Hofmann
  
  This file is part of ScaFaCoS.
  
  ScaFaCoS is free software: you can redistribute it and/or modify
  it under the terms of the GNU Lesser Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  ScaFaCoS is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU Lesser Public License for more details.
  
  You should have received a copy of the GNU Lesser Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/


#ifndef __GRIDSORT_BACK_TPROC_SL_H__
#define __GRIDSORT_BACK_TPROC_SL_H__


#define GRIDSORT_PREFIX
#define GRIDSORT_INDEX(_b_, _d_, _x_)  (_b_)->keys[_x_]

/* back_fp */
#define GRIDSORT_ELEM_BUF_T            fcs_back_fp_elements_t
#define GRIDSORT_ELEM_INDEX_T          fcs_back_fp_slint_t
#define GRIDSORT_INDEX(_b_, _d_, _x_)  (_b_)->keys[_x_]
/*#ifdef GRIDSORT_BACK_TPROC_EXDEF
# define DEFINE_GRIDSORT_BACK_TPROC_EXDEF(_xd_, _tp_, _s_...)  fcs_back_fp_SPEC_DEFINE_TPROC(_xd_, _tp_, _s_)
#endif*/
#define GRIDSORT_BACK_TPROC_NAME       gridsort_back_fp_tproc
#include "gridsort_back_tproc.h"

#undef DEFINE_GRIDSORT_BACK_TPROC_EXDEF
#undef GRIDSORT_ELEM_BUF_T
#undef GRIDSORT_ELEM_INDEX_T

/* back_f_ */
#define GRIDSORT_ELEM_BUF_T            fcs_back_f__elements_t
#define GRIDSORT_ELEM_INDEX_T          fcs_back_f__slint_t
/*#ifdef GRIDSORT_BACK_TPROC_EXDEF
# define DEFINE_GRIDSORT_BACK_TPROC_EXDEF(_xd_, _tp_, _s_...)  fcs_back_f__SPEC_DEFINE_TPROC(_xd_, _tp_, _s_)
#endif*/
#define GRIDSORT_BACK_TPROC_NAME       gridsort_back_f__tproc
#include "gridsort_back_tproc.h"

#undef DEFINE_GRIDSORT_BACK_TPROC_EXDEF
#undef GRIDSORT_ELEM_BUF_T
#undef GRIDSORT_ELEM_INDEX_T

/* back_f_ */
#define GRIDSORT_ELEM_BUF_T            fcs_back__p_elements_t
#define GRIDSORT_ELEM_INDEX_T          fcs_back__p_slint_t
/*#ifdef GRIDSORT_BACK_TPROC_EXDEF
# define DEFINE_GRIDSORT_BACK_TPROC_EXDEF(_xd_, _tp_, _s_...)  fcs_back__p_SPEC_DEFINE_TPROC(_xd_, _tp_, _s_)
#endif*/
#define GRIDSORT_BACK_TPROC_NAME       gridsort_back__p_tproc
#include "gridsort_back_tproc.h"

#undef DEFINE_GRIDSORT_BACK_TPROC_EXDEF
#undef GRIDSORT_ELEM_BUF_T
#undef GRIDSORT_ELEM_INDEX_T

#undef GRIDSORT_INDEX
#undef GRIDSORT_PREFIX


#endif /* __GRIDSORT_BACK_TPROC_SL_H__ */
