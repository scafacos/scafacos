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


#ifndef __GRIDSORT_FRONT_TPROC_SL_H__
#define __GRIDSORT_FRONT_TPROC_SL_H__


#define GRIDSORT_PREFIX
#define GRIDSORT_INT_T               fcs_forw_slint_t
#define GRIDSORT_ELEM_BUF_T          fcs_forw_elements_t
#define GRIDSORT_ELEM_INDEX_T        fcs_forw_slint_t
#define GRIDSORT_KEY(_b_, _x_)       (_b_)->keys[_x_]
#define GRIDSORT_DATA0(_b_, _x_)     (&(_b_)->data0[3 * (_x_)])
#define GRIDSORT_DATA1(_b_, _x_)     (_b_)->data1[_x_]
#define GRIDSORT_XYZ(_b_, _d_, _x_)  GRIDSORT_DATA0(_b_, _x_)
/*#ifdef GRIDSORT_FRONT_TPROC_EXDEF
# define DEFINE_GRIDSORT_FRONT_TPROC_EXDEF(_xd_, _tp_, _s_...)       fcs_forw_SPEC_DEFINE_TPROC(_xd_, _tp_, _s_)
# define DEFINE_GRIDSORT_FRONT_TPROCS_MOD_EXDEF(_xd_, _tp_, _s_...)  fcs_forw_SPEC_DEFINE_TPROCS_MOD(_xd_, _tp_, _s_)
#endif*/

#include "gridsort_front_tproc_combinations.h"

#undef DEFINE_GRIDSORT_FRONT_TPROC_EXDEF
#undef DEFINE_GRIDSORT_FRONT_TPROCS_MOD_EXDEF
#undef GRIDSORT_KEY
#undef GRIDSORT_DATA0
#undef GRIDSORT_DATA1
#undef GRIDSORT_XYZ
#undef GRIDSORT_INT_T
#undef GRIDSORT_ELEM_BUF_T
#undef GRIDSORT_ELEM_INDEX_T
#undef GRIDSORT_PREFIX


#endif /* __GRIDSORT_FRONT_TPROC_SL_H__ */
