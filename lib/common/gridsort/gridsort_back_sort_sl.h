/*
  Copyright (C) 2011, 2012, 2013, 2014, 2015 Michael Hofmann
  
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


#ifndef __GRIDSORT_BACK_SORT_SL_H__
#define __GRIDSORT_BACK_SORT_SL_H__


#define GRIDSORT_PREFIX

#define GRIDSORT_BACK_PREFIX      gridsort_back_

/* back_fp */
#define GRIDSORT_BACK_TPROC_NAME  fp_tproc
#define GRIDSORT_BACK_TLOC_NAME   fp_tloc
#define GRIDSORT_BACK_SORT_NAME   fp_sort
#define GRIDSORT_BACK_SL_PREFIX   fcs_back_fp_
#include "gridsort_back_sort.h"

#undef GRIDSORT_BACK_SL_PREFIX
#undef GRIDSORT_BACK_TPROC_NAME
#undef GRIDSORT_BACK_TLOC_NAME
#undef GRIDSORT_BACK_SORT_NAME

/* back_f_ */
#define GRIDSORT_BACK_TPROC_NAME  f__tproc
#define GRIDSORT_BACK_TLOC_NAME   f__tloc
#define GRIDSORT_BACK_SORT_NAME   f__sort
#define GRIDSORT_BACK_SL_PREFIX   fcs_back_f__
#include "gridsort_back_sort.h"

#undef GRIDSORT_BACK_SL_PREFIX
#undef GRIDSORT_BACK_TPROC_NAME
#undef GRIDSORT_BACK_TLOC_NAME
#undef GRIDSORT_BACK_SORT_NAME

/* back__p_ */
#define GRIDSORT_BACK_TPROC_NAME  _p_tproc
#define GRIDSORT_BACK_TLOC_NAME   _p_tloc
#define GRIDSORT_BACK_SORT_NAME   _p_sort
#define GRIDSORT_BACK_SL_PREFIX   fcs_back__p_
#include "gridsort_back_sort.h"

#undef GRIDSORT_BACK_SL_PREFIX
#undef GRIDSORT_BACK_TPROC_NAME
#undef GRIDSORT_BACK_TLOC_NAME
#undef GRIDSORT_BACK_SORT_NAME

#undef GRIDSORT_BACK_PREFIX


#if FCS_GRIDSORT_WITH_DIPOLES

#define GRIDSORT_BACK_PREFIX      gridsort_back_dipole_

/* dip_back_fp */
#define GRIDSORT_BACK_TPROC_NAME  fp_tproc
#define GRIDSORT_BACK_TLOC_NAME   fp_tloc
#define GRIDSORT_BACK_SORT_NAME   fp_sort
#define GRIDSORT_BACK_SL_PREFIX   fcs_dip_back_fp_
#include "gridsort_back_sort.h"

#undef GRIDSORT_BACK_SL_PREFIX
#undef GRIDSORT_BACK_TPROC_NAME
#undef GRIDSORT_BACK_TLOC_NAME
#undef GRIDSORT_BACK_SORT_NAME

/* dip_back_f_ */
#define GRIDSORT_BACK_TPROC_NAME  f__tproc
#define GRIDSORT_BACK_TLOC_NAME   f__tloc
#define GRIDSORT_BACK_SORT_NAME   f__sort
#define GRIDSORT_BACK_SL_PREFIX   fcs_dip_back_f__
#include "gridsort_back_sort.h"

#undef GRIDSORT_BACK_SL_PREFIX
#undef GRIDSORT_BACK_TPROC_NAME
#undef GRIDSORT_BACK_TLOC_NAME
#undef GRIDSORT_BACK_SORT_NAME

/* dip_back__p_ */
#define GRIDSORT_BACK_TPROC_NAME  _p_tproc
#define GRIDSORT_BACK_TLOC_NAME   _p_tloc
#define GRIDSORT_BACK_SORT_NAME   _p_sort
#define GRIDSORT_BACK_SL_PREFIX   fcs_dip_back__p_
#include "gridsort_back_sort.h"

#undef GRIDSORT_BACK_SL_PREFIX
#undef GRIDSORT_BACK_TPROC_NAME
#undef GRIDSORT_BACK_TLOC_NAME
#undef GRIDSORT_BACK_SORT_NAME

#undef GRIDSORT_BACK_PREFIX

#endif /* FCS_GRIDSORT_WITH_DIPOLES */


#undef GRIDSORT_PREFIX


#endif /* __GRIDSORT_BACK_SORT_SL_H__ */
