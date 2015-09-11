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

/* GRIDSORT_BACK_TPROC_NAME */


#define CONCAT_(_a_, _b_)  _a_##_b_
#define CONCAT(_a_, _b_)   CONCAT_(_a_, _b_)


static int CONCAT(GRIDSORT_PREFIX, CONCAT(GRIDSORT_BACK_PREFIX, GRIDSORT_BACK_TPROC_NAME))(GRIDSORT_ELEM_BUF_T *s, GRIDSORT_ELEM_INDEX_T x, void *data)
{
#ifdef GRIDSORT_INDEX_PTR
  fcs_gridsort_index_t *index_ptr = data;
#endif

  if (!GRIDSORT_INDEX_IS_VALID(GRIDSORT_INDEX(s, index_ptr, x))) return MPI_PROC_NULL;

  return GRIDSORT_INDEX_GET_PROC(GRIDSORT_INDEX(s, index_ptr, x));
}

#ifdef DEFINE_GRIDSORT_BACK_TPROC_EXDEF
DEFINE_GRIDSORT_BACK_TPROC_EXDEF(CONCAT(CONCAT(GRIDSORT_PREFIX, CONCAT(GRIDSORT_BACK_PREFIX, GRIDSORT_BACK_TPROC_NAME)), _exdef), CONCAT(GRIDSORT_PREFIX, CONCAT(GRIDSORT_BACK_PREFIX, GRIDSORT_BACK_TPROC_NAME)), static);
#endif


#undef CONCAT_
#undef CONCAT

#undef GRIDSORT_BACK_TPROC_NAME
