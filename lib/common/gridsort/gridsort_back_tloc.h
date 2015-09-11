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

/* GRIDSORT_BACK_TLOC_NAME */


#define CONCAT_(_a_, _b_)  _a_##_b_
#define CONCAT(_a_, _b_)   CONCAT_(_a_, _b_)


static GRIDSORT_ELEM_INDEX_T CONCAT(GRIDSORT_PREFIX, CONCAT(GRIDSORT_BACK_PREFIX, GRIDSORT_BACK_TLOC_NAME))(GRIDSORT_ELEM_BUF_T *b, GRIDSORT_ELEM_INDEX_T x, void *tloc_data)
{
  GRIDSORT_ELEM_INDEX_T r = GRIDSORT_INDEX_GET_POS(GRIDSORT_INDEX(b, tloc_data, x));

#if !SORT_BACKWARD_LOCAL_INPLACE
  GRIDSORT_INDEX(b, tloc_data, x) = GRIDSORT_INDEX(b, tloc_data, r);
#endif

  return r;
}


#undef CONCAT_
#undef CONCAT

#undef GRIDSORT_BACK_TLOC_NAME
