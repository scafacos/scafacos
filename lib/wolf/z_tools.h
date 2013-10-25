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


#ifndef __Z_TOOLS_H__
#define __Z_TOOLS_H__


#define Z_MOP(_mop_)  do { _mop_ } while (0)
#define Z_NOP()       Z_MOP()

#define z_max(_a_, _b_)           (((_a_)>(_b_))?(_a_):(_b_))
#define z_min(_a_, _b_)           (((_a_)<(_b_))?(_a_):(_b_))
#define z_max3(_a_, _b_, _c_)     z_max(_a_,z_max(_b_,_c_))
#define z_min3(_a_, _b_, _c_)     z_min(_a_,z_min(_b_,_c_))
#define z_minmax(_a_, _b_, _c_)   (((_b_)<(_a_))?(_a_):(((_b_)>(_c_))?(_c_):(_b_)))
#define z_abs(_a_)                (((_a_) >= 0)?(_a_):-(_a_))
#define z_swap(_a_, _b_, _t_)     Z_MOP((_t_) = (_a_); (_a_) = (_b_); (_b_) = (_t_);)
#define z_sqr(_a_)                ((_a_) * (_a_))

#define z_is_zero(_x_)            ((_x_) == 0)
#define z_is_not_zero(_x_)        ((_x_) != 0)

#define z_fp_is_zero(_x_)         ((_x_) <= 0 && (_x_) >= 0)
#define z_fp_is_not_zero(_x_)     ((_x_) < 0 || (_x_) > 0)

#define z_is_triclinic(_a_, _b_, _c_)  (z_fp_is_not_zero((_a_)[1]) || z_fp_is_not_zero((_a_)[2]) || z_fp_is_not_zero((_b_)[0]) || z_fp_is_not_zero((_b_)[2]) || z_fp_is_not_zero((_c_)[0]) || z_fp_is_not_zero((_c_)[1]))


#endif /* __Z_TOOLS_H__ */
