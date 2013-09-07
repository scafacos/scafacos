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


#ifndef __SL_ENVIRONMENT_H__
#define __SL_ENVIRONMENT_H__


#ifdef FCS_ENABLE_TIMING_SL_FMM
# define SL_TIMING  1
# define SL_TIMING_PRINT_PREFIX  "SL_FMM_TIMING: "
# define SL_TIMING_PRINT(_i_, _s_, _n_, _v_, _r_)  Z_MOP(if ((_r_) == 0) SL_TIMING_PRINT_DEFAULT(_i_, _s_, _n_, _v_, _r_);)
#endif

#ifdef FCS_ENABLE_DEBUG_SL_FMM
# define SLDEBUG_OUTPUT  5
#endif

#endif /* __SL_ENVIRONMENT_H__ */
