/*
 *  Copyright (C) 2011, 2012, 2013, 2014, 2015, 2016 Michael Hofmann
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


#ifndef __SL_CONFIG_GLOBAL_H__
#define __SL_CONFIG_GLOBAL_H__


#if defined(QWERTZ_MSEG_ROOT) && !defined(MSEG_ROOT)
# define MSEG_ROOT  QWERTZ_MSEG_ROOT
#endif

#if defined(QWERTZ_MSEG_BORDER_UPDATE_REDUCTION) && !defined(MSEG_BORDER_UPDATE_REDUCTION)
# define MSEG_BORDER_UPDATE_REDUCTION  QWERTZ_MSEG_BORDER_UPDATE_REDUCTION
#endif

#if defined(QWERTZ_MSEG_DISABLE_BEST_CHOICE) && !defined(MSEG_DISABLE_BEST_CHOICE)
# define MSEG_DISABLE_BEST_CHOICE  QWERTZ_MSEG_DISABLE_BEST_CHOICE
#endif

#if defined(QWERTZ_MSEG_DISABLE_MINMAX) && !defined(MSEG_DISABLE_MINMAX)
# define MSEG_DISABLE_MINMAX  QWERTZ_MSEG_DISABLE_MINMAX
#endif

#if defined(QWERTZ_MSEG_ENABLE_OPTIMZED_LOWHIGH) && !defined(MSEG_ENABLE_OPTIMZED_LOWHIGH)
# define MSEG_ENABLE_OPTIMZED_LOWHIGH  QWERTZ_MSEG_ENABLE_OPTIMZED_LOWHIGH
#endif

#if defined(QWERTZ_MSEG_FORWARD_ONLY) && !defined(MSEG_FORWARD_ONLY)
# define MSEG_FORWARD_ONLY  QWERTZ_MSEG_FORWARD_ONLY
#endif

#if defined(QWERTZ_MSEG_INFO) && !defined(MSEG_INFO)
# define MSEG_INFO  QWERTZ_MSEG_INFO
#endif

#if defined(QWERTZ_MSEG_TRACE_IF) && !defined(MSEG_TRACE_IF)
# define MSEG_TRACE_IF  QWERTZ_MSEG_TRACE_IF
#endif


#endif /* __SL_CONFIG_GLOBAL_H__ */
