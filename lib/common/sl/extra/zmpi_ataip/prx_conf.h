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


#ifndef __PRX_CONF_H__
#define __PRX_CONF_H__


#include "dash_conf.h"


#define PRX_RENAME

#ifdef ZMPI_PREFIX
# define PRX_PREFIX  PRX_CONCAT(ZMPI_PREFIX, zmpi_)
#else
# define PRX_PREFIX  zmpi_
#endif


typedef dsint_t prxint_t;
#define prxint_fmt  dsint_fmt


#endif /* __PRX_CONF_H__ */
