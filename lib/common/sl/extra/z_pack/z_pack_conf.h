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


#ifndef __Z_PACK_CONF_H__
#define __Z_PACK_CONF_H__


typedef long z_int_t;
#define z_int_fmt  "ld"


#ifdef Z_PREFIX
# define Z_PACK_RENAME
#endif


#define Z_PACK_NUMERIC


#define Z_PACK_ALLOC


#define Z_PACK_RANDOM


#define Z_PACK_DEBUG_OFF


#endif /* __Z_PACK_CONF_H__ */
