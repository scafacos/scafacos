/*
 * Copyright (C) 2013 Michael Pippig
 *
 * This file is part of ScaFaCoS.
 * 
 * ScaFaCoS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ScaFaCoS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *	
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _P2NFFT_PART_DERIVE_ONE_OVER_NORM_X_H_
#define _P2NFFT_PART_DERIVE_ONE_OVER_NORM_X_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "FCSCommon.h"
#include "FCSDefinitions.h"

fcs_float ifcs_p2nfft_part_derive_one_over_norm_x(
    fcs_int i, fcs_int j, fcs_int k,
    fcs_float x, fcs_float y, fcs_float z);

#endif
