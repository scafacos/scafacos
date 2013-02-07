/*
 * Copyright (c) 2012-2013 Michael Pippig
 *
 * This file is part of PNFFT.
 *
 * PNFFT is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * PNFFT is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PNFFT.  If not, see <http://www.gnu.org/licenses/>.
 *
 */


#include "ipnfft.h"

#ifndef __BSPLINE_H__
#define __BSPLINE_H__

R PNX(bspline)(
    int k, R x, R *scratch);
R PNX(fast_bspline)(
    int i, R x, int cao_value);
R PNX(fast_bspline_d)(
    int i, R x, int cao_value);


#endif
