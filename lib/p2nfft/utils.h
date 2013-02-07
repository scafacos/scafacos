/*
 * Copyright (C) 2011-2013 Michael Pippig
 * Copyright (C) 2011 Sebastian Banert
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

#ifndef _P2NFFT_UTILS_H
#define _P2NFFT_UTILS_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "FCSCommon.h"
#include "constants.h"

static inline fcs_float sinc(fcs_float d)
{
  const fcs_float epsi = 0.1;
  const fcs_float c2 = -0.1666666666667e-0;
  const fcs_float c4 = 0.8333333333333e-2;
  const fcs_float c6 = -0.1984126984127e-3;
  const fcs_float c8 = 0.2755731922399e-5;

  fcs_float PId = FCS_P2NFFT_PI*d, PId2;

  if (fcs_fabs(d)>epsi)
    return fcs_sin(PId)/PId;
  else {
    PId2 = FCS_P2NFFT_SQR(PId);
    return 1.0 + PId2*(c2+PId2*(c4+PId2*(c6+PId2*c8)));
  }
}

#endif
