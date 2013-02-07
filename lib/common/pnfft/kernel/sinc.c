/*
 * Copyright (c) 2002, 2009 Jens Keiner, Stefan Kunis, Daniel Potts
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

#include "pnfft.h"
#include "ipnfft.h"
#include "sinc.h"

R PNX(sinc)(
    const R x
    )
{
  /* Based on sinc function from Boost C++ library. */
  const R b =  PNFFT_EPSILON;
  const R bs = pnfft_sqrt(b);
  const R bs2 = pnfft_sqrt(bs);

  if (pnfft_fabs(x) >= bs2)
    return pnfft_sin(x)/x;
  else
  {
    R r = K(1.0);

    if (pnfft_fabs(x) >= b)
    {
      const R x2 = x * x;
      r -= x2 / K(6.0);

      if (pnfft_fabs(x) >= bs)
        r += (x2 * x2) / K(120.0);
    }

    return r;
  }
}

