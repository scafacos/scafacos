/*
  Copyright (C) 2011,2012 Olaf Lenz

  This file is part of ScaFaCoS.

  ScaFaCoS is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  ScaFaCoS is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>. 
 */
#include "CAF.hpp"

#include <cstdlib>
#include <cassert>
#include <sstream>
#include <stdexcept>
#include <limits>

#include "utils.hpp"

namespace P3M {
/** Factory method to generate the right CAF Cache. */
CAF* CAF::create(p3m_int cao,
		p3m_int n_interpol,
		bool derivative) {
	if (n_interpol > 0)
		return new InterpolatedCAF(cao, n_interpol, derivative);
	else
		return new CAF(cao, derivative);
}



/** Compute the P-th order charge assignment function from
 * Hockney/Eastwood 5-189 (or 8-61). The following charge fractions
 * are also tabulated in Deserno/Holm. */
p3m_float CAF::
compute(p3m_int i, p3m_float x, p3m_int cao) {
    switch (cao) {
    case 1: {
        return 1.0;
        break;
    }
    case 2: {
        switch (i) {
        case 0: return 0.5-x;
        case 1: return 0.5+x;
        }
        break;
    }
    case 3: {
        switch (i) {
        case 0: return 0.5*SQR(0.5 - x);
        case 1: return 0.75 - SQR(x);
        case 2: return 0.5*SQR(0.5 + x);
        }
        break;
    }
    case 4: {
        switch (i) {
        case 0: return ( 1.0+x*( -6.0+x*( 12.0-x* 8.0)))/48.0;
        case 1: return (23.0+x*(-30.0+x*(-12.0+x*24.0)))/48.0;
        case 2: return (23.0+x*( 30.0+x*(-12.0-x*24.0)))/48.0;
        case 3: return ( 1.0+x*(  6.0+x*( 12.0+x* 8.0)))/48.0;
        }
        break;
    }
    case 5: {
        switch (i) {
        case 0: return (  1.0+x*( -8.0+x*(  24.0+x*(-32.0+x*16.0))))/384.0;
        case 1: return ( 19.0+x*(-44.0+x*(  24.0+x*( 16.0-x*16.0))))/ 96.0;
        case 2: return (115.0+x*       x*(-120.0+x*       x*48.0))  /192.0;
        case 3: return ( 19.0+x*( 44.0+x*(  24.0+x*(-16.0-x*16.0))))/ 96.0;
        case 4: return (  1.0+x*(  8.0+x*(  24.0+x*( 32.0+x*16.0))))/384.0;
        }
        break;
    }
    case 6: {
        switch (i) {
        case 0: return (  1.0+x*( -10.0+x*(  40.0+x*( -80.0+x*(  80.0-x* 32.0)))))/3840.0;
        case 1: return (237.0+x*(-750.0+x*( 840.0+x*(-240.0+x*(-240.0+x*160.0)))))/3840.0;
        case 2: return (841.0+x*(-770.0+x*(-440.0+x*( 560.0+x*(  80.0-x*160.0)))))/1920.0;
        case 3: return (841.0+x*(+770.0+x*(-440.0+x*(-560.0+x*(  80.0+x*160.0)))))/1920.0;
        case 4: return (237.0+x*( 750.0+x*( 840.0+x*( 240.0+x*(-240.0-x*160.0)))))/3840.0;
        case 5: return (  1.0+x*(  10.0+x*(  40.0+x*(  80.0+x*(  80.0+x* 32.0)))))/3840.0;
        }
        break;
    }
    case 7: {
        switch (i) {
        case 0: return (    1.0+x*(   -12.0+x*(   60.0+x*( -160.0+x*(  240.0+x*(-192.0+x* 64.0))))))/46080.0;
        case 1: return (  361.0+x*( -1416.0+x*( 2220.0+x*(-1600.0+x*(  240.0+x*( 384.0-x*192.0))))))/23040.0;
        case 2: return (10543.0+x*(-17340.0+x*( 4740.0+x*( 6880.0+x*(-4080.0+x*(-960.0+x*960.0))))))/46080.0;
        case 3: return ( 5887.0+x*          x*(-4620.0+x*         x*( 1680.0-x*        x*320.0)))   /11520.0;
        case 4: return (10543.0+x*( 17340.0+x*( 4740.0+x*(-6880.0+x*(-4080.0+x*( 960.0+x*960.0))))))/46080.0;
        case 5: return (  361.0+x*(  1416.0+x*( 2220.0+x*( 1600.0+x*(  240.0+x*(-384.0-x*192.0))))))/23040.0;
        case 6: return (    1.0+x*(    12.0+x*(   60.0+x*(  160.0+x*(  240.0+x*( 192.0+x* 64.0))))))/46080.0;
        }
        break;
    }
    default: {
        std::ostringstream s;
        s << "Charge assignment order " << cao << " unknown.";
        throw std::logic_error(s.str());
    }
	}
	std::ostringstream s;
	s << "Tried to compute charge assignment function at "
			<< i << "'th point in scheme of order " << cao << ".";
	throw std::logic_error(s.str());
}

/** Computes the gradient of the charge assignment function of for
      the \a i'th degree at value \a x. */
p3m_float CAF::
computeDerivative(p3m_int i, p3m_float x, p3m_int cao) {
	p3m_int ip = cao - 1;
	switch (ip) {
	case 1: {
	    switch(i) {
	    case 0: return -1.0;
	    case 1: return 1.0;
	    }
	    break;
	}
	case 2: {
	    switch(i) {
	    case 0: return x-0.5;
	    case 1: return -2.0*x;
	    case 2: return x+0.5;
	    }
	    break;
	}
	case 3: {
	    switch(i) {
	    case 0: return (-1.0+x*(  4.0+x*( -4.0)))/8.0;
	    case 1: return (-5.0+x*( -4.0+x*( 12.0)))/8.0;
	    case 2: return ( 5.0+x*( -4.0+x*(-12.0)))/8.0;
	    case 3: return ( 1.0+x*(  4.0+x*(  4.0)))/8.0;
	    }
	    break;
	}
	case 4: {
	    switch(i) {
	    case 0: return ( -1.0+x*(  6.0+x*( -12.0+x*( 8.0))))/48.0;
	    case 1: return (-11.0+x*( 12.0+x*(  12.0+x*(-16.0))))/24.0;
	    case 2: return (      x*(-5.0+x*x*4.0))/4.0;
	    case 3: return ( 11.0+x*( 12.0+x*( -12.0+x*(-16.0))))/ 24.0;
	    case 4: return (  1.0+x*(  6.0+x*(  12.0+x*(  8.0))))/48.0;
	    }
	    break;
	}
	case 5: {
	    switch(i) {
	    case 0: return ( -1.0+x*(  8.0+x*( -24.0+x*(  32.0+x*(-16)))))/384.0;
	    case 1: return (-75.0+x*( 168.0+x*( -72.0+x*( -96.0+x*(80.0)))))/384.0;
	    case 2: return (-77.0+x*( -88.0+x*( 168.0+x*(  32.0+x*(-80.0)))))/192.0;
	    case 3: return ( 77.0+x*( -88.0+x*(-168.0+x*(  32.0+x*(80.0)))))/192.0;
	    case 4: return ( 75.0+x*( 168.0+x*(  72.0+x*( -96.0+x*(-80)))))/384.0;
	    case 5: return (  1.0+x*(   8.0+x*(  24.0+x*(  32.0+x*(16.0)))))/384.0;
	    }
	    break;
	}
	case 6: {
	    switch(i) {
	    case 0: return (  -1.0+x*( 10.0+x*( -40.0+x*(  80.0+x*(-80.0+x*32.0)))))/3840.0;
	    case 1: return ( -59.0+x*(185.0+x*(-200.0+x*(  40.0+x*( 80.0-x*48.0)))))/960.0;
	    case 2: return (-289.0+x*(158.0+x*( 344.0+x*(-272.0+x*(-80.0+x*96.0)))))/768.0;
	    case 3: return (       x*(-77.0+        x*x*(  56.0       -x*x*16.0) ) )/96.0;
	    case 4: return ( 289.0+x*(158.0+x*(-344.0+x*(-272.0+x*( 80.0+x*96.0)))))/768.0;
	    case 5: return (  59.0+x*(185.0+x*( 200.0+x*(  40.0+x*(-80.0-x*48.0)))))/960.0;
	    case 6: return (   1.0+x*( 10.0+x*(  40.0+x*(  80.0+x*( 80.0+x*32.0)))))/3840.0;
	    }
	    break;
	}
	default: {
	    std::ostringstream s;
	    s << "Charge assignment order " << cao << " unknown.";
	    throw std::logic_error(s.str());
	}
	}
	std::ostringstream s;
	s << "Tried to compute charge assignment function at "
			<< i << "'th point in scheme of order " << cao << ".";
	throw std::logic_error(s.str());
}

/* Cache for exact values */
CAF::DirectCache::DirectCache(CAF &_caf) : caf(_caf) {
	pbegin = new p3m_float[caf.cao];
	pend = pbegin + caf.cao;
}

CAF::DirectCache::~DirectCache() { delete[] pbegin; }

/* Fill the cache with the computed values. */
void CAF::DirectCache::update(p3m_float r0) {
#ifdef ADDITIONAL_CHECKS
	assert(r0 <= 0.5);
	assert(r0 >= -0.5);
#endif
	for (p3m_int i = 0; i < caf.cao; ++i)
		pbegin[i] = compute(i, r0, caf.cao);
}

/* Cache for exact values of the derivative */
CAF::DerivativeCache::DerivativeCache(CAF &_caf) : caf(_caf) {
	pbegin = new p3m_float[caf.cao];
	pend = pbegin + caf.cao;
}

CAF::DerivativeCache::~DerivativeCache() { delete[] pbegin; }

/* Fill the cache with the computed values. */
void CAF::DerivativeCache::update(p3m_float r0) {
#ifdef ADDITIONAL_CHECKS
	assert(r0 <= 0.5);
	assert(r0 >= -0.5);
#endif
	for (p3m_int i = 0; i < caf.cao; ++i)
		pbegin[i] = computeDerivative(i, r0, caf.cao);
}

CAF::CAF(p3m_int _cao, bool _derivative) : cao(_cao), derivative(_derivative) {
	if (cao > max_cao || cao < 1) {
		std::ostringstream s;
		s << "Charge assignment order " << cao << " unknown.";
		throw std::logic_error(s.str());
	}
}

CAF::Cache *CAF::createCache() {
	if (derivative)
		return new DerivativeCache(*this);
	else
		return new DirectCache(*this);
}

/* Cache for interpolated values. */
void InterpolatedCAF::Cache::update(p3m_float r0) {
	// get interpolated cache
	pbegin = caf.getBase(r0);
	pend = pbegin + caf.cao;
}

/* Interpolate the CAF. */
InterpolatedCAF::InterpolatedCAF(p3m_int cao, p3m_int _n_interpol,
		bool derivative)
: CAF(cao, derivative) {
	if (_n_interpol <= 0) {
		std::ostringstream s;
		s << "Bad interpolation order " << _n_interpol << ".";
		throw std::logic_error(s.str());
	}

	// interpolate CAF
	n_interpol = _n_interpol;
	data = new p3m_float[(2*n_interpol+1) * cao];

	const p3m_float dInterpol = 0.5 / n_interpol;

	P3M_DEBUG_LOCAL(printf("    InterpolatedCAF: interpolating caf of "       \
			"order %d with %d points\n",                       \
			cao,n_interpol));
	/* loop over all interpolation points */
	for (p3m_int i = -n_interpol; i <= n_interpol; i++)
		/* loop over all cao */
		for (p3m_int j = 0; j < cao; j++)
			data[(i+n_interpol)*cao + j] =
					(derivative ?
							P3M::CAF::computeDerivative(j, i*dInterpol, cao) :
							P3M::CAF::compute(j, i*dInterpol, cao));
}

InterpolatedCAF::~InterpolatedCAF() { delete[] data; }

p3m_float* InterpolatedCAF::getBase(p3m_float r0) {
#ifdef ADDITIONAL_CHECKS
	assert(r0 <= 0.5);
	assert(r0 >= -0.5);
#endif
	const p3m_int ix =
			static_cast<p3m_int>((2*n_interpol+1) * (r0 + 0.5)) * cao;
#ifdef ADDITIONAL_CHECKS
	assert(ix >= 0);
	assert(ix <= (2*n_interpol+1) * cao);
#endif
	return &data[ix];
}

CAF::Cache *InterpolatedCAF::createCache() { return new Cache(*this); }
}
