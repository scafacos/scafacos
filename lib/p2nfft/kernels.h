/*
 * Copyright (C) 2011-2013 Michael Pippig
 * Copyright (c) 2002, 2009 Jens Keiner, Stefan Kunis, Daniel Potts
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

/*! \file kernels.h
 *  \brief Header file with predefined kernels for the fast summation algorithm.
 */
#ifndef _P2NFFT_KERNELS_H
#define _P2NFFT_KERNELS_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "FCSCommon.h"
#include "FCSDefinitions.h"

//#ifdef HAVE_COMPLEX_H
#include <complex.h>
//#endif

/* At the moment we only need real valued kernel function.
 * Switch back to complex kernel if neccessary */
// #define FCS_P2NFFT_KERNEL_TYPE fcs_float _Complex
#define FCS_P2NFFT_KERNEL_TYPE fcs_float

typedef FCS_P2NFFT_KERNEL_TYPE (*ifcs_p2nfft_kernel)(fcs_float , fcs_int , const fcs_float *);

FCS_P2NFFT_KERNEL_TYPE ifcs_p2nfft_gaussian(fcs_float x, fcs_int der, const fcs_float *param);              /* K(x)=exp(-x^2/c^2) */
FCS_P2NFFT_KERNEL_TYPE ifcs_p2nfft_multiquadric(fcs_float x, fcs_int der, const fcs_float *param);          /* K(x)=sqrt(x^2+c^2) */
FCS_P2NFFT_KERNEL_TYPE ifcs_p2nfft_inverse_multiquadric(fcs_float x, fcs_int der, const fcs_float *param);  /* K(x)=1/sqrt(x^2+c^2) */
FCS_P2NFFT_KERNEL_TYPE ifcs_p2nfft_logarithm(fcs_float x, fcs_int der, const fcs_float *param);             /* K(x)=log |x| */
FCS_P2NFFT_KERNEL_TYPE ifcs_p2nfft_thinplate_spline(fcs_float x, fcs_int der, const fcs_float *param);      /* K(x) = x^2 log |x| */
FCS_P2NFFT_KERNEL_TYPE ifcs_p2nfft_one_over_square(fcs_float x, fcs_int der, const fcs_float *param);       /* K(x) = 1/x^2 */
FCS_P2NFFT_KERNEL_TYPE ifcs_p2nfft_one_over_modulus(fcs_float x, fcs_int der, const fcs_float *param);      /* K(x) = 1/|x| */
FCS_P2NFFT_KERNEL_TYPE ifcs_p2nfft_one_over_x(fcs_float x, fcs_int der, const fcs_float *param);            /* K(x) = 1/x */
FCS_P2NFFT_KERNEL_TYPE ifcs_p2nfft_inverse_multiquadric3(fcs_float x, fcs_int der, const fcs_float *param); /* K(x) = 1/sqrt(x^2+c^2)^3 */
FCS_P2NFFT_KERNEL_TYPE ifcs_p2nfft_sinc_kernel(fcs_float x, fcs_int der, const fcs_float *param);           /* K(x) = sin(cx)/x */
FCS_P2NFFT_KERNEL_TYPE ifcs_p2nfft_cosc(fcs_float x, fcs_int der, const fcs_float *param);                  /* K(x) = cos(cx)/x */
FCS_P2NFFT_KERNEL_TYPE ifcs_p2nfft_kcot(fcs_float x, fcs_int der, const fcs_float *param);                  /* K(x) = cot(cx) */
FCS_P2NFFT_KERNEL_TYPE ifcs_p2nfft_one_over_cube(fcs_float x, fcs_int der, const fcs_float *param);         /* K(x) = 1/x^3 */

FCS_P2NFFT_KERNEL_TYPE ifcs_p2nfft_ewald_2dp_kneq0(fcs_float x, fcs_int der, const fcs_float *param);       /* K(x) = exp(2*pi*k*x) * erf(pi*k/alpha + alpha*x) */
FCS_P2NFFT_KERNEL_TYPE ifcs_p2nfft_ewald_2dp_keq0(fcs_float x, fcs_int der, const fcs_float *param);        /* K(x) = 1/alpha * exp(alpha*x) + sqrt(pi)*x*erf(alpha*x) */
FCS_P2NFFT_KERNEL_TYPE ifcs_p2nfft_x_times_erf(fcs_float x, fcs_int der, const fcs_float *param);           /* K(x) = sqrt(PI) * x * erf(alpha*x) */

FCS_P2NFFT_KERNEL_TYPE ifcs_p2nfft_ewald_1dp_kneq0(fcs_float r, fcs_int der, const fcs_float *param);       /* K(x) = K_nu(pi^2*k^2/alpha^2, alpha^2*r^2) */
FCS_P2NFFT_KERNEL_TYPE ifcs_p2nfft_ewald_1dp_keq0(fcs_float r, fcs_int der, const fcs_float *param);        /* K(x) = gamma + Gamma(0,a^2*r^2) + ln(a^2*r^2) */


#endif
