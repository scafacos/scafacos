/*
  Copyright (C) 2011
  
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

/** \file specfunc.h
    This file contains implementations for some special functions which are needed by the MMM family of
    algorithms. This are the modified Hurwitz zeta function and the modified Bessel functions of first
    and second kind. The implementations are based on the GSL code (see \ref specfunc.c "specfunc.c"
    for the original GSL header).

    The Hurwitz zeta function is evaluated using the Euler-MacLaurin summation formula, the Bessel functions
    are evaluated using several different Chebychev expansions. Both achieve a precision of nearly machine
    precision, which is no problem for the Hurwitz zeta function, which is only used when determining the
    coefficients for the modified polygamma functions (see \ref mmm-common.h "mmm-common.h"). However, the
    Bessel functions are actually used in the near formula of MMM2D, which is therefore slightly slower than
    necessary. On the other hand, the number of terms in the Bessel sum is quite small normally, so that a less
    precise version will probably not generate a huge computational speed improvement.
*/
#ifndef SPECFUNC_H
#define SPECFUNC_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "FCSCommon.h"
#include "mmm-common.h"

/** Hurwitz zeta function. This function was taken from the GSL code. */
fcs_float mmm_hzeta(fcs_float order, fcs_float x);

/** Besselfunctions K0 and K1 at x.
    The implementation has an absolute precision of around 10^(-14), which is
    comparable to the relative precision sqrt implementation of current hardware.
*/
void mmm_LPK01(fcs_float x, fcs_float *K0, fcs_float *K1);

/** Modified Bessel function of second kind, order 0. This function was taken from
    the GSL code. Precise roughly up to machine precision. */
fcs_float mmm_K0(fcs_float x);

/** Modified Bessel function of second kind, order 1. This function was taken from
    the GSL code. Precise roughly up to machine precision. */
fcs_float mmm_K1(fcs_float x);

/** evaluate the polynomial interpreted as a Chebychev series. Requires a series with at least
    three coefficients, i.e. no linear approximations! */
fcs_float mmm_evaluateAsChebychevSeriesAt(SizedList *series, fcs_float x);

/** Hurwitz zeta function. This function was taken from the GSL code. */
fcs_float mmm_hzeta(fcs_float order, fcs_float x);

/** evaluate the polynomial interpreted as a Taylor series via the Horner scheme */
fcs_float mmm_evaluateAsTaylorSeriesAt(SizedList *series, fcs_float x);

/** modified polygamma for even order 2*n, n >= 0 */
fcs_float mmm_mod_psi_even(mmm_data_struct *polTaylor, fcs_int n, fcs_float x);

/** modified polygamma for odd order 2*n+1, n>= 0 */
fcs_float mmm_mod_psi_odd(mmm_data_struct *polTaylor, fcs_int n, fcs_float x);

/** create the both the even and odd polygamma functions up to order 2*n */
void mmm_create_mod_psi_up_to(mmm_data_struct *polTaylor, fcs_int new_n);
#endif
