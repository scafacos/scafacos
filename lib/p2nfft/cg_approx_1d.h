/*
 * Copyright (C) 2012,2013 Michael Pippig
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

#ifndef _P2NFFT_CG_APPROX_1D_H_
#define _P2NFFT_CG_APPROX_1D_H_

/* !!! Here N is the degree of the cosine polynomial.
 * !!! Therefore, it is only one half of the exponential degree. */

/* Use CG algorithm to calculate coefficients of 1D cosine approximation that is smooth at 1/2,
 * since we want to use it in 3D. Returns the error (maximum error) of the optimal approximation. */
double optimize_cos_coeff_via_cg_guru(
    int N, int M, int M_check, double eps_I, double eps_B, int p, int m, int iter, int verbose,
    double *f_hat_opt);

/* This wrapper is used to hide the parameters of the CG algorithm. */
double optimize_cos_coeff_via_cg(
    int N, double eps_I, double eps_B, int p,
    double *f_hat_opt);

#endif /* cg_approx_1d.h */

