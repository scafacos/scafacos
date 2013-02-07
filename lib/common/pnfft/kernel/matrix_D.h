/*
 * Copyright (c) 2010-2013 Michael Pippig
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

#ifndef __MATRIX_D_H__
#define __MATRIX_D_H__

void PNX(trafo_D)(
    PNX(plan) ths);
void PNX(adjoint_D)(
    PNX(plan) ths);

void PNX(precompute_inv_phi_hat_trafo)(
    PNX(plan) ths,
    C *pre_inv_phi_hat_trafo);
void PNX(precompute_inv_phi_hat_adj)(
    PNX(plan) ths,
    C *pre_inv_phi_hat_adj);

#endif /* __MATRIX_D_H__ */
