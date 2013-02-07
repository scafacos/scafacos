/*
 * Copyright (c) 2011-2013 Michael Pippig
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

#include <complex.h>
#include "pnfft.h"
#include "ipnfft.h"


void PNX(spread_f_c2c)(
    PNX(plan) ths, INT ind,
    C f, R *pre_psi, INT m0, INT *grid_size, int cutoff, int interlaced,
    C *grid
    )
{
  R* plan_pre_psi = (interlaced) ? ths->pre_psi_il : ths->pre_psi;

  if(ths->pnfft_flags & PNFFT_PRE_FULL_PSI)
    PNX(spread_f_c2c_pre_full_psi)(
        f, plan_pre_psi + ind*PNFFT_POW3(cutoff), m0, grid_size, cutoff, 
        grid);
  else if(ths->pnfft_flags & PNFFT_PRE_PSI)
    PNX(spread_f_c2c_pre_psi)(
        f, plan_pre_psi + ind*3*cutoff, m0, grid_size, cutoff, 
        grid);
  else
    PNX(spread_f_c2c_pre_psi)(
        f, pre_psi, m0, grid_size, cutoff, 
        grid);
}

void PNX(spread_f_r2r)(
    PNX(plan) ths, INT ind,
    R f, R *pre_psi, INT m0, INT *grid_size, int cutoff, INT ostride, int interlaced,
    R *grid
    )
{
  R* plan_pre_psi = (interlaced) ? ths->pre_psi_il : ths->pre_psi;

  if(ths->pnfft_flags & PNFFT_PRE_FULL_PSI)
    PNX(spread_f_r2r_pre_full_psi)(
        f, plan_pre_psi + ind*PNFFT_POW3(cutoff), m0, grid_size, cutoff, ostride,
        grid);
  else if(ths->pnfft_flags & PNFFT_PRE_PSI)
    PNX(spread_f_r2r_pre_psi)(
        f, plan_pre_psi + ind*3*cutoff, m0, grid_size, cutoff, ostride,
        grid);
  else
    PNX(spread_f_r2r_pre_psi)(
        f, pre_psi, m0, grid_size, cutoff, ostride,
        grid);
}

void PNX(assign_f_c2c)(
    PNX(plan) ths, INT ind,
    C *grid, R *pre_psi, INT m0, INT *grid_size, int cutoff, int interlaced,
    C *f
    )
{
  R* plan_pre_psi = (interlaced) ? ths->pre_psi_il : ths->pre_psi;

  if(ths->pnfft_flags & PNFFT_PRE_FULL_PSI)
    PNX(assign_f_c2c_pre_full_psi)(
        ths->g2, plan_pre_psi + ind*PNFFT_POW3(cutoff), m0, grid_size, cutoff,
        f);
  else if(ths->pnfft_flags & PNFFT_PRE_PSI)
    PNX(assign_f_c2c_pre_psi)(
        grid, plan_pre_psi + ind*3*cutoff, m0, grid_size, cutoff,
        f);
  else
    PNX(assign_f_c2c_pre_psi)(
        grid, pre_psi, m0, grid_size, cutoff,
        f);
}

void PNX(assign_f_r2r)(
    PNX(plan) ths, INT ind,
    R *grid, R *pre_psi, INT m0, INT *grid_size, int cutoff, INT istride, int interlaced,
    R *f
    )
{ 
  R* plan_pre_psi = (interlaced) ? ths->pre_psi_il : ths->pre_psi;

  if(ths->pnfft_flags & PNFFT_PRE_FULL_PSI)
    PNX(assign_f_r2r_pre_full_psi)(
        grid, plan_pre_psi + ind*PNFFT_POW3(cutoff), m0, grid_size, cutoff, istride,
        f);
  else if(ths->pnfft_flags & PNFFT_PRE_PSI)
    PNX(assign_f_r2r_pre_psi)(
        grid, plan_pre_psi + ind*3*cutoff, m0, grid_size, cutoff, istride,
        f);
  else
    PNX(assign_f_r2r_pre_psi)(
        grid, pre_psi, m0, grid_size, cutoff, istride,
        f);
}

void PNX(assign_grad_f_c2c)(
    PNX(plan) ths, INT ind,
    C *grid, R *pre_psi, R *pre_dpsi, INT m0, INT *grid_size, int cutoff, int interlaced,
    C *grad_f
    )
{
  R* plan_pre_psi  = (interlaced) ? ths->pre_psi_il  : ths->pre_psi;
  R* plan_pre_dpsi = (interlaced) ? ths->pre_dpsi_il : ths->pre_dpsi;

  if(ths->pnfft_flags & PNFFT_PRE_FULL_PSI)
    PNX(assign_grad_f_c2c_pre_full_psi)(
        grid, plan_pre_psi + ind*PNFFT_POW3(cutoff), plan_pre_dpsi + 3*ind*PNFFT_POW3(cutoff),
        m0, grid_size, cutoff,
        grad_f);
  else if(ths->pnfft_flags & PNFFT_PRE_PSI)
    PNX(assign_grad_f_c2c_pre_psi)(
        grid, plan_pre_psi + ind*3*cutoff, plan_pre_dpsi + ind*3*cutoff,
        m0, grid_size, cutoff,
        grad_f);
  else
    PNX(assign_grad_f_c2c_pre_psi)(
        grid, pre_psi, pre_dpsi,
        m0, grid_size, cutoff,
        grad_f);
}

void PNX(assign_grad_f_r2r)(
    PNX(plan) ths, INT ind,
    R *grid, R *pre_psi, R *pre_dpsi, INT m0, INT *grid_size, int cutoff,
    INT istride, INT ostride, int interlaced,
    R *grad_f
    )
{
  R* plan_pre_psi  = (interlaced) ? ths->pre_psi_il  : ths->pre_psi;
  R* plan_pre_dpsi = (interlaced) ? ths->pre_dpsi_il : ths->pre_dpsi;

  if(ths->pnfft_flags & PNFFT_PRE_FULL_PSI)
    PNX(assign_grad_f_r2r_pre_full_psi)(
        grid, plan_pre_psi + ind*PNFFT_POW3(cutoff), plan_pre_dpsi + 3*ind*PNFFT_POW3(cutoff),
        m0, grid_size, cutoff, istride, ostride,
        grad_f);
  else if(ths->pnfft_flags & PNFFT_PRE_PSI)
    PNX(assign_grad_f_r2r_pre_psi)(
        grid, plan_pre_psi + ind*3*cutoff, plan_pre_dpsi + ind*3*cutoff,
        m0, grid_size, cutoff, istride, ostride,
        grad_f);
  else
    PNX(assign_grad_f_r2r_pre_psi)(
        grid, pre_psi, pre_dpsi,
        m0, grid_size, cutoff, istride, ostride,
        grad_f);
}

void PNX(assign_f_and_grad_f_c2c)(
    PNX(plan) ths, INT ind,
    C *grid, R *pre_psi, R *pre_dpsi, INT m0, INT *grid_size, int cutoff, int interlaced,
    C *f, C *grad_f
    )
{ 
  R* plan_pre_psi  = (interlaced) ? ths->pre_psi_il  : ths->pre_psi;
  R* plan_pre_dpsi = (interlaced) ? ths->pre_dpsi_il : ths->pre_dpsi;

  if(ths->pnfft_flags & PNFFT_PRE_FULL_PSI)
    PNX(assign_f_and_grad_f_c2c_pre_full_psi)(
        grid, plan_pre_psi + ind*PNFFT_POW3(cutoff), plan_pre_dpsi + 3*ind*PNFFT_POW3(cutoff),
        m0, grid_size, cutoff,
        f, grad_f);
  else if(ths->pnfft_flags & PNFFT_PRE_PSI)
    PNX(assign_f_and_grad_f_c2c_pre_psi)(
        grid, plan_pre_psi + ind*3*cutoff, plan_pre_dpsi + ind*3*cutoff,
        m0, grid_size, cutoff,
        f, grad_f);
  else
    PNX(assign_f_and_grad_f_c2c_pre_psi)(
        grid, pre_psi, pre_dpsi,
        m0, grid_size, cutoff,
        f, grad_f);
}

void PNX(assign_f_and_grad_f_r2r)(
    PNX(plan) ths, INT ind,
    R *grid, R *pre_psi, R *pre_dpsi, INT m0, INT *grid_size, int cutoff,
    INT istride, INT ostride, int interlaced,
    R *f, R *grad_f
    )
{ 
  R* plan_pre_psi  = (interlaced) ? ths->pre_psi_il  : ths->pre_psi;
  R* plan_pre_dpsi = (interlaced) ? ths->pre_dpsi_il : ths->pre_dpsi;

  if(ths->pnfft_flags & PNFFT_PRE_FULL_PSI)
    PNX(assign_f_and_grad_f_r2r_pre_full_psi)(
        grid, plan_pre_psi + ind*PNFFT_POW3(cutoff), plan_pre_dpsi + 3*ind*PNFFT_POW3(cutoff),
        m0, grid_size, cutoff, istride, ostride,
        f, grad_f);
  else if(ths->pnfft_flags & PNFFT_PRE_PSI)
    PNX(assign_f_and_grad_f_r2r_pre_psi)(
        grid, plan_pre_psi + ind*3*cutoff, plan_pre_dpsi + ind*3*cutoff,
        m0, grid_size, cutoff, istride, ostride,
        f, grad_f);
  else
    PNX(assign_f_and_grad_f_r2r_pre_psi)(
        grid, pre_psi, pre_dpsi,
        m0, grid_size, cutoff, istride, ostride,
        f, grad_f);
}









void PNX(spread_f_c2c_pre_psi)(
    C f, R *pre_psi, INT m0, INT *grid_size, int cutoff,
    C *grid
    )
{ 
  INT m1, m2, l0, l1, l2;
  R *pre_psi_x = &pre_psi[0*cutoff];
  R *pre_psi_y = &pre_psi[1*cutoff];
  R *pre_psi_z = &pre_psi[2*cutoff];

  for(l0=0; l0<cutoff; l0++, m0 += grid_size[1]*grid_size[2]){
    for(l1=0, m1=m0; l1<cutoff; l1++, m1 += grid_size[2]){
      R psi_xy = pre_psi_x[l0] * pre_psi_y[l1];
      for(l2=0, m2 = m1; l2<cutoff; l2++, m2++ ){
        grid[m2] += psi_xy * pre_psi_z[l2] * f;
      }
    }
  }
}

void PNX(spread_f_c2c_pre_full_psi)(
    C f, R *pre_psi, INT m0, INT *grid_size, int cutoff,
    C *grid
    )
{ 
  INT m1, m2, l0, l1, l2, m=0;
  for(l0=0; l0<cutoff; l0++, m0 += grid_size[1]*grid_size[2])
    for(l1=0, m1=m0; l1<cutoff; l1++, m1 += grid_size[2])
      for(l2=0, m2 = m1; l2<cutoff; l2++, m2++, m++ )
        grid[m2] += pre_psi[m] * f;
}

void PNX(spread_f_r2r_pre_psi)(
    R f, R *pre_psi, INT m0, INT *grid_size, int cutoff, INT ostride,
    R *grid
    )
{ 
  INT m1, m2, l0, l1, l2;
  R *pre_psi_x = &pre_psi[0*cutoff];
  R *pre_psi_y = &pre_psi[1*cutoff];
  R *pre_psi_z = &pre_psi[2*cutoff];

  for(l0=0; l0<cutoff; l0++, m0 += grid_size[1]*grid_size[2]*ostride){
    for(l1=0, m1=m0; l1<cutoff; l1++, m1 += grid_size[2]*ostride){
      R psi_xy = pre_psi_x[l0] * pre_psi_y[l1];
      for(l2=0, m2 = m1; l2<cutoff; l2++, m2+=ostride ){
        grid[m2] += psi_xy * pre_psi_z[l2] * f;
      }
    }
  }
}

void PNX(spread_f_r2r_pre_full_psi)(
    R f, R *pre_psi, INT m0, INT *grid_size, int cutoff, INT ostride,
    R *grid
    )
{ 
  INT m1, m2, l0, l1, l2, m=0;
  for(l0=0; l0<cutoff; l0++, m0 += grid_size[1]*grid_size[2]*ostride)
    for(l1=0, m1=m0; l1<cutoff; l1++, m1 += grid_size[2]*ostride)
      for(l2=0, m2 = m1; l2<cutoff; l2++, m2+=ostride, m++ )
        grid[m2] += pre_psi[m] * f;
}


void PNX(assign_f_c2c_pre_psi)(
    C *grid, R *pre_psi, INT m0, INT *grid_size, int cutoff,
    C *fv
    )
{ 
  INT m1, m2, l0, l1, l2;
  R *pre_psi_x = &pre_psi[0*cutoff];
  R *pre_psi_y = &pre_psi[1*cutoff];
  R *pre_psi_z = &pre_psi[2*cutoff];
  C f=0;

  for(l0=0; l0<cutoff; l0++, m0 += grid_size[1]*grid_size[2]){
    for(l1=0, m1=m0; l1<cutoff; l1++, m1 += grid_size[2]){
      R psi_xy = pre_psi_x[l0] * pre_psi_y[l1];
      for(l2=0, m2 = m1; l2<cutoff; l2++, m2++ ){
        f += psi_xy * pre_psi_z[l2] * grid[m2];
      }
    }
  }
  *fv += f;
}

void PNX(assign_f_c2c_pre_full_psi)(
    C *grid, R *pre_psi, INT m0, INT *grid_size, int cutoff,
    C *fv
    )
{ 
  INT m1, m2, l0, l1, l2, m=0;
  C f=0;

  for(l0=0; l0<cutoff; l0++, m0 += grid_size[1]*grid_size[2]){
    for(l1=0, m1=m0; l1<cutoff; l1++, m1 += grid_size[2]){
      for(l2=0, m2 = m1; l2<cutoff; l2++, m2++, m++ ){
        f += pre_psi[m] * grid[m2];
      }
    }
  }
  *fv += f;
}

void PNX(assign_f_r2r_pre_psi)(
    R *grid, R *pre_psi, INT m0, INT *grid_size, int cutoff, INT istride,
    R *fv
    )
{ 
  INT m1, m2, l0, l1, l2;
  R *pre_psi_x = &pre_psi[0*cutoff];
  R *pre_psi_y = &pre_psi[1*cutoff];
  R *pre_psi_z = &pre_psi[2*cutoff];
  R f=0;

  for(l0=0; l0<cutoff; l0++, m0 += grid_size[1]*grid_size[2]*istride){
    for(l1=0, m1=m0; l1<cutoff; l1++, m1 += grid_size[2]*istride){
      R psi_xy = pre_psi_x[l0] * pre_psi_y[l1];
      for(l2=0, m2 = m1; l2<cutoff; l2++, m2+=istride ){
        f += psi_xy * pre_psi_z[l2] * grid[m2];
      }
    }
  }
  *fv += f;;
}

void PNX(assign_f_r2r_pre_full_psi)(
    R *grid, R *pre_psi, INT m0, INT *grid_size, int cutoff, int istride,
    R *fv
    )
{ 
// fprintf(stderr, "Hier bin ich, r2r\n");
  INT m1, m2, l0, l1, l2, m=0;
  R f=0;
  
  for(l0=0; l0<cutoff; l0++, m0 += grid_size[1]*grid_size[2]*istride){
    for(l1=0, m1=m0; l1<cutoff; l1++, m1 += grid_size[2]*istride){
      for(l2=0, m2 = m1; l2<cutoff; l2++, m2+=istride, m++ ){
        f += pre_psi[m] * grid[m2];
      }
    }
  }
  *fv += f;;
}

void PNX(assign_grad_f_c2c_pre_psi)(
    C *grid, R *pre_psi, R *pre_dpsi, INT m0, INT *grid_size, int cutoff,
    C *grad_f
    )
{ 
  INT m1, m2, l0, l1, l2;
  R *pre_psi_x = &pre_psi[0*cutoff], *pre_dpsi_x = &pre_dpsi[0*cutoff];
  R *pre_psi_y = &pre_psi[1*cutoff], *pre_dpsi_y = &pre_dpsi[1*cutoff];
  R *pre_psi_z = &pre_psi[2*cutoff], *pre_dpsi_z = &pre_dpsi[2*cutoff];
  C g0=0, g1=0, g2=0;

  for(l0=0; l0<cutoff; l0++, m0 += grid_size[1]*grid_size[2]){
    for(l1=0, m1=m0; l1<cutoff; l1++, m1 += grid_size[2]){
      R psi_xy  = pre_psi_x[l0]  * pre_psi_y[l1];
      R psi_dxy = pre_dpsi_x[l0] * pre_psi_y[l1]; 
      R psi_xdy = pre_psi_x[l0]  * pre_dpsi_y[l1];
      for(l2=0, m2 = m1; l2<cutoff; l2++, m2++ ){
        g0 += psi_dxy * pre_psi_z[l2]  * grid[m2];
        g1 += psi_xdy * pre_psi_z[l2]  * grid[m2];
        g2 += psi_xy  * pre_dpsi_z[l2] * grid[m2];
      }
    }
  }
  grad_f[0] += g0; grad_f[1] += g1; grad_f[2] += g2;
}

void PNX(assign_grad_f_c2c_pre_full_psi)(
    C *grid, R *pre_psi, R *pre_dpsi, INT m0, INT *grid_size, int cutoff,
    C *grad_f
    )
{ 
  INT m1, m2, l0, l1, l2, dm=0;
  C g0=0, g1=0, g2=0;

  for(l0=0; l0<cutoff; l0++, m0 += grid_size[1]*grid_size[2]){
    for(l1=0, m1=m0; l1<cutoff; l1++, m1 += grid_size[2]){
      for(l2=0, m2 = m1; l2<cutoff; l2++, m2++, dm+=3 ){
        g0 += pre_dpsi[dm+0] * grid[m2];
        g1 += pre_dpsi[dm+1] * grid[m2];
        g2 += pre_dpsi[dm+2] * grid[m2];
      }
    }
  }
  grad_f[0] += g0; grad_f[1] += g1; grad_f[2] += g2;
}

void PNX(assign_grad_f_r2r_pre_psi)(
    R *grid, R *pre_psi, R *pre_dpsi, INT m0, INT *grid_size, int cutoff, INT istride, INT ostride,
    R *grad_f
    )
{ 
  INT m1, m2, l0, l1, l2;
  R *pre_psi_x = &pre_psi[0*cutoff], *pre_dpsi_x = &pre_dpsi[0*cutoff];
  R *pre_psi_y = &pre_psi[1*cutoff], *pre_dpsi_y = &pre_dpsi[1*cutoff];
  R *pre_psi_z = &pre_psi[2*cutoff], *pre_dpsi_z = &pre_dpsi[2*cutoff];
  R g0=0, g1=0, g2=0;

  for(l0=0; l0<cutoff; l0++, m0 += grid_size[1]*grid_size[2]*istride){
    for(l1=0, m1=m0; l1<cutoff; l1++, m1 += grid_size[2]*istride){
      R psi_xy  = pre_psi_x[l0]  * pre_psi_y[l1];
      R psi_dxy = pre_dpsi_x[l0] * pre_psi_y[l1]; 
      R psi_xdy = pre_psi_x[l0]  * pre_dpsi_y[l1];
      for(l2=0, m2 = m1; l2<cutoff; l2++, m2+=istride ){
        g0 += psi_dxy * pre_psi_z[l2]  * grid[m2];
        g1 += psi_xdy * pre_psi_z[l2]  * grid[m2];
        g2 += psi_xy  * pre_dpsi_z[l2] * grid[m2];
      }
    }
  }
  grad_f[0*ostride] += g0; grad_f[1*ostride] += g1; grad_f[2*ostride] += g2;
}

void PNX(assign_grad_f_r2r_pre_full_psi)(
    R *grid, R *pre_psi, R *pre_dpsi, INT m0, INT *grid_size, int cutoff, INT istride, INT ostride,
    R *grad_f
    )
{ 
  INT m1, m2, l0, l1, l2, dm=0;
  R g0=0, g1=0, g2=0;

  for(l0=0; l0<cutoff; l0++, m0 += grid_size[1]*grid_size[2]*istride){
    for(l1=0, m1=m0; l1<cutoff; l1++, m1 += grid_size[2]*istride){
      for(l2=0, m2 = m1; l2<cutoff; l2++, m2+=istride, dm+=3 ){
        g0 += pre_dpsi[dm+0] * grid[m2];
        g1 += pre_dpsi[dm+1] * grid[m2];
        g2 += pre_dpsi[dm+2] * grid[m2];
      }
    }
  }
  grad_f[0*ostride] += g0; grad_f[1*ostride] += g1; grad_f[2*ostride] += g2;
}



void PNX(assign_f_and_grad_f_c2c_pre_psi)(
    C *grid, R *pre_psi, R *pre_dpsi, INT m0, INT *grid_size, int cutoff,
    C *fv, C *grad_f
    )
{ 
  INT m1, m2, l0, l1, l2;
  R *pre_psi_x = &pre_psi[0*cutoff], *pre_dpsi_x = &pre_dpsi[0*cutoff];
  R *pre_psi_y = &pre_psi[1*cutoff], *pre_dpsi_y = &pre_dpsi[1*cutoff];
  R *pre_psi_z = &pre_psi[2*cutoff], *pre_dpsi_z = &pre_dpsi[2*cutoff];
  C f=0, g0=0, g1=0, g2=0;

  for(l0=0; l0<cutoff; l0++, m0 += grid_size[1]*grid_size[2]){
    for(l1=0, m1=m0; l1<cutoff; l1++, m1 += grid_size[2]){
      R psi_xy  = pre_psi_x[l0] * pre_psi_y[l1];
      R psi_dxy = pre_dpsi_x[l0] * pre_psi_y[l1]; 
      R psi_xdy = pre_psi_x[l0] * pre_dpsi_y[l1];
      for(l2=0, m2 = m1; l2<cutoff; l2++, m2++ ){
        f  += psi_xy  * pre_psi_z[l2]  * grid[m2];
        g0 += psi_dxy * pre_psi_z[l2]  * grid[m2];
        g1 += psi_xdy * pre_psi_z[l2]  * grid[m2];
        g2 += psi_xy  * pre_dpsi_z[l2] * grid[m2];
      }
    }
  }
  *fv += f;
  grad_f[0] += g0; grad_f[1] += g1; grad_f[2] += g2;
}

void PNX(assign_f_and_grad_f_c2c_pre_full_psi)(
    C *grid, R *pre_psi, R *pre_dpsi, INT m0, INT *grid_size, int cutoff,
    C *fv, C *grad_f
    )
{ 
  INT m1, m2, l0, l1, l2, m=0, dm=0;
  C f=0, g0=0, g1=0, g2=0;

  for(l0=0; l0<cutoff; l0++, m0 += grid_size[1]*grid_size[2]){
    for(l1=0, m1=m0; l1<cutoff; l1++, m1 += grid_size[2]){
      for(l2=0, m2 = m1; l2<cutoff; l2++, m2++, m++, dm+=3 ){
        f  += pre_psi[m]  * grid[m2];
        g0 += pre_dpsi[dm+0] * grid[m2];
        g1 += pre_dpsi[dm+1] * grid[m2];
        g2 += pre_dpsi[dm+2] * grid[m2];
      }
    }
  }
  *fv += f;
  grad_f[0] += g0; grad_f[1] += g1; grad_f[2] += g2;
}

void PNX(assign_f_and_grad_f_r2r_pre_psi)(
    R *grid, R *pre_psi, R *pre_dpsi, INT m0, INT *grid_size, int cutoff, INT istride, INT ostride,
    R *fv, R *grad_f
    )
{ 
  INT m1, m2, l0, l1, l2;
  R *pre_psi_x = &pre_psi[0*cutoff], *pre_dpsi_x = &pre_dpsi[0*cutoff];
  R *pre_psi_y = &pre_psi[1*cutoff], *pre_dpsi_y = &pre_dpsi[1*cutoff];
  R *pre_psi_z = &pre_psi[2*cutoff], *pre_dpsi_z = &pre_dpsi[2*cutoff];
  R f=0, g0=0, g1=0, g2=0;

  for(l0=0; l0<cutoff; l0++, m0 += grid_size[1]*grid_size[2]*istride){
    for(l1=0, m1=m0; l1<cutoff; l1++, m1 += grid_size[2]*istride){
      R psi_xy  = pre_psi_x[l0] * pre_psi_y[l1];
      R psi_dxy = pre_dpsi_x[l0] * pre_psi_y[l1]; 
      R psi_xdy = pre_psi_x[l0] * pre_dpsi_y[l1];
      for(l2=0, m2 = m1; l2<cutoff; l2++, m2+=istride ){
        f  += psi_xy  * pre_psi_z[l2]  * grid[m2];
        g0 += psi_dxy * pre_psi_z[l2]  * grid[m2];
        g1 += psi_xdy * pre_psi_z[l2]  * grid[m2];
        g2 += psi_xy  * pre_dpsi_z[l2] * grid[m2];
      }
    }
  }
  *fv += f;
  grad_f[0*ostride] += g0; grad_f[1*ostride] += g1; grad_f[2*ostride] += g2;
}

void PNX(assign_f_and_grad_f_r2r_pre_full_psi)(
    R *grid, R *pre_psi, R *pre_dpsi, INT m0, INT *grid_size, int cutoff, INT istride, INT ostride,
    R *fv, R *grad_f
    )
{
  INT m1, m2, l0, l1, l2, m=0, dm=0;
  R f=0, g0=0, g1=0, g2=0;

  for(l0=0; l0<cutoff; l0++, m0 += grid_size[1]*grid_size[2]*istride){
    for(l1=0, m1=m0; l1<cutoff; l1++, m1 += grid_size[2]*istride){
      for(l2=0, m2 = m1; l2<cutoff; l2++, m2+=istride, m++, dm+=3 ){
        f  += pre_psi[m]  * grid[m2];
        g0 += pre_dpsi[dm+0] * grid[m2];
        g1 += pre_dpsi[dm+1] * grid[m2];
        g2 += pre_dpsi[dm+2] * grid[m2];
      }
    }
  }
  *fv += f;
  grad_f[0*ostride] += g0; grad_f[1*ostride] += g1; grad_f[2*ostride] += g2;
}


