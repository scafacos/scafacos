/*
  Copyright (C) 2011,2012,2013 Olaf Lenz
  
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
#include "error_estimate.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/** maximal precision */
#ifdef FCS_FLOAT_IS_DOUBLE
static const fcs_float ROUND_ERROR_PREC = 1.0e-14;
#else
static const fcs_float ROUND_ERROR_PREC = 1.0e-6;
#endif

static const fcs_float FULL_ESTIMATE_ALPHA_H_THRESHOLD = 0.5;

/***************************************************/
/* FORWARD DECLARATIONS OF INTERNAL FUNCTIONS */
/***************************************************/
static fcs_float 
ifcs_p3m_k_space_error_sum1(fcs_int n, fcs_float grid_i, 
                            fcs_int cao);

#ifdef P3M_AD
#ifdef P3M_INTERLACE
static void 
ifcs_p3m_k_space_error_sum2_adi(fcs_int nx, fcs_int ny, fcs_int nz, 
                                fcs_int grid[3], fcs_float grid_i[3], 
                                fcs_int cao, fcs_float alpha_L_i, 
                                fcs_float *alias1, fcs_float *alias2,
                                fcs_float *alias3, fcs_float *alias4,
                                fcs_float *alias5, fcs_float *alias6);
#endif
#endif

static void 
ifcs_p3m_k_space_error_sum2(fcs_int nx, fcs_int ny, fcs_int nz, 
                            fcs_int grid[3], fcs_float grid_i[3], 
                            fcs_int cao, fcs_float alpha_L_i, 
                            fcs_float *alias1, fcs_float *alias2);

static void 
ifcs_p3m_compute_error_estimate(ifcs_p3m_data_struct *d);

static void 
ifcs_p3m_determine_good_alpha(ifcs_p3m_data_struct *d);

static void
ifcs_p3m_ks_error_broadcast(ifcs_p3m_data_struct *d);

/***************************************************/
/* IMPLEMENTATION */
/***************************************************/
/** Calculates the SQuaRe of 'fcs_float' x, returning 'fcs_float'. */
static inline fcs_float SQR(fcs_float x) { return x*x; }

/** Calculates the sinc-function as sin(PI*x)/(PI*x).
 *
 * (same convention as in Hockney/Eastwood). In order to avoid
 * divisions by 0, arguments, whose modulus is smaller than epsi, will
 * be evaluated by an 8th order Taylor expansion of the sinc
 * function. Note that the difference between sinc(x) and this
 * expansion is smaller than 0.235e-12, if x is smaller than 0.1. (The
 * next term in the expansion is the 10th order contribution
 * PI^10/39916800 * x^10 = 0.2346...*x^12).  This expansion should
 * also save time, since it reduces the number of function calls to
 * sin().  
*/
static inline fcs_float sinc(fcs_float d)
{
  const fcs_float epsi = 0.1;
  const fcs_float c2 = -0.1666666666667e-0;
  const fcs_float c4 = 0.8333333333333e-2;
  const fcs_float c6 = -0.1984126984127e-3;
  const fcs_float c8 = 0.2755731922399e-5;

  fcs_float PId = FCS_PI*d, PId2;

  if (fabs(d)>epsi)
    return sin(PId)/PId;
  else {
    PId2 = SQR(PId);
    return 1.0 + PId2*(c2+PId2*(c4+PId2*(c6+PId2*c8)));
  }
}

/** Calculates the real space contribution to the rms error in the
    force (as described by Kolafa and Perram).
*/
void 
ifcs_p3m_real_space_error(ifcs_p3m_data_struct *d) {
  d->rs_error = (2.0*d->sum_q2*exp(-SQR(d->r_cut*d->alpha))) 
    / (sqrt((fcs_float)d->sum_qpart*d->r_cut*d->box_l[0]*d->box_l[1]*d->box_l[2]));
}

/** Calculates the reciprocal space contribution to the rms error in the
    force (as described in the book of Hockney and Eastwood
    (Eqn. 8.23) (for a system of N randomly distributed particles in a
    cubic box).
    \param N        	number of charged particles in the system
    \param sum_q2   	sum of square of charges in the system
    \param box_l    	system size
    \param grid     	number of grid points in the different directions
    \param alpha	ewald splitting parameter
    \param cao		charge assignment order
    \return 		reciprocal space error
*/
void 
ifcs_p3m_k_space_error(ifcs_p3m_data_struct *d) {
  int comm_rank, comm_size;
  MPI_Comm_rank(d->comm.mpicomm, &comm_rank);
  MPI_Comm_size(d->comm.mpicomm, &comm_size);

  // receiving function in ifcs_p3m_param_broadcast_slave
  if (comm_rank == 0) ifcs_p3m_ks_error_broadcast(d);

/* #ifdef FCS_ENABLE_DEBUG */
/*   printf(  */
/* 	  "        #%2d: N=%d Q2=%g box_l=(%g, %g, %g) grid=(%d, %d, %d) alpha=%g cao=%d\n",  */
/* 	  comm_rank, N, sum_q2, box_l[0], box_l[1], box_l[2],  */
/* 	  grid[0], grid[1], grid[2], alpha, cao); */
/* #endif */
      
  fcs_float local_he_q = 0.0;
  fcs_float grid_i[3] = 
    {1.0/d->grid[0], 1.0/d->grid[1], 1.0/d->grid[2]};
  /* @todo Handle non-cubic case? */
  fcs_float alpha_L_i = 1./(d->alpha*d->box_l[0]);

  // Distribute indices onto parallel tasks
  fcs_int num_ix = d->grid[0]*d->grid[1]*d->grid[2];
  fcs_int ix_per_task = num_ix / comm_size;
  fcs_int ix_rem = num_ix % comm_size;

  // Determine minimal and maximal index
  fcs_int min_ix = comm_rank * ix_per_task;
  fcs_int max_ix;
  if (comm_rank < ix_rem) {
    min_ix += comm_rank;
    max_ix = min_ix + ix_per_task + 1;
  } else {
    min_ix += ix_rem;
    max_ix = min_ix + ix_per_task;
  }

  // Loop over local indices
  for (fcs_int ix = min_ix; ix < max_ix; ix++) {
    fcs_int nx = ix / (d->grid[1]*d->grid[2]) - d->grid[0]/2;
    fcs_int ny = ix % (d->grid[1]*d->grid[2]) / d->grid[2] - d->grid[1]/2;
    fcs_int nz = ix % d->grid[2] - d->grid[2]/2;

/* #ifdef FCS_ENABLE_DEBUG */
/*     printf( "#%d: (%d, %d, %d)\n", comm_rank, nx, ny, nz); */
/* #endif */

#ifdef P3M_INTERLACE
#ifdef P3M_AD
    if (nx != 0 || ny != 0 || nz != 0) {
      fcs_float alias1, alias2, alias3, alias4, alias5, alias6;
      fcs_float D;

      ifcs_p3m_k_space_error_sum2_adi
        (nx,ny,nz, d->grid,grid_i, d->cao,alpha_L_i,
         &alias1,&alias2,&alias3,&alias4,&alias5,&alias6);

      D = alias1  -  SQR(alias2) / (0.5*(alias3*alias4 + alias5*alias6));

      if (D > 0.0) 
	local_he_q += D;
    }
#endif /* P3M_AD */
#else
    fcs_int old_nx = 1000000;
    fcs_int old_ny = 1000000;
    fcs_float ctan_x, ctan_y;
    if (ny != old_ny) {
      if (nx != old_nx) {
	old_nx = nx;
	ctan_x = ifcs_p3m_k_space_error_sum1(nx,grid_i[0],d->cao);
      }
      old_ny = ny;
      ctan_y = ifcs_p3m_k_space_error_sum1(ny,grid_i[1],d->cao);
    }

    if (nx != 0 || ny != 0 || nz != 0) {
      fcs_float n2 = nx*nx+ny*ny+nz*nz;
      fcs_float cs = ctan_x * ctan_y * ifcs_p3m_k_space_error_sum1(nz,grid_i[2],d->cao);
      fcs_float alias1, alias2;
      // TODO: sum2_ad for IKI?
      ifcs_p3m_k_space_error_sum2(nx,ny,nz,d->grid,grid_i,d->cao,alpha_L_i,
                                  &alias1,&alias2);
      fcs_float d = alias1  -  SQR(alias2/cs) / n2;
      /* at high precisions, d can become negative due to extinction;
	 also, don't take values that have no significant digits left*/
      if (d > 0.0 && (fabs(d/alias1) > ROUND_ERROR_PREC))
	local_he_q += d;
    }
#endif /* P3M_INTERLACE */
  }
  fcs_float he_q;
  MPI_Reduce(&local_he_q, &he_q, 1, FCS_MPI_FLOAT, MPI_SUM, 0, 
             d->comm.mpicomm);
#ifdef P3M_INTERLACE
#ifdef P3M_AD
  fcs_float ks_error = 
    2.0*d->sum_q2*sqrt( he_q / (fcs_float)d->sum_qpart ) 
    / (d->box_l[0] * d->box_l[0]);
#endif
#else
  fcs_float ks_error = 2.0*d->sum_q2*sqrt
    ( he_q / (fcs_float)d->sum_qpart / (d->box_l[0]*d->box_l[1]*d->box_l[2]) );
#endif

  d->ks_error = ks_error;
}

/** One of the aliasing sums used by \ref fcs_fftcommon_k_space_error. 
    (fortunately the one which is most important (because it converges
    most slowly, since it is not damped exponentially)) can be
    calculated analytically. The result (which depends on the order of
    the spline interpolation) can be written as an even trigonometric
    polynomial. The results are tabulated here (The employed formula
    is Eqn. 7.66 in the book of Hockney and Eastwood). */
fcs_float 
ifcs_p3m_k_space_error_sum1(fcs_int n, fcs_float grid_i, 
                            fcs_int cao) {
  fcs_float c, res=0.0;
  c = SQR(cos(FCS_PI*grid_i*(fcs_float)n));
  
  switch (cao) {
  case 1 : { 
    res = 1; 
    break; }
  case 2 : { 
    res = (1.0+c*2.0)/3.0; 
    break; }
  case 3 : { 
    res = (2.0+c*(11.0+c*2.0))/15.0; 
    break; }
  case 4 : { 
    res = (17.0+c*(180.0+c*(114.0+c*4.0)))/315.0; 
    break; }
  case 5 : { 
    res = (62.0+c*(1072.0+c*(1452.0+c*(247.0+c*2.0))))/2835.0; 
    break; }
  case 6 : { 
    res = (1382.0+c*(35396.0+c*(83021.0+c*(34096.0+c*(2026.0+c*4.0)))))/155925.0; 
    break; }
  case 7 : { 
    res = (21844.0+c*(776661.0+c*(2801040.0+c*(2123860.0+c*(349500.0+c*(8166.0+c*4.0))))))/6081075.0; 
    break; }
  default : {
    printf("INTERNAL_ERROR: The value %d for the interpolation order should not occur!\n", cao);
    exit(1);
  }
  }
  
  return res;
}

#ifdef P3M_AD
#ifdef P3M_INTERLACE
/** aliasing sum used by \ref ifcs_p3m_k_space_error. */
void 
ifcs_p3m_k_space_error_sum2_adi(fcs_int nx, fcs_int ny, fcs_int nz, 
                                fcs_int grid[3], fcs_float grid_i[3], 
                                fcs_int cao, fcs_float alpha_L_i, 
                                fcs_float *alias1, fcs_float *alias2,
                                fcs_float *alias3, fcs_float *alias4,
                                fcs_float *alias5, fcs_float *alias6)
{
  fcs_float prefactor = SQR(FCS_PI*alpha_L_i);

  *alias1 = *alias2 = *alias3 = *alias4 = *alias5 = *alias6 = 0.0;
  for (fcs_int mx=-P3M_BRILLOUIN; mx<=P3M_BRILLOUIN; mx++) {
    fcs_float nmx = nx + mx*grid[0];
    fcs_float fnmx = grid_i[0] * nmx;
    for (fcs_int my=-P3M_BRILLOUIN; my<=P3M_BRILLOUIN; my++) {
      fcs_float nmy = ny + my*grid[1];
      fcs_float fnmy = grid_i[1] * nmy;
      for (fcs_int mz=-P3M_BRILLOUIN; mz<=P3M_BRILLOUIN; mz++) {
	fcs_float nmz = nz + mz*grid[2];
	fcs_float fnmz = grid_i[2] * nmz;

	fcs_float nm2 = SQR(nmx) + SQR(nmy) + SQR(nmz);
	fcs_float ex = exp(-prefactor*nm2);
	
	fcs_float U2 = pow(sinc(fnmx)*sinc(fnmy)*sinc(fnmz), 2.0*cao);
	
	*alias1 += ex*ex / nm2;
	*alias2 += U2 * ex;
	*alias3 += U2 * nm2;
	*alias4 += U2;

        if (((mx+my+mz)%2)==0) {					//even term
	   *alias5 += U2 * nm2;
	   *alias6 += U2;
	 } else {						//odd term: minus sign!
	   *alias5 -= U2 * nm2;
	   *alias6 -= U2;
	 }
      }
    }
  }
}
#endif
#endif

/** aliasing sum used by \ref ifcs_p3m_k_space_error. */
void 
ifcs_p3m_k_space_error_sum2(fcs_int nx, fcs_int ny, fcs_int nz, 
                            fcs_int grid[3], fcs_float grid_i[3], 
                            fcs_int cao, fcs_float alpha_L_i, 
                            fcs_float *alias1, fcs_float *alias2)
{
  fcs_float prefactor = SQR(FCS_PI*alpha_L_i);

  *alias1 = *alias2 = 0.0;
  for (fcs_int mx=-P3M_BRILLOUIN; mx<=P3M_BRILLOUIN; mx++) {
    fcs_float nmx = nx + mx*grid[0];
    fcs_float fnmx = grid_i[0] * nmx;
    for (fcs_int my=-P3M_BRILLOUIN; my<=P3M_BRILLOUIN; my++) {
      fcs_float nmy = ny + my*grid[1];
      fcs_float fnmy = grid_i[1] * nmy;
      for (fcs_int mz=-P3M_BRILLOUIN; mz<=P3M_BRILLOUIN; mz++) {
	fcs_float nmz = nz + mz*grid[2];
	fcs_float fnmz = grid_i[2] * nmz;

	fcs_float nm2 = SQR(nmx) + SQR(nmy) + SQR(nmz);
	fcs_float ex = exp(-prefactor*nm2);
	
	fcs_float U2 = pow(sinc(fnmx)*sinc(fnmy)*sinc(fnmz), 2.0*cao);
	
	*alias1 += ex*ex / nm2;
	*alias2 += U2 * ex * (nx*nmx + ny*nmy + nz*nmz) / nm2;
      }
    }
  }
}


/** Calculate the analytical approximation for the k-space part of the
    error (Eq. 38 in Deserno, Holm; JCP 109,18; 1998). */
void
ifcs_p3m_k_space_error_approx(ifcs_p3m_data_struct *d) {
  /* grid spacing */
  /* TODO: non-cubic case*/
  fcs_float h = d->box_l[0]/d->grid[0];
  fcs_float ha = h*d->alpha;

  /* compute the sum in eq. 38 */
  fcs_float sum;
  switch (d->cao) {
  case 1:
    sum = 2./3.; 
    break;
  case 2:
    sum = 5./294.*pow(ha,2) + 1./50.; 
    break;
  case 3:
    sum = 
      21./3872.*pow(ha, 4) + 7./1440.*pow(ha,2) + 
      1./588.;
    break;
  case 4:
    sum = 
      143./28800.*pow(ha, 6) + 
      7601./2271360.*pow(ha, 4) + 3./1936.*pow(ha,2) + 1./4320.;
    break;
  case 5:
    sum = 
      106640677./11737571328.*pow(ha,8) + 517231./106536960.*pow(ha, 6) + 
      143./69120.*pow(ha, 4) + 7601./13628160.*pow(ha,2) + 
      1./23232.;
    break;
  case 6:
    sum = 
      326190917./11700633600.*pow(ha,10) + 
      733191589./59609088000.*pow(ha,8) + 9694607./2095994880.*pow(ha, 6) + 
      47021./35512320.*pow(ha, 4) + 13./57600.*pow(ha,2) + 
      691./68140800.;
    break;
  case 7:
    sum =
      4887769399./37838389248.*pow(ha,12) + 1755948832039./36229939200000.*pow(ha,10) + 
      25091609./1560084480.*pow(ha,8) + 56399353./12773376000.*pow(ha, 6) + 
      745739./838397952.*pow(ha, 4) + 3617./35512320.*pow(ha,2) + 
      1./345600.;
    break;
  default: 
    printf("INTERNAL_ERROR: ifcs_p3m_k_space_error_approx: "
           "Charge assignment order of %d should not occur!\n", d->cao);
    exit(1);
  };

  d->ks_error =
    d->sum_q2 / (d->box_l[0]*d->box_l[0]) *
    pow(h*d->alpha, d->cao) * 
    sqrt(d->alpha*d->box_l[0]/d->sum_qpart*sqrt(2.0*FCS_PI)*sum);
}

/** Calculates the rms error estimate in the force (as described in
    the book of Hockney and Eastwood (Eqn. 8.23) for a system of N
    randomly distributed particles.
*/
void
ifcs_p3m_compute_error_estimate(ifcs_p3m_data_struct *d) {
#ifdef FCS_ENABLE_DEBUG
  int comm_rank;
  MPI_Comm_rank(d->comm.mpicomm, &comm_rank);
#endif

  /* calculate real space and k space error */
  ifcs_p3m_real_space_error(d);

  fcs_int full_estimate = 0;
  // use the full estimate if alpha*h is larger than the threshold in any dimension
  for (fcs_int i = 0; i < 3 && !full_estimate; i++) {
    fcs_float alpha_h = d->alpha * d->box_l[i] / d->grid[i];
    full_estimate = alpha_h > FULL_ESTIMATE_ALPHA_H_THRESHOLD;
#ifdef FCS_ENABLE_DEBUG
    if (comm_rank == 0) {
      if (full_estimate)
	printf( "        alpha*h[%d]=%" FCS_LMOD_FLOAT "g > %" 
                FCS_LMOD_FLOAT "g => full estimate\n", 
                i, alpha_h, FULL_ESTIMATE_ALPHA_H_THRESHOLD);
    }
#endif
  }
#ifdef FCS_ENABLE_DEBUG
  if (comm_rank == 0) {
    if (!full_estimate)
      printf( "        alpha*h < %" FCS_LMOD_FLOAT "g => approximation\n", 
              FULL_ESTIMATE_ALPHA_H_THRESHOLD);
  }
#endif

  if (full_estimate)
    ifcs_p3m_k_space_error(d);
  else
    ifcs_p3m_k_space_error_approx(d);

  d->error = sqrt(SQR(d->rs_error)+SQR(d->ks_error));

#ifdef FCS_ENABLE_DEBUG
  if (comm_rank == 0)
    printf( "        error estimate: rs_err=%" FCS_LMOD_FLOAT 
            "e ks_err=%" FCS_LMOD_FLOAT "e err=%" FCS_LMOD_FLOAT "e\n",
	    d->rs_error, d->ks_error, d->error);
#endif
}

/** Determines a value for alpha that achieves the wanted_error, if at
    all possible. Also returns the achieved errors with these
    parameters. Check whether wanted_error > achieved_error to see
    whether the required error can actually be met.
*/
void
ifcs_p3m_determine_good_alpha(ifcs_p3m_data_struct *d) {
  /* Get the real space error for alpha=0 */
  d->alpha = 0.0;
  ifcs_p3m_real_space_error(d);
  fcs_float max_rs_err = d->rs_error;

  /* We know how the real space error behaves, so we can compute the
     alpha where the real space error is half of the wanted
     error. This is the alpha that we return. */
  if(FCS_SQRT2*max_rs_err > d->tolerance_field) {
    d->alpha = sqrt(log(FCS_SQRT2*max_rs_err/d->tolerance_field)) / d->r_cut;
  } else {
    /* if the error is small enough even for alpha=0 */
    d->alpha = 0.1 * d->box_l[0];
  }
  
#ifdef FCS_ENABLE_DEBUG
  int comm_rank;
  MPI_Comm_rank(d->comm.mpicomm, &comm_rank);
  if (comm_rank == 0)
    printf( "        determined alpha=%" FCS_LMOD_FLOAT "g\n", d->alpha);
#endif

  ifcs_p3m_compute_error_estimate(d);
}

void
ifcs_p3m_compute_error_and_tune_alpha(ifcs_p3m_data_struct *d) {
  if (d->tune_alpha)
    ifcs_p3m_determine_good_alpha(d);
  else
    ifcs_p3m_compute_error_estimate(d);
}

void
ifcs_p3m_ks_error_broadcast(ifcs_p3m_data_struct *d) {
  fcs_int int_buffer[5];
  
  // pack int data
  int_buffer[0] = 0;		/* param set for tuning */
  int_buffer[1] = d->grid[0];
  int_buffer[2] = d->grid[1];
  int_buffer[3] = d->grid[2];
  int_buffer[4] = d->cao;
  MPI_Bcast(int_buffer, 5, FCS_MPI_INT, 0, d->comm.mpicomm);
  
  // pack float data
  MPI_Bcast(&d->alpha, 1, FCS_MPI_FLOAT, 0, d->comm.mpicomm);
}

/* Broadcast final parameters */
void
ifcs_p3m_param_broadcast(ifcs_p3m_data_struct *d) {
  // broadcast parameters and mark as final
  fcs_int int_buffer[5];
  fcs_float float_buffer[5];

  // pack int data
  int_buffer[0] = 1;		/* final param set */
  int_buffer[1] = d->grid[0];
  int_buffer[2] = d->grid[1];
  int_buffer[3] = d->grid[2];
  int_buffer[4] = d->cao;
  MPI_Bcast(int_buffer, 5, FCS_MPI_INT, 0, d->comm.mpicomm);

  // pack float data
  float_buffer[0] = d->alpha;
  float_buffer[1] = d->r_cut;
  float_buffer[2] = d->error;
  float_buffer[3] = d->rs_error;
  float_buffer[4] = d->ks_error;
  MPI_Bcast(float_buffer, 5, FCS_MPI_FLOAT, 0, d->comm.mpicomm);
}

/* Slave error estimation loop. 
   Run this on a slave. When the function returns, r_cut, grid, alpha,
   cao, error, rs_error and ks_error are set the same on all nodes.
*/
void
ifcs_p3m_param_broadcast_slave(ifcs_p3m_data_struct *d) {
  fcs_int int_buffer[5];
  fcs_float float_buffer[5];
  fcs_int final;

  for (;;) {
    // unpack int data
    MPI_Bcast(int_buffer, 5, FCS_MPI_INT, 0, d->comm.mpicomm);
    final = int_buffer[0];
    d->grid[0] = int_buffer[1];
    d->grid[1] = int_buffer[2];
    d->grid[2] = int_buffer[3];
    d->cao = int_buffer[4];

    if (final) {
      // unpack float data
      MPI_Bcast(float_buffer, 5, FCS_MPI_FLOAT, 0, d->comm.mpicomm);
      d->alpha = float_buffer[0];
      d->r_cut = float_buffer[1];
      d->error = float_buffer[2];
      d->rs_error = float_buffer[3];
      d->ks_error = float_buffer[4];

      break;
    } else {
      MPI_Bcast(float_buffer, 1, FCS_MPI_FLOAT, 0, d->comm.mpicomm);
      d->alpha = float_buffer[0];

      // do slave portion of k_space error estimation
      ifcs_p3m_k_space_error(d);
    }
  }
}

