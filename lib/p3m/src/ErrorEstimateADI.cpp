/*
 Copyright (C) 2014 Olaf Lenz, Gabriel Sichardt
 Copyright (C) 2010,2011 The ESPResSo project
 Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
 Max-Planck-Institute for Polymer Research, Theory Group

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

#include "utils.hpp"
#include "ErrorEstimateADI.hpp"

namespace P3M {

p3m_float ErrorEstimateADI::compute_ks_error(Parameters& p, p3m_int num_charges,
		p3m_float sum_q2, p3m_float box_l[3]) {
	/* #ifdef P3M_ENABLE_DEBUG */
	/*   printf(  */
	/* 	  "        #%2d: N=%d Q2=%g box_l=(%g, %g, %g) grid=(%d, %d, %d) alpha=%g cao=%d\n",  */
	/* 	  comm_rank, N, sum_q2, box_l[0], box_l[1], box_l[2],  */
	/* 	  grid[0], grid[1], grid[2], alpha, cao); */
	/* #endif */
           
	p3m_float local_he_q = 0.0;
	p3m_float grid_i[3] =
			{ 1.0 / p.grid[0], 1.0 / p.grid[1], 1.0 / p.grid[2] };
	/* @todo Handle non-cubic case? */
	p3m_float alpha_L_i = 1. / (p.alpha * box_l[0]);

	// Distribute indices onto parallel tasks
	p3m_int num_ix = p.grid[0] * p.grid[1] * p.grid[2];
	p3m_int ix_per_task = num_ix / comm.size;
	p3m_int ix_rem = num_ix % comm.size;

	// Determine minimal and maximal index
	p3m_int min_ix = comm.rank * ix_per_task;
	p3m_int max_ix;
	if (comm.rank < ix_rem) {
		min_ix += comm.rank;
		max_ix = min_ix + ix_per_task + 1;
	} else {
		min_ix += ix_rem;
		max_ix = min_ix + ix_per_task;
	}

	// Loop over local indices
	for (p3m_int ix = min_ix; ix < max_ix; ix++) {
		p3m_int nx = ix / (p.grid[1] * p.grid[2]) - p.grid[0] / 2;
		p3m_int ny = ix % (p.grid[1] * p.grid[2]) / p.grid[2]
				- p.grid[1] / 2;
		p3m_int nz = ix % p.grid[2] - p.grid[2] / 2;

		/* #ifdef P3M_ENABLE_DEBUG */
		/*     printf( "#%d: (%d, %d, %d)\n", comm_rank, nx, ny, nz); */
		/* #endif */

		if (nx != 0 || ny != 0 || nz != 0) {
			p3m_float alias1, alias2, alias3, alias4, alias5, alias6;
			p3m_float D;

			k_space_error_sum2
			(nx,ny,nz, p.grid,grid_i, p.cao,alpha_L_i,
					&alias1,&alias2,&alias3,&alias4,&alias5,&alias6);

			D = alias1 - SQR(alias2) / (0.5*(alias3*alias4 + alias5*alias6));

			if (D > 0.0)
			local_he_q += D;
		}
	}
	p3m_float he_q;
	MPI_Reduce(&local_he_q, &he_q, 1, P3M_MPI_FLOAT, MPI_SUM, 0,
			comm.mpicomm);
        
	return 2.0 * sum_q2 *
			sqrt(he_q / num_charges) / (box_l[0] * box_l[1]);
}


p3m_float ErrorEstimateADI::compute_ks_error_triclinic(Parameters& p, p3m_int num_charges,
		p3m_float sum_q2, p3m_float  box_vectors[3][3], bool isTriclinic) {
    /* #ifdef P3M_ENABLE_DEBUG */
	/*   prp3m_intf(  */
	/* 	  "        #%2d: N=%d Q2=%g box_l=(%g, %g, %g) grid=(%d, %d, %d) alpha=%g cao=%d\n",  */
	/* 	  comm_rank, N, sum_q2, box_l[0], box_l[1], box_l[2],  */
	/* 	  grid[0], grid[1], grid[2], alpha, cao); */
	/* #endif */
    p3m_int comm_rank = 0; p3m_int comm_size=1;
	p3m_float local_he_q = 0.0;
	p3m_float grid_i[3] = { 1.0 / p.grid[0], 1.0 / p.grid[1], 1.0 / p.grid[2] };

	// Distribute indices onto parallel tasks
	p3m_int num_ix = p.grid[0] * p.grid[1] * p.grid[2];
	p3m_int ix_per_task = num_ix / comm_size;
	p3m_int ix_rem = num_ix % comm_size;

	// Determine minimal and maximal index
	p3m_int min_ix = comm_rank * ix_per_task;
	p3m_int max_ix;
	if (comm_rank < ix_rem) {
		min_ix += comm_rank;
		max_ix = min_ix + ix_per_task + 1;
	} else {
		min_ix += ix_rem;
		max_ix = min_ix + ix_per_task;
	}
        
	// Loop over local indices
	for (p3m_int ix = min_ix; ix < max_ix; ix++) {
		p3m_int nx = ix / (p.grid[1] * p.grid[2]) - p.grid[0] / 2;
		p3m_int ny = ix % (p.grid[1] * p.grid[2]) / p.grid[2]
				- p.grid[1] / 2;
		p3m_int nz = ix % p.grid[2] - p.grid[2] / 2;

		/* #ifdef P3M_ENABLE_DEBUG */
		/*     prp3m_intf( "#%d: (%d, %d, %d)\n", comm_rank, nx, ny, nz); */
		/* #endif */
                
		if (nx != 0 || ny != 0 || nz != 0) {
			p3m_float sum_Fref2, sqrt_nominator, sum_U2, sum_U2k2, alias5, alias6;
			p3m_float D;

			k_space_error_sum2_triclinic(nx,ny,nz, p.grid,grid_i, p.cao,p.alpha,
					&sum_Fref2,&sqrt_nominator,&sum_U2k2,&sum_U2,&alias5,&alias6, box_vectors, isTriclinic);

			D = sum_Fref2 - SQR(sqrt_nominator) / (0.5*(sum_U2*sum_U2k2 + alias5*alias6));


			if (D > 0.0)
			local_he_q += D;
		}
	}
        p3m_float he_q;
	MPI_Reduce(&local_he_q, &he_q, 1, P3M_MPI_FLOAT, MPI_SUM, 0,
			comm.mpicomm);

	return 2.0 * sum_q2 *
			sqrt(he_q / num_charges)/(box_vectors[0][0]*box_vectors[1][1]*box_vectors[2][2]);
                
    
    
}

/** aliasing sum used by \ref k_space_error. */
void ErrorEstimateADI::k_space_error_sum2(p3m_int nx, p3m_int ny, p3m_int nz,
		p3m_int grid[3], p3m_float grid_i[3], p3m_int cao, p3m_float alpha_L_i,
		p3m_float *alias1, p3m_float *alias2, p3m_float *alias3,
		p3m_float *alias4, p3m_float *alias5, p3m_float *alias6) {
	p3m_float prefactor = SQR(M_PI * alpha_L_i);

	*alias1 = *alias2 = *alias3 = *alias4 = *alias5 = *alias6 = 0.0;
	for (p3m_int mx = -P3M_BRILLOUIN; mx <= P3M_BRILLOUIN; mx++) {
		p3m_float nmx = nx + mx * grid[0];
		p3m_float fnmx = grid_i[0] * nmx;
		for (p3m_int my = -P3M_BRILLOUIN; my <= P3M_BRILLOUIN; my++) {
			p3m_float nmy = ny + my * grid[1];
			p3m_float fnmy = grid_i[1] * nmy;
			for (p3m_int mz = -P3M_BRILLOUIN; mz <= P3M_BRILLOUIN; mz++) {
				p3m_float nmz = nz + mz * grid[2];
				p3m_float fnmz = grid_i[2] * nmz;

				p3m_float nm2 = SQR(nmx) + SQR(nmy) + SQR(nmz);
				p3m_float ex = exp(-prefactor * nm2);

				p3m_float U2 = pow(sinc(fnmx) * sinc(fnmy) * sinc(fnmz),
						2.0 * cao);

				*alias1 += ex * ex / nm2;
				*alias2 += U2 * ex;
				*alias3 += U2 * nm2;
				*alias4 += U2;

				if (((mx + my + mz) % 2) == 0) {					//even term
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
/** aliasing sum used by \ref k_space_error. */
void ErrorEstimateADI::k_space_error_sum2_triclinic(p3m_int nx, p3m_int ny, p3m_int nz,
		p3m_int grid[3], p3m_float grid_i[3], p3m_int cao, p3m_float alpha,
		p3m_float *sum_Fref2, p3m_float *sqrt_nominator, p3m_float *sum_U2k2,
		p3m_float *sum_U2, p3m_float *alias5, p3m_float *alias6, p3m_float box_vectors[3][3], bool isTriclinic) {
	p3m_float pi2_alpha2 = SQR(M_PI / alpha);

	*sum_Fref2 = *sqrt_nominator = *sum_U2k2 = *sum_U2 = *alias5 = *alias6 = 0.0;
	for (p3m_int mx = -P3M_BRILLOUIN; mx <= P3M_BRILLOUIN; mx++) {
		p3m_int nmx = nx + mx * grid[0];
		p3m_float fnmx = grid_i[0] * nmx;
		for (p3m_int my = -P3M_BRILLOUIN; my <= P3M_BRILLOUIN; my++) {
			p3m_int nmy = ny + my * grid[1];
			p3m_float fnmy = grid_i[1] * nmy;
			for (p3m_int mz = -P3M_BRILLOUIN; mz <= P3M_BRILLOUIN; mz++) {
				p3m_int nmz = nz + mz * grid[2];
				p3m_float fnmz = grid_i[2] * nmz;
                                p3m_float k2_4pi2 = 0.0;
                                if (!isTriclinic) {
                                    k2_4pi2 =
                                            SQR(nmx / box_vectors[RX][RX]) +
                                            SQR(nmy / box_vectors[RY][RY]) +
                                            SQR(nmz / box_vectors[RZ][RZ]);
                                } else {
                                    p3m_int i;
                                    const p3m_int nNm[3] = {nmx, nmy, nmz};

                                    for (i = 0; i < 3; i++) {
                                        p3m_int j = (i + 1) % 3;
                                        p3m_int k = (i + 2) % 3;

                                        k2_4pi2 += SQR((nNm[i]*(box_vectors[j][j] * box_vectors[k][k]\
                                 - box_vectors[j][k] * box_vectors[k][j]) + nNm[j]\
                                 *(box_vectors[k][j] * box_vectors[i][k] - box_vectors[k][k]\
                                 * box_vectors[i][j]) + nNm[k]*(box_vectors[i][j]\
                                 * box_vectors[j][k] - box_vectors[i][k] * box_vectors[j][j])));
                                    }

                                    k2_4pi2 *= 1 / SQR(box_vectors[0][0]*box_vectors[1][1]*box_vectors[2][2]);

                                }
                                
				p3m_float ex = exp(-pi2_alpha2 * k2_4pi2); // exp(-k_{n+Nm}^2 / (4*alpha^2))

				p3m_float U2 = pow(sinc(fnmx) * sinc(fnmy) * sinc(fnmz),
						2.0 * cao);             // U(k_{n+Nm})^2

				*sum_Fref2 += ex * ex / k2_4pi2;        // 1/4 * Sum_m [ |R(k_{n+Nm})|^2 ]
				*sqrt_nominator += U2 * ex;             // Sum_m [ U(k_{n+Nm})^2 * exp(-k_{n+Nm}^2 / (4*alpha^2)) ]
				*sum_U2k2 += U2 * k2_4pi2 ;             // 1/(4*alpha^2) * Sum_m [ U(k_{n+Nm})^2 * k_{n+Nm}^2 ]
				*sum_U2 += U2;                          // Sum_m [ U(k_{n+Nm})^2 ]
                                
				if (((mx + my + mz) % 2) == 0) {	//even term
					*alias5 += U2 * k2_4pi2;
					*alias6 += U2;
				} else {				//odd term: minus sign!
					*alias5 -= U2 * k2_4pi2;
					*alias6 -= U2;
				}
			}
		}
	}
}
}



