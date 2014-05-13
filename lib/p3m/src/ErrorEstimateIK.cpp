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
#include "ErrorEstimateIK.hpp"

namespace P3M {
static const p3m_float FULL_ESTIMATE_ALPHA_H_THRESHOLD = 0.5;

p3m_float ErrorEstimateIK::compute_ks_error(Parameters& p, p3m_int num_charges,
		p3m_float sum_q2, p3m_float box_l[3]) {
	bool full_estimate = false;
	// use the full estimate if alpha*h is larger than the threshold in any dimension
	for (int i = 0; i < 3; i++) {
	    p3m_float alpha_h = p.alpha * box_l[i] / p.grid[i];
	    full_estimate = alpha_h > FULL_ESTIMATE_ALPHA_H_THRESHOLD;
	    if (full_estimate) {
	        P3M_DEBUG(printf("        alpha*h[%d]=" FFLOAT " > "        \
	                FFLOAT " => full estimate\n", i, alpha_h,   \
	                FULL_ESTIMATE_ALPHA_H_THRESHOLD));
	        break;
	    }
	}

#ifdef P3M_ENABLE_DEBUG
	if (!full_estimate)
	    printf("        alpha*h < " FFLOAT " => approximation\n",
	            FULL_ESTIMATE_ALPHA_H_THRESHOLD);
#endif

	if (full_estimate)
	    return compute_ks_error_full(p, num_charges, sum_q2, box_l);
	else
	    return compute_ks_error_approx(p, num_charges, sum_q2, box_l);
}
p3m_float ErrorEstimateIK::compute_ks_error_triclinic(Parameters& p, p3m_int num_charges,
		p3m_float sum_q2, p3m_float box_vectors[3][3], bool isTriclinic) {
//	bool full_estimate = false;
	// use the full estimate if alpha*h is larger than the threshold in any dimension
//	for (int i = 0; i < 3; i++) {
//	    p3m_float alpha_h = p.alpha * box_vectors[i][i] / p.grid[i];
//	    full_estimate = alpha_h > FULL_ESTIMATE_ALPHA_H_THRESHOLD;
//	    if (full_estimate) {
//	        P3M_DEBUG(printf("        alpha*h[%d]=" FFLOAT " > "        \
//	                FFLOAT " => full estimate\n", i, alpha_h,   \
//	                FULL_ESTIMATE_ALPHA_H_THRESHOLD));
//	        break;
//	    }
//	}

//#ifdef P3M_ENABLE_DEBUG
//	if (!full_estimate)
//	    printf("        alpha*h < " FFLOAT " => approximation\n",
//	            FULL_ESTIMATE_ALPHA_H_THRESHOLD);
//#endif

	//if (full_estimate)
    fcs_float box_l[3]={box_vectors[0][0],box_vectors[1][1],box_vectors[2][2]};
    return compute_ks_error_full(p, num_charges, sum_q2, box_l);//todo: remove this and get the triclinic estimate right.
//	    return compute_ks_error_full_triclinic(p, num_charges, sum_q2, box_vectors, isTriclinic);
	//else
	//    return compute_ks_error_approx_triclinic(p, num_charges, sum_q2, box_vectors, isTriclinic);
}

p3m_float ErrorEstimateIK::compute_ks_error_approx(Parameters& p, p3m_int num_charges,
		p3m_float sum_q2, p3m_float box_l[3]) {

	p3m_float h = box_l[0] / p.grid[0];
	p3m_float ha = h * p.alpha;

	/* compute the sum in eq. 38 */
	p3m_float sum;
	switch (p.cao) {
	case 1:
		sum = 2. / 3.;
		break;
	case 2:
		sum = 5. / 294. * pow(ha, 2) + 1. / 50.;
		break;
	case 3:
		sum = 21. / 3872. * pow(ha, 4) + 7. / 1440. * pow(ha, 2) + 1. / 588.;
		break;
	case 4:
		sum = 143. / 28800. * pow(ha, 6) + 7601. / 2271360. * pow(ha, 4)
		+ 3. / 1936. * pow(ha, 2) + 1. / 4320.;
		break;
	case 5:
		sum = 106640677. / 11737571328. * pow(ha, 8)
		+ 517231. / 106536960. * pow(ha, 6) + 143. / 69120. * pow(ha, 4)
		+ 7601. / 13628160. * pow(ha, 2) + 1. / 23232.;
		break;
	case 6:
		sum = 326190917. / 11700633600. * pow(ha, 10)
		+ 733191589. / 59609088000. * pow(ha, 8)
		+ 9694607. / 2095994880. * pow(ha, 6)
		+ 47021. / 35512320. * pow(ha, 4) + 13. / 57600. * pow(ha, 2)
		+ 691. / 68140800.;
		break;
	case 7:
		sum = 4887769399. / 37838389248. * pow(ha, 12)
		+ 1755948832039. / 36229939200000. * pow(ha, 10)
		+ 25091609. / 1560084480. * pow(ha, 8)
		+ 56399353. / 12773376000. * pow(ha, 6)
		+ 745739. / 838397952. * pow(ha, 4)
		+ 3617. / 35512320. * pow(ha, 2) + 1. / 345600.;
		break;
	default:
		throw std::logic_error("INTERNAL_ERROR: k_space_error_approx: "
				"Charge assignment order should not occur!\n");
	};

	return sum_q2 / (box_l[0] * box_l[0])
					* pow(h * p.alpha, p.cao)
					* sqrt(p.alpha * box_l[0] / num_charges * sqrt(2.0 * M_PI) * sum);
}

/** Calculates the reciprocal space contribution to the rms error in the
 force (as described in the book of Hockney and Eastwood
 (Eqn. 8.23) (for a system of N randomly distributed particles in a
 cubic box).
 */
p3m_float ErrorEstimateIK::compute_ks_error_full(Parameters& p, p3m_int num_charges,
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

		p3m_int old_nx = 1000000;
		p3m_int old_ny = 1000000;
		p3m_float ctan_x, ctan_y;
		if (ny != old_ny) {
			if (nx != old_nx) {
				old_nx = nx;
				ctan_x = k_space_error_sum1(nx, grid_i[0], p.cao);
			}
			old_ny = ny;
			ctan_y = k_space_error_sum1(ny, grid_i[1], p.cao);
		}

		if (nx != 0 || ny != 0 || nz != 0) {
			p3m_float n2 = nx * nx + ny * ny + nz * nz;
			p3m_float cs = ctan_x * ctan_y
					* k_space_error_sum1(nz, grid_i[2], p.cao);
			p3m_float alias1, alias2;
			// TODO: sum2_ad for IKI?
			k_space_error_sum2(nx, ny, nz, p.grid, grid_i, p.cao, alpha_L_i,
					&alias1, &alias2);
			p3m_float d = alias1 - SQR(alias2 / cs) / n2;
			/* at high precisions, d can become negative due to extinction;
			 also, don't take values that have no significant digits left*/
			if (d > 0.0 && (fabs(d / alias1) > ROUND_ERROR_PREC))
				local_he_q += d;
		}
	}
	p3m_float he_q;
	MPI_Reduce(&local_he_q, &he_q, 1, P3M_MPI_FLOAT, MPI_SUM, 0, comm.mpicomm);

	p3m_float ks_error = 2.0 * sum_q2 *
			sqrt(he_q / (num_charges * box_l[0] * box_l[1] * box_l[2]));

	return ks_error;
}

/** One of the aliasing sums used by \ref p3m_fftcommon_k_space_error.
 (fortunately the one which is most important (because it converges
 most slowly, since it is not damped exponentially)) can be
 calculated analytically. The result (which depends on the order of
 the spline interpolation) can be written as an even trigonometric
 polynomial. The results are tabulated here (The employed formula
 is Eqn. 7.66 in the book of Hockney and Eastwood). */
p3m_float ErrorEstimateIK::k_space_error_sum1(p3m_int n, p3m_float grid_i, p3m_int cao) {
	p3m_float c, res = 0.0;
	c = SQR(cos(M_PI * grid_i * (p3m_float) n));

	switch (cao) {
	case 1: {
		res = 1;
		break;
	}
	case 2: {
		res = (1.0 + c * 2.0) / 3.0;
		break;
	}
	case 3: {
		res = (2.0 + c * (11.0 + c * 2.0)) / 15.0;
		break;
	}
	case 4: {
		res = (17.0 + c * (180.0 + c * (114.0 + c * 4.0))) / 315.0;
		break;
	}
	case 5: {
		res = (62.0 + c * (1072.0 + c * (1452.0 + c * (247.0 + c * 2.0))))
				/ 2835.0;
		break;
	}
	case 6: {
		res =
				(1382.0
						+ c
								* (35396.0
										+ c
												* (83021.0
														+ c
																* (34096.0
																		+ c
																				* (2026.0
																						+ c
																								* 4.0)))))
						/ 155925.0;
		break;
	}
	case 7: {
		res =
				(21844.0
						+ c
								* (776661.0
										+ c
												* (2801040.0
														+ c
																* (2123860.0
																		+ c
																				* (349500.0
																						+ c
																								* (8166.0
																										+ c
																												* 4.0))))))
						/ 6081075.0;
		break;
	}
	default:
		throw std::logic_error("INTERNAL_ERROR: This value for the interpolation order should not occur!");
	}

	return res;
}

/** aliasing sum used by \ref k_space_error. */
void
ErrorEstimateIK::k_space_error_sum2(p3m_int nx, p3m_int ny, p3m_int nz, p3m_int grid[3],
		p3m_float grid_i[3], p3m_int cao, p3m_float alpha_L_i,
		p3m_float *alias1, p3m_float *alias2) {
	p3m_float prefactor = SQR(M_PI * alpha_L_i);

	*alias1 = *alias2 = 0.0;
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
				*alias2 += U2 * ex * (nx * nmx + ny * nmy + nz * nmz) / nm2;
			}
		}
	}
}

//p3m_float ErrorEstimateIK::compute_ks_error_approx_triclinic(Parameters& p, p3m_int num_charges,
//		p3m_float sum_q2, p3m_float box_vectors[3][3], bool isTriclinic) {
//
//	p3m_float h = box_vectors[0][0] / p.grid[0];
//	p3m_float ha = h * p.alpha;
//
//	/* compute the sum in eq. 38 */
//	p3m_float sum;
//	switch (p.cao) {
//	case 1:
//		sum = 2. / 3.;
//		break;
//	case 2:
//		sum = 5. / 294. * pow(ha, 2) + 1. / 50.;
//		break;
//	case 3:
//		sum = 21. / 3872. * pow(ha, 4) + 7. / 1440. * pow(ha, 2) + 1. / 588.;
//		break;
//	case 4:
//		sum = 143. / 28800. * pow(ha, 6) + 7601. / 2271360. * pow(ha, 4)
//		+ 3. / 1936. * pow(ha, 2) + 1. / 4320.;
//		break;
//	case 5:
//		sum = 106640677. / 11737571328. * pow(ha, 8)
//		+ 517231. / 106536960. * pow(ha, 6) + 143. / 69120. * pow(ha, 4)
//		+ 7601. / 13628160. * pow(ha, 2) + 1. / 23232.;
//		break;
//	case 6:
//		sum = 326190917. / 11700633600. * pow(ha, 10)
//		+ 733191589. / 59609088000. * pow(ha, 8)
//		+ 9694607. / 2095994880. * pow(ha, 6)
//		+ 47021. / 35512320. * pow(ha, 4) + 13. / 57600. * pow(ha, 2)
//		+ 691. / 68140800.;
//		break;
//	case 7:
//		sum = 4887769399. / 37838389248. * pow(ha, 12)
//		+ 1755948832039. / 36229939200000. * pow(ha, 10)
//		+ 25091609. / 1560084480. * pow(ha, 8)
//		+ 56399353. / 12773376000. * pow(ha, 6)
//		+ 745739. / 838397952. * pow(ha, 4)
//		+ 3617. / 35512320. * pow(ha, 2) + 1. / 345600.;
//		break;
//	default:
//		throw std::logic_error("INTERNAL_ERROR: k_space_error_approx: "
//				"Charge assignment order should not occur!\n");
//	};
//
//	return sum_q2 / (box_vectors[0][0] * box_vectors[0][0])
//					* pow(h * p.alpha, p.cao)
//					* sqrt(p.alpha * box_vectors[0][0] / num_charges * sqrt(2.0 * M_PI) * sum);
//}

/** Calculates the reciprocal space contribution to the rms error in the
 force (as described in the book of Hockney and Eastwood
 (Eqn. 8.23) (for a system of N randomly distributed particles in a
 cubic box).
 */
p3m_float ErrorEstimateIK::compute_ks_error_full_triclinic(Parameters& p, p3m_int num_charges,
		p3m_float sum_q2,  p3m_float box_vectors[3][3], bool isTriclinic) {
	/* #ifdef P3M_ENABLE_DEBUG */
	/*   printf(  */
	/* 	  "        #%2d: N=%d Q2=%g box_l=(%g, %g, %g) grid=(%d, %d, %d) alpha=%g cao=%d\n",  */
	/* 	  comm_rank, N, sum_q2, box_l[0], box_l[1], box_l[2],  */
	/* 	  grid[0], grid[1], grid[2], alpha, cao); */
	/* #endif */

	p3m_float local_he_q = 0.0;
	p3m_float grid_i[3] =
			{ 1.0 / p.grid[0], 1.0 / p.grid[1], 1.0 / p.grid[2] };
	
	//p3m_float alpha_L_i = 1. / (p.alpha * box_vectors[0][0]);

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

		p3m_int old_nx = 1000000;
		p3m_int old_ny = 1000000;
		p3m_float ctan_x, ctan_y;
		if (ny != old_ny) {
			if (nx != old_nx) {
				old_nx = nx;
				ctan_x = k_space_error_sum1_triclinic(nx, grid_i[0], p.cao);
			}
			old_ny = ny;
			ctan_y = k_space_error_sum1_triclinic(ny, grid_i[1], p.cao);
		}

		if (nx != 0 || ny != 0 || nz != 0) {
			p3m_float n2 = nx * nx + ny * ny + nz * nz;
			p3m_float cs = ctan_x * ctan_y
					* k_space_error_sum1_triclinic(nz, grid_i[2], p.cao); //Prod_a k_space_error_sum1(n_a,1/grid_a, cao)
			p3m_float sum_Fref2, alias2;
                        
			k_space_error_sum2_triclinic(nx, ny, nz, p.grid, grid_i, p.cao, p.alpha,
					&sum_Fref2, &alias2, box_vectors, isTriclinic);
			p3m_float d = sum_Fref2 - SQR(alias2 / cs) / n2;
			/* at high precisions, d can become negative due to extinction;
			 also, don't take values that have no significant digits left*/
			if (d > 0.0 && (fabs(d / sum_Fref2) > ROUND_ERROR_PREC))
				local_he_q += d;
		}
	}
	p3m_float he_q;
	MPI_Reduce(&local_he_q, &he_q, 1, P3M_MPI_FLOAT, MPI_SUM, 0, comm.mpicomm);

	p3m_float ks_error = 2.0 * sum_q2 *
			sqrt(he_q / num_charges)/ (box_vectors[0][0] * box_vectors[1][1] * box_vectors[2][2]);

	return ks_error;
}

/** One of the aliasing sums used by \ref p3m_fftcommon_k_space_error.
 (fortunately the one which is most important (because it converges
 most slowly, since it is not damped exponentially)) can be
 calculated analytically. The result (which depends on the order of
 the spline interpolation) can be written as an even trigonometric
 polynomial. The results are tabulated here (The employed formula
 is Eqn. 7.66 in the book of Hockney and Eastwood). */ // (-1)^s/s! * d^s/dx^s cot(x) = sum_-inf to inf (x-pi n)^-(s+1)
p3m_float ErrorEstimateIK::k_space_error_sum1_triclinic(p3m_int n, p3m_float grid_i, p3m_int cao) {
	p3m_float c, res = 0.0;
	c = SQR(cos(M_PI * grid_i * (p3m_float) n)); // cos(pi/grid*n)^2

	switch (cao) {
	case 1: {
		res = 1;
		break;
	}
	case 2: {
		res = (1.0 + c * 2.0) / 3.0;
		break;
	}
	case 3: {
		res = (2.0 + c * (11.0 + c * 2.0)) / 15.0;
		break;
	}
	case 4: {
		res = (17.0 + c * (180.0 + c * (114.0 + c * 4.0))) / 315.0;
		break;
	}
	case 5: {
		res = (62.0 + c * (1072.0 + c * (1452.0 + c * (247.0 + c * 2.0))))
				/ 2835.0;
		break;
	}
	case 6: {
		res =
				(1382.0
						+ c
								* (35396.0
										+ c
												* (83021.0
														+ c
																* (34096.0
																		+ c
																				* (2026.0
																						+ c
																								* 4.0)))))
						/ 155925.0;
		break;
	}
	case 7: {
		res =
				(21844.0
						+ c
								* (776661.0
										+ c
												* (2801040.0
														+ c
																* (2123860.0
																		+ c
																				* (349500.0
																						+ c
																								* (8166.0
																										+ c
																												* 4.0))))))
						/ 6081075.0;
		break;
	}
	default:
		throw std::logic_error("INTERNAL_ERROR: This value for the interpolation order should not occur!");
	}

	return res;
}

/** aliasing sum used by \ref k_space_error. */
void
ErrorEstimateIK::k_space_error_sum2_triclinic(p3m_int nx, p3m_int ny, p3m_int nz, p3m_int grid[3],
		p3m_float grid_i[3], p3m_int cao, p3m_float alpha,
		p3m_float *sum_Fref2, p3m_float *alias2, p3m_float box_vectors[3][3], bool isTriclinic) {
	p3m_float prefactor = SQR(M_PI / alpha);

	*sum_Fref2 = *alias2 = 0.0;
	for (p3m_int mx = -P3M_BRILLOUIN; mx <= P3M_BRILLOUIN; mx++) {
		p3m_float nmx = nx + mx * grid[0];
		p3m_float fnmx = grid_i[0] * nmx;
		for (p3m_int my = -P3M_BRILLOUIN; my <= P3M_BRILLOUIN; my++) {
			p3m_float nmy = ny + my * grid[1];
			p3m_float fnmy = grid_i[1] * nmy;
			for (p3m_int mz = -P3M_BRILLOUIN; mz <= P3M_BRILLOUIN; mz++) {
				p3m_float nmz = nz + mz * grid[2];
				p3m_float fnmz = grid_i[2] * nmz;

				p3m_float k2_4pi2 =
                                            SQR(nmx / box_vectors[RX][RX]) +
                                            SQR(nmy / box_vectors[RY][RY]) +
                                            SQR(nmz / box_vectors[RZ][RZ]);
				p3m_float ex = exp(-prefactor * k2_4pi2); //exp(-pi^2/alpha^2 * nm2)

				p3m_float U2 = pow(sinc(fnmx) * sinc(fnmy) * sinc(fnmz),
						2.0 * cao);

				*sum_Fref2 += ex * ex / k2_4pi2;
				*alias2 += U2 * ex * (nx * nmx + ny * nmy + nz * nmz) / k2_4pi2;
			}
		}
	}
}

}


