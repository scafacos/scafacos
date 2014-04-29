/*
 Copyright (C) 2014 Olaf Lenz, Gabriel Sichardt
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
#include "ErrorEstimate.hpp"
#include "ErrorEstimateIK.hpp"
#include "ErrorEstimateADI.hpp"
#include "utils.hpp"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

namespace P3M {

/** Factory method to create ErrorEstimates */
ErrorEstimate*
ErrorEstimate::create(Communication &comm) {
#if defined(P3M_AD) && defined(P3M_INTERLACE)
	return new ErrorEstimateADI(comm);
#elif defined(P3M_IK) && !defined(P3M_INTERLACE)
	return new ErrorEstimateIK(comm);
#endif
}

ErrorEstimate::CantGetRequiredAccuracy::CantGetRequiredAccuracy() :
	std::logic_error("Cannot achieve the required accuracy.") {}

void ErrorEstimate::compute_alpha(p3m_float required_accuracy, Parameters& p,
		p3m_int num_charges, p3m_float sum_q2, p3m_float box_l[3]) {
	/* Get the real space error for alpha=0 */
	p.alpha = 0.0;
	p3m_float max_rs_error = compute_rs_error(p, num_charges, sum_q2, box_l);

	/* We know how the real space error behaves, so we can compute the
	 alpha where the real space error is half of the wanted
	 error. This is the alpha that we return. */
	if (M_SQRT2 * max_rs_error > required_accuracy) {
		p.alpha = sqrt(log(M_SQRT2 * max_rs_error / required_accuracy)) / p.r_cut;
	} else {
		/* if the error is small enough even for alpha=0 */
		p.alpha = 0.1 * box_l[0];
	}
}

void ErrorEstimate::compute(Parameters& p, p3m_int num_charges,
		p3m_float sum_q2, p3m_float box_l[3],
		p3m_float &error, p3m_float &rs_error, p3m_float &ks_error, p3m_float box_vectors[3][3], bool isTriclinic) {
	rs_error = compute_rs_error(p, num_charges, sum_q2, box_l);
	//ks_error = compute_ks_error(p, num_charges, sum_q2, box_l);
        //printf("ks error-orig %e\n",ks_error);
        ks_error = compute_ks_error_triclinic(p, num_charges, sum_q2, box_vectors, isTriclinic);
       // printf("ks error-tric %e\n",ks_error);
	error = sqrt(SQR(rs_error) + SQR(ks_error));

#ifdef P3M_ENABLE_DEBUG
	if (comm.onMaster())
		printf("        error estimate: rs_err=" FFLOATE ", "
		"ks_err=" FFLOATE ", err=" FFLOATE "\n", rs_error, ks_error, error);
#endif
}

p3m_float ErrorEstimate::compute(Parameters& p, p3m_int num_charges,
		p3m_float sum_q2, p3m_float box_l[3], p3m_float box_vectors[3][3], bool isTriclinic) {
	p3m_float ks_error, rs_error, error;
	this->compute(p, num_charges, sum_q2, box_l, error, rs_error, ks_error,box_vectors, isTriclinic);
	return error;
}

p3m_float ErrorEstimate::compute_master(Parameters &p,
        p3m_int num_charges, p3m_float sum_q2, p3m_float box_l[3], p3m_float box_vectors[3][3], bool isTriclinic) {
    p3m_float ks_error, rs_error, error;
    this->compute_master(p, num_charges, sum_q2, box_l, error, rs_error, ks_error, box_vectors, isTriclinic);
    return error;
}

void
ErrorEstimate::compute_master(Parameters &p,
        p3m_int num_charges, p3m_float sum_q2, p3m_float box_l[3],
        p3m_float &error, p3m_float &rs_error, p3m_float &ks_error, p3m_float box_vectors[3][3], bool isTriclinic) {
    if (!comm.onMaster())
        throw std::logic_error("Do not call ErrorEstimate::compute_master() on slave.");

    // broadcast parameters
    p3m_int int_buffer[6];
    p3m_float float_buffer[5];

    // pack int data
    int_buffer[0] = p.cao;
    int_buffer[1] = p.grid[0];
    int_buffer[2] = p.grid[1];
    int_buffer[3] = p.grid[2];
    int_buffer[4] = num_charges;
    int_buffer[5] = isTriclinic?1:0;
    MPI_Bcast(int_buffer, 6, P3M_MPI_INT,
            Communication::MPI_MASTER, comm.mpicomm);

    // pack float data
    float_buffer[0] = p.alpha;
    float_buffer[1] = sum_q2;
    float_buffer[2] = box_l[0];
    float_buffer[3] = box_l[1];
    float_buffer[4] = box_l[2];
    float_buffer[5] = box_vectors[2][1];
    float_buffer[6] = box_vectors[2][0];
    float_buffer[7] = box_vectors[1][0];
    MPI_Bcast(float_buffer, 8, P3M_MPI_FLOAT,
            Communication::MPI_MASTER, comm.mpicomm);
    
    // run master job
    this->compute(p, num_charges, sum_q2, box_l, error, rs_error, ks_error,box_vectors, isTriclinic);
}

void ErrorEstimate::compute_slave() {
    if (comm.onMaster())
        throw std::logic_error("Do not call ErrorEstimate::compute_slave() on master.");
    // receive parameters
    p3m_int int_buffer[6];
    p3m_float float_buffer[5];

    Parameters p;
    p3m_float box_l[3];
    p3m_int num_charges;
    p3m_float sum_q2;
    p3m_float box_vectors[3][3];
    
    bool isTriclinic;

    MPI_Bcast(int_buffer, 6, P3M_MPI_INT,
            Communication::MPI_MASTER, comm.mpicomm);
    p.cao = int_buffer[0];
    p.grid[0] = int_buffer[1];
    p.grid[1] = int_buffer[2];
    p.grid[2] = int_buffer[3];
    num_charges = int_buffer[4];
    isTriclinic = (int_buffer[5]==1)?true:false;

    // unpack float data
    MPI_Bcast(float_buffer, 8, P3M_MPI_FLOAT,
            Communication::MPI_MASTER, comm.mpicomm);
    p.alpha = float_buffer[0];
    sum_q2 = float_buffer[1];
    box_l[0] = float_buffer[2];
    box_l[1] = float_buffer[3];
    box_l[2] = float_buffer[4];
    box_vectors[0][0]=float_buffer[2];box_vectors[0][1]=0.0;box_vectors[0][2]=0.0;
    box_vectors[1][0]=float_buffer[7];box_vectors[1][1]=float_buffer[3];box_vectors[1][2]=0.0;
    box_vectors[2][0]=float_buffer[6];box_vectors[2][1]=float_buffer[5];box_vectors[2][2]=float_buffer[4];
    // run slave job
    this->compute(p, num_charges, sum_q2, box_l,box_vectors, isTriclinic);
}

p3m_float ErrorEstimate::compute_rs_error(Parameters& p, p3m_int num_charges,
        p3m_float sum_q2, p3m_float box_l[3]) {
    return (2.0 * sum_q2 * exp(-SQR(p.r_cut * p.alpha)))
            / (sqrt(num_charges * p.r_cut * box_l[0] * box_l[1] * box_l[2]));
}

}
