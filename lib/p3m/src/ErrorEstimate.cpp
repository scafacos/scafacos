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
#include "IK/ErrorEstimate.hpp"
#include "ADI/ErrorEstimate.hpp"
#include "utils.hpp"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

namespace P3M {

/** Factory method to create ErrorEstimates */
ErrorEstimate*
ErrorEstimate::create(Communication &comm) {
#if defined(P3M_AD) && defined(P3M_INTERLACE)
	return new ADI::ErrorEstimate(comm);
#elif defined(P3M_IK) && !defined(P3M_INTERLACE)
	return new IK::ErrorEstimate(comm);
#endif
}

ErrorEstimate::CantGetRequiredAccuracy::CantGetRequiredAccuracy() :
	std::logic_error("Cannot achieve the required accuracy.") {}

void ErrorEstimate::computeAlpha(p3m_float required_accuracy, TuneParameters& p,
		p3m_int num_charges, p3m_float sum_q2, p3m_float box_l[3]) {
	/* Get the real space error for alpha=0 */
	p.alpha = 0.0;
	computeRSError(p, num_charges, sum_q2, box_l);
	p3m_float max_rs_error = p.rs_error;

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
		p3m_float sum_q2, p3m_float box_l[3], p3m_float box_vectors[3][3], bool isTriclinic) {
	if (comm.onMaster()) {
        computeRSError(p, num_charges, sum_q2, box_l);
        computeKSError(p, num_charges, sum_q2, box_l);
        //printf("ks error-orig %e\n",p.ks_error);
        //computeKSError_triclinic(p, num_charges, sum_q2, box_vectors, isTriclinic);
        // printf("ks error-tric %e\n",p.ks_error);
        p.error = sqrt(SQR(p.rs_error) + SQR(p.ks_error));

#ifdef P3M_ENABLE_DEBUG
        printf("        error estimate: rs_err=" FFLOATE ", "
                "ks_err=" FFLOATE ", err=" FFLOATE "\n", p.rs_error, p.ks_error, p.error);
#endif
    } else
        computeKSError(p, num_charges, sum_q2, box_l);
}

void
ErrorEstimate::computeMaster(TuneParameters &p,
        p3m_int num_charges, p3m_float sum_q2, p3m_float box_l[3],
        p3m_float box_vectors[3][3], bool isTriclinic) {
    if (!comm.onMaster())
        throw std::logic_error("Do not call ErrorEstimate::computeMaster() on slave.");

    bool finished = false;
    MPI_Bcast(&finished, 1, MPI_C_BOOL,
            Communication::MPI_MASTER, comm.mpicomm);

    // broadcast parameters
    p3m_int int_buffer[6];
    p3m_float float_buffer[8];

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
    this->compute(p, num_charges, sum_q2, box_l, box_vectors, isTriclinic);
}

void
ErrorEstimate::endLoop() {
    if (!comm.onMaster())
        throw std::logic_error("Do not call ErrorEstimate::endLoop() on slave.");
    bool finished = true;
    MPI_Bcast(&finished, 1, MPI_C_BOOL, Communication::MPI_MASTER, comm.mpicomm);
}

void ErrorEstimate::loopSlave() {
    if (comm.onMaster())
        throw std::logic_error("Do not call ErrorEstimate::loopSlave() on master.");

    bool finished;
    while (true) {
        P3M_DEBUG_LOCAL(printf( \
                "      %2d: ErrorEstimate::loopSlave() waiting.\n", \
                comm.rank));
        // check whether to continue in slave loop
        MPI_Bcast(&finished, 1, MPI_C_BOOL,
                Communication::MPI_MASTER, comm.mpicomm);
        P3M_DEBUG_LOCAL(\
                printf("      %2d: ErrorEstimate::loopSlave() got %s.\n", \
                comm.rank, (finished ? "FINISH":"ERROR ESTIMATE")));
        if (finished) break;

        // receive parameters
        p3m_int int_buffer[6];
        p3m_float float_buffer[8];

        TuneParameters p;
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
        box_vectors[0][0] = float_buffer[2];
        box_vectors[0][1] = 0.0;
        box_vectors[0][2] = 0.0;
        box_vectors[1][0] = float_buffer[7];
        box_vectors[1][1] = float_buffer[3];
        box_vectors[1][2] = 0.0;
        box_vectors[2][0] = float_buffer[6];
        box_vectors[2][1] = float_buffer[5];
        box_vectors[2][2] = float_buffer[4];

        // run slave job
        this->compute(p, num_charges, sum_q2, box_l,box_vectors, isTriclinic);
    }
}

void ErrorEstimate::computeRSError(TuneParameters& p, p3m_int num_charges,
        p3m_float sum_q2, p3m_float box_l[3]) {
    p.rs_error = (2.0 * sum_q2 * exp(-SQR(p.r_cut * p.alpha)))
            / (sqrt(num_charges * p.r_cut * box_l[0] * box_l[1] * box_l[2]));
}

}
