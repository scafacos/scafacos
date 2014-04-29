/*
 Copyright (C) 2014 Olaf Lenz, Gabriel Sichardt

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
#ifndef _P3M_ERRORESTIMATEADI_HPP
#define _P3M_ERRORESTIMATEADI_HPP

#include "ErrorEstimate.hpp"

namespace P3M {
/** Estimate the errors in the IK differentiated P3M Algorithm */
class ErrorEstimateADI : public ErrorEstimate {
public:
    ErrorEstimateADI(Communication &comm) : ErrorEstimate(comm) {}

	virtual p3m_float compute_ks_error(Parameters &p,
			p3m_int num_charges, p3m_float sum_q2, p3m_float box_l[3]);
        virtual p3m_float compute_ks_error_triclinic(Parameters& p, p3m_int num_charges,
		p3m_float sum_q2, p3m_float box_vectors[3][3], bool isTriclinic);

protected:
    void
	k_space_error_sum2(p3m_int nx, p3m_int ny, p3m_int nz, p3m_int grid[3],
			p3m_float grid_i[3], p3m_int cao, p3m_float alpha_L_i,
			p3m_float *alias1, p3m_float *alias2, p3m_float *alias3,
			p3m_float *alias4, p3m_float *alias5, p3m_float *alias6);
    void k_space_error_sum2_triclinic(p3m_int nx, p3m_int ny, p3m_int nz, p3m_int grid[3],
			p3m_float grid_i[3], p3m_int cao, p3m_float alpha,p3m_float *alias1, p3m_float *alias2, p3m_float *alias3,
			p3m_float *alias4, p3m_float *alias5, p3m_float *alias6, p3m_float box_vectors[3][3], bool isTriclinic);
};
}

#endif
