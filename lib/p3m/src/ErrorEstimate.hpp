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
#ifndef _P3M_ERRORESTIMATE_HPP
#define _P3M_ERRORESTIMATE_HPP

#include "p3mconfig.hpp"
#include "types.hpp"
#include "Communication.hpp"
#include <stdexcept>

namespace P3M {
/* Base class of the P3M error estimates. */
class ErrorEstimate {
public:
    static ErrorEstimate *create(Communication &comm);

    class CantGetRequiredAccuracy: public std::logic_error {
    public:
        CantGetRequiredAccuracy();
    };

    ErrorEstimate(Communication &comm)
            : comm(comm) {
    }
    virtual ~ErrorEstimate() {
    }

    /** Determines a value for alpha that achieves the required_accuracy in
     * the near field. It still needs to be checked whether the far field
     * can achieve the required accuracy.
     */
    virtual void
    computeAlpha(p3m_float required_accuracy, TuneParameters &p,
            p3m_int num_charges, p3m_float sum_q2, p3m_float box_l[3]);

    /** Master variant of the error computation. It first broadcasts the
     * parameters to all tasks, then runs compute. */
    virtual void
    computeMaster(TuneParameters &p,
            p3m_int num_charges, p3m_float sum_q2, p3m_float box_l[3]);

    /* Call this on the master to end the slave loop. */
    void endLoop();

    /** Slave variant of the error computation. It runs compute until
     * finishedMaster is called. */
    virtual void loopSlave();

    /** Computes the error estimate. When called in parallel, the result is
     * undefined on the slaves. */
    virtual void
    compute(TuneParameters &p,
            p3m_int num_charges, p3m_float sum_q2, p3m_float box_l[3],
            p3m_float box_vectors[3][3], bool isTriclinic);

    /** Calculates the real space contribution to the rms error in the
     * force (as described by Kolafa and Perram). When called in parallel, the
     * result is undefined on the slaves.
     */
    virtual void computeRSError(TuneParameters &p,
            p3m_int num_charges, p3m_float sum_q2, p3m_float box_l[3]);

    /** Calculates the reciprocal space contribution to the rms error in the
     * force. When called in parallel, the result is undefined on the slaves. */
    virtual void computeKSError(TuneParameters &p,
            p3m_int num_charges, p3m_float sum_q2, p3m_float box_l[3]) = 0;
/** Calculates the reciprocal space contribution to the rms error in the
     force.
     */
    virtual p3m_float computeKSError_triclinic(TuneParameters &p,
            p3m_int num_charges, p3m_float sum_q2, p3m_float box_vectors[3][3], bool isTriclinic) = 0;
	protected:
	    Communication &comm;

};
}
#endif
