/*
 Copyright (C) 2011,2012,2013,2014 Olaf Lenz
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

#ifndef _P3M_FARSOLVER_HPP
#define _P3M_FARSOLVER_HPP

#include "types.hpp"
#include "Communication.hpp"
#include "Parallel3DFFT.hpp"
#include "ErrorEstimate.hpp"
#include "CAF.hpp"

namespace P3M {

class FarSolver {
public:
    enum TimingComponent {
        TOTAL, CA, GATHER, FORWARD, BACK, INFLUENCE, SPREAD, POTENTIALS, FIELDS
    };
    static const int NUM_TIMINGS = 9;

    /* enumeration to specify the type of timings */
    enum TimingType {
        NONE, // do not measure timings
        ESTIMATE_ALL, // estimate all components that you can
        ESTIMATE_FFT, // estimate the timing of the FFT
        ESTIMATE_ASSIGNMENT, // estimate the timing of the assignment
        FULL // time all components
    };

    FarSolver(Communication &comm, p3m_float box_l[3],
            p3m_float r_cut, p3m_float alpha, p3m_int grid[3], p3m_int cao, p3m_float box_vectors[3][3], p3m_float volume, bool isTriclinic);
    virtual ~FarSolver();

    void runADI(p3m_int num_charges, p3m_float *positions, p3m_float *charges,
                p3m_float *fields, p3m_float *potentials);
    void runIK(p3m_int num_charges, p3m_float *positions, p3m_float *charges,
                p3m_float *fields, p3m_float *potentials);

    void setRequireTotalEnergy(bool flag = true);
    fcs_float getTotalEnergy();

    void setRequireTimings(TimingType type = NONE);
    /** Test run the method with the current parameters.
     * Return the total run time. */
    const double* measureTimings(p3m_int num_particles,
            p3m_float *positions, p3m_float *charges);
    /** Fetch the timings. */
    const double* getTimings();
protected:
    Communication &comm;
    Parallel3DFFT fft;
    ErrorEstimate *errorEstimate;

    /****************************************************
     * SYSTEM PARAMETERS
     ****************************************************/
    /* System size in x,y,z */
    p3m_float box_l[3];
    /* the complete box vectors */
    p3m_float box_vectors[3][3];
    /* Volume of the box */
    p3m_float volume;
    /* whether the box is triclinic */
    bool isTriclinic; 
    

    /****************************************************
     * PARAMETERS OF THE METHOD
     ****************************************************/
    /* Skin */
    p3m_float skin;
    /** number of interpolation points for charge assignment function */
    p3m_int n_interpol;

    /** cutoff radius */
    p3m_float r_cut;
    /** Ewald splitting parameter */
    p3m_float alpha;
    /** number of grid points per coordinate direction (>0). */
    p3m_int grid[3];
    /** charge assignment order ([0,P3M_MAX_CAO]). */
    p3m_int cao;

    /****************************************************
     * FLAGS TO TURN ON/OFF COMPUTATION OF DIFFERENT COMPONENTS
     ****************************************************/
    /** Whether or not the total energy is to be computed. */
    bool require_total_energy;
    /** The total energy. */
    p3m_float total_energy;
    /** Whether and how timings are to be taken. */
    TimingType require_timings;

    double timings[NUM_TIMINGS];

    /****************************************************
     * METHOD DATA
     ****************************************************/
    /** The number of nodes in each spatial dimension. */
    p3m_int node_grid[3];
    /** position of this node in the node grid */
    p3m_int node_pos[3];
    /** the six nearest neighbors of a node in the node grid. */
    p3m_int node_neighbors[6];

    /** Size of the local box. */
    p3m_float local_box_l[3];
    /** Left (bottom, front) corner of this nodes local box. */
    p3m_float my_left[3];
    /** Right (top, back) corner of this nodes local box. */
    p3m_float my_right[3];

    /** offset of the first grid point (lower left corner) from the
     coordinate origin ([0,1[). */
    p3m_float grid_off[3];

    /** Cutoff for charge assignment. */
    p3m_float cao_cut[3];
    /** grid constant. */
    p3m_float a[3];
    /** inverse grid constant. */
    p3m_float ai[3];
    /** additional points around the charge assignment grid, for method like dielectric ELC
     creating virtual charges. */
    p3m_float additional_grid[3];

    /** local grid. */
    local_grid_t local_grid;
    /** real space grid (local) for CA/FFT. */
    p3m_float *rs_grid;
    /** k space grid (local) for k space calculation and FFT. */
    p3m_float *ks_grid;
    /** additional grid for buffering. */
    p3m_float *buffer;

    /** number of charged particles */
    p3m_int sum_qpart;
    /** Sum of square of charges */
    p3m_float sum_q2;
    /** square of sum of charges */
    p3m_float square_sum_q;

    /** charge assignment function. */
    CAF *caf;
    CAF::Cache *cafx;
    CAF::Cache *cafy;
    CAF::Cache *cafz;
    /** gradient of charge assignment function */
    CAF *caf_d;
    CAF::Cache *cafx_d;
    CAF::Cache *cafy_d;
    CAF::Cache *cafz_d;

    /** position shift for calc. of first assignment grid point. */
    p3m_float pos_shift;

    /** Spatial differential operator in k-space. We use an i*k differentiation. */
    p3m_int *d_op[3];
    /** Force optimized influence function (k-space) */
    p3m_float *g_force;
    /** Energy optimized influence function (k-space) */
    p3m_float *g_energy;

    /** number of permutations in k_space */
    p3m_int ks_pnum;

    /** send/recv grid sizes */
    send_grid_t sm;

    /** Field to store grid points to send. */
    p3m_float *send_grid;
    /** Field to store grid points to recv */
    p3m_float *recv_grid;

    void prepareSendGrid();
    void prepareLocalCAGrid();
    void printLocalGrid();
    void printSendGrid();

    void computeInfluenceFunctionIK();
    void computeInfluenceFunctionADI();
    void computeInfluenceFunctionIKI();
    void performAliasingSumsIK(
            p3m_int nmx0, p3m_int nmy0, p3m_int nmz0,
            p3m_float numerator_force[3],
            p3m_float &numerator_energy,
            p3m_float &denominator);
    void performAliasingSumsADI(
            p3m_int nmx0, p3m_int nmy0, p3m_int nmz0,
            p3m_float &numerator_force,
            p3m_float &numerator_energy, p3m_float denominator[4]);
    void performAliasingSumsIKI(
            p3m_int nmx0, p3m_int nmy0, p3m_int nmz0,
            p3m_float numerator_force[3],
            p3m_float &numerator_energy,
            p3m_float denominator[2]);
    void cartesianizeFields(p3m_float *fields, p3m_int num_particles);
    p3m_int *computeGridShift(int dir, p3m_int size);

    /* charge assignment */
    p3m_int getCAPoints(p3m_float real_pos[3], p3m_int shifted);
    void assignCharges(p3m_float *data,
            p3m_int num_charges, p3m_float *positions, p3m_float *charges, p3m_int shifted);

    /* collect grid from neighbor processes */
    void gatherGrid(p3m_float* rs_grid);
    /* spread grid to neighbor processors */
    void spreadGrid(p3m_float* rs_grid);

    /* apply influence function g to the array in, yielding out */
    void applyInfluenceFunction(p3m_float *in, p3m_float *out, p3m_float *g);

    /* differentiate in in direction dim, yielding out*/
    void differentiateIK(int dim, p3m_float* in, p3m_float* out);

    /* compute the total energy (in k-space, so no backtransform) */
    p3m_float computeTotalEnergy();
    /* assign the potentials to the positions */
    void
    assignPotentials(p3m_float *data,
            p3m_int num_particles, p3m_float* positions, p3m_float* charges,
            p3m_int shifted, p3m_float* potentials);

    /* Assign the fields to the positions in dimension dim [IK] */
    void
    assignFieldsIK(p3m_float *data,
            p3m_int dim, p3m_int num_particles, p3m_float* positions,
            p3m_int shifted, p3m_float* fields);

    /* Backinterpolate the forces obtained from k-space to the positions [AD]*/
    void
    assignFieldsAD(p3m_float *data,
            p3m_int num_particles, p3m_float* positions,
            p3m_int shifted, p3m_float* fields);

    void countCharges(p3m_int num_particles, p3m_float *charges);

    void gatherTimings();

    void resetTimers() {
        for (int i=0; i < NUM_TIMINGS; i++)
            timings[i] = 0.0;
    }

    void startTimer(TimingComponent comp) {
        if (require_timings != NONE)
            timings[comp] += -MPI_Wtime();
    }

    void stopTimer(TimingComponent comp) {
        if (require_timings != NONE)
            timings[comp] += MPI_Wtime();
    }

    void switchTimer(TimingComponent comp1,
            TimingComponent comp2) {
        if (require_timings != NONE) {
            timings[comp1] += MPI_Wtime();
            timings[comp2] += -MPI_Wtime();
        }
    }
};


}



#endif /* FARSOLVER_HPP_ */
