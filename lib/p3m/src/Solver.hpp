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
#ifndef _P3M_SOLVER_HPP
#define _P3M_SOLVER_HPP

#include "types.hpp"
#include "Communication.hpp"
#include "Parallel3DFFT.hpp"
#include "ErrorEstimate.hpp"
#include "CAF.hpp"
#include "common/gridsort/gridsort.h"
#include <list>


namespace P3M {

/** Structure that holds all data of the P3M algorithm */
class Solver {
public:
    Solver(MPI_Comm mpicomm);
    ~Solver();

    void prepare();

    void tune(p3m_int num_particles, p3m_float *positions, p3m_float *charges);

    void run(p3m_int num_particles, p3m_float *positions, p3m_float *charges,
            p3m_float *fields, p3m_float *potentials);

    /** Test run the method with the current parameters. Afterwards, the
     timing variables in the data struct are set. */
    void timing(p3m_int num_particles, p3m_float *positions, p3m_float *charges);

    Communication comm;
    Parallel3DFFT fft;
    ErrorEstimate *errorEstimate;

    /****************************************************
     * SYSTEM PARAMETERS
     ****************************************************/
    /* System size in x,y,z */
    p3m_float box_l[3];
    /* Skin */
    p3m_float skin;
    /* the complete box vectors */
    p3m_float box_vectors[3][3];
    /* Volume of the box */
    p3m_float volume;
    /****************************************************
     * PARAMETERS OF THE METHOD
     ****************************************************/
    /** Tolerance in the field rms. */
    p3m_float tolerance_field;
    /** number of interpolation points for charge assignment function */
    p3m_int n_interpol;
    /** whether to compute the near field in the method */
    bool near_field_flag;
    /* whether the box is triclinic */
    bool isTriclinic; 
    /** flag that determines if the Gaussian potentials are shifted */
    bool shiftGaussians;
    
    /* TUNABLE PARAMETERS */
    /** cutoff radius */
    p3m_float r_cut;
    /** Ewald splitting parameter */
    p3m_float alpha;
    /** number of grid points per coordinate direction (>0). */
    p3m_int grid[3];
    /** charge assignment order ([0,P3M_MAX_CAO]). */
    p3m_int cao;

    /* Whether or not it is necessary to retune the method before running it. */
    bool needs_retune;
    /** Whether or not rcut is to be automatically tuned. */
    bool tune_r_cut;
    /** Whether or not alpha is to be automatically tuned. */
    bool tune_alpha;
    /** Whether or not the grid is to be automatically tuned. */
    bool tune_grid;
    /** Whether or not the charge assignment order is to be automatically tuned. */
    bool tune_cao;

    /****************************************************
     * FLAGS TO TURN ON/OFF COMPUTATION OF DIFFERENT COMPONENTS
     ****************************************************/
    /** Whether or not the total energy is to be computed. */
    bool require_total_energy;
    /** The total energy. */
    p3m_float total_energy;
    /** Whether and how timings are to be taken.
     * 0: run without timings
     * 1: partial timing to estimate without results
     * 2: full means all timings and correct results*/
    timingEnum require_timings;

    static const int TIMING = 0;
    static const int TIMING_NEAR = 1;
    static const int TIMING_FAR = 2;
    static const int TIMING_CA = 3;
    static const int TIMING_GATHER = 4;
    static const int TIMING_FORWARD = 5;
    static const int TIMING_BACK = 6;
    static const int TIMING_INFLUENCE = 7;
    static const int TIMING_SPREAD = 8;
    static const int TIMING_POTENTIALS = 9;
    static const int TIMING_FIELDS = 10;
    static const int TIMING_DECOMP = 11;
    static const int TIMING_COMP = 12;
    static const int NUM_TIMINGS = 13;
    double timings[NUM_TIMINGS];

    /****************************************************
     * DERIVED PARAMETERS
     ****************************************************/
    /** The errors of the method according to the error formula. */
    p3m_float error, ks_error, rs_error;

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

    /****************************************************
     * METHOD DATA
     ****************************************************/

    /** local grid. */
    local_grid_t local_grid;
    /** real space grid (local) for CA/FFT.*/
    p3m_float *rs_grid;
    /** k space grid (local) for k space calculation and FFT.*/
    p3m_float *ks_grid;

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
    
    /** potential shift required to make real space potentials continuous */
    p3m_float potentialOffset;

private:
    // submethods of prepare()
    void prepareSendGrid();
    void prepareLocalCAGrid();
    void prepareDifferentialOperator();
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
    p3m_int *computeGridShift(int dir, p3m_int size);

    // submethods of run()
    
    /* conversion of cartesian positions to triclinic positions */
    void
    calculateTriclinicPositions(p3m_float *positions,
            p3m_float *triclinic_positions, p3m_int number);  
    
    /* domain decomposition */
    void
    decompose(fcs_gridsort_t *gridsort,
            p3m_int _num_particles,
            p3m_float *_positions, p3m_float *_charges,
            p3m_int *num_real_particles,
            p3m_float **positions, p3m_float **charges,
            fcs_gridsort_index_t **indices,
            p3m_int *num_ghost_particles,
            p3m_float **ghost_positions, p3m_float **ghost_charges,
            fcs_gridsort_index_t **ghost_indices);

    void computeFarADI(p3m_int num_charges, p3m_float* positions, p3m_float* charges,
            p3m_float* fields, p3m_float* potentials);
    void computeFarIK(p3m_int num_charges, p3m_float* positions, p3m_float* charges,
            p3m_float* fields, p3m_float* potentials);

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

    void collectPrintTimings();

    /* Assign the fields to the positions in dimension dim [IK] */
    void
    assignFieldsIK(p3m_float *data,
            p3m_int dim, p3m_int num_particles, p3m_float* positions,
            p3m_int shifted, p3m_float* fields);

    /* Backinterpolate the forces obtained from k-space to the positions [AD]*/
    void
    assignFieldsAD(p3m_float *data,
            p3m_int num_real_particles, p3m_float* positions,
            p3m_int shifted, p3m_float* fields);
    
    void
    cartesianizeFields(p3m_float *fields, p3m_int num_particles);

    // submethods of tune()
    struct TuneParameters {
        Parameters p;
        /** Errors */
        p3m_float rs_error, ks_error, error;
        /** Timings */
        p3m_float timing, timing_near, timing_far;
    };
    typedef std::list<TuneParameters*> TuneParametersList;

    TuneParameters* tuneFar(p3m_int num_particles, p3m_float *positions,
            p3m_float *charges, p3m_float r_cut);

    /** Slave variant of tune_far. */
    void tuneFar(p3m_int num_particles,
            p3m_float *positions, p3m_float *charges);

    void tuneAlpha(TuneParametersList &params_to_try);
    void tuneCAO(TuneParametersList &params_to_try);
    void tuneGrid(TuneParametersList &params_to_try);

    TuneParameters* timeParams(p3m_int num_particles, p3m_float *positions,
            p3m_float *charges, TuneParametersList &params_to_try);

    void countCharges(p3m_int num_particles, p3m_float *charges);

    /* Events during tuning */
    static const int CMD_FAILED = -1;
    static const int CMD_FINISHED = 0;
    static const int CMD_COMPUTE_ERROR_ESTIMATE = 1;
    static const int CMD_TIMING = 2;

    void tuneBroadcastParams(Parameters &p);
    void tuneReceiveParams();

    void tuneBroadcastCommand(p3m_int command);
    void tuneBroadcastCommand(p3m_int command, Parameters &p);

    void tuneBroadcastSlave(p3m_int num_particles,
            p3m_float *positions, p3m_float *charges);

};


}

#endif
