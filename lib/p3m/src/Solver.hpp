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
#include "ErrorEstimate.hpp"
#include "FarSolver.hpp"
#include "CAF.hpp"
#include "common/gridsort/gridsort.h"
#include <list>


namespace P3M {

  /** Structure that holds all data of the P3M algorithm */
  class Solver {
  public:
    enum TimingComponent {
      TOTAL, COMP, DECOMP, NEAR,
      FAR, CA, GATHER, FORWARD, BACK, INFLUENCE, SPREAD, POTENTIALS, FIELDS
    };
    static const int NUM_TIMINGS_NOTFAR = 4;
    static const int NUM_TIMINGS = 13;
    enum TimingType {
      NONE, // do not time anything
      NOTFAR, // do not time the far field
      FULL // time everything
    };
    
    Solver(MPI_Comm mpicomm);
    ~Solver();

    void prepare();

    void tune(p3m_int num_particles, p3m_float *positions, p3m_float *charges);

    void run(p3m_int num_particles, p3m_float *positions, p3m_float *charges,
             p3m_float *fields, p3m_float *potentials);
    
    void setRequireTotalEnergy(bool flag = true);
    fcs_float getTotalEnergy();
    
    void setRequireTimings(TimingType type = NONE);
    TimingType getRequireTimings();
    
    /** Test run the method with the current parameters.
     * Put the different run times in the TuneParameters. */
    const double* measureTimings(p3m_int num_particles,
                                 p3m_float *positions, p3m_float *charges);
    /** Fetch the detailed timings. */
    const double* getTimings();
    
    Communication comm;
    ErrorEstimate *errorEstimate;
    FarSolver *farSolver;

    /****************************************************
     * SYSTEM PARAMETERS
     ****************************************************/
    /* System size in x,y,z */
    p3m_float box_l[3];
    p3m_int sum_qpart;
    p3m_float sum_q2, square_sum_q;
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
    /** Tolerance in the field rms. */
    p3m_float tolerance_field;
    /** whether to compute the near field in the method */
    bool near_field_flag;
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

private:
    /****************************************************
     * FLAGS TO TURN ON/OFF COMPUTATION OF DIFFERENT COMPONENTS
     ****************************************************/
    /** Whether or not the total energy is to be computed. */
    bool require_total_energy;
    double total_energy;

    /** Whether and how timings are to be taken. */
    TimingType require_timings;
    double timings[NUM_TIMINGS];

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

    // submethods of tune()

    typedef std::list<TuneParameters> TuneParameterList;

    /* Get a list of parameters that achieve the required tolerance.
     * Does NOT time. */
    TuneParameterList tuneFar(Parameters p);
    /** Slave variant of tune_far. */
    void tuneFarSlave();

    void tuneAlpha(TuneParameterList &params_to_try);
    void tuneCAO(TuneParameterList &params_to_try);
    void tuneGrid(TuneParameterList &params_to_try);

    TuneParameters timeParams(p3m_int num_particles, p3m_float *positions,
            p3m_float *charges, TuneParameterList &params_to_try);

    void countCharges(p3m_int num_particles, p3m_float *charges);

    // tuning master functions
    void tuneBroadcastSendParams(TuneParameters p);
    void tuneBroadcastFail();
    void tuneBroadcastTiming(TuneParameters &p, p3m_int num_particles,
            p3m_float *positions, p3m_float *charges);
    void tuneBroadcastFinish(TuneParameters p);
    void tuneBroadcastNoTune();
    TuneParameterList tuneBroadcastTuneFar(TuneParameters p);
    
    void
    cartesianizeFields(p3m_float *fields, p3m_int num_particles);

    // tuning slave functions
    void tuneBroadcastReceiveParams();
    void tuneLoopSlave(p3m_int num_particles,
            p3m_float *positions, p3m_float *charges);

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

    void gatherTimings();
};


}

#endif
