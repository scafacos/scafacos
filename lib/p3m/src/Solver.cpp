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
#include "Solver.hpp"
#include "FarSolver.hpp"
#include "utils.hpp"
#include "FCSCommon.h"
#include "common/near/near.h"
#include <stdexcept>

namespace P3M {

Solver::Solver(MPI_Comm mpicomm) : comm(mpicomm), errorEstimate(NULL) {
    P3M_DEBUG(printf( "P3M::Solver() started...\n"));

    errorEstimate = ErrorEstimate::create(comm);
    farSolver = NULL;

    /* SYSTEM PARAMETERS */
    box_l[0] = 1.0;
    box_l[1] = 1.0;
    box_l[2] = 1.0;
    sum_qpart = 0;
    sum_q2 = 0.0;
    square_sum_q = 0.0;

    /* P3M PARAMETERS */
    skin = 0.0;
    tolerance_field = P3M_DEFAULT_TOLERANCE_FIELD;

    /* tunable */
    r_cut = 0.0;
    alpha = 0.0;
    grid[0] = 0;
    grid[1] = 0;
    grid[2] = 0;
    cao = 0;

    /* Everything needs to be retuned at the beginning */
    needs_retune = true;
    tune_r_cut = true;
    tune_alpha = true;
    tune_grid = true;
    tune_cao = true;

    /* Which components to compute? */
    require_total_energy = false;
    near_field_flag = false;

#ifdef P3M_PRINT_TIMINGS
    require_timings = FULL;
#else
    require_timings = NONE;
#endif

    resetTimers();

    P3M_DEBUG(printf( "P3M::Solver() finished.\n"));
}

Solver::~Solver() {
    if (errorEstimate != NULL) delete errorEstimate;
    if (farSolver != NULL) delete farSolver;
}

void Solver::prepare() {
    if (farSolver != NULL) delete farSolver;
    comm.prepare(box_l);
    farSolver = new FarSolver(comm, box_l, r_cut, alpha, grid, cao);
}

/* callback function for near field computations */
inline void
compute_near(const void *param, p3m_float dist, p3m_float *field, p3m_float *potential)
{
    const p3m_float alpha = *(static_cast<const p3m_float*>(param));
    const p3m_float adist = alpha * dist;

#ifdef P3M_USE_ERFC_APPROXIMATION

    /* approximate \f$ \exp(d^2) \mathrm{erfc}(d)\f$ by applying a formula from:
       Abramowitz/Stegun: Handbook of Mathematical Functions, Dover
       (9. ed.), chapter 7 */
    p3m_float t = 1.0 / (1.0 + 0.3275911 * adist);
    p3m_float erfc_part_ri = exp(-adist*adist) *
            (t * (0.254829592 +
                    t * (-0.284496736 +
                            t * (1.421413741 +
                                    t * (-1.453152027 +
                                            t * 1.061405429)))))
                                            / dist;

    *potential = erfc_part_ri;
    *field = -(erfc_part_ri + 2.0*alpha*0.56418958354775627928034964498*exp(-adist*adist))
              / dist;

#else

    p3m_float erfc_part_ri = (1.0 - erf(adist)) / dist; /* use erf instead of erfc to fix ICC performance problems */
    *potential = erfc_part_ri;
    *field = -(erfc_part_ri + 2.0*alpha*0.56418958354775627928034964498*exp(-adist*adist))
              / dist;
#endif

}

/* callback function for performing a whole loop of near field computations (using compute_near) */
FCS_NEAR_LOOP_FP(compute_near_loop, compute_near);

/* domain decomposition */
void Solver::decompose(fcs_gridsort_t *gridsort,
        p3m_int _num_particles,
        p3m_float *_positions, p3m_float *_charges,
        p3m_int *num_real_particles,
        p3m_float **positions, p3m_float **charges,
        fcs_gridsort_index_t **indices,
        p3m_int *num_ghost_particles,
        p3m_float **ghost_positions, p3m_float **ghost_charges,
        fcs_gridsort_index_t **ghost_indices
) {
    p3m_float box_base[3] = {0.0, 0.0, 0.0 };
    p3m_float box_a[3] = {box_l[0], 0.0, 0.0 };
    p3m_float box_b[3] = {0.0, box_l[1], 0.0 };
    p3m_float box_c[3] = {0.0, 0.0, box_l[2] };
    p3m_int num_particles;

    fcs_gridsort_create(gridsort);

    fcs_gridsort_set_system(gridsort, box_base, box_a, box_b, box_c, NULL);
    fcs_gridsort_set_particles(gridsort, _num_particles, _num_particles,
            _positions, _charges);

    P3M_DEBUG(printf( "  calling fcs_gridsort_sort_forward()...\n"));
    /* @todo: Set skin to r_cut only, when near field is wanted! */
    fcs_gridsort_sort_forward(gridsort,
            (near_field_flag ? r_cut : 0.0),
            comm.mpicomm);
    P3M_DEBUG(printf( "  returning from fcs_gridsort_sort_forward().\n"));
    fcs_gridsort_separate_ghosts(gridsort,
            num_real_particles,
            num_ghost_particles);

    fcs_gridsort_get_sorted_particles(gridsort,
            &num_particles, NULL, NULL, NULL, NULL);

    fcs_gridsort_get_real_particles(gridsort,
            num_real_particles,
            positions, charges,
            indices);

    fcs_gridsort_get_ghost_particles(gridsort, num_ghost_particles,
            ghost_positions, ghost_charges,
            ghost_indices);

    P3M_DEBUG_LOCAL(MPI_Barrier(comm.mpicomm));
    P3M_DEBUG_LOCAL(printf(                                               \
            "    %d: num_particles=%d"                     \
            " num_real_particles=%d"                       \
            " num_ghost_particles=%d\n",                   \
            comm.rank, num_particles,                   \
            *num_real_particles, *num_ghost_particles));
}

void Solver::run(
        p3m_int _num_particles, p3m_float *_positions, p3m_float *_charges,
        p3m_float *_fields, p3m_float *_potentials) {
    P3M_INFO(printf( "P3M::Solver::run() started...\n"));
    if (farSolver == NULL)
        throw std::logic_error("FarSolver is not initialized.");

    P3M_INFO(printf("    system parameters: box_l=" F3FLOAT "\n", \
            box_l[0], box_l[1], box_l[2]));
    P3M_DEBUG_LOCAL(MPI_Barrier(comm.mpicomm));
    P3M_DEBUG_LOCAL(printf("    %d: num_particles=%d\n",    \
            comm.rank, _num_particles));

    P3M_DEBUG(printf("  type of timing: %d\n", require_timings));
    /* reset all timers */
    if (require_timings != NONE) resetTimers();

    startTimer(TOTAL);
    startTimer(DECOMP);

    /* decompose system */
    p3m_int num_real_particles;
    p3m_int num_ghost_particles;
    p3m_float *positions, *ghost_positions;
    p3m_float *charges, *ghost_charges;
    fcs_gridsort_index_t *indices, *ghost_indices;
    fcs_gridsort_t gridsort;
    this->decompose(&gridsort,
            _num_particles, _positions, _charges,
            &num_real_particles,
            &positions, &charges, &indices,
            &num_ghost_particles,
            &ghost_positions, &ghost_charges, &ghost_indices);

    /* allocate local fields and potentials */
    p3m_float *fields = NULL;
    p3m_float *potentials = NULL;
    if (_fields != NULL)
        fields = new p3m_float[3*num_real_particles];
    if (_potentials != NULL || require_total_energy)
        potentials = new p3m_float[num_real_particles];

    stopTimer(DECOMP);

    if (require_timings != NOTFAR) {
#if defined(P3M_INTERLACE) && defined(P3M_AD)
        farSolver->runADI(num_real_particles, positions, charges, fields, potentials);
#else
        farSolver->runIK(num_real_particles, positions, charges, fields, potentials);
#endif
    }

    if (near_field_flag) {
        /* start near timer */
        startTimer(NEAR);

        /* compute near field */
        fcs_near_t near;

        fcs_near_create(&near);
        /*  fcs_near_set_field_potential(&near, compute_near);*/
        fcs_near_set_loop(&near, compute_near_loop);

        p3m_float box_base[3] = {0.0, 0.0, 0.0 };
        p3m_float box_a[3] = {box_l[0], 0.0, 0.0 };
        p3m_float box_b[3] = {0.0, box_l[1], 0.0 };
        p3m_float box_c[3] = {0.0, 0.0, box_l[2] };
        fcs_near_set_system(&near, box_base, box_a, box_b, box_c, NULL);

        fcs_near_set_particles(&near, num_real_particles, num_real_particles,
                positions, charges, indices,
                (_fields != NULL) ? fields : NULL,
                        (_potentials != NULL) ? potentials : NULL);

        fcs_near_set_ghosts(&near, num_ghost_particles,
                ghost_positions, ghost_charges, ghost_indices);

        P3M_DEBUG(printf( "  calling fcs_near_compute()...\n"));
        fcs_near_compute(&near, r_cut, &alpha, comm.mpicomm);
        P3M_DEBUG(printf( "  returning from fcs_near_compute().\n"));

        fcs_near_destroy(&near);

        stopTimer(NEAR);
    }

    startTimer(COMP);
    /* sort particles back */
    P3M_DEBUG(printf( "  calling fcs_gridsort_sort_backward()...\n"));
    fcs_gridsort_sort_backward(&gridsort,
            fields, potentials,
            _fields, _potentials, 1,
            comm.mpicomm);
    P3M_DEBUG(printf( "  returning from fcs_gridsort_sort_backward().\n"));

    fcs_gridsort_free(&gridsort);
    fcs_gridsort_destroy(&gridsort);

    stopTimer(COMP);
    stopTimer(TOTAL);

    // gather timings
    if (comm.onMaster())
        MPI_Reduce(MPI_IN_PLACE, timings,
                NUM_TIMINGS_NOTFAR, MPI_DOUBLE, MPI_MAX,
                0, comm.mpicomm);
    else
        MPI_Reduce(timings, 0,
                NUM_TIMINGS_NOTFAR, MPI_DOUBLE, MPI_MAX,
                0, comm.mpicomm);

    // copy the far field timings to the end of the timings
    const double *farTimings = farSolver->getTimings();
    memcpy(timings+NUM_TIMINGS_NOTFAR, farTimings, sizeof(double)*FarSolver::NUM_TIMINGS);

#ifdef P3M_PRINT_TIMINGS
#define PRINT(s, ID) printf("%10s=%le (%lf)\n", s, \
    timings[ID], (timings[ID]/timings[TOTAL]))

    if (require_timings != NONE && comm.onMaster()) {
        printf("  P3M TIMINGS:\n");
        printf("%10s=%le\n", "total", timings[TOTAL]);
        PRINT("far", FAR);
        PRINT("ca", CA);
        PRINT("potentials", POTENTIALS);
        PRINT("fields", FIELDS);
        //if(require_timings == ESTIMATE_ALL //not yet implemented
        //|| require_timings == ESTIMATE_ASSIGNMENT)
        //    printf(" (empirical estimate)");
        //if(require_timings == ESTIMATE_ALL //not yet implemented
        //|| require_timings == ESTIMATE_ASSIGNMENT)
        //    printf(" (empirical estimate)");
        PRINT("gather", GATHER);
        PRINT("spread", SPREAD);
        PRINT("forward", FORWARD);
        PRINT("back", BACK);
        printf("\n");
        PRINT("near", NEAR);
        PRINT("comp", COMP);
        PRINT("decomp", DECOMP);
//        if(require_timings == ESTIMATE_ALL
//                || require_timings == ESTIMATE_FFT)
//            printf(" (empirical estimate)");
    }
#endif

    sdelete(fields);
    sdelete(potentials);

    P3M_INFO(printf( "P3M::Solver::run() finished.\n"));
}



/***************************************************/
/******************** TUNING ***********************/
/***************************************************/
/** Good mesh sizes for fftw
 */
const p3m_int good_gridsize[] = {
        0, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32,
        36, 40, 42, 44, 48, 50, 52, 54, 56, 60, 64,
        66, 70, 72, 78, 80, 84, 88, 90, 96, 98, 100, 104, 108, 110, 112, 120, 126, 128,
        130, 132, 140, 144, 150, 154, 156, 160, 162, 168, 176, 180, 182, 192,
        196, 198, 200, 208, 210, 216, 220, 224, 234, 240, 242, 250, 252, 256,
        260, 264, 270, 280, 288, 294,
        300, 360, 400, 480, 512,
        576, 648, 729, 768, 864, 972, 1024,
        1152, 1296, 1458, 1536,
        1728, 1944, 2048,
        2187, 2304, 2592,
        2916, 3072,
        3456,
        3888};

/** Steps to do when trying to determine smallest possible grid size. */
const p3m_int step_good_gridsize[] =
{ 0, 15, 26, 44, 58, 72, 78, 83, 90, 94, 98, 101, 103, 104 } ;
const p3m_int num_steps_good_gridsize = 14;


void
Solver::tune(p3m_int num_particles, p3m_float *positions, p3m_float *charges) {
    /* Prepare the communicator before tuning */
    comm.prepare(box_l);

    /* Count the charges */
    p3m_float sum_q2_before = sum_q2;
    this->countCharges(num_particles, charges);

    if (!comm.onMaster()) {
        this->tuneLoopSlave(num_particles, positions, charges);
    } else {
        /* Retune if the number of charges has changed */
        if (!(float_is_equal(sum_q2_before, sum_q2))) {
            P3M_INFO(printf( "  Number of charges changed, retuning is required.\n"));
            needs_retune = true;
        }

        /* Do not retune if there are no charges */
        if (float_is_zero(sum_q2)) {
            P3M_INFO(printf( "  No charges in the system.\n"));
            needs_retune = false;
        }

        /* Exit if retuning is unnecessary */
        if (!needs_retune) {
            P3M_INFO(printf( "  Retuning is not required.\n"));
            P3M_INFO(printf("    Finished tuning.\n"));
            tuneBroadcastNoTune();
            return;
        }

        P3M_INFO(printf( "  Retuning is required.\n"));

        if (!tune_r_cut) {
            P3M_INFO(printf( "    r_cut=" FFLOAT " (fixed)\n", r_cut));

            TuneParameters p;
            p.cao = cao;
            p.alpha = alpha;
            p.grid[0] = grid[0];
            p.grid[1] = grid[1];
            p.grid[2] = grid[2];
            p.r_cut = r_cut;
            TuneParameterList params_to_try =
                    this->tuneBroadcastTuneFar(p);
            TuneParameters best =
                    this->timeParams(num_particles, positions, charges,
                            params_to_try);
            this->tuneBroadcastFinish(best);
        } else {
            TuneParameters p;
            p.cao = cao;
            p.alpha = alpha;
            p.grid[0] = grid[0];
            p.grid[1] = grid[1];
            p.grid[2] = grid[2];

            /* compute the average distance between two charges  */
            p3m_float avg_dist =
                    pow((box_l[0]*box_l[1]*box_l[2])
                            / sum_qpart, 0.33333);

            /* FIRST ESTIMATE */
            /* set the initial r_cut to 3 times the average distance
             * between charges */
            p.r_cut = 3.0 * avg_dist;

            /* tune r_cut to half the box length */
            if (0.5*box_l[1]-skin < p.r_cut)
                p.r_cut = 0.5*box_l[0] - skin;
            if (0.5*box_l[1]-skin < p.r_cut)
                p.r_cut = 0.5*box_l[1] - skin;
            if (0.5*box_l[2]-skin < p.r_cut)
                p.r_cut = 0.5*box_l[2] - skin;

            P3M_INFO(printf( "    r_cut=" FFLOAT
                    " (first estimate)\n", p.r_cut));

            TuneParameterList params_to_try1 =
                    this->tuneBroadcastTuneFar(p);
            TuneParameters best1 =
                    this->timeParams(num_particles, positions, charges,
                            params_to_try1);

            /* SECOND ESTIMATE */
            /* use the fact that we know that timing_near scales like
             * r_cut**3 and timing_far like r_cut**(-3) to get the second
             * estimate */

            p3m_float rel_timing_diff =
                    fabs(best1.timing_near - best1.timing_far) /
                    (best1.timing_near + best1.timing_far);
            P3M_INFO(printf( "    rel_timing_diff=" FFLOAT "\n", rel_timing_diff));

            p3m_float rcut3 = pow(best1.r_cut, 3);
            p3m_float c_near = best1.timing_near/rcut3;
            p3m_float c_far = best1.timing_far*rcut3;
            p.r_cut = pow(c_far/c_near, 1./6.);

            P3M_INFO(printf( "    r_cut=" FFLOAT " (second estimate)\n", \
                    p.r_cut));

            TuneParameterList params_to_try2 =
                    this->tuneBroadcastTuneFar(p);
            TuneParameters best2 =
                    this->timeParams(num_particles, positions, charges,
                            params_to_try2);

            rel_timing_diff =
                    fabs(best2.timing_near - best2.timing_far) /
                    (best2.timing_near + best2.timing_far);
            P3M_INFO(printf("    rel_timing_diff=" FFLOAT "\n", rel_timing_diff));
            P3M_INFO(printf("    Finished tuning.\n"));

            this->tuneBroadcastFinish(best2);
        }
    }
}

/** Slave variant of tuneFar. */
void
Solver::tuneFarSlave() {
    if (comm.onMaster())
        throw std::logic_error("tuneFarSlave cannot be called on master.");
    /* wait for calls to error estimate until finished */
    errorEstimate->loopSlave();
}

Solver::TuneParameterList
Solver::tuneFar(Parameters p) {
    P3M_INFO(printf( "P3M::Solver::tuneFar() started...\n"));

    if (!comm.onMaster())
        throw std::logic_error("Internal error: Do not call tuneFar() on slave!");

    TuneParameterList params;
    try {
        /* check whether the input parameters are sane */
        if (p.r_cut < 0.0)
            throw std::domain_error("r_cut is negative!");

        if (float_is_zero(p.r_cut))
            throw std::domain_error("r_cut is too small!");

        /* check whether cutoff is larger than half a box length */
        if ((p.r_cut > 0.5 * box_l[0]) ||
                (p.r_cut > 0.5 * box_l[1]) ||
                (p.r_cut > 0.5 * box_l[2]))
            throw std::domain_error(
                    "r_cut is larger than half a system box length.");

        /* check whether cutoff is larger than domain size */

        params.push_back(TuneParameters(p));
        this->tuneAlpha(params);
        this->tuneCAO(params);
        this->tuneGrid(params);
        errorEstimate->endLoop();
    } catch (...) {
        P3M_INFO(printf("  Tuning failed.\n"));
        errorEstimate->endLoop();
        this->tuneBroadcastFail();
        P3M_INFO(printf("P3M::Solver::tuneFar() finished.\n"));
        throw;
    }

    P3M_INFO(printf( "  Tuning was successful.\n"));

    P3M_INFO(printf( "P3M::Solver::tuneFar() finished.\n"));

    return params;
}


/* Tune alpha */
void
Solver::tuneAlpha(TuneParameterList &params_to_try) {
    if (tune_alpha) {
        // compute the alpha for all params to try
        for (TuneParameterList::iterator pit = params_to_try.begin();
                pit != params_to_try.end(); ++pit) {
            errorEstimate->computeAlpha(tolerance_field,
                    *pit, sum_qpart, sum_q2, box_l);
            P3M_INFO(printf("    => alpha=" FFLOAT "\n", pit->alpha));
        }
    } else
        P3M_INFO(printf("    alpha=" FFLOAT " (fixed)\n", \
                params_to_try.front().alpha));

}

/* Tune cao */
void
Solver::tuneCAO(TuneParameterList &params_to_try) {
    if (tune_cao) {
#ifdef P3M_AD
        const p3m_int min_cao = 2;
#else
        const p3m_int min_cao = 1;
#endif
        P3M_INFO(printf("    Testing cao={ "));
        for (TuneParameterList::iterator pit = params_to_try.begin();
                pit != params_to_try.end(); ++pit) {
            for (p3m_int cao = CAF::MAX_CAO; cao > min_cao; cao--) {
                TuneParameters p = *pit;
                p.cao = cao;
                params_to_try.insert(pit, p);
                P3M_INFO(printf(FINT " ", cao));
            }
            pit->cao = min_cao;
        }
        P3M_INFO(printf(FINT " }\n", min_cao));
    } else {
        P3M_INFO(printf( "    cao=" FINT " (fixed)\n", \
                params_to_try.front().cao));
    }
}

/* params should have decreasing cao */
void
Solver::tuneGrid(TuneParameterList &params_to_try) {
    if (tune_grid) {
        p3m_int step_ix = 0;
        // store the minimal grid size seen so far
        p3m_int min_grid1d = good_gridsize[step_good_gridsize[num_steps_good_gridsize-1]];

        for (TuneParameterList::iterator pit = params_to_try.begin();
                pit != params_to_try.end(); ++pit) {
            // Find smallest possible grid for this set of parameters
            // Test whether accuracy can be achieved with given parameters
            p3m_int upper_ix;
            // test this step
            P3M_DEBUG(printf("    Trying to find grid for r_cut=" FFLOAT ", " \
                    "alpha=" FFLOAT ", "                          \
                    "cao=" FINT "\n",                             \
                    pit->r_cut, pit->alpha, pit->cao));

            TuneParameters p = *pit;

            do {
                step_ix++;
                if (step_ix >= num_steps_good_gridsize) break;
                upper_ix = step_good_gridsize[step_ix];
                p3m_int grid1d = good_gridsize[upper_ix];
                p.grid[0] = p.grid[1] = p.grid[2] = grid1d;
                P3M_DEBUG(printf("      rough grid=" FINT "\n", grid1d));
                errorEstimate->computeMaster(p, sum_qpart, sum_q2, box_l);
            } while (p.error > tolerance_field);
            // store the (working) rough grid param set
            *pit = p;

            // reached largest possible grid, remove the rest of the parameter sets
            if (step_ix >= num_steps_good_gridsize) {
                P3M_INFO(printf("    Too large grid size, skipping rest of parameter sets.\n"));
                TuneParameterList::iterator next = pit;
                ++next;
                params_to_try.erase(next, params_to_try.end());
                break;
            }

            // error would be small enough at upper_ix, but not at lower_ix,
            // so bisect to find optimal gridsize
            p3m_int lower_ix = step_good_gridsize[step_ix-1];
            while (lower_ix+1 < upper_ix) {
                p3m_int test_ix = (lower_ix+upper_ix)/2;
                p3m_int grid1d = good_gridsize[test_ix];
                p.grid[0] = p.grid[1] = p.grid[2] = grid1d;
                P3M_DEBUG(printf("      fine grid=" FINT "\n", grid1d));
                errorEstimate->computeMaster(p, sum_qpart, sum_q2, box_l);
                if (p.error < tolerance_field) {
                    // parameters achieve error
                    upper_ix = test_ix;
                    // save the current param set
                    *pit = p;
                } else {
                    // parameters do not achieve error
                    lower_ix = test_ix;
                }
            }

            P3M_INFO(printf( "      => "                          \
                    "r_cut=" FFLOAT ", "                          \
                    "alpha=" FFLOAT ", "                          \
                    "cao=" FINT "\n"                              \
                    "grid=" F3INT ", "                            \
                    "error=" FFLOATE "\n",                        \
                    pit->r_cut, pit->alpha, pit->cao,             \
                    pit->grid[0], pit->grid[1], pit->grid[2],     \
                    pit->error));

            // decrease step_ix so that the same step_ix is tested for the
            // next param set
            step_ix--;

            // compare grid size to previous minimal grid size
            // if it is larger than any previous size + P3M_MAX_GRID_DIFF,
            // skip it
            if (min_grid1d > pit->grid[0]) {
                min_grid1d = pit->grid[0];
            } else if (min_grid1d + P3M_MAX_GRID_DIFF < pit->grid[0]) {
                P3M_INFO(printf("      grid too large => removing data set\n"));
                // remove the rest of the params
                TuneParameterList::iterator next = pit;
                ++next;
                params_to_try.erase(next, params_to_try.end());
                break;
            }

            if (pit != params_to_try.begin()) {
                TuneParameterList::iterator prev = pit;
                prev--;
                if (prev->grid[0] >= pit->grid[0]) {
                    P3M_INFO(printf("      better than previous => removing previous data set.\n"));
                    params_to_try.erase(prev);
                }
            }
        }
    } else {
        // !tune_grid
        TuneParameterList::iterator pit = params_to_try.begin();
        while (pit != params_to_try.end()) {
            // test whether accuracy can be achieved with given parameters
            P3M_INFO(printf("    grid=" F3INT " (fixed)\n",                \
                    pit->grid[0], pit->grid[1], pit->grid[2]));
            errorEstimate->computeMaster(*pit, sum_qpart, sum_q2, box_l);

            if (pit->error < tolerance_field) {
                // error is small enough for this parameter set, so keep it
                P3M_INFO(printf("    error (%le) < tolerance (%le), keeping params\n", \
                        pit->error, tolerance_field));
                ++pit;
            } else {
                // otherwise remove this parameter set
                P3M_INFO(printf("    error (%le) > tolerance (%le), removing params\n", \
                        pit->error, tolerance_field));
                TuneParameterList::iterator badpit = pit;
                ++pit;
                params_to_try.erase(badpit);
            }
        }
    }
}

TuneParameters
Solver::timeParams(
        p3m_int num_particles, p3m_float *positions, p3m_float *charges,
        TuneParameterList &params_to_try) {

    if (params_to_try.empty())
        throw std::logic_error("No parameter set left to time.");

    /* Now time the different parameter sets */
    TuneParameters best_params;
    double best_timing = 1.e100;

#ifdef P3M_ENABLE_INFO
    printf("Timing %ld param sets...\n", params_to_try.size());
    for (TuneParameterList::iterator pit = params_to_try.begin();
            pit != params_to_try.end();
            ++pit) {
        printf( "  r_cut=" FFLOAT ", "
                "alpha=" FFLOAT ", "
                "grid=" F3INT ", "
                "cao=" FINT "\n",
                pit->r_cut, pit->alpha,
                pit->grid[0], pit->grid[1], pit->grid[2],
                pit->cao);
    }
#endif

    for (TuneParameterList::iterator pit = params_to_try.begin();
            pit != params_to_try.end();
            ++pit) {
        this->tuneBroadcastTiming(*pit, num_particles, positions, charges);
        P3M_INFO(printf( "  Timing r_cut=" FFLOAT ", "         \
                "alpha=" FFLOAT ", "                           \
                "grid=" F3INT ", "                             \
                "cao=" FINT "\n    "                           \
                "=> timing=" FFLOAT " "                        \
                "(" FFLOAT " near, " FFLOAT " far)\n",         \
                r_cut, alpha,                                  \
                grid[0], grid[1], grid[2],                     \
                cao, pit->timing,                              \
                pit->timing_near, pit->timing_far));

        if (pit->timing < best_timing) {
            best_timing = pit->timing;
            best_params = *pit;
        }
    }

    return best_params;
}

void P3M::Solver::setRequireTimings(TimingType type) {
    require_timings = type;
}

/** Test run the method with the current parameters.
 * Return the total run time. */
const double* Solver::measureTimings(p3m_int num_particles,
        p3m_float *positions, p3m_float *charges) {
    p3m_float *fields = new p3m_float[3*num_particles];
    p3m_float *potentials = new p3m_float[num_particles];

    /* store require_timings */
    TimingType require_timings_before = require_timings;
    if (require_timings_before == NONE)
        require_timings = FULL;
    this->run(num_particles, positions, charges, fields, potentials);

    /* restore require_timings */
    require_timings = require_timings_before;

    delete fields;
    delete potentials;

    return timings;
}

const double* Solver::getTimings() {
    if (require_timings == NONE)
        throw std::logic_error("Wanting to get timings, but timings were not measured.");
    return timings;
}

/** Calculate number of charged particles, the sum of the squared
      charges and the squared sum of the charges. Called in parallel at
      the beginning of tuning. */
void Solver::countCharges(p3m_int num_particles, p3m_float *charges) {
    p3m_float node_sums[3], tot_sums[3];

    for (int i=0; i<3; i++) {
        node_sums[i] = 0.0;
        tot_sums[i] = 0.0;
    }

    for (int i = 0; i < num_particles; i++) {
        if (!float_is_zero(charges[i])) {
            node_sums[0] += 1.0;
            node_sums[1] += SQR(charges[i]);
            node_sums[2] += charges[i];
        }
    }

    MPI_Allreduce(node_sums, tot_sums, 3, P3M_MPI_FLOAT, MPI_SUM, comm.mpicomm);
    sum_qpart    = (p3m_int)(tot_sums[0]+0.1);
    sum_q2       = tot_sums[1];
    square_sum_q = SQR(tot_sums[2]);

    P3M_DEBUG(printf("  countCharges(): "          \
            "num_charged_particles=" FINT ", "                   \
            "sum_squared_charges=" FFLOAT ", "                   \
            "sum_charges=" FFLOAT "\n",                          \
            sum_qpart, sum_q2, sqrt(square_sum_q)));
}

void Solver::setRequireTotalEnergy(bool flag) {
    require_total_energy = flag;
}

p3m_float Solver::getTotalEnergy() {
    if (!require_total_energy)
        throw std::logic_error("");
    return total_energy;
}


/* Events during tuning */
static const int CMD_FINISHED = 0;
static const int CMD_FAILED = 1;
static const int CMD_NO_TUNE = 2;
static const int CMD_TUNE_FAR = 3;
static const int CMD_TIMING = 4;
static const char* CMD_NAMES[5] =
{ "FINISHED", "FAILED", "NO_TUNE", "TUNE_FAR", "TIMING" };


static void tuneBroadcastCommand(Communication &comm, p3m_int command) {
    /* Send the command */
    P3M_DEBUG_LOCAL(printf("       %2d: Broadcasting command %s.\n", \
            comm.rank, CMD_NAMES[command]));
    MPI_Bcast(&command, 1, P3M_MPI_INT, 0, comm.mpicomm);
}

void Solver::tuneBroadcastSendParams(TuneParameters p) {
    // broadcast the parameters
    p3m_int int_buffer[4];
    p3m_float float_buffer[2];

    // pack int data
    int_buffer[0] = p.grid[0];
    int_buffer[1] = p.grid[1];
    int_buffer[2] = p.grid[2];
    int_buffer[3] = p.cao;
    MPI_Bcast(int_buffer, 4, P3M_MPI_INT, 0, comm.mpicomm);

    // pack float data
    float_buffer[0] = p.alpha;
    float_buffer[1] = p.r_cut;
    MPI_Bcast(float_buffer, 2, P3M_MPI_FLOAT, 0, comm.mpicomm);

    r_cut = p.r_cut;
    alpha = p.alpha;
    grid[0] = p.grid[0];
    grid[1] = p.grid[1];
    grid[2] = p.grid[2];
    cao = p.cao;
}

void Solver::tuneBroadcastReceiveParams() {
    // receive parameters
    p3m_int int_buffer[4];
    p3m_float float_buffer[2];

    MPI_Bcast(int_buffer, 4, P3M_MPI_INT, 0, comm.mpicomm);
    grid[0] = int_buffer[0];
    grid[1] = int_buffer[1];
    grid[2] = int_buffer[2];
    cao = int_buffer[3];

    // unpack float data
    MPI_Bcast(float_buffer, 2, P3M_MPI_FLOAT, 0,  comm.mpicomm);
    alpha = float_buffer[0];
    r_cut = float_buffer[1];
}


void Solver::tuneBroadcastFail() {
    tuneBroadcastCommand(comm, CMD_FAILED);
    needs_retune = true;
}

void Solver::tuneBroadcastFinish(TuneParameters p) {
    tuneBroadcastCommand(comm, CMD_FINISHED);
    this->tuneBroadcastSendParams(p);
    needs_retune = false;
    this->prepare();
}

void Solver::tuneBroadcastNoTune() {
    tuneBroadcastCommand(comm, CMD_NO_TUNE);
    needs_retune = false;
}

Solver::TuneParameterList
Solver::tuneBroadcastTuneFar(TuneParameters p) {
    tuneBroadcastCommand(comm, CMD_TUNE_FAR);
    return this->tuneFar(p);
}

void Solver::tuneBroadcastTiming(TuneParameters &p, p3m_int num_particles,
        p3m_float *positions, p3m_float *charges) {
    tuneBroadcastCommand(comm, CMD_TIMING);
    this->tuneBroadcastSendParams(p);
    this->prepare();

    const double *timings = this->measureTimings(num_particles, positions, charges);
    p.timing = timings[TOTAL];
    p.timing_far = timings[FAR];
    p.timing_near = timings[NEAR];

}

void Solver::tuneLoopSlave(p3m_int num_particles, p3m_float* positions,
        p3m_float* charges) {
    P3M_DEBUG_LOCAL(printf( "      %2d: P3M::Solver::tuneLoopSlave() " \
            "started...\n", comm.rank));
    if (comm.onMaster())
        throw std::logic_error("Internal error: tuneLoopSlave "
                "should not be called on master!");

    for (;;) {
        /* Receive the command */
        p3m_int command;
        P3M_DEBUG_LOCAL(printf("      %2d: Solver::tuneLoopSlave(): " \
                "Waiting to receive command.\n",                      \
                comm.rank));
        MPI_Bcast(&command, 1, P3M_MPI_INT, 0, comm.mpicomm);
        P3M_DEBUG_LOCAL(printf("      %2d: Solver::tuneLoopSlave(): " \
                "Received command %s.\n", \
                comm.rank, CMD_NAMES[command]));

        switch (command) {
        case CMD_TUNE_FAR:
            this->tuneFarSlave();
            break;
        case CMD_TIMING:
            this->tuneBroadcastReceiveParams();
            this->prepare();
            this->measureTimings(num_particles, positions, charges);
            break;
        case CMD_FINISHED:
            this->tuneBroadcastReceiveParams();
            P3M_DEBUG_LOCAL(printf( "P3M::Solver::tuneLoopSlave() " \
                    "finished.\n"));
            this->prepare();
            needs_retune = false;
            break;
        case CMD_NO_TUNE:
            P3M_DEBUG_LOCAL(printf( "P3M::Solver::tuneLoopSlave() " \
                    "finished.\n"));
            needs_retune = false;
            break;
        case CMD_FAILED: {
            needs_retune = true;
            char msg[255];
            sprintf(msg,
                    "Cannot achieve required accuracy (p3m_tolerance_field="
                    FFLOATE ") for given parameters.",
                    tolerance_field);
            throw std::logic_error(msg);
        }
        }
        if (command == CMD_FAILED
                || command == CMD_FINISHED
                || command == CMD_NO_TUNE)
            break;
    }
}

}
