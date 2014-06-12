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
#include "FarSolver.hpp"
#include <stdexcept>

P3M::FarSolver::FarSolver(Communication &comm, p3m_float box_l[3],
        p3m_float r_cut, p3m_float alpha, p3m_int grid[3], p3m_int cao, p3m_float box_vectors[3][3], p3m_float volume, bool isTriclinic)
: comm(comm), fft(comm), errorEstimate(NULL) {
    P3M_DEBUG(printf( "P3M::FarSolver() started...\n"));

    errorEstimate = ErrorEstimate::create(comm);

    this->g_force = NULL;
    this->g_energy = NULL;
    this->d_op[0] = NULL;
    this->d_op[1] = NULL;
    this->d_op[2] = NULL;

    for(int i = 0; i < 3 ; ++i){
    this->box_l[i] = box_l[i];
    for(int j = 0; j < 3 ; ++j)
    this->box_vectors[i][j]=box_vectors[i][j];
    }
    this->isTriclinic=isTriclinic;
    this->volume = volume;
    
    this->sum_qpart = 0;
    this->sum_q2 = 0.0;
    this->square_sum_q = 0.0;
    P3M_INFO(printf("    box_l=" F3FLOAT "\n", box_l[0], box_l[1], box_l[2]));

    this->skin = 0.0;
    this->n_interpol = P3M_DEFAULT_N_INTERPOL;
    for (int i = 0; i < 3; i++) {
        this->additional_grid[i] = 0;
        this->grid_off[i] = P3M_DEFAULT_GRIDOFF;
    }

    this->r_cut = r_cut;
    this->alpha = alpha;
    this->grid[0] = grid[0];
    this->grid[1] = grid[1];
    this->grid[2] = grid[2];
    this->cao = cao;

    P3M_INFO(printf(                                              \
            "    p3m params: "                                    \
            "r_cut=" FFLOAT ", grid=" F3INT ", cao=" FINT ", "    \
            "alpha=" FFLOAT "\n",                                 \
            r_cut, grid[0], grid[1], grid[2], cao, alpha));

    /* Which components to compute? */
    require_total_energy = false;
    total_energy = 0.0;

#ifdef P3M_PRINT_TIMINGS
    require_timings = FULL;
#else
    require_timings = NONE;
#endif
    for (int i=0; i < NUM_TIMINGS; i++)
        timings[i] = 0.0;

    /* Initializes the (inverse) grid constant ai and the cutoff for charge
     * assignment cao_cut. */
    for (p3m_int i=0; i<3; i++) {
        ai[i]      = grid[i]/box_l[i];
        a[i]       = 1.0/ai[i];
        cao_cut[i] = 0.5*a[i]*cao;
    }
    this->prepareLocalCAGrid();
    P3M_DEBUG(this->printLocalGrid());
    this->prepareSendGrid();
    P3M_DEBUG(this->printSendGrid());
    send_grid = new p3m_float[sm.max];
    recv_grid = new p3m_float[sm.max];

    P3M_DEBUG(printf("    Interpolating charge assignment function...\n"));
    caf = P3M::CAF::create(cao, n_interpol);
    cafx = caf->createCache();
    cafy = caf->createCache();
    cafz = caf->createCache();
#ifdef P3M_AD
    caf_d = P3M::CAF::create(cao, n_interpol, true);
    cafx_d = caf_d->createCache();
    cafy_d = caf_d->createCache();
    cafz_d = caf_d->createCache();
#else
    caf_d = NULL;
    cafx_d = cafy_d = cafz_d = NULL;
#endif

    /* position offset for calc. of first gridpoint */
    pos_shift = (p3m_float)((cao-1)/2) - (cao%2)/2.0;
    P3M_DEBUG(printf("    pos_shift=" FFLOAT "\n",pos_shift));

    /* FFT */
    P3M_INFO(printf("    Preparing FFTs...\n"));
    fft.prepare(local_grid.dim, local_grid.margin, grid, grid_off, &ks_pnum);
    rs_grid = fft.malloc_data();
    ks_grid = fft.malloc_data();
    buffer = fft.malloc_data();

    /* k-space part */
    /* Calculates the Fourier transformed differential operator.
     *  Remark: This is done on the level of n-vectors and not k-vectors,
     *           i.e. the prefactor i*2*PI/L is missing! */
    for (int i = 0; i < 3; i++) {
        d_op[i] = new p3m_int[grid[i]];
        d_op[i][0] = 0;
        d_op[i][grid[i]/2] = 0;

        for (int j = 1; j < grid[i]/2; j++) {
            d_op[i][j] = j;
            d_op[i][grid[i] - j] = -j;
        }
    }

    P3M_INFO(printf("    Calculating influence function...\n"));
#if !defined(P3M_INTERLACE) && defined(P3M_IK)
    computeInfluenceFunctionIK();
#elif defined(P3M_INTERLACE) && defined(P3M_IK)
    computeInfluenceFunctionIKI();
#else
    computeInfluenceFunctionADI();
#endif

    P3M_DEBUG(printf( "P3M::FarSolver() finished.\n"));
}

P3M::FarSolver::~FarSolver() {
    delete errorEstimate;
    fft.free_data(rs_grid);
    fft.free_data(ks_grid);
    fft.free_data(buffer);
    delete send_grid;
    delete recv_grid;
    delete caf;
    delete cafx;
    delete cafy;
    delete cafz;
    delete caf_d;
    delete cafx_d;
    delete cafy_d;
    delete cafz_d;
    delete g_force;
    delete g_energy;
    delete d_op[0];
    delete d_op[1];
    delete d_op[2];
}

void P3M::FarSolver::runADI(p3m_int num_charges, p3m_float* positions,
        p3m_float* charges, p3m_float* fields, p3m_float* potentials) {
    P3M_INFO(printf( "P3M::Solver::runADI() started...\n"));

    countCharges(num_charges, charges);

    startTimer(CA);

    /* charge assignment */
    this->assignCharges(buffer, num_charges, positions, charges, 0);
    switchTimer(CA, GATHER);
    /* gather the ca grid */
    this->gatherGrid(buffer);

    // Complexify
    for (p3m_int i = local_grid.size-1; i >= 0; i--)
        rs_grid[2*i] = buffer[i];

    switchTimer(GATHER, CA);

    // Second (shifted) run
    /* charge assignment */
    this->assignCharges(buffer, num_charges, positions, charges, 1);

    switchTimer(CA, GATHER);

    /* gather the ca grid */
    this->gatherGrid(buffer);
    /* now rs_grid should contain the local ca grid */

    // Complexify
    for (p3m_int i = local_grid.size-1; i >= 0; i--)
        rs_grid[2*i+1] = buffer[i];
    stopTimer(GATHER);

    /* forward transform */
    P3M_DEBUG(printf("  calling fft_perform_forw()...\n"));
    startTimer(FORWARD);
    fft.forward(rs_grid, buffer);
    stopTimer(FORWARD);
    P3M_DEBUG(printf("  returned from fft_perform_forw().\n"));

    /********************************************/
    /* POTENTIAL COMPUTATION */
    /********************************************/
    if (require_total_energy || potentials != NULL) {
        /* apply energy optimized influence function */
        startTimer(INFLUENCE);
        this->applyInfluenceFunction(rs_grid, ks_grid, g_energy);
        /* result is in ks_grid */
        stopTimer(INFLUENCE);

        /* compute total energy, but not potentials */
        if (require_total_energy && potentials == NULL) {
            startTimer(POTENTIALS);
            total_energy = this->computeTotalEnergy();
            stopTimer(POTENTIALS);
        }

        if (potentials != NULL) {
            /* backtransform the grid */

            if( require_timings != ESTIMATE_ALL
                    && require_timings != ESTIMATE_FFT ){
                P3M_DEBUG(printf( "  calling fft_perform_back (potentials)...\n"));
                startTimer(BACK);
                fft.backward(ks_grid, buffer);
                stopTimer(BACK);
                P3M_DEBUG(printf( "  returned from fft_perform_back.\n"));
            }
            /** First (unshifted) run */
            startTimer(SPREAD);
            for (p3m_int i = 0; i < local_grid.size; i++)
                buffer[i] = ks_grid[2*i];
            this->spreadGrid(buffer);
            switchTimer(SPREAD, POTENTIALS);

            this->assignPotentials(buffer, num_charges, positions,
                    charges, 0, potentials);

            switchTimer(POTENTIALS, SPREAD);

            /** Second (shifted) run */
            for (p3m_int i = 0; i < local_grid.size; i++)
                buffer[i] = ks_grid[2*i+1];
            this->spreadGrid(buffer);

            switchTimer(SPREAD, POTENTIALS);

            this->assignPotentials(buffer, num_charges, positions,
                    charges, 1, potentials);

            stopTimer(POTENTIALS);
        }
    }

    /********************************************/
    /* FIELD COMPUTATION */
    /********************************************/
    if (fields != NULL) {
        /* apply force optimized influence function */
        startTimer(INFLUENCE);
        this->applyInfluenceFunction(rs_grid, ks_grid, g_force);
        stopTimer(INFLUENCE);

        /* backtransform the grid */
        if(require_timings != ESTIMATE_ALL
                && require_timings != ESTIMATE_FFT){
            P3M_DEBUG(printf("  calling fft_perform_back...\n"));
            startTimer(BACK);
            fft.backward(ks_grid, buffer);
            stopTimer(BACK);
            P3M_DEBUG(printf("  returned from fft_perform_back.\n"));
        }

        startTimer(SPREAD);
        /* First (unshifted) run */
        P3M_INFO(printf("  computing unshifted grid\n"));
        for (p3m_int i=0; i<local_grid.size; i++) {
            buffer[i] = ks_grid[2*i];
        }

        this->spreadGrid(buffer);

        switchTimer(SPREAD, FIELDS);

        this->assignFieldsAD(buffer, num_charges, positions, 0, fields);

        switchTimer(FIELDS, SPREAD);

        /* Second (shifted) run */
        P3M_INFO(printf("  computing shifted grid\n"));
        for (p3m_int i=0; i<local_grid.size; i++) {
            buffer[i] = ks_grid[2*i+1];
        }

        this->spreadGrid(buffer);

        switchTimer(SPREAD, FIELDS);

        this->assignFieldsAD(buffer, num_charges, positions, 1, fields);
        
        if(isTriclinic) this->cartesianizeFields(fields, num_charges);
        
        stopTimer(FIELDS);
    }

    /* estimate FFT back timing */
    if (require_timings == ESTIMATE_ALL || require_timings == ESTIMATE_FFT) {
        if (potentials != NULL)
            timings[BACK] = 2 * timings[FORWARD];
        else
            timings[BACK] = timings[FORWARD];
    }

    this->gatherTimings();

    P3M_INFO(printf( "P3M::FarSolver::runADI() finished.\n"));
}

void P3M::FarSolver::runIK(p3m_int num_charges, p3m_float* positions,
        p3m_float* charges, p3m_float* fields, p3m_float* potentials) {
    P3M_INFO(printf( "P3M::FarSolver::runIK() started...\n"));

    countCharges(num_charges, charges);

    startTimer(CA);

    /* charge assignment */
    this->assignCharges(rs_grid, num_charges, positions, charges, 0);
    switchTimer(CA, GATHER);
    /* gather the ca grid */
    this->gatherGrid(rs_grid);
    /* now rs_grid should contain the local ca grid */
    stopTimer(GATHER);

    /* forward transform */

    P3M_DEBUG(printf( "  calling fft.forward()...\n"));
    startTimer(FORWARD);
    fft.forward(rs_grid, buffer);
    stopTimer(FORWARD);
    P3M_DEBUG(printf("  returned from fft.forward().\n"));

    /********************************************/
    /* POTENTIAL COMPUTATION */
    /********************************************/
    if (require_total_energy || potentials != NULL) {
        /* apply energy optimized influence function */
        startTimer(INFLUENCE);
        this->applyInfluenceFunction(rs_grid, ks_grid, g_energy);
        /* result is in ks_grid */
        stopTimer(INFLUENCE);

        /* compute total energy, but not potentials */
        if (require_total_energy && potentials == NULL) {
            startTimer(POTENTIALS);
            total_energy = this->computeTotalEnergy();
            stopTimer(POTENTIALS);
        }

        if (potentials != NULL) {
            /* backtransform the grid */

            if (require_timings != ESTIMATE_ALL
                    && require_timings != ESTIMATE_FFT){
                P3M_DEBUG(printf( "  calling fft.backward (potentials)...\n"));
                startTimer(BACK);
                fft.backward(ks_grid, buffer);
                stopTimer(BACK);
                P3M_DEBUG(printf( "  returned from fft.backward.\n"));
            }
            /* redistribute energy grid */
            startTimer(SPREAD);
            this->spreadGrid(ks_grid);

            switchTimer(SPREAD, POTENTIALS);

            /* compute potentials */
            this->assignPotentials(ks_grid,
                    num_charges, positions, charges, 0, potentials);

            stopTimer(POTENTIALS);
        }
    }

    /********************************************/
    /* FIELD COMPUTATION */
    /********************************************/
    if (fields != NULL) {
        /* apply force optimized influence function */
        startTimer(INFLUENCE);
        this->applyInfluenceFunction(rs_grid, ks_grid, g_force);
        stopTimer(INFLUENCE);

        /* result is in ks_grid */
        for (int dim = 0; dim < 3; dim++) {
            /* differentiate in direction dim */
            /* result is stored in rs_grid */
            startTimer(FIELDS);
            this->differentiateIK(dim, ks_grid, rs_grid);
            stopTimer(FIELDS);

            if (require_timings != ESTIMATE_ALL
                    && require_timings != ESTIMATE_FFT) {

                /* backtransform the grid */
                P3M_DEBUG(printf( "  calling fft.backward (field dim=%d)...\n", dim));
                startTimer(BACK);
                fft.backward(rs_grid, buffer);
                stopTimer(BACK);
                P3M_DEBUG(printf("  returned from fft.backward.\n"));
            }

            /* redistribute force grid */
            startTimer(SPREAD);
            this->spreadGrid(rs_grid);
            switchTimer(SPREAD, FIELDS);
            this->assignFieldsIK(rs_grid, dim, num_charges, positions, 0, fields);
            stopTimer(FIELDS);
        }
    }

    /* estimate FFT back timing*/
    if (require_timings == ESTIMATE_ALL || require_timings == ESTIMATE_FFT) {
        if (potentials != NULL)
            timings[BACK] = 4 * timings[FORWARD];
        else
            timings[BACK] = 3 * timings[FORWARD];
    }

    this->gatherTimings();

    P3M_INFO(printf( "P3M::FarSolver::runIK() finished.\n"));
}

/** Calculates the properties of the send/recv sub-grides of the local
 *  FFT grid.  In order to calculate the recv sub-grides there is a
 *  communication of the margins between neighbouring nodes. */
void P3M::FarSolver::prepareSendGrid() {
    p3m_int done[3]={0,0,0};

    /* send grids */
    for (int i=0; i<3; i++) {
        for (int j=0; j<3; j++) {
            /* left */
            sm.s_ld[i*2][j] = 0 + done[j]*local_grid.margin[j*2];
            if (j==i) sm.s_ur[i*2][j] = local_grid.margin[j*2];
            else     sm.s_ur[i*2][j] = local_grid.dim[j]-done[j]*local_grid.margin[(j*2)+1];
            /* right */
            if (j==i) sm.s_ld[(i*2)+1][j] = local_grid.in_ur[j];
            else     sm.s_ld[(i*2)+1][j] = 0 + done[j]*local_grid.margin[j*2];
            sm.s_ur[(i*2)+1][j] = local_grid.dim[j] - done[j]*local_grid.margin[(j*2)+1];
        }
        done[i]=1;
    }

    sm.max=0;
    for (int i=0; i<6; i++) {
        sm.s_size[i] = 1;
        for (int j=0; j<3; j++) {
            sm.s_dim[i][j] = sm.s_ur[i][j]-sm.s_ld[i][j];
            sm.s_size[i] *= sm.s_dim[i][j];
        }
        if (sm.s_size[i]>sm.max) sm.max=sm.s_size[i];
    }

    /* communication */
    for (int i=0; i<6; i++) {
        int j = ((i%2==0) ? i+1 : i-1);
        if (comm.node_neighbors[i] != comm.rank) {
            P3M_DEBUG(printf("      %d: sending local_grid.margin to %d, receiving from %d\n", \
                    comm.rank, comm.node_neighbors[i], comm.node_neighbors[j]));
            /* two step communication: first all even positions than all odd */
            MPI_Sendrecv(&(local_grid.margin[i]), 1, P3M_MPI_INT,
                    comm.node_neighbors[i], 0,
                    &(local_grid.r_margin[j]), 1, P3M_MPI_INT,
                    comm.node_neighbors[j], 0, comm.mpicomm, MPI_STATUS_IGNORE);
        } else {
            local_grid.r_margin[j] = local_grid.margin[i];
        }
    }

    /* recv grids */
    for (int i=0; i<3; i++)
        for (int j=0; j<3; j++) {
            if (j==i) {
                sm.r_ld[ i*2   ][j] = sm.s_ld[ i*2   ][j] + local_grid.margin[2*j];
                sm.r_ur[ i*2   ][j] = sm.s_ur[ i*2   ][j] + local_grid.r_margin[2*j];
                sm.r_ld[(i*2)+1][j] = sm.s_ld[(i*2)+1][j] - local_grid.r_margin[(2*j)+1];
                sm.r_ur[(i*2)+1][j] = sm.s_ur[(i*2)+1][j] - local_grid.margin[(2*j)+1];
            } else {
                sm.r_ld[ i*2   ][j] = sm.s_ld[ i*2   ][j];
                sm.r_ur[ i*2   ][j] = sm.s_ur[ i*2   ][j];
                sm.r_ld[(i*2)+1][j] = sm.s_ld[(i*2)+1][j];
                sm.r_ur[(i*2)+1][j] = sm.s_ur[(i*2)+1][j];
            }
        }
    for (int i=0; i<6; i++) {
        sm.r_size[i] = 1;
        for (int j=0; j<3; j++) {
            sm.r_dim[i][j] = sm.r_ur[i][j]-sm.r_ld[i][j];
            sm.r_size[i] *= sm.r_dim[i][j];
        }
        if (sm.r_size[i]>sm.max) sm.max=sm.r_size[i];
    }
}

/** Calculates properties of the local FFT grid for the
         charge assignment process. */
void P3M::FarSolver::prepareLocalCAGrid() {
    p3m_int ind[3];
    /* total skin size */
    p3m_float full_skin[3];

    for (int i=0;i<3;i++)
        full_skin[i]= cao_cut[i]+skin+additional_grid[i];

    /* inner bottom left grid point (global index) */
    for (int i=0;i<3;i++)
        local_grid.in_ld[i] =
                (p3m_int)ceil(comm.my_left[i]*ai[i]-grid_off[i]);

    /* inner upper right grid point (global index) */
    for (int i=0;i<3;i++)
        local_grid.in_ur[i] =
                (p3m_int)floor(comm.my_right[i]*ai[i]-grid_off[i]);

    /* correct for roundoff errors at boundary */
    for (int i=0;i<3;i++) {
        if (float_is_zero((comm.my_right[i] * ai[i] - grid_off[i])
                - local_grid.in_ur[i]))
            local_grid.in_ur[i]--;
        if (float_is_zero(1.0+(comm.my_left[i] * ai[i] - grid_off[i])
                - local_grid.in_ld[i]))
            local_grid.in_ld[i]--;
    }

    /* inner grid dimensions */
    for (int i=0; i<3; i++)
        local_grid.inner[i] = local_grid.in_ur[i] - local_grid.in_ld[i] + 1;

    /* index of the bottom left grid point. */
    for (int i=0; i<3; i++)
        local_grid.ld_ind[i] = static_cast<p3m_int>(
                ceil((comm.my_left[i]-full_skin[i])*ai[i]-grid_off[i])
        );

    /* spatial position of the bottom left grid point. */
    for (int i=0;i<3;i++)
        local_grid.ld_pos[i] =
                (local_grid.ld_ind[i]+ grid_off[i])*a[i];

    /* bottom left margin */
    for (int i=0;i<3;i++)
        local_grid.margin[i*2] = local_grid.in_ld[i]-local_grid.ld_ind[i];

    /* upper right grid point */
    for (int i=0;i<3;i++)
        ind[i] = (p3m_int)floor((comm.my_right[i]+full_skin[i])*ai[i]-grid_off[i]);

    /* correct roundoff errors at upper right boundary */
    for (int i=0;i<3;i++)
        if (((comm.my_right[i]+full_skin[i])*ai[i]-grid_off[i])-ind[i]==0)
            ind[i]--;

    /* upper right margin */
    for (int i=0;i<3;i++)
        local_grid.margin[(i*2)+1] = ind[i] - local_grid.in_ur[i];

    /* grid dimension */
    local_grid.size=1;
    for (int i=0;i<3;i++) {
        local_grid.dim[i] = ind[i] - local_grid.ld_ind[i] + 1;
        local_grid.size *= local_grid.dim[i];
    }

    /* reduce inner grid indices from global to local */
    for (int i=0;i<3;i++) {
        local_grid.in_ld[i] = local_grid.margin[i*2];
        local_grid.in_ur[i] = local_grid.margin[i*2]+local_grid.inner[i];
    }

    local_grid.q_2_off  = local_grid.dim[2] - cao;
    local_grid.q_21_off = local_grid.dim[2] * (local_grid.dim[1] - cao);}

/** Debug function printing p3m structures */
void P3M::FarSolver::printLocalGrid() {
    printf( "    local_grid:\n");
    printf( "      dim=" F3INT ", size=" FINT "\n",
            local_grid.dim[0], local_grid.dim[1], local_grid.dim[2], local_grid.size);
    printf("      ld_ind=" F3INT ", ld_pos=" F3FLOAT "\n",
            local_grid.ld_ind[0],local_grid.ld_ind[1],local_grid.ld_ind[2],
            local_grid.ld_pos[0],local_grid.ld_pos[1],local_grid.ld_pos[2]);
    printf("      inner=" F3INT "[" F3INT "-" F3INT "]\n",
            local_grid.inner[0],local_grid.inner[1],local_grid.inner[2],
            local_grid.in_ld[0],local_grid.in_ld[1],local_grid.in_ld[2],
            local_grid.in_ur[0],local_grid.in_ur[1],local_grid.in_ur[2]);
    printf("      margin=(" FINT "," FINT " ," FINT "," FINT " ," FINT "," FINT ")\n",
            local_grid.margin[0],local_grid.margin[1],local_grid.margin[2],
            local_grid.margin[3],local_grid.margin[4],local_grid.margin[5]);
    printf("      r_margin=(" FINT "," FINT " ," FINT "," FINT " ," FINT "," FINT ")\n",
            local_grid.r_margin[0],local_grid.r_margin[1],
            local_grid.r_margin[2],local_grid.r_margin[3],
            local_grid.r_margin[4],local_grid.r_margin[5]);
}

/** Debug function printing p3m structures */
void P3M::FarSolver::printSendGrid() {
    printf( "    send_grid:\n");
    printf( "      max=%d\n",sm.max);
    for (int i=0;i<6;i++) {
        printf("      dir=%d: s_dim (%d,%d,%d)  s_ld (%d,%d,%d) s_ur (%d,%d,%d) s_size=%d\n",
                i, sm.s_dim[i][0], sm.s_dim[i][1], sm.s_dim[i][2],
                sm.s_ld[i][0], sm.s_ld[i][1], sm.s_ld[i][2],
                sm.s_ur[i][0], sm.s_ur[i][1], sm.s_ur[i][2], sm.s_size[i]);
        printf("             r_dim (%d,%d,%d)  r_ld (%d,%d,%d) r_ur (%d,%d,%d) r_size=%d\n",
                sm.r_dim[i][0], sm.r_dim[i][1], sm.r_dim[i][2],
                sm.r_ld[i][0], sm.r_ld[i][1], sm.r_ld[i][2],
                sm.r_ur[i][0], sm.r_ur[i][1], sm.r_ur[i][2], sm.r_size[i]);
    }

}

/** Calculates the optimal influence function of Hockney and Eastwood
 * for force calculations and energy calculations.
 *
 *  Each node calculates only the values for its domain in k-space
 *  (see fft.plan[3].grid and fft.plan[3].start).
 *
 *  See also: Hockney/Eastwood 8-22 (p275). Note the somewhat
 *  different convention for the prefactors, which is described in
 *  Deserno/Holm. */
/** Calculates the optimal influence function of Hockney and Eastwood
 * for force calculations and energy calculations.
 *
 *  Each node calculates only the values for its domain in k-space
 *  (see fft.plan[3].grid and fft.plan[3].start).
 *
 *  See also: Hockney/Eastwood 8-22 (p275). Note the somewhat
 *  different convention for the prefactors, which is described in
 *  Deserno/Holm. */
void P3M::FarSolver::computeInfluenceFunctionIK() {
    P3M_DEBUG(printf("  FarSolver::computeInfluenceFunctionIK() started...\n"));

    int size = 1;
    int end[3];
    const p3m_int *extent;
    const p3m_int *offset;
    fft.getKSExtent(offset, extent);
    for (p3m_int i=0;i<3;i++) {
        size *= extent[i];
        end[i] = offset[i] + extent[i];
    }

    g_force = new p3m_float[size];
    g_energy = new p3m_float[size];

    p3m_int *gridshift_x = this->computeGridShift(RX, grid[0]);
    p3m_int *gridshift_y = this->computeGridShift(RY, grid[1]);
    p3m_int *gridshift_z = this->computeGridShift(RZ, grid[2]);

    int n[3];
    for (n[0]=offset[0]; n[0]<end[0]; n[0]++) {
        for (n[1]=offset[1]; n[1]<end[1]; n[1]++) {
            for (n[2]=offset[2]; n[2]<end[2]; n[2]++) {
                const p3m_int ind = (n[2]-offset[2]) + extent[2] *
                        ((n[1]-offset[1]) + extent[1]*(n[0]-offset[0]));
                if ((n[KX]%(grid[RX]/2)==0) &&
                        (n[KY]%(grid[RY]/2)==0) &&
                        (n[KZ]%(grid[RZ]/2)==0) ) {
                    g_force[ind] = 0.0;
                    g_energy[ind] = 0.0;
                } else {
                    p3m_float numerator_force[3];
                    p3m_float numerator_energy;
                    p3m_float sumU2;// up to here comparable to AD
                    this->performAliasingSumsIK(
                            gridshift_x[n[KX]],
                            gridshift_y[n[KY]],
                            gridshift_z[n[KZ]],
                            numerator_force,
                            numerator_energy,
                            sumU2);

                    p3m_float fak1; // k scalar numerator force (disregarding prefactors)
                    p3m_float fak2; // k squared (disregarding prefactors)
                    if (!isTriclinic) {
                        fak1 = d_op[RX][n[KX]] * numerator_force[RX] / box_l[RX] +
                                d_op[RY][n[KY]] * numerator_force[RY] / box_l[RY] +
                                d_op[RZ][n[KZ]] * numerator_force[RZ] / box_l[RZ];
                        fak2 = SQR(d_op[RX][n[KX]] / box_l[RX]) +
                                SQR(d_op[RY][n[KY]] / box_l[RY]) +
                                SQR(d_op[RZ][n[KZ]] / box_l[RZ]);
                        } else{
                        fak1=numerator_force[RX]*(d_op[RX][n[KX]]*(box_vectors[1][1]*box_vectors[2][2]))
                                +numerator_force[RY]*(d_op[RX][n[KX]]*(-box_vectors[1][0]*box_vectors[2][2])+d_op[RY][n[KY]]*(box_vectors[2][2]*box_vectors[0][0]))
                                +numerator_force[RZ]*(d_op[RX][n[KX]]*(box_vectors[1][0]*box_vectors[2][1]-box_vectors[1][1]*box_vectors[2][0])-d_op[RY][n[KY]]*box_vectors[2][1]*box_vectors[0][0]+d_op[RZ][n[KZ]]*(box_vectors[0][0]*box_vectors[1][1]));//correct k
                        fak1*=1/volume;
                        fak2 = SQR(d_op[RX][n[KX]]*(box_vectors[1][1]*box_vectors[2][2]))
                                +SQR(d_op[RX][n[KX]]*(box_vectors[1][0]*box_vectors[2][2])+d_op[RY][n[KY]]*(box_vectors[2][2]*box_vectors[0][0]))
                                +SQR(d_op[RX][n[KX]]*(box_vectors[1][0]*box_vectors[2][1]-box_vectors[1][1]*box_vectors[2][0])-d_op[RY][n[KY]]*box_vectors[2][1]*box_vectors[0][0]+d_op[RZ][n[KZ]]*(box_vectors[0][0]*box_vectors[1][1]));//correct k
                        fak2*=1/SQR(volume);
                        }
                    const p3m_float fak3 = fak1/(fak2 * SQR(sumU2));
                    g_force[ind] = M_2_PI*fak3;
                    g_energy[ind] = 0.5 * g_force[ind];
                    /* g_energy[ind] = M_1_PI*numerator_energy/SQR(denominator); */
                    /* if (fabs(g_energy[ind]) > 1.e-5) */
                    /*   printf("%d: %15.10e\n", ind, g_energy[ind]/g_force[ind]); */
                }
            }
        }
    }

    delete gridshift_x;
    delete gridshift_y;
    delete gridshift_z;

    P3M_DEBUG(printf("  FarSolver::computeInfluenceFunctionIK() finished.\n"));
}

void P3M::FarSolver::computeInfluenceFunctionADI() {
    P3M_DEBUG(printf("  FarSolver::computeInfluenceFunctionADI() started...\n"));
    int size = 1;
    p3m_int end[3];
    const p3m_int *start, *extent;
    fft.getKSExtent(start, extent);
    for (p3m_int i=0;i<3;i++) {
      size *= extent[i];
      end[i] = start[i] + extent[i];
    }

    g_force = new p3m_float[size];
    g_energy = new p3m_float[size];

    p3m_int *gridshift_x = this->computeGridShift(RX, grid[0]);
    p3m_int *gridshift_y = this->computeGridShift(RY, grid[1]);
    p3m_int *gridshift_z = this->computeGridShift(RZ, grid[2]);

    int n[3];
    for (n[0]=start[0]; n[0]<end[0]; n[0]++) {
      for (n[1]=start[1]; n[1]<end[1]; n[1]++) {
        for (n[2]=start[2]; n[2]<end[2]; n[2]++) {
          p3m_int ind = (n[2]-start[2]) + extent[2] *
            ((n[1]-start[1]) + extent[1]*(n[0]-start[0]));
          if ((n[KX]%(grid[RX]/2)==0) &&
              (n[KY]%(grid[RY]/2)==0) &&
              (n[KZ]%(grid[RZ]/2)==0) ) {
            g_force[ind] = 0.0;
            g_energy[ind] = 0.0;
          } else {
            p3m_float numerator_force;
            p3m_float numerator_energy;
            p3m_float denominator[4];
            this->performAliasingSumsADI(
                    gridshift_x[n[KX]],
                    gridshift_y[n[KY]],
                    gridshift_z[n[KZ]],
                    numerator_force, numerator_energy,
                    denominator);

            g_force[ind] = numerator_force /
              (0.5 * M_PI * (denominator[0] * denominator[1] + denominator[2] * denominator[3] )) ;
            g_energy[ind] = M_1_PI*numerator_energy/SQR(denominator[0]);
            /* g_energy[ind] = 0.5 * g_force[ind]; */
          }
        }
      }
    }

    delete gridshift_x;
    delete gridshift_y;
    delete gridshift_z;

    P3M_DEBUG(printf("  FarSolver::computeInfluenceFunctionADI() finished.\n"));
}

void P3M::FarSolver::computeInfluenceFunctionIKI() {
    P3M_DEBUG(printf("  FarSolver::computeInfluenceFunctionIKI() started...\n"));
    int size = 1;
    p3m_int end[3];
    const p3m_int *start, *extent;
    fft.getKSExtent(start, extent);
    for (p3m_int i=0;i<3;i++) {
      size *= extent[i];
      end[i] = start[i] + extent[i];
    }

    g_force = new p3m_float[size];
    g_energy = new p3m_float[size];

    p3m_int *gridshift_x = this->computeGridShift(RX, grid[0]);
    p3m_int *gridshift_y = this->computeGridShift(RY, grid[1]);
    p3m_int *gridshift_z = this->computeGridShift(RZ, grid[2]);

    int n[3];
    for (n[0]=start[0]; n[0]<end[0]; n[0]++) {
      for (n[1]=start[1]; n[1]<end[1]; n[1]++) {
        for (n[2]=start[2]; n[2]<end[2]; n[2]++) {
          p3m_int ind = (n[2]-start[2]) + extent[2] *
            ((n[1]-start[1]) + extent[1]*(n[0]-start[0]));
          if ((n[KX]%(grid[RX]/2)==0) &&
              (n[KY]%(grid[RY]/2)==0) &&
              (n[KZ]%(grid[RZ]/2)==0) ) {
            g_force[ind] = 0.0;
            g_energy[ind] = 0.0;
          } else {
            p3m_float numerator_force[3];
            p3m_float numerator_energy;
            p3m_float denominator[2];
            this->performAliasingSumsIKI(
                    gridshift_x[n[KX]],
                    gridshift_y[n[KY]],
                    gridshift_z[n[KZ]],
                    numerator_force, numerator_energy,
                    denominator);

            p3m_float fak1 =
              d_op[RX][n[KX]]*numerator_force[RX]/box_l[RX] +
              d_op[RY][n[KY]]*numerator_force[RY]/box_l[RY] +
              d_op[RZ][n[KZ]]*numerator_force[RZ]/box_l[RZ];
            p3m_float fak2 =
              SQR(d_op[RX][n[KX]]/box_l[RX]) +
              SQR(d_op[RY][n[KY]]/box_l[RY]) +
              SQR(d_op[RZ][n[KZ]]/box_l[RZ]);
            p3m_float fak3 = fak1/(fak2 * 0.5 * (SQR(denominator[0]) + SQR(denominator[1])) );
            g_force[ind] = M_2_PI*fak3;
            g_energy[ind] = M_1_PI*numerator_energy/SQR(denominator[0]);
          }
        }
      }
    }

    delete gridshift_x;
    delete gridshift_y;
    delete gridshift_z;
    P3M_DEBUG(printf("  FarSolver::computeInfluenceFunctionIKI() finished.\n"));
}

void P3M::FarSolver::performAliasingSumsIK(p3m_int nmx0, p3m_int nmy0,
        p3m_int nmz0, p3m_float numerator_force[3], p3m_float& numerator_energy,
        p3m_float& sumU2)
        {
    sumU2 = 0.0;
    numerator_energy = 0.0;
    for (p3m_int i = 0; i < 3; i++)
        numerator_force[i] = 0.0;

    p3m_float prefactor = SQR(M_PI/alpha);

    for (p3m_int mx = -P3M_BRILLOUIN; mx <= P3M_BRILLOUIN; mx++) {
        const p3m_int nmx = nmx0 + grid[RX]*mx;
        const p3m_float sx  = pow(sinc(nmx/(p3m_float)grid[RX]),2.0*cao);
        for (p3m_int my = -P3M_BRILLOUIN; my <= P3M_BRILLOUIN; my++) {
            const p3m_int nmy = nmy0 + grid[RY]*my;
            const p3m_float sy  = sx*pow(sinc(nmy/(p3m_float)grid[RY]),2.0*cao);
            for (p3m_int mz = -P3M_BRILLOUIN; mz <= P3M_BRILLOUIN; mz++) {
                const p3m_int nmz = nmz0 + grid[RZ]*mz;
                const p3m_float U2 = sy*pow(sinc(nmz/(p3m_float)grid[RZ]),2.0*cao);

                p3m_float nm2 = 0.0;
         if (!isTriclinic) {
            nm2 =
                    SQR(nmx / box_l[RX]) +
                    SQR(nmy / box_l[RY]) +
                    SQR(nmz / box_l[RZ]);

            } else {
                p3m_int i;
                const p3m_int nNm[3] = {nmx, nmy, nmz};

            for (i = 0; i < 3; i++) {
                int j = (i + 1) % 3;
                int k = (i + 2) % 3;

            nm2 += SQR((nNm[i]*(box_vectors[j][j] * box_vectors[k][k]\
             - box_vectors[j][k] * box_vectors[k][j]) + nNm[j]\
             *(box_vectors[k][j] * box_vectors[i][k] - box_vectors[k][k]\
             * box_vectors[i][j]) + nNm[k]*(box_vectors[i][j]\
             * box_vectors[j][k] - box_vectors[i][k] * box_vectors[j][j])));

                        }
                        nm2 *= 1 / SQR(volume);
                    }
                const p3m_float prefactor2 = U2*exp(-prefactor*nm2)/(nm2)/(isTriclinic?volume:1.0);

                numerator_energy += prefactor2;
                if(!isTriclinic){//numerator_force = SUM prefactor2 * k/2pi =SUM U^2 exp(-pi^2/alpha^2 k^2 / 4 pi^2)) / (k^2 / (4 pi^2)) * k/2pi;
                numerator_force[RX] += prefactor2*nmx/box_l[RX];
                numerator_force[RY] += prefactor2*nmy/box_l[RY];
                numerator_force[RZ] += prefactor2*nmz/box_l[RZ];
                }else{
                numerator_force[RX] += prefactor2*(nmx*(box_vectors[1][1]*box_vectors[2][2]))/volume;
                numerator_force[RY] += prefactor2*(nmx*(-box_vectors[1][0]*box_vectors[2][2])+nmy*(box_vectors[2][2]*box_vectors[0][0]))/volume;
                numerator_force[RZ] += prefactor2*(nmx*(box_vectors[1][0]*box_vectors[2][1]-box_vectors[1][1]*box_vectors[2][0])+nmy*(-box_vectors[2][1]*box_vectors[0][0])+nmz*(box_vectors[0][0]*box_vectors[1][1]))/volume;//correct k
                }
                sumU2 += U2; // denominator = SUM U^2
            }
        }
    }
}

void P3M::FarSolver::performAliasingSumsADI(
        p3m_int nmx0, p3m_int nmy0, p3m_int nmz0,
        p3m_float &numerator_force,
        p3m_float &numerator_energy, p3m_float denominator[4]) {
   denominator[0] = denominator[1] = denominator[2] = denominator[3] = 0.0;
   numerator_energy = 0.0;
   numerator_force = 0.0;

   p3m_float prefactor = SQR(M_PI/alpha);

   for (p3m_int mx = -P3M_BRILLOUIN; mx <= P3M_BRILLOUIN; mx++) {
     const p3m_int nmx = nmx0 + grid[RX]*mx;
     const p3m_float sx = pow(sinc(nmx/(p3m_float)grid[RX]), 2.0*cao);
     for (p3m_int my = -P3M_BRILLOUIN; my <= P3M_BRILLOUIN; my++) {
       const p3m_int nmy = nmy0 + grid[RY]*my;
       const p3m_float sy = sx*pow(sinc(nmy/(p3m_float)grid[RY]), 2.0*cao);
       for (p3m_int mz = -P3M_BRILLOUIN; mz <= P3M_BRILLOUIN; mz++) {
         const p3m_int nmz = nmz0 + grid[RZ]*mz;
         const p3m_float U2 = sy*pow(sinc(nmz/(p3m_float)grid[RZ]), 2.0*cao);
         
         p3m_float nm2 = 0.0;
         if (!isTriclinic) {
             nm2 =
                     SQR(nmx / box_l[RX]) +
                     SQR(nmy / box_l[RY]) +
                     SQR(nmz / box_l[RZ]);
         } else {
             p3m_int i;
             const p3m_int nNm[3] = {nmx, nmy, nmz};

             for (i = 0; i < 3; i++) {
                 int j = (i + 1) % 3;
                 int k = (i + 2) % 3;

                 nm2 += SQR((nNm[i]*(box_vectors[j][j] * box_vectors[k][k]\
                 - box_vectors[j][k] * box_vectors[k][j]) + nNm[j]\
                 *(box_vectors[k][j] * box_vectors[i][k] - box_vectors[k][k]\
                 * box_vectors[i][j]) + nNm[k]*(box_vectors[i][j]\
                 * box_vectors[j][k] - box_vectors[i][k] * box_vectors[j][j])));
            }
             
            nm2 *= 1 / SQR(volume);

        }
        const p3m_float prefactor2 = U2 * exp(-prefactor * nm2)
                / (isTriclinic ? volume:1.0);

         numerator_energy += prefactor2/nm2;
         numerator_force += prefactor2;

         denominator[0] += U2;
         denominator[1] += U2 * nm2;

         if(((mx+my+mz) % 2) == 0) {
           denominator[2] += U2;
           denominator[3] += U2 * nm2;
         } else {
           denominator[2] -= U2;
           denominator[3] -= U2 * nm2;
         }
       }
     }
   }
 }

void P3M::FarSolver::performAliasingSumsIKI(p3m_int nmx0, p3m_int nmy0,
        p3m_int nmz0, p3m_float numerator_force[3], p3m_float& numerator_energy,
        p3m_float denominator[2]) {
    denominator[0] = denominator[1] = 0.0;
    numerator_energy = 0.0;
    for (p3m_int i = 0; i < 3; i++)
        numerator_force[i] = 0.0;

    const p3m_float prefactor = SQR(M_PI/alpha);

    for (p3m_int mx = -P3M_BRILLOUIN; mx <= P3M_BRILLOUIN; mx++) {
        const p3m_int nmx = nmx0 + grid[RX]*mx;
        const p3m_float sx  = pow(sinc(nmx/(p3m_float)grid[RX]),2.0*cao);
        for (p3m_int my = -P3M_BRILLOUIN; my <= P3M_BRILLOUIN; my++) {
            const p3m_int nmy = nmy0 + grid[RY]*my;
            const p3m_float sy  = sx*pow(sinc(nmy/(p3m_float)grid[RY]),2.0*cao);
            for (p3m_int mz = -P3M_BRILLOUIN; mz <= P3M_BRILLOUIN; mz++) {
                const p3m_int nmz = nmz0 + grid[RZ]*mz;
                const p3m_float sz  = sy*pow(sinc(nmz/(p3m_float)grid[RZ]),2.0*cao);
                const p3m_float nm2 =
                        SQR(nmx/box_l[RX]) +
                        SQR(nmy/box_l[RY]) +
                        SQR(nmz/box_l[RZ]);
                const p3m_float prefactor2 = sz*exp(-prefactor*nm2)/nm2;

                numerator_energy += prefactor2;
                numerator_force[RX] += prefactor2*nmx/box_l[RX];
                numerator_force[RY] += prefactor2*nmy/box_l[RY];
                numerator_force[RZ] += prefactor2*nmz/box_l[RZ];

                denominator[0]  += sz;

                if(((mx+my+mz) % 2) == 0)
                    denominator[1] += sz;
                else
                    denominator[1] -= sz;

            }
        }
    }
}

/*transforms the field computation results to cartesian fields in ADI*/
void P3M::FarSolver::cartesianizeFields(p3m_float *fields, p3m_int num_particles){
    p3m_int part_no;
    
    for (part_no = 0; part_no < num_particles; part_no++) {
        fields[3 * part_no + 2] =
            (box_vectors[1][0] * box_vectors[2][1] - box_vectors[1][1] * box_vectors[2][0])
            / (box_vectors[0][0] * box_vectors[1][1] * box_vectors[2][2])
            * fields[3 * part_no]
            -(box_vectors[2][1]) / (box_vectors[1][1] * box_vectors[2][2])
            * fields[3 * part_no + 1]
            + 1 / box_vectors[2][2] * fields[3 * part_no + 2];
   
        fields[3 * part_no + 1] =
            -(box_vectors[1][0]) / (box_vectors[0][0] * box_vectors[1][1])
            * fields[3 * part_no]
            + 1 / box_vectors[1][1] * fields[3 * part_no + 1];
        
        fields[3 * part_no] = 1 / (box_vectors[0][0]) * fields[3 * part_no];
    }
}

/** shifts the grid points by grid/2 */
p3m_int*
P3M::FarSolver::computeGridShift(int dir, p3m_int size) {
    p3m_int *gridshift = new p3m_int[size];
    gridshift[0] = 0;
    for (int i = 1; i <= grid[dir]/2; i++) {
        gridshift[i] = i;
        gridshift[grid[dir] - i] = -i;
    }
    return gridshift;
}

/** Compute the data of the charge assignment grid points.

      The function returns the linear index of the top left grid point
      in the charge assignment grid that corresponds to real_pos. When
      "shifted" is set, it uses the shifted position for interlacing.
      After the call, caf_cache contains a cache of the values of the
      charge assignment fraction (caf) for x,y,z.
 */
p3m_int P3M::FarSolver::getCAPoints(p3m_float real_pos[3], p3m_int shifted) {
    /* linear index of the grid point */
    p3m_int linind = 0;

    for (p3m_int dim=0; dim<3; dim++) {
        /* position in normalized coordinates in [0,1] */
        p3m_float pos = (real_pos[dim] - local_grid.ld_pos[dim]) * ai[dim];
        /* shift position to the corner of the charge assignment area */
        pos -= pos_shift;
        /* if using the interlaced grid, shift it more */
        if (shifted) pos -= 0.5;
        /* nearest grid point in the ca grid */
        p3m_int grid_ind  = (p3m_int)floor(pos);
#ifdef ADDITIONAL_CHECKS
        if (grid_ind < 0) {
            printf("grid_ind[%d]=%d < 0\n", dim, grid_ind);
            printf("pos=%lf, pos_shift=%lf, real_pos=%lf, ld_pos=%lf, ai=%lf, shifted=%d\n",
                    pos, pos_shift, real_pos[dim], local_grid.ld_pos[dim],  ai[dim],
                    shifted);
        } else if (grid_ind > local_grid.dim[dim]) {
            printf("grid_ind[%d]=%d > %d\n", dim, grid_ind, local_grid.dim[dim]);
            printf("pos=%lf, pos_shift=%lf, real_pos=%lf, ld_pos=%lf, ai=%lf, shifted=%d\n",
                    pos, pos_shift, real_pos[dim], local_grid.ld_pos[dim],  ai[dim],
                    shifted);
        }
#endif
        /* linear index of grid point */
        linind = grid_ind + linind*local_grid.dim[dim];
        /* normalized distance to grid point */
        p3m_float dist = (pos-grid_ind)-0.5;

        switch (dim) {
        case 0: cafx->update(dist); break;
        case 1: cafy->update(dist); break;
        case 2: cafz->update(dist); break;
        }

#ifdef P3M_AD
        switch (dim) {
        case 0: cafx_d->update(dist); break;
        case 1: cafy_d->update(dist); break;
        case 2: cafz_d->update(dist); break;
        }
#endif

#ifdef ADDITIONAL_CHECKS
        if (real_pos[dim] < comm.my_left[dim]
                                            || real_pos[dim] > comm.my_right[dim]) {
            printf("%d: dim %d: position not in domain! " F3FLOAT "\n",
                    comm.rank, dim, real_pos[0], real_pos[1], real_pos[2]);
        }
#endif
    }

#ifdef ADDITIONAL_CHECKS
    if (linind < 0) {
        printf("ERROR: %d: linind %d < 0\n", comm.rank, linind);
    } else if (linind >= local_grid.size) {
        printf("ERROR: %d: linind %d > %d\n", comm.rank, linind, local_grid.size);
    }
#endif

    return linind;
}

/** Assign the charges to the grid */
void P3M::FarSolver::assignCharges(p3m_float* data, p3m_int num_charges,
        p3m_float* positions, p3m_float* charges, p3m_int shifted) {
    P3M_DEBUG(printf( "  P3M::FarSolver::assignCharges() started...\n"));

    const p3m_int q2off = local_grid.q_2_off;
    const p3m_int q21off = local_grid.q_21_off;

    /* init local charge grid */
    for (p3m_int i=0; i<local_grid.size; i++) data[i] = 0.0;

    /* now assign the charges */
    for (p3m_int pid = 0; pid < num_charges; pid++) {
        const p3m_float q = charges[pid];
        p3m_int linind_grid = this->getCAPoints(&positions[pid*3], shifted);

        /* Loop over all ca grid points nearby and compute charge assignment fraction */
        for (p3m_float *caf_x = cafx->begin(); caf_x < cafx->end(); caf_x++) {
            for (p3m_float *caf_y = cafy->begin(); caf_y < cafy->end(); caf_y++) {
                p3m_float caf_xy = *caf_x * *caf_y;
                for (p3m_float *caf_z = cafz->begin(); caf_z < cafz->end(); caf_z++) {
                    /* add it to the grid */
                    data[linind_grid] += q * caf_xy * *caf_z;
                    linind_grid++;
                }
                linind_grid += q2off;
            }
            linind_grid += q21off;
        }
    }

    P3M_DEBUG(printf( "  P3M::FarSolver::assignCharges() finished...\n"));

}

/* Gather information for FFT grid inside the nodes domain (inner local grid) */
void P3M::FarSolver::gatherGrid(p3m_float* rs_grid) {
    P3M_DEBUG(printf( "  P3M::FarSolver::gatherGrid() started...\n"));

    /* direction loop */
    for (p3m_int s_dir=0; s_dir<6; s_dir++) {
        p3m_int r_dir;
        if(s_dir%2==0) r_dir = s_dir+1;
        else           r_dir = s_dir-1;
        /* pack send block */
        if (sm.s_size[s_dir]>0)
            Parallel3DFFT::pack_block(rs_grid, send_grid, sm.s_ld[s_dir],
                    sm.s_dim[s_dir], local_grid.dim, 1);

        /* communication */
        if (comm.node_neighbors[s_dir] != comm.rank)
            MPI_Sendrecv(send_grid, sm.s_size[s_dir], P3M_MPI_FLOAT,
                    comm.node_neighbors[s_dir], REQ_P3M_GATHER,
                    recv_grid, sm.r_size[r_dir], P3M_MPI_FLOAT,
                    comm.node_neighbors[r_dir], REQ_P3M_GATHER,
                    comm.mpicomm, MPI_STATUS_IGNORE);
        else std::swap(recv_grid, send_grid);

        /* add recv block */
        if(sm.r_size[r_dir]>0) {
            Parallel3DFFT::add_block(recv_grid, rs_grid, sm.r_ld[r_dir],
                    sm.r_dim[r_dir], local_grid.dim);
        }
    }

    P3M_DEBUG(printf( "  P3M::FarSolver::gatherGrid() finished.\n"));
}

void P3M::FarSolver::spreadGrid(p3m_float* rs_grid) {
    P3M_DEBUG(printf( "  P3M::FarSolver::spreadGrid() started...\n"));

    /* direction loop */
    for (int s_dir=5; s_dir>=0; s_dir--) {
        int r_dir;
        if (s_dir%2==0) r_dir = s_dir+1;
        else           r_dir = s_dir-1;
        /* pack send block */
        if(sm.s_size[s_dir]>0)
            Parallel3DFFT::pack_block(rs_grid, send_grid, sm.r_ld[r_dir],
                    sm.r_dim[r_dir], local_grid.dim, 1);
        /* communication */
        if (comm.node_neighbors[r_dir] != comm.rank)
            MPI_Sendrecv(send_grid, sm.r_size[r_dir], P3M_MPI_FLOAT,
                    comm.node_neighbors[r_dir], REQ_P3M_SPREAD,
                    recv_grid, sm.s_size[s_dir], P3M_MPI_FLOAT,
                    comm.node_neighbors[s_dir], REQ_P3M_SPREAD,
                    comm.mpicomm, MPI_STATUS_IGNORE);
        else std::swap(recv_grid, send_grid);

        /* unpack recv block */
        if (sm.s_size[s_dir]>0)
            Parallel3DFFT::unpack_block(recv_grid, rs_grid, sm.s_ld[s_dir],
                    sm.s_dim[s_dir], local_grid.dim, 1);
    }

    P3M_DEBUG(printf( "  P3M::FarSolver::spreadGrid() finished.\n"));

}

/* apply the influence function */
void P3M::FarSolver::applyInfluenceFunction(p3m_float* in, p3m_float* out,
        p3m_float* g) {
    const p3m_int size = fft.getKSSize();
    for (p3m_int i=0; i < size; i++) {
        out[2*i] = g[i] * in[2*i];
        out[2*i+1] = g[i] * in[2*i+1];
    }
}

void P3M::FarSolver::differentiateIK(int dim, p3m_float* in, p3m_float* out) {
    /* direction in k space: */
    p3m_int dim_rs = (dim+ks_pnum)%3;

    p3m_int* d_operator = NULL;
    switch (dim) {
    case KX:
        d_operator = d_op[RX];
        break;
    case KY:
        d_operator = d_op[RY];
        break;
    case KZ:
        d_operator = d_op[RZ];
    }

    /* srqt(-1)*k differentiation */
    p3m_int ind=0;
    p3m_int j[3];
    const p3m_int *size, *start;
    fft.getKSExtent(start, size);
    for (j[0]=0; j[0] < size[0]; j[0]++) {
        for (j[1]=0; j[1] < size[1]; j[1]++) {
            for (j[2]=0; j[2]< size[2]; j[2]++) {
                /* i*k*(Re+i*Im) = - Im*k + i*Re*k     (i=sqrt(-1)) */
                if (!isTriclinic) { //orthorhombic case
                out[ind] =
                        -2.0*M_PI*in[ind+1] *
                        d_operator[ j[dim] + start[dim] ] /
                        box_l[dim_rs];
                out[ind+1] =
                        2.0*M_PI*in[ind] *
                        d_operator[ j[dim] + start[dim] ] /
                        box_l[dim_rs];
                    } else { //triclinic case //k correct
                    //todo: differentiation wrong: mixing of directions is wrong-->completely wrong.
                        out [ind] = -2.0 * M_PI * in[ind + 1] * (d_operator[ j[dim]+start[dim] ]*(box_vectors[(dim+1)%3][(dim+1)%3]*box_vectors[(dim+2)%3][(dim+2)%3]-box_vectors[(dim+1)%3][(dim+2)%3]*box_vectors[(dim+2)%3][(dim+1)%3])+d_operator[ j[(dim+1)%3]+start[(dim+1)%3] ]*(box_vectors[(dim+2)%3][(dim+1)%3]*box_vectors[dim][(dim+2)%3]-box_vectors[(dim+2)%3][(dim+2)%3]*box_vectors[dim][(dim+1)%3])+d_operator[ j[(dim+2)%3]+start[(dim+2)%3] ]*(box_vectors[dim][(dim+1)%3]*box_vectors[(dim+1)%3][(dim+2)%3]-box_vectors[dim][(dim+2)%3]*box_vectors[(dim+1)%3][(dim+1)%3]))
                                /volume;
                     out[ind+1] =//k correct
                                     2.0*M_PI*in[ind] * (d_operator[ j[dim]+start[dim] ]*(box_vectors[(dim+1)%3][(dim+1)%3]*box_vectors[(dim+2)%3][(dim+2)%3]-box_vectors[(dim+1)%3][(dim+2)%3]*box_vectors[(dim+2)%3][(dim+1)%3])+d_operator[ j[(dim+1)%3]+start[(dim+1)%3] ]*(box_vectors[(dim+2)%3][(dim+1)%3]*box_vectors[dim][(dim+2)%3]-box_vectors[(dim+2)%3][(dim+2)%3]*box_vectors[dim][(dim+1)%3])+d_operator[ j[(dim+2)%3]+start[(dim+2)%3] ]*(box_vectors[dim][(dim+1)%3]*box_vectors[(dim+1)%3][(dim+2)%3]-box_vectors[dim][(dim+2)%3]*box_vectors[(dim+1)%3][(dim+1)%3]))
                                / volume;
                    }
                    ind += 2;
            }
        }
    }
}

/** Compute the total energy of the system in kspace. No need to
      backtransform the FFT grid in this case! */
p3m_float P3M::FarSolver::computeTotalEnergy() {
    p3m_float k_space_energy = 0.0;

    P3M_DEBUG(printf( "  P3M::FarSolver::computeTotalEnergy() started...\n"));

    p3m_int size = fft.getKSSize();
    for (p3m_int i=0; i < size; i++)
        /* Use the energy optimized influence function */
        k_space_energy += g_energy[i] * ( SQR(rs_grid[2*i]) +
                SQR(rs_grid[2*i+1]) );

    MPI_Reduce(MPI_IN_PLACE, &k_space_energy, 1, P3M_MPI_FLOAT,
            MPI_SUM, 0, comm.mpicomm);
    p3m_float prefactor = 1.0 / (2.0 * box_l[0] * box_l[1] * box_l[2]);
    k_space_energy *= prefactor;

#ifdef P3M_INTERLACE
    /* In the case of interlacing we have calculated the sum of the
       shifted and unshifted charges, we have to take the average. */
    k_space_energy *= 0.5;
#endif

    /* self energy correction */
    k_space_energy -= sum_q2 * alpha * 0.5*M_2_SQRTPI;
    /* net charge correction */
    k_space_energy -= square_sum_q * M_PI * prefactor / SQR(alpha);

    P3M_DEBUG(printf( "  P3M::FarSolver::computeTotalEnergy() finished.\n"));
    return k_space_energy;
}

/* Backinterpolate the potentials obtained from k-space to the positions */
void P3M::FarSolver::assignPotentials(p3m_float* data, p3m_int num_particles,
        p3m_float* positions, p3m_float* charges, p3m_int shifted,
        p3m_float* potentials) {
    const p3m_int q2off = local_grid.q_2_off;
    const p3m_int q21off = local_grid.q_21_off;
    const p3m_float prefactor = 1.0 / (box_l[0] * box_l[1] * box_l[2]);

    P3M_DEBUG(printf( "  P3M::FarSolver::assignPotentials() started...\n"));
    /* Loop over all particles */
    for (p3m_int pid=0; pid < num_particles; pid++) {
        p3m_float potential = 0.0;
        p3m_int linind_grid =
                this->getCAPoints(&positions[pid*3], shifted);

        /* Loop over all ca grid points nearby and compute charge assignment fraction */
        for (p3m_float *caf_x = cafx->begin(); caf_x < cafx->end(); caf_x++) {
            for (p3m_float *caf_y = cafy->begin(); caf_y < cafy->end(); caf_y++) {
                p3m_float caf_xy = *caf_x * *caf_y;
                for (p3m_float *caf_z = cafz->begin(); caf_z < cafz->end(); caf_z++) {
                    potential += *caf_z * caf_xy * data[linind_grid];
                    linind_grid++;
                }
                linind_grid += q2off;
            }
            linind_grid += q21off;
        }

        potential *= prefactor;
        /* self energy correction */
        potential -= charges[pid] * M_2_SQRTPI * alpha;
        /* net charge correction */
        /* potential -= fabs(charges[pid]) * PI * prefactor / SQR(alpha); */

        /* store the result */
        if (!shifted) {
            potentials[pid] = potential;
        } else {
            potentials[pid] = 0.5*(potentials[pid] + potential);
        }

    }
    P3M_DEBUG(printf( "  P3M::FarSolver::assignPotentials() finished.\n"));

}

/* Backinterpolate the forces obtained from k-space to the positions */
void P3M::FarSolver::assignFieldsIK(p3m_float* data, p3m_int dim,
        p3m_int num_particles, p3m_float* positions, p3m_int shifted,
        p3m_float* fields){
    const p3m_int q2off = local_grid.q_2_off;
    const p3m_int q21off = local_grid.q_21_off;
    const p3m_float prefactor = 1.0 / (2.0 * box_l[0] * box_l[1] * box_l[2]);
    const p3m_int dim_rs = (dim+ks_pnum) % 3;

    P3M_DEBUG(printf( "  P3M::FarSolver::assignFieldsIK() started...\n"));
    /* Loop over all particles */
    for (p3m_int pid=0; pid < num_particles; pid++) {
        p3m_float field = 0.0;
        p3m_int linind_grid =
                this->getCAPoints(&positions[3*pid], shifted);

        /* loop over the local grid, compute the field */
        for (p3m_float *caf_x = cafx->begin(); caf_x < cafx->end(); caf_x++) {
            for (p3m_float *caf_y = cafy->begin(); caf_y < cafy->end(); caf_y++) {
                p3m_float caf_xy = *caf_x * *caf_y;
                for (p3m_float *caf_z = cafz->begin(); caf_z < cafz->end(); caf_z++) {
                    field -= *caf_z * caf_xy * data[linind_grid];
                    linind_grid++;
                }
                linind_grid += q2off;
            }
            linind_grid += q21off;
        }

        field *= prefactor;

        /* store the result */
        if (!shifted)
            fields[3*pid + dim_rs] = field;
        else
            fields[3*pid + dim_rs] = 0.5 * (fields[3*pid + dim_rs] + field);
    }
    P3M_DEBUG(printf( "  P3M::FarSolver::assignFieldsIK() finished.\n"));

}

/* Backinterpolate the forces obtained from k-space to the positions */
void P3M::FarSolver::assignFieldsAD(p3m_float* data, p3m_int num_particles,
        p3m_float* positions, p3m_int shifted, p3m_float* fields) {
    const p3m_int q2off = local_grid.q_2_off;
    const p3m_int q21off = local_grid.q_21_off;
    const p3m_float prefactor = 1.0 / (box_l[0] * box_l[1] * box_l[2]);
    const p3m_float l_x_inv = 1.0/box_l[0];
    const p3m_float l_y_inv = 1.0/box_l[1];
    const p3m_float l_z_inv = 1.0/box_l[2];

    P3M_DEBUG(printf( "  P3M::Solver::assign_fields_ad() started...\n"));
    /* Loop over all particles */
    for (p3m_int pid = 0; pid < num_particles; pid++) {
        p3m_float field[3] = { 0.0, 0.0, 0.0 };
        p3m_int linind_grid =
                this->getCAPoints(&positions[pid*3], shifted);

        p3m_float *caf_x_d = cafx_d->begin();
        for (p3m_float *caf_x = cafx->begin(); caf_x < cafx->end(); caf_x++) {
            p3m_float *caf_y_d = cafy_d->begin();
            for (p3m_float *caf_y = cafy->begin(); caf_y < cafy->end(); caf_y++) {
                p3m_float *caf_z_d = cafz_d->begin();
                for (p3m_float *caf_z = cafz->begin(); caf_z < cafz->end(); caf_z++) {
                    field[0] -= *caf_x_d * *caf_y * *caf_z * l_x_inv
                            * data[linind_grid] * grid[0];
                    field[1] -= *caf_x * *caf_y_d * *caf_z * l_y_inv
                            * data[linind_grid] * grid[1];
                    field[2] -= *caf_x * *caf_y * *caf_z_d * l_z_inv
                            * data[linind_grid] * grid[2];
                    linind_grid++;
                    caf_z_d++;
                }
                linind_grid += q2off;
                caf_y_d++;
            }
            linind_grid += q21off;
            caf_x_d++;
        }
        field[0] *= prefactor;
        field[1] *= prefactor;
        field[2] *= prefactor;

        if (!shifted) {
            fields[3*pid + 0] = field[0];
            fields[3*pid + 1] = field[1];
            fields[3*pid + 2] = field[2];
        } else {
            fields[3*pid + 0] = 0.5*(fields[3*pid + 0] + field[0]);
            fields[3*pid + 1] = 0.5*(fields[3*pid + 1] + field[1]);
            fields[3*pid + 2] = 0.5*(fields[3*pid + 2] + field[2]);
        }
    }
    P3M_DEBUG(printf( "  assign_fields_ad() finished.\n"));
}

/** Calculate number of charged particles, the sum of the squared
      charges and the squared sum of the charges. Called in parallel at
      the beginning of tuning. */
void P3M::FarSolver::countCharges(p3m_int num_particles, p3m_float *charges) {
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

void P3M::FarSolver::setRequireTotalEnergy(bool flag) {
    require_total_energy = flag;
}

p3m_float P3M::FarSolver::getTotalEnergy() {
    if (!require_total_energy)
        throw std::logic_error("");
    return total_energy;
}

void P3M::FarSolver::setRequireTimings(TimingType type) {
    require_timings = type;
}

const double* P3M::FarSolver::measureTimings(p3m_int num_particles,
            p3m_float *positions, p3m_float *charges) {
    p3m_float *fields = new p3m_float[3*num_particles];
    p3m_float *potentials = new p3m_float[num_particles];

    /* store require_timings */
    TimingType require_timings_before = require_timings;
    if (require_timings == NONE || require_timings == FULL)
        require_timings = ESTIMATE_ALL;
#if defined(P3M_INTERLACE) && defined(P3M_AD)
    //todo here we need the triclinic conversion!
    this->runADI(num_particles, positions, charges, fields, potentials);
#else
    this->runIK(num_particles, positions, charges, fields, potentials);
#endif

    /* restore require_timings */
    require_timings = require_timings_before;

    delete fields;
    delete potentials;

    return timings;
}

void P3M::FarSolver::gatherTimings() {
    // gather timings
    if (comm.onMaster())
        MPI_Reduce(MPI_IN_PLACE, timings,
                NUM_TIMINGS, MPI_DOUBLE, MPI_MAX,
                0, comm.mpicomm);
    else
        MPI_Reduce(timings, 0,
                NUM_TIMINGS, MPI_DOUBLE, MPI_MAX,
                0, comm.mpicomm);

    for (int i = 1; i < NUM_TIMINGS; i++)
        timings[TOTAL] += timings[i];
}

/* Fetch the last timings */
const double* P3M::FarSolver::getTimings() {
    if (require_timings == NONE)
        throw std::logic_error("Wanting to get timings, but timings were not measured.");
    return timings;
}


