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
#include "utils.hpp"
#include "FCSCommon.h"
#include "common/near/near.h"

namespace P3M {

#define START(ID)                                               \
        if (require_timings!=NONE) timings[ID] += -MPI_Wtime();
#define STOP(ID)                                                \
        if (require_timings!=NONE) timings[ID] += MPI_Wtime();
#define STOPSTART(ID1, ID2)                     \
        if (require_timings!=NONE) {               \
            timings[ID1] += MPI_Wtime();             \
            timings[ID2] += -MPI_Wtime();            \
        }                                             \

Solver::Solver(MPI_Comm mpicomm) :
                    comm(mpicomm), fft(comm), errorEstimate(NULL) {
    P3M_DEBUG(printf( "P3M::Solver() started...\n"));

    errorEstimate = ErrorEstimate::create(comm);

    /* Init the P3M parameters */
    box_l[0] = 1.0;
    box_l[1] = 1.0;
    box_l[2] = 1.0;

    box_vectors[0][0] = box_vectors[1][1] = box_vectors[2][2] = 1.0;
    box_vectors[0][1] = box_vectors[0][2] = box_vectors[1][0] = box_vectors[2][0]
            = box_vectors[2][1] = box_vectors[1][2] = 0.0;
    isTriclinic = false;
    
    shiftGaussians=false;

    skin = 0.0;
    tolerance_field = P3M_DEFAULT_TOLERANCE_FIELD;
    n_interpol = P3M_DEFAULT_N_INTERPOL;

    /* Tunable parameters */
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
    total_energy = 0.0;
    near_field_flag = false;

#ifdef P3M_PRINT_TIMINGS
    require_timings = FULL;
#else
    require_timings = NONE;
#endif
    for (int i=0; i < NUM_TIMINGS; i++)
        timings[i] = 0.0;

    /* Init the derived params */
    grid_off[0] = P3M_DEFAULT_GRIDOFF;
    grid_off[1] = P3M_DEFAULT_GRIDOFF;
    grid_off[2] = P3M_DEFAULT_GRIDOFF;
    cao_cut[0] = 0.0;
    cao_cut[1] = 0.0;
    cao_cut[2] = 0.0;
    a[0] = 0.0;
    a[1] = 0.0;
    a[2] = 0.0;
    ai[0] = 0.0;
    ai[1] = 0.0;
    ai[2] = 0.0;
    additional_grid[0] = 0.0;
    additional_grid[1] = 0.0;
    additional_grid[2] = 0.0;

    /* init the P3M data */
    sum_qpart = 0;
    sum_q2 = 0.0;
    square_sum_q = 0.0;

    caf = NULL;
    cafx = NULL;
    cafy = NULL;
    cafz = NULL;
    caf_d = NULL;
    cafx_d = NULL;
    cafy_d = NULL;
    cafz_d = NULL;

    pos_shift = 0.0;

    d_op[0] = NULL;
    d_op[1] = NULL;
    d_op[2] = NULL;
    g_force = NULL;
    g_energy = NULL;

    ks_pnum = 0;

    send_grid = NULL;
    recv_grid = NULL;

    rs_grid = NULL;
    ks_grid = NULL;

    P3M_DEBUG(printf( "P3M::Solver() finished.\n"));
}

Solver::~Solver() {
    delete errorEstimate;
    fft.free_data(rs_grid);
    fft.free_data(ks_grid);
    sfree(send_grid);
    sfree(recv_grid);
    delete caf;
    delete cafx;
    delete cafy;
    delete cafz;
    delete caf_d;
    delete cafx_d;
    delete cafy_d;
    delete cafz_d;
    if (g_force != NULL) delete g_force;
    if (g_energy != NULL) delete g_energy;
    for (p3m_int i=0; i<3; i++)
        sfree(d_op[i]);
}

/***************************************************/
/* FORWARD DECLARATIONS OF INTERNAL FUNCTIONS */
/***************************************************/
#ifdef P3M_ENABLE_DEBUG
void printLocalGrid(local_grid_t l);
void printSendGrid(send_grid_t sm);
#endif

/** Prepare the data structures and constants of the P3M algorithm.
     All parameters have to be set. */
void Solver::prepare() {
    P3M_DEBUG(printf("  P3M::Solver::prepare() started... \n"));
    comm.prepare(box_l);
    /* Initializes the (inverse) grid constant ai and the cutoff for charge
     * assignment cao_cut. */
    for (p3m_int i=0; i<3; i++) {
        ai[i]      = grid[i]/(isTriclinic?1.0:box_l[i]);
        a[i]       = 1.0/ai[i];
        cao_cut[i] = 0.5*a[i]*cao;
    }
    this->prepareLocalCAGrid();
    P3M_DEBUG(this->printLocalGrid());
    this->prepareSendGrid();
    P3M_DEBUG(this->printSendGrid());
    send_grid = (p3m_float *) realloc(send_grid, sizeof(p3m_float)*sm.max);
    recv_grid = (p3m_float *) realloc(recv_grid, sizeof(p3m_float)*sm.max);

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
#endif

    /* position offset for calc. of first gridpoint */
    pos_shift = (p3m_float)((cao-1)/2) - (cao%2)/2.0;
    P3M_DEBUG(printf("    pos_shift=" FFLOAT "\n",pos_shift));

    /* FFT */
    P3M_INFO(printf("    Preparing FFTs...\n"));
    fft.prepare(local_grid.dim, local_grid.margin, grid, grid_off, &ks_pnum);
    rs_grid = fft.malloc_data();
    ks_grid = fft.malloc_data();

    /* k-space part */
#if defined(P3M_IK) 
    /* Calculates the Fourier transformed differential operator.
     *  Remark: This is done on the level of n-vectors and not k-vectors,
     *           i.e. the prefactor i*2*PI/L is missing! */
    for (int i = 0; i < 3; i++) {
        d_op[i] = static_cast<p3m_int *>(realloc(d_op[i], grid[i]*sizeof(p3m_int)));
        d_op[i][0] = 0;
        d_op[i][grid[i]/2] = 0;

        for (p3m_int j = 1; j < grid[i]/2; j++) {
            d_op[i][j] = j;
            d_op[i][grid[i] - j] = -j;
        }
    }
#endif    

    P3M_INFO(printf("    Calculating influence function...\n"));
#if !defined(P3M_INTERLACE) && defined(P3M_IK)
    computeInfluenceFunctionIK();
#elif defined(P3M_INTERLACE) && defined(P3M_IK)
    computeInfluenceFunctionIKI();
#else
    computeInfluenceFunctionADI();
#endif

    P3M_DEBUG(printf("  P3m::Solver::prepare() finished.\n"));
}


/** Calculates properties of the local FFT grid for the
         charge assignment process. */
void Solver::prepareLocalCAGrid() {
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
    local_grid.q_21_off = local_grid.dim[2] * (local_grid.dim[1] - cao);
}

/** Calculates the properties of the send/recv sub-grides of the local
 *  FFT grid.  In order to calculate the recv sub-grides there is a
 *  communication of the margins between neighbouring nodes. */
void Solver::prepareSendGrid() {
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

/** Debug function printing p3m structures */
void Solver::printLocalGrid() {
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
void Solver::printSendGrid() {
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
void Solver::computeInfluenceFunctionIK() {
    P3M_DEBUG(printf("  Solver::computeInfluenceFunctionIK() started...\n"));

    int size = 1;
    int end[3];
    int *new_grid = fft.plan[3].new_grid;;
    int *start = fft.plan[3].start;
    for (p3m_int i=0;i<3;i++) {
        size *= new_grid[i];
        end[i] = fft.plan[3].start[i] + new_grid[i];
    }

    if (g_force != NULL) delete g_force;
    if (g_energy != NULL) delete g_energy;
    g_force = new p3m_float[size];
    g_energy = new p3m_float[size];

    p3m_int *gridshift_x = this->computeGridShift(RX, grid[0]);
    p3m_int *gridshift_y = this->computeGridShift(RY, grid[1]);
    p3m_int *gridshift_z = this->computeGridShift(RZ, grid[2]);

    int n[3];
    for (n[0]=start[0]; n[0]<end[0]; n[0]++) {
        for (n[1]=start[1]; n[1]<end[1]; n[1]++) {
            for (n[2]=start[2]; n[2]<end[2]; n[2]++) {
                const p3m_int ind = (n[2]-start[2]) + new_grid[2] *
                        ((n[1]-start[1]) + new_grid[1]*(n[0]-start[0]));
                if ((n[KX]%(grid[RX]/2)==0) &&
                        (n[KY]%(grid[RY]/2)==0) &&
                        (n[KZ]%(grid[RZ]/2)==0) ) {
                    g_force[ind] = 0.0;
                    g_energy[ind] = 0.0;
                } else {
                    p3m_float numerator_force[3];
                    p3m_float numerator_energy;
                    p3m_float sumU2;
                    this->performAliasingSumsIK(
                            gridshift_x[n[KX]],
                            gridshift_y[n[KY]],
                            gridshift_z[n[KZ]],
                            numerator_force,
                            numerator_energy,
                            sumU2);

                    p3m_float fak1; //  k scalar numerator force (besides prefactors)
                    p3m_float fak2; //  k squared (besides prefactors)
                    if (!isTriclinic) {
                        fak1 =  d_op[RX][n[KX]] * numerator_force[RX] / box_l[RX] +
                                d_op[RY][n[KY]] * numerator_force[RY] / box_l[RY] +
                                d_op[RZ][n[KZ]] * numerator_force[RZ] / box_l[RZ];
                        fak2 =  SQR(d_op[RX][n[KX]] / box_l[RX]) +
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
    P3M_DEBUG(printf("  Solver::computeInfluenceFunctionIK() finished.\n"));
}

void Solver::performAliasingSumsIK(
        p3m_int nmx0, p3m_int nmy0, p3m_int nmz0,
        p3m_float numerator_force[3],
        p3m_float &numerator_energy,
        p3m_float &sumU2) {
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
                const p3m_float U2  = sy*pow(sinc(nmz/(p3m_float)grid[RZ]),2.0*cao);

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
                sumU2  += U2; // denominator = SUM U^2
            }
        }
    }
}

void Solver::computeInfluenceFunctionADI() {
   P3M_DEBUG(printf("  Solver::computeInfluenceFunctionADI() started...\n"));
   int size = 1;
   int end[3];
   int *new_grid = fft.plan[3].new_grid;;
   int *start = fft.plan[3].start;
   for (p3m_int i=0;i<3;i++) {
     size *= new_grid[i];
     end[i] = fft.plan[3].start[i] + new_grid[i];
   }

   if (g_force != NULL) delete g_force;
   if (g_energy != NULL) delete g_energy;
   g_force = new p3m_float[size];
   g_energy = new p3m_float[size];

   p3m_int *gridshift_x = this->computeGridShift(RX, grid[0]);
   p3m_int *gridshift_y = this->computeGridShift(RY, grid[1]);
   p3m_int *gridshift_z = this->computeGridShift(RZ, grid[2]);

   int n[3];
   for (n[0]=start[0]; n[0]<end[0]; n[0]++) {
     for (n[1]=start[1]; n[1]<end[1]; n[1]++) {
       for (n[2]=start[2]; n[2]<end[2]; n[2]++) {
         p3m_int ind = (n[2]-start[2]) + new_grid[2] *
           ((n[1]-start[1]) + new_grid[1]*(n[0]-start[0]));
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
   P3M_DEBUG(printf("  Solver::computeInfluenceFunctionADI() finished.\n"));
 }

void Solver::performAliasingSumsADI(
        p3m_int nmx0, p3m_int nmy0, p3m_int nmz0,
        p3m_float &numerator_force,
        p3m_float &numerator_energy, p3m_float denominator[4]) {
   denominator[0] = denominator[1] = denominator[2] = denominator[3] = 0.0;
   numerator_energy = 0.0;
   numerator_force = 0.0;

   p3m_float prefactor = SQR(M_PI/alpha);

   for (p3m_int mx = -P3M_BRILLOUIN; mx <= P3M_BRILLOUIN; mx++) {
     const p3m_int nmx = nmx0 + grid[RX]*mx;
     const p3m_float sx  = pow(sinc(nmx/(p3m_float)grid[RX]), 2.0*cao);
     for (p3m_int my = -P3M_BRILLOUIN; my <= P3M_BRILLOUIN; my++) {
       const p3m_int nmy = nmy0 + grid[RY]*my;
       const p3m_float sy  = sx*pow(sinc(nmy/(p3m_float)grid[RY]), 2.0*cao);
       for (p3m_int mz = -P3M_BRILLOUIN; mz <= P3M_BRILLOUIN; mz++) {
         const p3m_int nmz = nmz0 + grid[RZ]*mz;
         const p3m_float U2  = sy*pow(sinc(nmz/(p3m_float)grid[RZ]), 2.0*cao);
         
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
         numerator_force  += prefactor2;

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

void Solver::computeInfluenceFunctionIKI() {
  P3M_DEBUG(printf("  Solver::computeInfluenceFunctionIKI() started...\n"));
  int size = 1;
  int end[3];
  int *new_grid = fft.plan[3].new_grid;;
  int *start = fft.plan[3].start;
  for (p3m_int i=0;i<3;i++) {
    size *= new_grid[i];
    end[i] = fft.plan[3].start[i] + new_grid[i];
  }

  if (g_force != NULL) delete g_force;
  if (g_energy != NULL) delete g_energy;
  g_force = new p3m_float[size];
  g_energy = new p3m_float[size];

  p3m_int *gridshift_x = this->computeGridShift(RX, grid[0]);
  p3m_int *gridshift_y = this->computeGridShift(RY, grid[1]);
  p3m_int *gridshift_z = this->computeGridShift(RZ, grid[2]);

  int n[3];
  for (n[0]=start[0]; n[0]<end[0]; n[0]++) {
    for (n[1]=start[1]; n[1]<end[1]; n[1]++) {
      for (n[2]=start[2]; n[2]<end[2]; n[2]++) {
        p3m_int ind = (n[2]-start[2]) + new_grid[2] *
          ((n[1]-start[1]) + new_grid[1]*(n[0]-start[0]));
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
  P3M_DEBUG(printf("  Solver::computeInfluenceFunctionIKI() finished.\n"));
}

void Solver::performAliasingSumsIKI(
        p3m_int nmx0, p3m_int nmy0, p3m_int nmz0,
        p3m_float numerator_force[3],
        p3m_float &numerator_energy,
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


/** shifts the grid points by grid/2 */
p3m_int *Solver::computeGridShift(int dir, int size) {
    p3m_int *gridshift = new p3m_int[size];
    gridshift[0] = 0;

    for (int i=1; i<=grid[dir]/2; i++) {
        gridshift[i] = i;
        gridshift[grid[dir] - i] = -i;
    }

    return gridshift;
}


/* callback function for near field computations */
inline void
compute_near(const void *param, p3m_float dist, p3m_float *field, p3m_float *potential)
{
const near_params_t params = *(static_cast<const near_params_t*>(param));
//    const p3m_float alpha = *(static_cast<const p3m_float*>(param));
    const p3m_float adist = params.alpha * dist;
    
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

    *potential = erfc_part_ri-params.potentialOffset;
    //printf("_______near loop offset: %e \n ",params.potentialOffset);
   // printf("_______near loop alpha: %e \n ",params.alpha);
    *field = -(erfc_part_ri + 2.0*params.alpha*0.56418958354775627928034964498*exp(-adist*adist))
              / dist;

#else

    p3m_float erfc_part_ri = (1.0 - erf(adist)) / dist; /* use erf instead of erfc to fix ICC performance problems */
    *potential = erfc_part_ri-params.potentialOffset;
     //   printf("erfc part ri: %f \n",erfc_part_ri);
#ifdef ADDITIONAL_CHECKS
        if(*potential<0.0) printf("ERROR: The calculated real space potential contribution is smaller than the potential offset.\n");
#endif
        
    *field = -(erfc_part_ri + 2.0*params.alpha*0.56418958354775627928034964498*exp(-adist*adist))
              / dist;
#endif

}

/* callback function for performing a whole loop of near field computations (using compute_near) */
FCS_NEAR_LOOP_FP(compute_near_loop, compute_near);

void Solver::run(
        p3m_int _num_particles, p3m_float *_positions, p3m_float *_charges,
        p3m_float *_fields, p3m_float *_potentials) {
    P3M_INFO(printf( "P3M::Solver::run() started...\n"));
    P3M_DEBUG(printf("  type of timing: %d\n", require_timings));
    /* reset all timers */
    if (require_timings!=NONE) {
        for (int i = 0; i < NUM_TIMINGS; i++)
            timings[i] = 0.0;
    }

        if (isTriclinic) {
            box_l[0] = box_l[1] = box_l[2] = 1.0;
            this->prepare();
        }
    
    P3M_INFO(printf("    "
            "system parameters: box_l=" F3FLOAT "\n", \
            box_l[0], box_l[1], box_l[2]));
    P3M_INFO(printf(                                                      \
            "    p3m params: "                                    \
            "r_cut=" FFLOAT ", grid=" F3INT ", cao=" FINT ", "    \
            "alpha=" FFLOAT ", grid_off=" F3FLOAT "\n",           \
            r_cut, grid[0], grid[1], grid[2], cao, \
            alpha, grid_off[0], grid_off[1], grid_off[2]));

    P3M_DEBUG_LOCAL(MPI_Barrier(comm.mpicomm));
    P3M_DEBUG_LOCAL(printf("    %d: num_particles=%d\n",    \
            comm.rank, _num_particles));


    /* decompose system */
    p3m_int num_real_particles;
    p3m_int num_ghost_particles;
    p3m_float *positions, *ghost_positions;
    p3m_float *charges, *ghost_charges;
    fcs_gridsort_index_t *indices, *ghost_indices;
    fcs_gridsort_t gridsort;

    START(TIMING_DECOMP);
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

    STOP(TIMING_DECOMP);

#if defined(P3M_INTERLACE) && defined(P3M_AD)
    if(!isTriclinic){ //orthorhombic
    this->computeFarADI(num_real_particles, positions, charges, fields, potentials);
    } else {//triclinic
      p3m_float *positions_triclinic= new p3m_float[num_real_particles*3];
      this->calculateTriclinicPositions(positions, positions_triclinic,num_real_particles);
      this->computeFarADI(num_real_particles, positions_triclinic, charges, fields, potentials);
    }
    
#else
    if(!isTriclinic){ //orthorhombic
    this->computeFarIK(num_real_particles, positions, charges, fields, potentials);
    } else {//triclinic
      p3m_float *positions_triclinic= new p3m_float[num_real_particles*3];
      this->calculateTriclinicPositions(positions, positions_triclinic,num_real_particles);
      this->computeFarIK(num_real_particles, positions_triclinic, charges, fields, potentials);
    }
    
#endif

    if(isTriclinic){
        for(p3m_int i=0;i<3;i++) box_l[i] = box_vectors[i][i];
    }
    
    if (near_field_flag) {
        /* start near timer */
        START(TIMING_NEAR);

        /* compute near field */
        fcs_near_t near;

        fcs_near_create(&near);
        /*  fcs_near_set_field_potential(&near, compute_near);*/
        fcs_near_set_loop(&near, compute_near_loop);

        p3m_float box_base[3] = {0.0, 0.0, 0.0 };
        fcs_near_set_system(&near, box_base, box_vectors[0], box_vectors[1], box_vectors[2], NULL);
        
        fcs_near_set_particles(&near, num_real_particles, num_real_particles,
                positions, charges, indices,
                (_fields != NULL) ? fields : NULL,
                        (_potentials != NULL) ? potentials : NULL);

        fcs_near_set_ghosts(&near, num_ghost_particles,
                ghost_positions, ghost_charges, ghost_indices);
        
        near_params_t params;
      //  std::cout<<"in Solver run: flag "<<shiftGaussians<<std::endl;
        params.alpha=alpha; params.potentialOffset=(shiftGaussians?(erfc(alpha*r_cut))/r_cut:0.0);
      //  std::cout<<"_____erfc(ar) in Solver run: "<<erfc(alpha*r_cut)<<std::endl;
      //  std::cout<<"_____offset in Solver run: "<<params.potentialOffset<<std::endl;
      //   printf("\033[31;1m solver run: off %e \033[0m \n", params.potentialOffset);
        P3M_DEBUG(printf( "  calling fcs_near_compute()...\n"));
        fcs_near_compute(&near, r_cut, &params, comm.mpicomm);
        P3M_DEBUG(printf( "  returning from fcs_near_compute().\n"));

        fcs_near_destroy(&near);

        STOP(TIMING_NEAR);
    }

    START(TIMING_COMP);
    /* sort particles back */
    P3M_DEBUG(printf( "  calling fcs_gridsort_sort_backward()...\n"));
    fcs_gridsort_sort_backward(&gridsort,
            fields, potentials,
            _fields, _potentials, 1,
            comm.mpicomm);
    P3M_DEBUG(printf( "  returning from fcs_gridsort_sort_backward().\n"));

    fcs_gridsort_free(&gridsort);
    fcs_gridsort_destroy(&gridsort);

    sdelete(fields);
    sdelete(potentials);

    STOP(TIMING_COMP);

    /* collect timings from the different nodes */
    if (require_timings != NONE) this->collectPrintTimings();

    P3M_INFO(printf( "P3M::Solver::run() finished.\n"));
}

void Solver::computeFarADI(
        p3m_int num_charges, p3m_float* positions, p3m_float* charges,
        p3m_float* fields, p3m_float* potentials) {
    P3M_INFO(printf( "P3M::Solver::compute_far_adi() started...\n"));

    START(TIMING_CA);

    /* charge assignment */
    this->assignCharges(fft.data_buf, num_charges, positions, charges, 0);
    STOPSTART(TIMING_CA, TIMING_GATHER);
    /* gather the ca grid */
    this->gatherGrid(fft.data_buf);

    // Complexify
    for (p3m_int i = local_grid.size-1; i >= 0; i--)
        rs_grid[2*i] = fft.data_buf[i];

    STOPSTART(TIMING_GATHER, TIMING_CA);

    // Second (shifted) run
    /* charge assignment */
    this->assignCharges(fft.data_buf, num_charges, positions, charges, 1);

    STOPSTART(TIMING_CA, TIMING_GATHER);

    /* gather the ca grid */
    this->gatherGrid(fft.data_buf);
    /* now rs_grid should contain the local ca grid */

    // Complexify
    for (p3m_int i = local_grid.size-1; i >= 0; i--)
        rs_grid[2*i+1] = fft.data_buf[i];
    STOP(TIMING_GATHER);

    /* forward transform */
    P3M_DEBUG(printf("  calling fft_perform_forw()...\n"));
    START(TIMING_FORWARD);
    fft.forward(rs_grid);
    STOP(TIMING_FORWARD);
    P3M_DEBUG(printf("  returned from fft_perform_forw().\n"));

    /********************************************/
    /* POTENTIAL COMPUTATION */
    /********************************************/
    if (require_total_energy || potentials != NULL) {
        /* apply energy optimized influence function */
        START(TIMING_INFLUENCE);
        this->applyInfluenceFunction(rs_grid, ks_grid, g_energy);
        /* result is in ks_grid */
        STOP(TIMING_INFLUENCE);

        /* compute total energy, but not potentials */
        if (require_total_energy && potentials == NULL) {
            START(TIMING_POTENTIALS);
            total_energy = this->computeTotalEnergy();
            STOP(TIMING_POTENTIALS);
        }

        if (potentials != NULL) {
            /* backtransform the grid */

            if( require_timings != ESTIMATE_ALL
                    && require_timings != ESTIMATE_FFT ){
                P3M_DEBUG(printf( "  calling fft_perform_back (potentials)...\n"));
                START(TIMING_BACK);
                fft.backward(ks_grid);
                STOP(TIMING_BACK);
                P3M_DEBUG(printf( "  returned from fft_perform_back.\n"));
            }
            /** First (unshifted) run */
            START(TIMING_SPREAD)
            for (p3m_int i=0; i<local_grid.size; i++) {
                fft.data_buf[i] = ks_grid[2*i];
            }

            this->spreadGrid(fft.data_buf);

            STOPSTART(TIMING_SPREAD, TIMING_POTENTIALS);

            this->assignPotentials(fft.data_buf,
                    num_charges, positions,
                    charges, 0, potentials);

            STOPSTART(TIMING_POTENTIALS, TIMING_SPREAD);

            /** Second (shifted) run */
            for (p3m_int i=0; i<local_grid.size; i++) {
                fft.data_buf[i] = ks_grid[2*i+1];
            }
            this->spreadGrid(fft.data_buf);

            STOPSTART(TIMING_SPREAD, TIMING_POTENTIALS);

            this->assignPotentials(fft.data_buf,
                    num_charges, positions,
                    charges, 1, potentials);

            STOP(TIMING_POTENTIALS);
        }
    }

    /********************************************/
    /* FIELD COMPUTATION */
    /********************************************/
    if (fields != NULL) {
        /* apply force optimized influence function */
        START(TIMING_INFLUENCE);
        this->applyInfluenceFunction(rs_grid, ks_grid, g_force);
        STOP(TIMING_INFLUENCE);

        /* backtransform the grid */
        if(require_timings != ESTIMATE_ALL
                && require_timings != ESTIMATE_FFT){
            P3M_DEBUG(printf("  calling fft_perform_back...\n"));
            START(TIMING_BACK);
            fft.backward(ks_grid);
            STOP(TIMING_BACK);
            P3M_DEBUG(printf("  returned from fft_perform_back.\n"));
        }

        START(TIMING_SPREAD);
        /* First (unshifted) run */
        P3M_INFO(printf("  computing unshifted grid\n"));
        for (p3m_int i=0; i<local_grid.size; i++) {
            fft.data_buf[i] = ks_grid[2*i];
        }

        this->spreadGrid(fft.data_buf);

        STOPSTART(TIMING_SPREAD, TIMING_FIELDS);

        this->assignFieldsAD(fft.data_buf, num_charges, positions, 0, fields);

        STOPSTART(TIMING_FIELDS, TIMING_SPREAD);

        /* Second (shifted) run */
        P3M_INFO(printf("  computing shifted grid\n"));
        for (p3m_int i=0; i<local_grid.size; i++) {
            fft.data_buf[i] = ks_grid[2*i+1];
        }

        this->spreadGrid(fft.data_buf);

        STOPSTART(TIMING_SPREAD, TIMING_FIELDS);

        this->assignFieldsAD(fft.data_buf, num_charges, positions, 1, fields);

        if(isTriclinic) this->cartesianizeFields(fields, num_charges);
                
        STOP(TIMING_FIELDS);
    }


    /* estimate FFT back timing*/
    if(require_timings == ESTIMATE_ALL
            || require_timings == ESTIMATE_FFT) {
        if (potentials != NULL)
            timings[TIMING_BACK] = 2 * timings[TIMING_FORWARD];
        else
            timings[TIMING_BACK] = timings[TIMING_FORWARD];
    }

    P3M_INFO(printf( "P3M::Solver::compute_far_adi() finished.\n"));
}

void Solver::computeFarIK(
        p3m_int num_charges, p3m_float* positions, p3m_float* charges,
        p3m_float* fields, p3m_float* potentials) {
    P3M_INFO(printf( "P3M::Solver::compute_far_ik() started...\n"));

    START(TIMING_CA);

    /* charge assignment */
    this->assignCharges(rs_grid, num_charges, positions, charges, 0);
    STOPSTART(TIMING_CA, TIMING_GATHER);
    /* gather the ca grid */
    this->gatherGrid(rs_grid);
    /* now rs_grid should contain the local ca grid */
    STOP(TIMING_GATHER);

    /* forward transform */

    P3M_DEBUG(printf( "  calling fft.forward()...\n"));
    START(TIMING_FORWARD);
    fft.forward(rs_grid);
    STOP(TIMING_FORWARD);
    P3M_DEBUG(printf("  returned from fft.forward().\n"));

    /********************************************/
    /* POTENTIAL COMPUTATION */
    /********************************************/
    if (require_total_energy || potentials != NULL) {
        /* apply energy optimized influence function */
        START(TIMING_INFLUENCE);
        this->applyInfluenceFunction(rs_grid, ks_grid, g_energy);
        /* result is in ks_grid */
        STOP(TIMING_INFLUENCE);

        /* compute total energy, but not potentials */
        if (require_total_energy && potentials == NULL) {
            START(TIMING_POTENTIALS);
            total_energy = this->computeTotalEnergy();
            STOP(TIMING_POTENTIALS);
        }

        if (potentials != NULL) {
            /* backtransform the grid */

            if (require_timings != ESTIMATE_ALL
                    && require_timings != ESTIMATE_FFT){
                P3M_DEBUG(printf( "  calling fft.backward (potentials)...\n"));
                START(TIMING_BACK);
                fft.backward(ks_grid);
                STOP(TIMING_BACK);
                P3M_DEBUG(printf( "  returned from fft.backward.\n"));
            }
            /* redistribute energy grid */
            START(TIMING_SPREAD);
            this->spreadGrid(ks_grid);

            STOPSTART(TIMING_SPREAD, TIMING_POTENTIALS);

            /* compute potentials */
            this->assignPotentials(ks_grid,
                    num_charges, positions, charges, 0, potentials);

            STOP(TIMING_POTENTIALS);
        }
    }

    /********************************************/
    /* FIELD COMPUTATION */
    /********************************************/
    if (fields != NULL) {
        /* apply force optimized influence function */
        START(TIMING_INFLUENCE);
        this->applyInfluenceFunction(rs_grid, ks_grid, g_force);
        STOP(TIMING_INFLUENCE);

        /* result is in ks_grid */
        for (int dim = 0; dim < 3; dim++) {
            /* differentiate in direction dim */
            /* result is stored in rs_grid */
            START(TIMING_FIELDS);
            this->differentiateIK(dim, ks_grid, rs_grid);
            STOP(TIMING_FIELDS);

            if (require_timings != ESTIMATE_ALL
                    && require_timings != ESTIMATE_FFT) {

                /* backtransform the grid */
                P3M_DEBUG(printf( "  calling fft.backward (field dim=%d)...\n", dim));
                START(TIMING_BACK);
                fft.backward(rs_grid);
                STOP(TIMING_BACK);
                P3M_DEBUG(printf("  returned from fft.backward.\n"));
            }

            /* redistribute force grid */
            START(TIMING_SPREAD);
            this->spreadGrid(rs_grid);
            STOPSTART(TIMING_SPREAD, TIMING_FIELDS);
            this->assignFieldsIK(rs_grid, dim, num_charges, positions, 0, fields);
            STOP(TIMING_FIELDS);
        }
    }

    /* estimate FFT back timing*/
    if (require_timings == ESTIMATE_ALL || require_timings == ESTIMATE_FFT) {
        if (potentials != NULL)
            timings[TIMING_BACK] = 4 * timings[TIMING_FORWARD];
        else
            timings[TIMING_BACK] = 3 * timings[TIMING_FORWARD];
    }

    P3M_INFO(printf( "P3M::Solver::compute_far_ik() finished.\n"));
}


/***************************************************/
/* RUN COMPONENTS */

    void Solver::calculateTriclinicPositions(p3m_float *positions,
            p3m_float *triclinic_positions, p3m_int number) {

        for (p3m_int i = 0; i < number; i++) {
            triclinic_positions[3 * i] = 1 / box_vectors[0][0] * positions[3 * i]
                    - box_vectors[1][0] / (box_vectors[0][0] * box_vectors[1][1])
                    * positions[3 * i + 1] + (box_vectors[1][0]
                    * box_vectors[2][1] - box_vectors[1][1] * box_vectors[2][0])
                    / (box_vectors[0][0] * box_vectors[1][1] * box_vectors[2][2])
                    * positions[3 * i + 2];
            triclinic_positions[3 * i + 1] = 1 / box_vectors[1][1]
                    * positions[3 * i + 1] - box_vectors[2][1]
                    / (box_vectors[1][1] * box_vectors[2][2])
                    * positions[3 * i + 2];
            triclinic_positions[3 * i + 2] = 1 / box_vectors[2][2]
                    * positions[3 * i + 2];
            
            P3M_DEBUG_LOCAL(if(i<5)printf("triclinic position %d: %f %f %f \n", i,
                    triclinic_positions[3*i],triclinic_positions[3*i+1],
                    triclinic_positions[3*i+2]));
        }
    }

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

    p3m_int num_particles;

    fcs_gridsort_create(gridsort);

            fcs_gridsort_set_system(gridsort, box_base,
                    box_vectors[0], box_vectors[1], box_vectors[2], NULL);         
    
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

/* CHARGE ASSIGNMENT */
/** Compute the data of the charge assignment grid points.

      The function returns the linear index of the top left grid point
      in the charge assignment grid that corresponds to real_pos. When
      "shifted" is set, it uses the shifted position for interlacing.
      After the call, caf_cache contains a cache of the values of the
      charge assignment fraction (caf) for x,y,z.
 */
p3m_int
Solver::getCAPoints(p3m_float real_pos[3], p3m_int shifted) {
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
            printf("pos=%lf, pos_shift=%lf, real_pos=%lf, ld_pos=%lf, ai=%lf, shifted=%d",
                    pos, pos_shift, real_pos[dim], local_grid.ld_pos[dim],  ai[dim],
                    shifted);
        } else if (grid_ind > local_grid.dim[dim]) {
            printf("grid_ind[%d]=%d > %d\n", dim, grid_ind, local_grid.dim[dim]);
            printf("pos=%lf, pos_shift=%lf, real_pos=%lf, ld_pos=%lf, ai=%lf, shifted=%d",
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
void Solver::assignCharges(p3m_float *data,
        p3m_int num_real_particles, p3m_float *positions, p3m_float *charges,
        p3m_int shifted) {
    P3M_DEBUG(printf( "  P3M::Solver::assign_charges() started...\n"));

    const p3m_int q2off = local_grid.q_2_off;
    const p3m_int q21off = local_grid.q_21_off;

    /* init local charge grid */
    for (p3m_int i=0; i<local_grid.size; i++) data[i] = 0.0;

    /* now assign the charges */
    for (p3m_int pid=0; pid < num_real_particles; pid++) {
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

    P3M_DEBUG(printf( "  P3M::Solver::assign_charges() finished...\n"));
}

/* Gather information for FFT grid inside the nodes domain (inner local grid) */
void Solver::gatherGrid(p3m_float* rs_grid) {
    P3M_DEBUG(printf( "  P3M::Solver::gather_grid() started...\n"));

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
//            Parallel3DFFT::unpack_block(recv_grid, rs_grid, sm.r_ld[r_dir],
//                    sm.r_dim[r_dir], local_grid.dim, 1);
        }
    }

    P3M_DEBUG(printf( "  P3M::Solver::gather_grid() finished.\n"));
}

void Solver::spreadGrid(p3m_float* rs_grid) {
    P3M_DEBUG(printf( "  P3M::Solver::spread_grid() started...\n"));

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

    P3M_DEBUG(printf( "  P3M::Solver::spread_grid() finished.\n"));
}

/* apply the influence function */
void Solver::applyInfluenceFunction(p3m_float *in, p3m_float *out, p3m_float *g) {
    const p3m_int size = fft.plan[3].new_size;
    for (p3m_int i=0; i < size; i++) {
        out[2*i] = g[i] * in[2*i];
        out[2*i+1] = g[i] * in[2*i+1];
    }
}

/* Add up the measured timings and collect information from all nodes.
 * Print the timings if requested so.
 */
void Solver::collectPrintTimings() {
    timings[TIMING_FAR] += timings[TIMING_CA];
    timings[TIMING_FAR] += timings[TIMING_GATHER];
    timings[TIMING_FAR] += timings[TIMING_FORWARD];
    timings[TIMING_FAR] += timings[TIMING_BACK];
    timings[TIMING_FAR] += timings[TIMING_INFLUENCE];
    timings[TIMING_FAR] += timings[TIMING_SPREAD];
    timings[TIMING_FAR] += timings[TIMING_POTENTIALS];
    timings[TIMING_FAR] += timings[TIMING_FIELDS];

    timings[TIMING] += timings[TIMING_DECOMP];
    timings[TIMING] += timings[TIMING_FAR];
    timings[TIMING] += timings[TIMING_NEAR];
    timings[TIMING] += timings[TIMING_COMP];

    if (comm.onMaster())
        MPI_Reduce(MPI_IN_PLACE, timings,
                NUM_TIMINGS, MPI_DOUBLE, MPI_MAX,
                0, comm.mpicomm);
    else
        MPI_Reduce(timings, 0,
                NUM_TIMINGS, MPI_DOUBLE, MPI_MAX,
                0, comm.mpicomm);

#ifdef P3M_PRINT_TIMINGS
    if (comm.onMaster()) {
        printf("  P3M TIMINGS:\n");
        printf("    total=%le (%lf)\n", timings[TIMING], 1.0);
        printf("      far=%le (%lf)\n", timings[TIMING_FAR],
                timings[TIMING_FAR]/timings[TIMING]);
        printf("     near=%le (%lf)\n", timings[TIMING_NEAR],
                timings[TIMING_NEAR]/timings[TIMING]);
        printf("       ca=%le (%lf)\n", timings[TIMING_CA],
                timings[TIMING_CA]/timings[TIMING]);
        printf("      pot=%le (%lf)", timings[TIMING_POTENTIALS],
                timings[TIMING_POTENTIALS]/timings[TIMING]);
        //if(require_timings == ESTIMATE_ALL //not yet implemented
        //|| require_timings == ESTIMATE_ASSIGNMENT)
        //    printf(" (estimate)");
        printf("\n");
        printf("   fields=%le (%lf)", timings[TIMING_FIELDS],
                timings[TIMING_FIELDS]/timings[TIMING]);
        //if(require_timings == ESTIMATE_ALL //not yet implemented
        //|| require_timings == ESTIMATE_ASSIGNMENT)
        //    printf(" (estimate)");
        printf("\n");
        printf("   gather=%le (%lf)\n", timings[TIMING_GATHER],
                timings[TIMING_GATHER]/timings[TIMING]);
        printf("   spread=%le (%lf)\n", timings[TIMING_SPREAD],
                timings[TIMING_SPREAD]/timings[TIMING]);
        printf("  forward=%le (%lf)\n", timings[TIMING_FORWARD],
                timings[TIMING_FORWARD]/timings[TIMING]);
        printf("     back=%le (%lf)", timings[TIMING_BACK],
                timings[TIMING_BACK]/timings[TIMING]);
        if(require_timings == ESTIMATE_ALL
                || require_timings == ESTIMATE_FFT)
            printf(" (estimate)");
        printf("\n");
        printf("   decomp=%le (%lf)\n", timings[TIMING_DECOMP],
                timings[TIMING_DECOMP]/timings[TIMING]);
        printf("     comp=%le (%lf)\n", timings[TIMING_COMP],
                timings[TIMING_COMP]/timings[TIMING]);
    }
#endif

}

void Solver::differentiateIK(int dim, p3m_float* in, p3m_float* out) {
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
    for (j[0]=0; j[0]<fft.plan[3].new_grid[0]; j[0]++) {
        for (j[1]=0; j[1]<fft.plan[3].new_grid[1]; j[1]++) {
            for (j[2]=0; j[2]<fft.plan[3].new_grid[2]; j[2]++) {
                /* i*k*(Re+i*Im) = - Im*k + i*Re*k     (i=sqrt(-1)) */
                if (!isTriclinic) { //orthorhombic case
                        out[ind] =
                        -2.0*M_PI*(in[ind+1] * d_operator[ j[dim]+fft.plan[3].start[dim] ])
                                / box_l[dim_rs];
                        out[ind+1] =
                        2.0*M_PI*in[ind] * d_operator[ j[dim]+fft.plan[3].start[dim] ]
                                / box_l[dim_rs];
                    } else { //triclinic case //k correct
                        out [ind] = -2.0 * M_PI * in[ind + 1] * (d_operator[ j[dim]+fft.plan[3].start[dim] ]*(box_vectors[(dim+1)%3][(dim+1)%3]*box_vectors[(dim+2)%3][(dim+2)%3]-box_vectors[(dim+1)%3][(dim+2)%3]*box_vectors[(dim+2)%3][(dim+1)%3])+d_operator[ j[(dim+1)%3]+fft.plan[3].start[(dim+1)%3] ]*(box_vectors[(dim+2)%3][(dim+1)%3]*box_vectors[dim][(dim+2)%3]-box_vectors[(dim+2)%3][(dim+2)%3]*box_vectors[dim][(dim+1)%3])+d_operator[ j[(dim+2)%3]+fft.plan[3].start[(dim+2)%3] ]*(box_vectors[dim][(dim+1)%3]*box_vectors[(dim+1)%3][(dim+2)%3]-box_vectors[dim][(dim+2)%3]*box_vectors[(dim+1)%3][(dim+1)%3]))
                                /volume;
                     out[ind+1] =//k correct
                                     2.0*M_PI*in[ind] * (d_operator[ j[dim]+fft.plan[3].start[dim] ]*(box_vectors[(dim+1)%3][(dim+1)%3]*box_vectors[(dim+2)%3][(dim+2)%3]-box_vectors[(dim+1)%3][(dim+2)%3]*box_vectors[(dim+2)%3][(dim+1)%3])+d_operator[ j[(dim+1)%3]+fft.plan[3].start[(dim+1)%3] ]*(box_vectors[(dim+2)%3][(dim+1)%3]*box_vectors[dim][(dim+2)%3]-box_vectors[(dim+2)%3][(dim+2)%3]*box_vectors[dim][(dim+1)%3])+d_operator[ j[(dim+2)%3]+fft.plan[3].start[(dim+2)%3] ]*(box_vectors[dim][(dim+1)%3]*box_vectors[(dim+1)%3][(dim+2)%3]-box_vectors[dim][(dim+2)%3]*box_vectors[(dim+1)%3][(dim+1)%3]))
                                / volume;
                    }
                    ind += 2;
            }
        }
    }

    /* store the result in rs_grid */
}

/** Compute the total energy of the system in kspace. No need to
      backtransform the FFT grid in this case! */
p3m_float Solver::computeTotalEnergy() {
    p3m_float local_k_space_energy;
    p3m_float k_space_energy;

    P3M_DEBUG(printf( "  P3M::Solver::compute_total_energy() started...\n"));

    local_k_space_energy = 0.0;
    for (p3m_int i=0; i < fft.plan[3].new_size; i++)
        /* Use the energy optimized influence function */
        local_k_space_energy += g_energy[i] * ( SQR(rs_grid[2*i]) +
                SQR(rs_grid[2*i+1]) );

    MPI_Reduce(&local_k_space_energy, &k_space_energy, 1, P3M_MPI_FLOAT,
            MPI_SUM, 0, comm.mpicomm);
    p3m_float prefactor = 1.0 / (2.0 * volume);
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

    P3M_DEBUG(printf( "  P3M::Solver::compute_total_energy() finished.\n"));
    return k_space_energy;
}

/* Backinterpolate the potentials obtained from k-space to the positions */
void
Solver::assignPotentials(p3m_float *data,
        p3m_int num_real_particles, p3m_float* positions, p3m_float* charges,
        p3m_int shifted, p3m_float* potentials) {
    const p3m_int q2off = local_grid.q_2_off;
    const p3m_int q21off = local_grid.q_21_off;
    const p3m_float V_inverse = 1.0 / (box_l[0] * box_l[1] * box_l[2]);

    P3M_DEBUG(printf( "  P3M::Solver::assign_potentials() started...\n"));
    /* Loop over all particles */
    for (p3m_int pid=0; pid < num_real_particles; pid++) {
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

        potential *= V_inverse;
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
    P3M_DEBUG(printf( "  P3M::Solver::assign_potentials() finished.\n"));
}

/* Backinterpolate the forces obtained from k-space to the positions */
void
Solver::assignFieldsIK(p3m_float *data, p3m_int dim,
        p3m_int num_real_particles, p3m_float* positions,
        p3m_int shifted, p3m_float* fields) {
    const p3m_int q2off = local_grid.q_2_off;
    const p3m_int q21off = local_grid.q_21_off;
    const p3m_float prefactor = 1.0 / (2.0 * box_l[0] * box_l[1] * box_l[2]);
    const p3m_int dim_rs = (dim+ks_pnum) % 3;

    P3M_DEBUG(printf( "  P3M::Solver::assign_fields_ik() started...\n"));
    /* Loop over all particles */
    for (p3m_int pid=0; pid < num_real_particles; pid++) {
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
    P3M_DEBUG(printf( "  P3M::Solver::assign_fields_ik() finished.\n"));
}

/* Backinterpolate the forces obtained from k-space to the positions */
void
Solver::assignFieldsAD(p3m_float *data,
        p3m_int num_real_particles, p3m_float* positions,
        p3m_int shifted, p3m_float* fields) {
    const p3m_int q2off = local_grid.q_2_off;
    const p3m_int q21off = local_grid.q_21_off;
    const p3m_float prefactor = 1.0 / (box_l[0] * box_l[1] * box_l[2]);
    const p3m_float l_x_inv = 1.0/box_l[0];
    const p3m_float l_y_inv = 1.0/box_l[1];
    const p3m_float l_z_inv = 1.0/box_l[2];

    P3M_DEBUG(printf( "  P3M::Solver::assign_fields_ad() started...\n"));
    /* Loop over all particles */
    for (p3m_int pid = 0; pid < num_real_particles; pid++) {
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

void
Solver::cartesianizeFields(p3m_float *fields, p3m_int num_particles){
    p3m_int part_no;   
    
    for (part_no = 0; part_no < num_particles; part_no++) {
        fields[3 * part_no + 2] = 
            (box_vectors[1][0] * box_vectors[2][1] - box_vectors[1][1] * box_vectors[2][0]) 
            / (box_vectors[0][0] * box_vectors[1][1] * box_vectors[2][2]) 
            * fields[3 * part_no]
            -(box_vectors[2][1]) / (box_vectors[1][1] * box_vectors[2][2]) 
            * fields[3 * part_no + 1]
            + 1 / box_vectors[2][2] * fields[3 * part_no + 2];
   //todo: -boxbx oder boxby?     
        fields[3 * part_no + 1] =
            -(box_vectors[1][0]) / (box_vectors[0][0] * box_vectors[1][1])
            * fields[3 * part_no]
            + 1 / box_vectors[1][1] * fields[3 * part_no + 1];
        
        fields[3 * part_no] = 1 / (box_vectors[0][0]) * fields[3 * part_no];
    }
}


/***************************************************/
/******************** TUNING ***********************/
/***************************************************/
/** Good mesh sizes for fftw
 */
static const p3m_int good_gridsize[] = {
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
static const p3m_int step_good_gridsize[] =
{ 0, 15, 26, 44, 58, 72, 78, 83, 90, 94, 98, 101, 103, 104 } ;
static const p3m_int num_steps_good_gridsize = 14;


void
Solver::tune(p3m_int num_particles, p3m_float *positions, p3m_float *charges) {

    /* Count the charges */
    p3m_float sum_q2_before = sum_q2;
    this->countCharges(num_particles, charges);

    if (!comm.onMaster()) {
        int howoften;
        P3M_DEBUG_LOCAL(printf("  %d: How often to run tune_far?\n", comm.rank));
        MPI_Bcast(&howoften, 1, MPI_INT, Communication::MPI_MASTER, comm.mpicomm);
        P3M_DEBUG_LOCAL(printf("  %d: Running tune_far %d times.\n", comm.rank, howoften));
        for (; howoften > 0; howoften--)
            this->tuneFar(num_particles, positions, charges);

        needs_retune = false;
        return;
    }

    /* Retune if the number of charges has changed */
    if (!(float_is_equal(sum_q2_before, sum_q2))) {
        P3M_INFO(printf( "  Number of charges changed, retuning is needed.\n"));
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
        int howoften = 0;
        MPI_Bcast(&howoften, 1, MPI_INT, Communication::MPI_MASTER, comm.mpicomm);
        P3M_INFO(printf( "tune() finished.\n"));
        return;
    }

    P3M_INFO(printf( "  Retuning is required.\n"));

    if (!tune_r_cut) {
        P3M_INFO(printf( "    r_cut=" FFLOAT " (fixed)\n", r_cut));
        int howoften = 1;
        MPI_Bcast(&howoften, 1, MPI_INT, Communication::MPI_MASTER, comm.mpicomm);
        this->tuneFar(num_particles, positions, charges, r_cut);
    } else {
        int howoften = 2;
        MPI_Bcast(&howoften, 1, MPI_INT, Communication::MPI_MASTER, comm.mpicomm);

        /* compute the average distance between two charges  */
        p3m_float avg_dist =
                pow((box_l[0]*box_l[1]*box_l[2])
                        / sum_qpart, 0.33333);

        /* FIRST ESTIMATE */
        /* set the initial r_cut to 3 times the average distance between
         charges */
        p3m_float r_cut = 3.0 * avg_dist;

        /* tune r_cut to half the box length */
        if (0.5*box_l[1]-skin < r_cut)
            r_cut = 0.5*box_l[0] - skin;
        if (0.5*box_l[1]-skin < r_cut)
            r_cut = 0.5*box_l[1] - skin;
        if (0.5*box_l[2]-skin < r_cut)
            r_cut = 0.5*box_l[2] - skin;

        P3M_INFO(printf( "    r_cut=" FFLOAT " (first estimate)\n", r_cut));

        // @todo get near timing
        TuneParameters *p =
                this->tuneFar(num_particles, positions, charges, r_cut);

        /* SECOND ESTIMATE */
        /* use the fact that we know that timing_near scales like r_cut**3
         and timing_far like r_cut**(-3) to get the second estimate */

        p3m_float rel_timing_diff =
                fabs(p->timing_near - p->timing_far) /
                (p->timing_near + p->timing_far);
        P3M_INFO(printf( "    rel_timing_diff=" FFLOAT "\n", rel_timing_diff));

        p3m_float rcut3 = pow(r_cut, 3);
        p3m_float c_near = p->timing_near/rcut3;
        p3m_float c_far = p->timing_far*rcut3;
        p3m_float rcut_new = pow(c_far/c_near, 1./6.);

        r_cut = rcut_new;
        P3M_INFO(printf( "    r_cut=" FFLOAT " (second estimate)\n", r_cut));

        // @todo get near timing
        // second far tuning
        p = this->tuneFar(num_particles, positions, charges, r_cut);

        rel_timing_diff =
                fabs(p->timing_near - p->timing_far) /
                (p->timing_near + p->timing_far);
        P3M_INFO(printf("    rel_timing_diff=" FFLOAT "\n", rel_timing_diff));
        P3M_INFO(printf("    Finished tuning.\n"));
    }

    /* mark that the method was retuned */
    needs_retune = false;
}

/** Slave variant of tune_far. */
void
Solver::tuneFar(p3m_int num_particles, p3m_float *positions, p3m_float *charges) {
    if (comm.onMaster()) throw std::runtime_error("tuneFar without r_cut cannot be called on master node.");
    this->tuneFar(num_particles, positions, charges, 0.0);
}

Solver::TuneParameters*
Solver::tuneFar(p3m_int num_particles, p3m_float *positions, p3m_float *charges,
        p3m_float _r_cut) {
    /* r_cut is ignored on the slaves */
    /* Distinguish between two types of parameters:
       Input params:
     * box length
     * #charges
     * cutoff
     * skin (only if comm is cartesian)
     * error (default=1e-4)
     * component flags

       Performance related params:
     * caf interpolation

       Tuned params:
     * alpha
     * cao
     * grid
     */
    P3M_INFO(printf( "P3M::Solver::tuneFar() started...\n"));

    TuneParameters *final_params = NULL;

    if (!comm.onMaster()) {
        this->tuneBroadcastSlave(num_particles, positions, charges);
    } else {

        /* check whether the input parameters are sane */
        if (_r_cut < 0.0)
            throw std::domain_error("r_cut is negative!");

        if (float_is_zero(_r_cut))
            throw std::domain_error("r_cut is too small!");

        /* check whether cutoff is larger than half a box length */
//        if ((_r_cut > 0.5*box_l[0]) ||
//                (_r_cut > 0.5*box_l[1]) ||
//                (_r_cut > 0.5*box_l[2]))
//            throw std::domain_error("r_cut is larger than half a system box length.");

        /* check whether cutoff is larger than domain size */

        try {
            TuneParameters *p = new TuneParameters();
            p->p.r_cut = _r_cut;
            p->p.alpha = alpha;
            p->p.cao = cao;
            p->p.grid[0] = grid[0];
            p->p.grid[1] = grid[1];
            p->p.grid[2] = grid[2];
            TuneParametersList params_to_try;
            params_to_try.push_back(p);
            this->tuneAlpha(params_to_try);
            this->tuneCAO(params_to_try);
            this->tuneGrid(params_to_try);
            if (params_to_try.size() == 0)
                throw std::logic_error("No parameter set left to time.");
            final_params =
                    this->timeParams(num_particles, positions, charges,
                            params_to_try);

            // @todo: move to tune function
            /* Set and broadcast the final parameters. */
            this->r_cut = final_params->p.r_cut;
            alpha = final_params->p.alpha;
            grid[0] = final_params->p.grid[0];
            grid[1] = final_params->p.grid[1];
            grid[2] = final_params->p.grid[2];
            cao = final_params->p.cao;
            this->tuneBroadcastCommand(CMD_FINISHED, final_params->p);
        } catch (...) {
            P3M_INFO(printf( "  Tuning failed.\n"));
            this->tuneBroadcastCommand(CMD_FAILED);
            P3M_INFO(printf( "P3M::Solver::tune_far() finished.\n"));
            throw;
        }
    }

    P3M_INFO(printf( "  Tuning was successful.\n"));

    /* At the end of retuning, prepare the method */
    this->prepare();

    P3M_INFO(printf( "P3M::Solver::tuneFar() finished.\n"));

    return final_params;
}


/* Tune alpha */
void
Solver::tuneAlpha(TuneParametersList &params_to_try) {
    if (tune_alpha) {
        // compute the alpha for all params to try
        for (TuneParametersList::iterator p = params_to_try.begin();
                p != params_to_try.end(); ++p) {
            errorEstimate->compute_alpha(tolerance_field,
                    (*p)->p, sum_qpart, sum_q2, box_l);
            P3M_INFO(printf("    => alpha=" FFLOAT "\n", (*p)->p.alpha));
        }
    } else {
        P3M_INFO(printf("    alpha=" FFLOAT " (fixed)\n", params_to_try.front()->p.alpha));
    }
}

/* Tune cao */
void
Solver::tuneCAO(TuneParametersList &params_to_try) {
#ifdef P3M_AD
    const p3m_int min_cao = 2;
#else
    const p3m_int min_cao = 1;
#endif

    if (tune_cao) {
        P3M_INFO(printf("    Testing cao={ "));
        for (TuneParametersList::iterator p = params_to_try.begin();
                p != params_to_try.end(); ++p) {
            (*p)->p.cao = min_cao;
            for (p3m_int cao = CAF::MAX_CAO; cao > min_cao; cao--) {
                TuneParameters *tp_new = new TuneParameters();
                tp_new->p = (*p)->p;
                tp_new->p.cao = cao;
                params_to_try.insert(p, tp_new);
                P3M_INFO(printf(FINT " ", cao));
            }
        }
        P3M_INFO(printf(FINT " }\n", min_cao));
    } else {
        P3M_INFO(printf( "    cao=" FINT " (fixed)\n", params_to_try.front()->p.cao));
    }
}

/* params should have decreasing cao */
void
Solver::tuneGrid(TuneParametersList &params_to_try) {
    if (tune_grid) {
        p3m_int step_ix = 0;
        // store the minimal grid size seen so far
        p3m_int min_grid1d = good_gridsize[step_good_gridsize[num_steps_good_gridsize-1]];

        for (TuneParametersList::iterator p = params_to_try.begin();
                p != params_to_try.end(); ++p) {
            // Find smallest possible grid for this set of parameters
            // Test whether accuracy can be achieved with given parameters
            p3m_int upper_ix;
            // test this step
            P3M_INFO(printf("    Trying to find grid for r_cut=" FFLOAT ", " \
                    "alpha=" FFLOAT ", "                          \
                    "cao=" FINT "\n",                             \
                    (*p)->p.r_cut, (*p)->p.alpha, (*p)->p.cao));
            do {
                step_ix++;
                if (step_ix >= num_steps_good_gridsize) break;
                upper_ix = step_good_gridsize[step_ix];
                p3m_int grid1d = good_gridsize[upper_ix];
                (*p)->p.grid[0] = (*p)->p.grid[1] = (*p)->p.grid[2] = grid1d;
                P3M_DEBUG(printf("      rough grid=" FINT "\n", grid1d));
                this->tuneBroadcastCommand(CMD_COMPUTE_ERROR_ESTIMATE);
                errorEstimate->compute_master((*p)->p, sum_qpart, sum_q2, box_l, (*p)->error, (*p)->rs_error, (*p)->ks_error, box_vectors, isTriclinic);
                P3M_DEBUG(printf("        => error=" FFLOATE "\n", (*p)->error));
            } while ((*p)->error > tolerance_field);

            // reached largest possible grid, remove the rest of the parameter sets
            if (step_ix >= num_steps_good_gridsize) {
                P3M_INFO(printf("    Too large grid size, skipping rest of parameter sets.\n"));
                for (TuneParametersList::iterator pit = params_to_try.end();
                        pit != p; --pit) {
                    delete[] *pit;
                    params_to_try.erase(pit);
                }
                break;
            }

            // error would be small enough at upper_ix, but not at lower_ix,
            // so bisect to find optimal gridsize
            p3m_int lower_ix = step_good_gridsize[step_ix-1];
            while (lower_ix+1 < upper_ix) {
                p3m_int test_ix = (lower_ix+upper_ix)/2;
                p3m_int grid1d = good_gridsize[test_ix];
                (*p)->p.grid[0] = (*p)->p.grid[1] = (*p)->p.grid[2] = grid1d;
                P3M_DEBUG(printf("      fine grid=" FINT "\n", grid1d));
                this->tuneBroadcastCommand(CMD_COMPUTE_ERROR_ESTIMATE);
                errorEstimate->compute_master((*p)->p, sum_qpart, sum_q2, box_l, (*p)->error, (*p)->rs_error, (*p)->ks_error, box_vectors, isTriclinic);
                P3M_DEBUG(printf("          => error=" FFLOATE "\n", (*p)->error));
                if ((*p)->error < tolerance_field) {
                    // parameters achieve error
                    upper_ix = test_ix;
                } else {
                    // parameters do not achieve error
                    lower_ix = test_ix;
                }
            }

            // now the right size is at upper_ix
            p3m_int grid1d = good_gridsize[upper_ix];

            // store the new grid size and alpha
            if (min_grid1d > grid1d) min_grid1d = grid1d;
            (*p)->p.grid[0] = (*p)->p.grid[1] = (*p)->p.grid[2] = grid1d;
            P3M_INFO(printf( "      => grid=" F3INT ", "                      \
                    "error=" FFLOATE "\n",                           \
                    (*p)->p.grid[0], (*p)->p.grid[1], (*p)->p.grid[2],     \
                    (*p)->error));

            // decrease step_ix so that the same step_ix is tested for the
            // next param set
            step_ix--;

            // compare grid size to previous data set
            // if it is larger than any previous size + P3M_MAX_GRID_DIFF, skip it
            if (min_grid1d + P3M_MAX_GRID_DIFF < grid1d) {
                P3M_INFO(printf("      grid too large => removing data set\n"));
                // remove the rest of the params
                for (TuneParametersList::iterator pit = p;
                        pit != params_to_try.end(); ++pit) {
                    delete[] *pit;
                    *pit = NULL;
                }
                break;
            }

            if (p != params_to_try.begin()) {
                TuneParametersList::iterator prev = p;
                prev--;
                if ((*prev)->p.grid[0] >= (*p)->p.grid[0]) {
                    P3M_INFO(printf("      better than previous => removing previous data set.\n"));
                    delete[] *prev;
                    *prev = NULL;
                    params_to_try.erase(prev);
                }
            }
        }

        // remove empty param sets
        params_to_try.remove(NULL);


    } else {
        // !tune_grid
        for (TuneParametersList::iterator p = params_to_try.begin();
                p != params_to_try.end(); ++p) {
            // test whether accuracy can be achieved with given parameters
            P3M_INFO(printf("    grid=" F3INT " (fixed)\n",                \
                    (*p)->p.grid[0], (*p)->p.grid[1], (*p)->p.grid[2]));
            this->tuneBroadcastCommand(CMD_COMPUTE_ERROR_ESTIMATE);
            errorEstimate->compute_master((*p)->p, sum_qpart, sum_q2, box_l,
                    (*p)->error, (*p)->rs_error, (*p)->ks_error,box_vectors,isTriclinic);

            if ((*p)->error < tolerance_field) {
                // error is small enough for this parameter set, so keep it
                P3M_INFO(printf("    error (%le) < tolerance (%le), keeping params\n", \
                        (*p)->error, tolerance_field));
            } else {
                // otherwise remove this parameter set
                P3M_INFO(printf("    error (%le) > tolerance (%le), removing params\n", \
                        (*p)->error, tolerance_field));
                delete[] *p;
                *p = NULL;
            }
        }

        // really remove the removed param sets
        params_to_try.remove(NULL);
    }
}

Solver::TuneParameters*
Solver::timeParams(
        p3m_int num_particles, p3m_float *positions, p3m_float *charges,
        TuneParametersList &params_to_try) {
    /* Now time the different parameter sets */
    TuneParameters *best_params = NULL;
    double best_timing = 1.e100;

#ifdef P3M_ENABLE_INFO
    printf("Timing %ld param sets...\n", params_to_try.size());
    for (TuneParametersList::iterator p = params_to_try.begin();
            p != params_to_try.end();
            p++) {
        printf( "  r_cut=" FFLOAT ", "
                "alpha=" FFLOAT ", "
                "grid=" F3INT ", "
                "cao=" FINT "\n",
                (*p)->p.r_cut, (*p)->p.alpha, (*p)->p.grid[0],
                (*p)->p.grid[1], (*p)->p.grid[2], (*p)->p.cao);
    }
#endif

    for (TuneParametersList::iterator it = params_to_try.begin();
            it != params_to_try.end();
            it++) {
        TuneParameters &p = **it;
        /* use the parameters */
        alpha = p.p.alpha;
        grid[0] = p.p.grid[0];
        grid[1] = p.p.grid[1];
        grid[2] = p.p.grid[2];
        cao = p.p.cao;
        r_cut = p.p.r_cut;

        this->timing(num_particles, positions, charges);
        p.timing = timings[TIMING];
        p.timing_near = timings[TIMING_NEAR];
        p.timing_far = timings[TIMING_FAR];
        P3M_INFO(printf( "  Timing r_cut=" FFLOAT ", "                  \
                "alpha=" FFLOAT ", "                           \
                "grid=" F3INT ", "                             \
                "cao=" FINT " "                                \
                "=> timing=" FFLOAT " "                        \
                "(" FFLOAT " near, " FFLOAT " far)\n",         \
                r_cut, alpha,                            \
                grid[0], grid[1], grid[2],            \
                cao, p.timing,                \
                p.timing_near, p.timing_far));

        if (p.timing < best_timing) {
            best_timing = p.timing;
            best_params = &p;
        }
    }

    if (best_params == NULL)
        throw std::runtime_error("Internal error: Tuning could not get a best timing.");

    for (TuneParametersList::iterator p = params_to_try.begin();
            p != params_to_try.end(); ++p) {
        if (*p != best_params)
            delete[] *p;
    }

    return best_params;
}

void Solver::timing(p3m_int num_particles,
        p3m_float *positions, p3m_float *charges) {
    if (comm.rank == 0) {
        Parameters p;
        p.cao = cao;
        p.alpha = alpha;
        p.r_cut = r_cut;
        p.grid[0] = grid[0];
        p.grid[1] = grid[1];
        p.grid[2] = grid[2];
        this->tuneBroadcastCommand(CMD_TIMING, p);

    }

    p3m_float *fields = new p3m_float[3*num_particles];
    p3m_float *potentials = new p3m_float[num_particles];

    this->prepare();

    /* store require_timings */
    timingEnum require_timings_before = require_timings;
    if(require_timings == NONE || require_timings == FULL)
        require_timings = ESTIMATE_ALL;
    this->run(num_particles, positions, charges, fields, potentials);
    /* Afterwards, d->timings is set */
    /* restore require_timings */
    require_timings = require_timings_before;

    delete[] fields;
    delete[] potentials;
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

void Solver::tuneBroadcastParams(Parameters &p) {
    // broadcast parameters and mark as final
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
}

void Solver::tuneReceiveParams()
{
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

void Solver::tuneBroadcastCommand(p3m_int command){
    /* Send the command */
    P3M_DEBUG_LOCAL(printf("       %2d: Broadcasting command %d.\n", \
            comm.rank, command));
    MPI_Bcast(&command, 1, P3M_MPI_INT, 0, comm.mpicomm);

    /* Now send the parameters, depending on the command */
    switch (command) {
    case CMD_FINISHED:
    case CMD_TIMING:
        throw std::logic_error("Internal Error: Bad call to tuneBroadcastCommand()");
    case CMD_COMPUTE_ERROR_ESTIMATE:
    case CMD_FAILED:
        return;
    }
}

void Solver::tuneBroadcastCommand(p3m_int command, Parameters &p){
    /* Send the command */
    P3M_DEBUG_LOCAL(printf("       %2d: Broadcasting command %d.\n", \
            comm.rank, command));
    MPI_Bcast(&command, 1, P3M_MPI_INT, 0, comm.mpicomm);

    /* Now send the parameters, depending on the command */
    switch (command) {
    case CMD_FINISHED:
    case CMD_TIMING:
        this->tuneBroadcastParams(p);
        return;
    case CMD_COMPUTE_ERROR_ESTIMATE:
    case CMD_FAILED:
        throw std::logic_error("Internal Error: Bad call to tuneBroadcastCommand()");
        return;
    }
}

void Solver::tuneBroadcastSlave(p3m_int num_particles, p3m_float* positions,
        p3m_float* charges) {
    P3M_DEBUG(printf( "P3M::Solver::tuneBroadcastSlave() started...\n"));
    if (comm.onMaster())
        throw std::logic_error("Internal error: tuneBroadcastSlave "
                "should not be called on master!");

    for (;;) {
        /* Receive the command */
        p3m_int command;
        P3M_DEBUG_LOCAL(printf("      %2d: Waiting to receive command.\n", \
                comm.rank));
        MPI_Bcast(&command, 1, P3M_MPI_INT, 0, comm.mpicomm);
        P3M_DEBUG_LOCAL(printf("      %2d: Received command %d.\n", \
                comm.rank, command));

        switch (command) {
        case CMD_COMPUTE_ERROR_ESTIMATE:
            errorEstimate->compute_slave();
            break;
        case CMD_TIMING:
            this->tuneReceiveParams();
            timing(num_particles, positions, charges);
            break;
        case CMD_FINISHED:
            this->tuneReceiveParams();
            return;
        case CMD_FAILED: {
            char msg[255];
            sprintf(msg,
                    "Cannot achieve required accuracy (p3m_tolerance_field=" FFLOATE \
                    ") for given parameters.",
                    tolerance_field);
            throw std::logic_error(msg);
        }
        }
    }
    P3M_DEBUG(printf( "P3M::Solver::tuneBroadcastSlave() finished.\n"));
}

}
