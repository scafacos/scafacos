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
    P3M_DEBUG(printf( "P3M::P3M() started...\n"));

    errorEstimate = ErrorEstimate::create(comm);

    /* Init the P3M parameters */
    box_l[0] = 1.0;
    box_l[1] = 1.0;
    box_l[2] = 1.0;
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
    needs_retune = 1;
    tune_r_cut = 1;
    tune_alpha = 1;
    tune_grid = 1;
    tune_cao = 1;

    /* Which components to compute? */
    require_total_energy = 0;
    total_energy = 0.0;

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

    P3M_DEBUG(printf( "P3M::P3M() finished.\n"));
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
void print_local_grid(local_grid_t l);
void print_send_grid(send_grid_t sm);
#endif

/** Prepare the data structures and constants of the P3M algorithm.
     All parameters have to be set. */
void Solver::prepare() {
    P3M_DEBUG(printf("  P3M::Solver::prepare() started... \n"));

    /* initializes the (inverse) grid constant a
          (ai) and the cutoff for charge assignment
          cao_cut */
    this->prepare_a_ai_cao_cut();
    this->calc_local_ca_grid();
    this->calc_send_grid();
    P3M_DEBUG(this->print_local_grid());
    P3M_DEBUG(this->print_send_grid());
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
    this->calc_differential_operator();
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

/** Initializes the (inverse) grid constant \ref a (\ref
         ai) and the cutoff for charge assignment \ref
         cao_cut, which has to be done by \ref init_charges
         once and by \ref scaleby_box_l whenever the \ref box_l
         changed.  */
void Solver::prepare_a_ai_cao_cut() {
    P3M_DEBUG(printf("    Solver::prepare_a_ai_cao_cut() started... \n"));
    for (p3m_int i=0; i<3; i++) {
        ai[i]      = (p3m_float)grid[i]/box_l[i];
        a[i]       = 1.0/ai[i];
        cao_cut[i] = 0.5*a[i]*cao;
    }
    P3M_DEBUG(printf("    Solver::prepare_a_ai_cao_cut() finished. \n"));
}

/** Calculate the spacial position of the left down grid point of the
         local grid, to be stored in \ref local_grid::ld_pos; function
         called by \ref calc_local_ca_grid once and by \ref
         scaleby_box_l whenever the \ref box_l changed. */
void Solver::calc_lm_ld_pos() {
    /* spacial position of left bottom grid point */
    for (int i=0;i<3;i++) {
        local_grid.ld_pos[i] =
                (local_grid.ld_ind[i]+ grid_off[i])*a[i];
    }
}

/** Calculates properties of the local FFT grid for the
         charge assignment process. */
void Solver::calc_local_ca_grid() {
    p3m_int ind[3];
    /* total skin size */
    p3m_float full_skin[3];

    P3M_DEBUG(printf("    calc_local_ca_grid() started... \n"));
    for (int i=0;i<3;i++)
        full_skin[i]= cao_cut[i]+skin+additional_grid[i];

    /* inner left down grid point (global index) */
    for (int i=0;i<3;i++)
        local_grid.in_ld[i] =
                (p3m_int)ceil(comm.my_left[i]*ai[i]-grid_off[i]);
    /* inner up right grid point (global index) */
    for (int i=0;i<3;i++)
        local_grid.in_ur[i] =
                (p3m_int)floor(comm.my_right[i]*ai[i]-grid_off[i]);

    /* correct roundof errors at boundary */
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
    /* index of left down grid point in global grid */
    for (int i=0; i<3; i++)
        local_grid.ld_ind[i] =
                (p3m_int)ceil((comm.my_left[i]-full_skin[i])*ai[i]-grid_off[i]);
    /* spatial position of left down grid point */
    this->calc_lm_ld_pos();

    /* left down margin */
    for (int i=0;i<3;i++)
        local_grid.margin[i*2] = local_grid.in_ld[i]-local_grid.ld_ind[i];
    /* up right grid point */
    for (int i=0;i<3;i++)
        ind[i] = (p3m_int)floor((comm.my_right[i]+full_skin[i])*ai[i]-grid_off[i]);
    /* correct roundof errors at up right boundary */
    for (int i=0;i<3;i++)
        if (((comm.my_right[i]+full_skin[i])*ai[i]-grid_off[i])-ind[i]==0)
            ind[i]--;
    /* up right margin */
    for (int i=0;i<3;i++) local_grid.margin[(i*2)+1] = ind[i] - local_grid.in_ur[i];

    /* grid dimension */
    local_grid.size=1;
    for (int i=0;i<3;i++) {
        local_grid.dim[i] = ind[i] - local_grid.ld_ind[i] + 1;
        local_grid.size *= local_grid.dim[i];
    }
    /* reduce inner grid indices from global to local */
    for (int i=0;i<3;i++)
        local_grid.in_ld[i] = local_grid.margin[i*2];
    for (int i=0;i<3;i++)
        local_grid.in_ur[i] = local_grid.margin[i*2]+local_grid.inner[i];

    local_grid.q_2_off  = local_grid.dim[2] - cao;
    local_grid.q_21_off = local_grid.dim[2] * (local_grid.dim[1] - cao);

    P3M_DEBUG(printf("    calc_local_ca_grid() finished. \n"));
}

/** Calculates the properties of the send/recv sub-grides of the local
 *  FFT grid.  In order to calculate the recv sub-grides there is a
 *  communication of the margins between neighbouring nodes. */
void Solver::calc_send_grid() {
    p3m_int evenodd;
    p3m_int done[3]={0,0,0};

    P3M_DEBUG(printf("    P3M::Solver::calc_send_grid() started... \n"));
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

    // fixme
    /* communication */
    for (int i=0; i<6; i++) {
        int j = ((i%2==0) ? i+1 : i-1);
        if (comm.node_neighbors[i] != comm.rank) {
            /* two step communication: first all even positions than all odd */
            for (evenodd=0; evenodd<2; evenodd++) {
                if ((comm.node_pos[i/2]+evenodd)%2 == 0) {
                    P3M_DEBUG(printf("      %d: sending local_grid.margin to %d\n", \
                            comm.rank, comm.node_neighbors[i]));
                    MPI_Send(&(local_grid.margin[i]), 1, P3M_MPI_INT,
                            comm.node_neighbors[i], 0, comm.mpicomm);
                } else {
                    P3M_DEBUG(printf("      %d: receiving local_grid.margin from %d\n", \
                            comm.rank, comm.node_neighbors[j]));
                    MPI_Recv(&(local_grid.r_margin[j]), 1, P3M_MPI_INT,
                            comm.node_neighbors[j], 0, comm.mpicomm, MPI_STATUS_IGNORE);
                }
            }
        } else {
            local_grid.r_margin[j] = local_grid.margin[i];
        }
    }

    /* /\* communication *\/ */
    /* for (i = 0; i < 3; i++) { */
    /*   /\* upshift *\/ */
    /*   MPI_Sendrecv(&(local_grid.margin[2*i]), 1, P3M_MPI_INT, */
    /*       comm.node_neighbors[2*i+1], 0, */
    /*       &(local_grid.r_margin[2*i]), 1, P3M_MPI_INT, */
    /*       comm.node_neighbors[2*i], 0, */
    /*       comm.mpicomm, &status); */
    /*   /\* downshift *\/ */
    /*   MPI_Sendrecv(&(local_grid.margin[2*i+1]), 1, P3M_MPI_INT, */
    /*       comm.node_neighbors[2*i], 0, */
    /*       &(local_grid.r_margin[2*i+1]), 1, P3M_MPI_INT, */
    /*       comm.node_neighbors[2*i+1], 0, */
    /*       comm.mpicomm, &status); */
    /* } */

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
    P3M_DEBUG(printf("    P3M::Solver::calc_send_grid() finished. \n"));
}

/** Calculates the Fourier transformed differential operator.
 *  Remark: This is done on the level of n-vectors and not k-vectors,
 *           i.e. the prefactor i*2*PI/L is missing! */
void Solver::calc_differential_operator() {
    for (int i = 0; i < 3; i++) {
        d_op[i] = static_cast<p3m_int *>(realloc(d_op[i], grid[i]*sizeof(p3m_int)));
        d_op[i][0] = 0;
        d_op[i][grid[i]/2] = 0;

        for (p3m_int j = 1; j < grid[i]/2; j++) {
            d_op[i][j] = j;
            d_op[i][grid[i] - j] = -j;
        }
    }
}

/** Debug function printing p3m structures */
void Solver::print_local_grid() {
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
void Solver::print_send_grid() {
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
                    p3m_float denominator;
                    this->performAliasingSumsIK(
                            gridshift_x[n[KX]],
                            gridshift_y[n[KY]],
                            gridshift_z[n[KZ]],
                            numerator_force,
                            numerator_energy,
                            denominator);

                    const p3m_float fak1 =
                            d_op[RX][n[KX]]*numerator_force[RX]/box_l[RX] +
                            d_op[RY][n[KY]]*numerator_force[RY]/box_l[RY] +
                            d_op[RZ][n[KZ]]*numerator_force[RZ]/box_l[RZ];
                    const p3m_float fak2 =
                            SQR(d_op[RX][n[KX]]/box_l[RX]) +
                            SQR(d_op[RY][n[KY]]/box_l[RY]) +
                            SQR(d_op[RZ][n[KZ]]/box_l[RZ]);
                    const p3m_float fak3 = fak1/(fak2 * SQR(denominator));
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
        p3m_float &denominator) {
    denominator = 0.0;
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

                denominator  += sz;
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
         const p3m_float sz  = sy*pow(sinc(nmz/(p3m_float)grid[RZ]), 2.0*cao);
         const p3m_float nm2 =
           SQR(nmx/box_l[RX]) +
           SQR(nmy/box_l[RY]) +
           SQR(nmz/box_l[RZ]);
         const p3m_float prefactor2 = sz*exp(-prefactor*nm2);

         numerator_energy += prefactor2/nm2;
         numerator_force  += prefactor2;

         denominator[0] += sz;
         denominator[1] += sz * nm2;

         if(((mx+my+mz) % 2) == 0) {
           denominator[2] += sz;
           denominator[3] += sz * nm2;
         } else {
           denominator[2] -= sz;
           denominator[3] -= sz * nm2;
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

    P3M_INFO(printf("    system parameters: box_l=" F3FLOAT "\n", \
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
    this->domain_decompose(&gridsort,
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
    this->compute_far_adi(num_real_particles, positions, charges, fields, potentials);
#else
    this->compute_far_ik(num_real_particles, positions, charges, fields, potentials);
#endif

    if (near_field_flag) {
        /* start near timer */
        START(TIMING_NEAR);

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
    P3M_INFO(printf( "P3M::Solver::run() finished.\n"));
}

void Solver::compute_far_adi(
        p3m_int num_charges, p3m_float* positions, p3m_float* charges,
        p3m_float* fields, p3m_float* potentials) {
    P3M_INFO(printf( "P3M::Solver::compute_far_adi() started...\n"));

    START(TIMING_CA);

    /* charge assignment */
    this->assign_charges(fft.data_buf, num_charges, positions, charges, 0);
    STOPSTART(TIMING_CA, TIMING_GATHER);
    /* gather the ca grid */
    this->gather_grid(fft.data_buf);

    // Complexify
    for (p3m_int i = local_grid.size-1; i >= 0; i--)
        rs_grid[2*i] = fft.data_buf[i];

    STOPSTART(TIMING_GATHER, TIMING_CA);

    // Second (shifted) run
    /* charge assignment */
    this->assign_charges(fft.data_buf, num_charges, positions, charges, 1);

    STOPSTART(TIMING_CA, TIMING_GATHER);

    /* gather the ca grid */
    this->gather_grid(fft.data_buf);
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
        this->apply_energy_influence_function();
        /* result is in ks_grid */
        STOP(TIMING_INFLUENCE);

        /* compute total energy, but not potentials */
        if (require_total_energy && potentials == NULL) {
            START(TIMING_POTENTIALS);
            total_energy = this->compute_total_energy();
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

            this->spread_grid(fft.data_buf);

            STOPSTART(TIMING_SPREAD, TIMING_POTENTIALS);

            this->assign_potentials(fft.data_buf,
                    num_charges, positions,
                    charges, 0, potentials);

            STOPSTART(TIMING_POTENTIALS, TIMING_SPREAD);

            /** Second (shifted) run */
            for (p3m_int i=0; i<local_grid.size; i++) {
                fft.data_buf[i] = ks_grid[2*i+1];
            }
            this->spread_grid(fft.data_buf);

            STOPSTART(TIMING_SPREAD, TIMING_POTENTIALS);

            this->assign_potentials(fft.data_buf,
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
        this->apply_force_influence_function();
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

        this->spread_grid(fft.data_buf);

        STOPSTART(TIMING_SPREAD, TIMING_FIELDS);

        this->assign_fields_ad(fft.data_buf, num_charges, positions, 0, fields);

        STOPSTART(TIMING_FIELDS, TIMING_SPREAD);

        /* Second (shifted) run */
        P3M_INFO(printf("  computing shifted grid\n"));
        for (p3m_int i=0; i<local_grid.size; i++) {
            fft.data_buf[i] = ks_grid[2*i+1];
        }

        this->spread_grid(fft.data_buf);

        STOPSTART(TIMING_SPREAD, TIMING_FIELDS);

        this->assign_fields_ad(fft.data_buf, num_charges, positions, 1, fields);

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

    /* collect timings from the different nodes */
    if (require_timings!=NONE) this->collect_print_timings();

    P3M_INFO(printf( "P3M::Solver::compute_far_adi() finished.\n"));
}

void Solver::compute_far_ik(
        p3m_int num_charges, p3m_float* positions, p3m_float* charges,
        p3m_float* fields, p3m_float* potentials) {
    P3M_INFO(printf( "P3M::Solver::compute_far_ik() started...\n"));

    START(TIMING_CA);

    /* charge assignment */
    this->assign_charges(rs_grid, num_charges, positions, charges, 0);
    STOPSTART(TIMING_CA, TIMING_GATHER);
    /* gather the ca grid */
    this->gather_grid(rs_grid);
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
        this->apply_energy_influence_function();
        /* result is in ks_grid */
        STOP(TIMING_INFLUENCE);

        /* compute total energy, but not potentials */
        if (require_total_energy && potentials == NULL) {
            START(TIMING_POTENTIALS);
            total_energy = this->compute_total_energy();
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
            this->spread_grid(ks_grid);

            STOPSTART(TIMING_SPREAD, TIMING_POTENTIALS);

            /* compute potentials */
            this->assign_potentials(ks_grid,
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
        this->apply_force_influence_function();
        STOP(TIMING_INFLUENCE);

        /* result is in ks_grid */
        for (int dim = 0; dim < 3; dim++) {
            /* differentiate in direction dim */
            /* result is stored in rs_grid */
            START(TIMING_FIELDS);
            this->ik_diff(dim);
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
            this->spread_grid(rs_grid);
            STOPSTART(TIMING_SPREAD, TIMING_FIELDS);
            this->assign_fields_ik(rs_grid, dim, num_charges, positions, 0, fields);
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

    /* collect timings from the different nodes */
    if (require_timings != NONE) this->collect_print_timings();

    P3M_INFO(printf( "P3M::Solver::compute_far_ik() finished.\n"));
}


/***************************************************/
/* RUN COMPONENTS */

/* domain decomposition */
void Solver::domain_decompose(fcs_gridsort_t *gridsort,
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

/* CHARGE ASSIGNMENT */
/** Compute the data of the charge assignment grid points.

      The function returns the linear index of the top left grid point
      in the charge assignment grid that corresponds to real_pos. When
      "shifted" is set, it uses the shifted position for interlacing.
      After the call, caf_cache contains a cache of the values of the
      charge assignment fraction (caf) for x,y,z.
 */
p3m_int
Solver::get_ca_points(p3m_float real_pos[3], p3m_int shifted) {
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
void Solver::assign_charges(p3m_float *data,
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
        p3m_int linind_grid = this->get_ca_points(&positions[pid*3], shifted);

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
void Solver::gather_grid(p3m_float* rs_grid) {
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

void Solver::spread_grid(p3m_float* rs_grid) {
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
void Solver::apply_energy_influence_function() {
    P3M_DEBUG(printf( "  apply_energy_influence_function() started...\n"));
    const p3m_int size = fft.plan[3].new_size;
    for (p3m_int i=0; i < size; i++) {
        ks_grid[2*i] = g_energy[i] * rs_grid[2*i];
        ks_grid[2*i+1] = g_energy[i] * rs_grid[2*i+1];
    }
    P3M_DEBUG(printf( "  apply_energy_influence_function() finished.\n"));
}

/* apply the influence function */
void Solver::apply_force_influence_function() {
    P3M_DEBUG(printf( "  apply_force_influence_function() started...\n"));
    const p3m_int size = fft.plan[3].new_size;
    for (p3m_int i=0; i < size; i++) {
        ks_grid[2*i] = g_force[i] * rs_grid[2*i];
        ks_grid[2*i+1] = g_force[i] * rs_grid[2*i+1];
    }
    P3M_DEBUG(printf( "  apply_force_influence_function() finished.\n"));
}

/* Add up the measured timings and collect information from all nodes.
 * Print the timings if requested so.
 */
void Solver::collect_print_timings() {
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
        //    printf(" (empirical estimate)");
        printf("\n");
        printf("   fields=%le (%lf)", timings[TIMING_FIELDS],
                timings[TIMING_FIELDS]/timings[TIMING]);
        //if(require_timings == ESTIMATE_ALL //not yet implemented
        //|| require_timings == ESTIMATE_ASSIGNMENT)
        //    printf(" (empirical estimate)");
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
            printf(" (theoretical estimate)");
        printf("\n");
        printf("   decomp=%le (%lf)\n", timings[TIMING_DECOMP],
                timings[TIMING_DECOMP]/timings[TIMING]);
        printf("     comp=%le (%lf)\n", timings[TIMING_COMP],
                timings[TIMING_COMP]/timings[TIMING]);
    }
#endif

}

void Solver::ik_diff(int dim) {
    /* direction in k space: */
    p3m_int dim_rs = (dim+ks_pnum)%3;

    P3M_DEBUG(printf( "  P3M::Solver::ik_diff() started...\n"));
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
                rs_grid[ind] =
                        -2.0*M_PI*(ks_grid[ind+1] * d_operator[ j[dim]+fft.plan[3].start[dim] ])
                        / box_l[dim_rs];
                rs_grid[ind+1] =
                        2.0*M_PI*ks_grid[ind] * d_operator[ j[dim]+fft.plan[3].start[dim] ]
                                                               / box_l[dim_rs];
                ind+=2;
            }
        }
    }

    P3M_DEBUG(printf( "  P3M::Solver::ik_diff() finished.\n"));
    /* store the result in rs_grid */
}


/** Compute the total energy of the system in kspace. No need to
      backtransform the FFT grid in this case! */
p3m_float Solver::compute_total_energy() {
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

    P3M_DEBUG(printf( "  P3M::Solver::compute_total_energy() finished.\n"));
    return k_space_energy;
}

/* Backinterpolate the potentials obtained from k-space to the positions */
void
Solver::assign_potentials(p3m_float *data,
        p3m_int num_real_particles, p3m_float* positions, p3m_float* charges,
        p3m_int shifted, p3m_float* potentials) {
    const p3m_int q2off = local_grid.q_2_off;
    const p3m_int q21off = local_grid.q_21_off;
    const p3m_float prefactor = 1.0 / (box_l[0] * box_l[1] * box_l[2]);

    P3M_DEBUG(printf( "  P3M::Solver::assign_potentials() started...\n"));
    /* Loop over all particles */
    for (p3m_int pid=0; pid < num_real_particles; pid++) {
        p3m_float potential = 0.0;
        p3m_int linind_grid =
                this->get_ca_points(&positions[pid*3], shifted);

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
    P3M_DEBUG(printf( "  P3M::Solver::assign_potentials() finished.\n"));
}

/* Backinterpolate the forces obtained from k-space to the positions */
void
Solver::assign_fields_ik(p3m_float *data, p3m_int dim,
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
                this->get_ca_points(&positions[3*pid], shifted);

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
Solver::assign_fields_ad(p3m_float *data,
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
                this->get_ca_points(&positions[pid*3], shifted);

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
    /* Prepare the communicator before tuning */
    comm.prepare(box_l);

    /* Count the charges */
    p3m_float sum_q2_before = sum_q2;
    this->count_charges(num_particles, charges);

    if (!comm.onMaster()) {
        int howoften;
        P3M_DEBUG_LOCAL(printf("  %d: How often to run tune_far?\n", comm.rank));
        MPI_Bcast(&howoften, 1, MPI_INT, Communication::MPI_MASTER, comm.mpicomm);
        P3M_DEBUG_LOCAL(printf("  %d: Running tune_far %d times.\n", comm.rank, howoften));
        for (; howoften > 0; howoften--)
            this->tune_far(num_particles, positions, charges);

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
        this->tune_far(num_particles, positions, charges, r_cut);
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
        tune_params *p =
                this->tune_far(num_particles, positions, charges, r_cut);

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
        p = this->tune_far(num_particles, positions, charges, r_cut);

        rel_timing_diff =
                fabs(p->timing_near - p->timing_far) /
                (p->timing_near + p->timing_far);
        P3M_INFO(printf("    rel_timing_diff=" FFLOAT "\n", rel_timing_diff));
        P3M_INFO(printf("    Finished tuning.\n"));
    }

    /* mark that the method was retuned */
    needs_retune = false;
}

Solver::tune_params*
Solver::tune_far(p3m_int num_particles, p3m_float *positions, p3m_float *charges,
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
    P3M_INFO(printf( "P3M::Solver::tune_far() started...\n"));

    tune_params *final_params = NULL;

    if (!comm.onMaster()) {
        tune_broadcast_slave(this, num_particles, positions, charges);
    } else {

        /* check whether the input parameters are sane */
        if (_r_cut < 0.0)
            throw std::domain_error("r_cut is negative!");

        if (float_is_zero(_r_cut))
            throw std::domain_error("r_cut is too small!");

        /* check whether cutoff is larger than half a box length */
        if ((_r_cut > 0.5*box_l[0]) ||
                (_r_cut > 0.5*box_l[1]) ||
                (_r_cut > 0.5*box_l[2]))
            throw std::domain_error("r_cut is larger than half a system box length.");

        /* check whether cutoff is larger than domain size */

        try {
            this->r_cut = _r_cut;
            final_params =
                    this->tune_alpha_cao_grid(num_particles, positions,
                            charges, _r_cut);

            // @todo: move to tune function
            /* Set and broadcast the final parameters. */
            this->r_cut = _r_cut;
            alpha = final_params->alpha;
            grid[0] = final_params->grid[0];
            grid[1] = final_params->grid[1];
            grid[2] = final_params->grid[2];
            cao = final_params->cao;
            tune_broadcast_command(this, CMD_FINISHED);
        } catch (...) {
            P3M_INFO(printf( "  Tuning failed.\n"));
            tune_broadcast_command(this, CMD_FAILED);
            P3M_INFO(printf( "P3M::Solver::tune_far() finished.\n"));
            throw;
        }
    }

    P3M_INFO(printf( "  Tuning was successful.\n"));

    /* At the end of retuning, prepare the method */
    this->prepare();

    P3M_INFO(printf( "P3M::Solver::tune_far() finished.\n"));

    return final_params;
}

/** Slave variant of tune_far. */
void
Solver::tune_far(p3m_int num_particles, p3m_float *positions, p3m_float *charges) {
    if (comm.onMaster()) throw std::runtime_error("tune_far without r_cut cannot be called on master node.");
    this->tune_far(num_particles, positions, charges, 0.0);
}


/* Tune alpha */
Solver::tune_params*
Solver::tune_alpha_cao_grid(p3m_int num_particles, p3m_float *positions,
        p3m_float *charges, p3m_float r_cut) {
    if (tune_alpha) {
        // fixme
        Parameters p;
        p.cao = cao;
        p.alpha = alpha;
        p.r_cut = r_cut;
        p.grid[0] = grid[0];
        p.grid[1] = grid[1];
        p.grid[2] = grid[2];
        errorEstimate->compute_alpha(tolerance_field,
                p, sum_qpart, sum_q2, box_l);
        alpha = p.alpha;
        P3M_INFO(printf("    => alpha=" FFLOAT "\n", alpha));
    } else {
        P3M_INFO(printf("    alpha=" FFLOAT " (fixed)\n", alpha));
    }

    return this->tune_cao_grid(num_particles, positions, charges, r_cut, alpha);
}

/* Tune cao */
Solver::tune_params*
Solver::tune_cao_grid(p3m_int num_particles, p3m_float *positions, p3m_float *charges,
        p3m_float r_cut, p3m_float alpha) {
#ifdef P3M_AD
    const p3m_int min_cao = 2;
#else
    const p3m_int min_cao = 1;
#endif

    tune_params_l params_to_try;

    if (tune_cao) {
        P3M_INFO(printf("    Testing cao={ "));
        for (p3m_int cao = CAF::max_cao; cao >= min_cao; --cao) {
            tune_params *p = new tune_params;
            p->cao = cao;
            p->alpha = alpha;
            params_to_try.push_back(p);
            P3M_INFO(printf(FINT " ", cao));
        }
        P3M_INFO(printf("}\n"));
    } else {
        tune_params *p = new tune_params;
        p->cao = cao;
        p->alpha = alpha;
        params_to_try.push_back(p);
        P3M_INFO(printf( "    cao=" FINT " (fixed)\n", cao));
    }

    return this->tune_grid_(num_particles, positions, charges, params_to_try);
}

/* params should have decreasing cao */
Solver::tune_params*
Solver::tune_grid_(p3m_int num_particles,
        p3m_float *positions, p3m_float *charges,
        tune_params_l &params_to_try) {
    if (tune_grid) {
        p3m_int step_ix = 0;
        // store the minimal grid size seen so far
        p3m_int min_grid1d = good_gridsize[step_good_gridsize[num_steps_good_gridsize-1]];

        for (tune_params_l::iterator p = params_to_try.begin();
                p != params_to_try.end(); ++p) {
            // Find smallest possible grid for this set of parameters
            // Test whether accuracy can be achieved with given parameters
            cao = (*p)->cao;

            p3m_int upper_ix;
            // test this step
            P3M_INFO(printf("    Trying to find grid for r_cut=" FFLOAT ", " \
                    "alpha=" FFLOAT ", "                          \
                    "cao=" FINT "\n",                             \
                    r_cut, alpha, cao));
            p3m_float error, rs_error, ks_error;
            do {
                step_ix++;
                if (step_ix >= num_steps_good_gridsize) break;
                upper_ix = step_good_gridsize[step_ix];
                p3m_int grid1d = good_gridsize[upper_ix];
                P3M_DEBUG(printf("      rough grid=" FINT "\n", grid1d));
                // fixme
                Parameters p2;
                p2.grid[0] = grid1d;
                p2.grid[1] = grid1d;
                p2.grid[2] = grid1d;
                p2.cao = (*p)->cao;
                p2.alpha = (*p)->alpha;
                p2.r_cut = r_cut;
                tune_broadcast_command(this, CMD_COMPUTE_ERROR_ESTIMATE);
                errorEstimate->compute_master(p2, sum_qpart, sum_q2, box_l, error, rs_error, ks_error);
                (*p)->grid[0] = grid1d;
                (*p)->grid[1] = grid1d;
                (*p)->grid[2] = grid1d;
                (*p)->error = error;
                (*p)->rs_error = rs_error;
                (*p)->ks_error = ks_error;
                P3M_DEBUG(printf("        => error=" FFLOATE "\n", error));
            } while (error > tolerance_field);

            // reached largest possible grid, remove the rest of the parameter sets
            if (step_ix >= num_steps_good_gridsize) {
                P3M_INFO(printf("    Too large grid size, skipping rest of parameter sets.\n"));
                for (tune_params_l::iterator pit = params_to_try.end();
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
                P3M_DEBUG(printf("      fine grid=" FINT "\n", grid1d));
                // fixme
                Parameters p2;
                p2.grid[0] = grid1d;
                p2.grid[1] = grid1d;
                p2.grid[2] = grid1d;
                p2.cao = (*p)->cao;
                p2.alpha = (*p)->alpha;
                p2.r_cut = r_cut;
                tune_broadcast_command(this, CMD_COMPUTE_ERROR_ESTIMATE);
                errorEstimate->compute_master(p2, sum_qpart, sum_q2, box_l, error, rs_error, ks_error);
                (*p)->grid[0] = grid1d;
                (*p)->grid[1] = grid1d;
                (*p)->grid[2] = grid1d;
                (*p)->error = error;
                (*p)->rs_error = rs_error;
                (*p)->ks_error = ks_error;
                P3M_DEBUG(printf("          => error=" FFLOATE "\n", error));
                if (error < tolerance_field) {
                    // parameters achieve error
                    upper_ix = test_ix;
                    (*p)->error = error;
                    (*p)->rs_error = rs_error;
                    (*p)->ks_error = ks_error;
                } else {
                    // parameters do not achieve error
                    lower_ix = test_ix;
                }
            }

            // now the right size is at upper_ix
            p3m_int grid1d = good_gridsize[upper_ix];

            // store the new grid size and alpha
            if (min_grid1d > grid1d) min_grid1d = grid1d;

            (*p)->grid[0] = grid1d;
            (*p)->grid[1] = grid1d;
            (*p)->grid[2] = grid1d;
            P3M_INFO(printf( "      => grid=" F3INT ", "                      \
                    "error=" FFLOATE "\n",                           \
                    (*p)->grid[0], (*p)->grid[1], (*p)->grid[2],     \
                    (*p)->error));

            // decrease step_ix so that the same step_ix is tested for the
            // next param set
            step_ix--;

            // compare grid size to previous data set
            // if it is larger than any previous size + P3M_MAX_GRID_DIFF, skip it
            if (min_grid1d + P3M_MAX_GRID_DIFF < grid1d) {
                P3M_INFO(printf("      grid too large => removing data set\n"));
                // remove the rest of the params
                for (tune_params_l::iterator pit = p;
                        pit != params_to_try.end(); ++pit) {
                    delete[] *pit;
                    *pit = NULL;
                }
                break;
            }

            if (p != params_to_try.begin()) {
                tune_params_l::iterator prev = p;
                prev--;
                if ((*prev)->grid[0] >= (*p)->grid[0]) {
                    P3M_INFO(printf("      better than previous => removing previous data set.\n"));
                    delete[] *prev;
                    params_to_try.erase(prev);
                }
            }
        }

        // remove empty param sets
        params_to_try.remove(NULL);


    } else {
        // !tune_grid
        for (tune_params_l::iterator p = params_to_try.begin();
                p != params_to_try.end(); ++p) {
            // test whether accuracy can be achieved with given parameters
            // fixme
            Parameters p2;
            p3m_float error, rs_error, ks_error;
            p2.grid[0] = grid[0];
            p2.grid[1] = grid[1];
            p2.grid[2] = grid[2];
            p2.cao = (*p)->cao;
            p2.alpha = (*p)->alpha;
            p2.r_cut = r_cut;
            P3M_INFO(printf("    r_cut=" FFLOAT ", "                   \
                    "cao=" FINT ", "                           \
                    "grid=" F3INT " (fixed)\n",                \
                    p2.r_cut, p2.cao,                          \
                    p2.grid[0], p2.grid[1], p2.grid[2]));
            tune_broadcast_command(this, CMD_COMPUTE_ERROR_ESTIMATE);
            errorEstimate->compute_master(p2, sum_qpart, sum_q2, box_l, error, rs_error, ks_error);

            if (error < tolerance_field) {
                // error is small enough for this parameter set, so keep it
                P3M_INFO(printf("    error (%le) < tolerance (%le), keeping params\n", \
                        error, tolerance_field));
                (*p)->grid[0] = grid[0];
                (*p)->grid[1] = grid[1];
                (*p)->grid[2] = grid[2];
                (*p)->alpha = alpha;
                (*p)->error = error;
                (*p)->rs_error = rs_error;
                (*p)->ks_error = ks_error;
            } else {
                // otherwise remove this parameter set
                P3M_INFO(printf("    error (%le) > tolerance (%le), removing params\n", \
                        error, tolerance_field));
                delete[] *p;
                *p = NULL;
            }
        }

        // really remove the removed param sets
        params_to_try.remove(NULL);
    }

    return this->time_params(num_particles, positions, charges, params_to_try);
}

Solver::tune_params*
Solver::time_params(
        p3m_int num_particles, p3m_float *positions, p3m_float *charges,
        tune_params_l &params_to_try) {
    /* Now time the different parameter sets */
    tune_params *best_params = NULL;
    double best_timing = 1.e100;

#ifdef P3M_ENABLE_INFO
    printf("Timing %ld param sets...\n", params_to_try.size());
#endif

    for (tune_params_l::iterator it = params_to_try.begin();
            it != params_to_try.end();
            it++) {
        tune_params &p = **it;
        /* use the parameters */
        alpha = p.alpha;
        grid[0] = p.grid[0];
        grid[1] = p.grid[1];
        grid[2] = p.grid[2];
        cao = p.cao;

        timing(this, num_particles, positions, charges);
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

    for (tune_params_l::iterator p = params_to_try.begin();
            p != params_to_try.end(); ++p) {
        if (*p != best_params)
            delete[] *p;
    }

    return best_params;
}

/** Calculate number of charged particles, the sum of the squared
      charges and the squared sum of the charges. Called in parallel at
      the beginning of tuning. */
void Solver::count_charges(p3m_int num_particles, p3m_float *charges) {
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

    P3M_DEBUG(printf("  count_charges(): "          \
            "num_charged_particles=" FINT ", "                   \
            "sum_squared_charges=" FFLOAT ", "                   \
            "sum_charges=" FFLOAT "\n",                          \
            sum_qpart, sum_q2, sqrt(square_sum_q)));
}

}


