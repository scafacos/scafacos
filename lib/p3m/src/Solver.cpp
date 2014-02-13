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

namespace P3M {
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
    meshift_x = NULL;
    meshift_y = NULL;
    meshift_z = NULL;

    d_op[0] = NULL;
    d_op[1] = NULL;
    d_op[2] = NULL;
    g_force = NULL;
    g_energy = NULL;

    ks_pnum = 0;

    send_grid = NULL;
    recv_grid = NULL;


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
    sfree(g_energy);
    sfree(g_force);
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
    Communication &comm = comm;

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
    cafx_d = caf_createCache();
    cafy_d = caf_createCache();
    cafz_d = caf_createCache();
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
    calc_influence_function_ik(this);
#elif defined(P3M_INTERLACE) && defined(P3M_IK)
    calc_influence_function_iki(this);
#else
    calc_influence_function_adi(this);
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

}
