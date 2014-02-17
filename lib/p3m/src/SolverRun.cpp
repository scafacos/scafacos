/*
  Copyright (C) 2013 Olaf Lenz, Florian Weik
  Copyright (C) 2011,2012 Olaf Lenz
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
#include "utils.hpp"
#include "types.hpp"
#include "Solver.hpp"
#include <sys/time.h>
#include <sys/resource.h>
#include <cstdlib>
#include <cstdio>

#include "fcs_p3m_p.h"
#include "FCSCommon.h"
#include "common/near/near.h"

namespace P3M {

/***************************************************/
/* IMPLEMENTATION */
/***************************************************/
#define START(ID)                                               \
		if (require_timings!=NONE) timings[ID] += -MPI_Wtime();
#define STOP(ID)                                                \
		if (require_timings!=NONE) timings[ID] += MPI_Wtime();
#define STOPSTART(ID1, ID2)                     \
		if (require_timings!=NONE) {               \
			timings[ID1] += MPI_Wtime();             \
			timings[ID2] += -MPI_Wtime();            \
		}                                             \


/* callback function for near field computations */
inline void
compute_near(const void *param, p3m_float dist, p3m_float *f, p3m_float *p)
{
    p3m_float alpha = *((p3m_float *) param);
    fcs_p3m_compute_near(alpha, dist, p, f);
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
    P3M_DEBUG_LOCAL(printf("    %d: num_particles=%d\n",	\
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

}
