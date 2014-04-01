/*
 Copyright (C) 2014 Olaf Lenz
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
/** \file fft.c
 *
 *  Routines, row decomposition, data structures and communication for the 3D-FFT. 
 *
 */
#include "types.hpp"
#include "Parallel3DFFT.hpp"
#include "utils.hpp"
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <cassert>
#include <stdexcept>

#define fftw_complex  FFTW_MANGLE(complex)
#define fftw_free  FFTW_MANGLE(free)
#define fftw_malloc  FFTW_MANGLE(malloc)
#define fftw_plan_many_dft  FFTW_MANGLE(plan_many_dft)
#define fftw_destroy_plan  FFTW_MANGLE(destroy_plan)
#define fftw_execute  FFTW_MANGLE(execute)
#define fftw_execute_dft  FFTW_MANGLE(execute_dft)
#define fftw_import_wisdom_from_file  FFTW_MANGLE(import_wisdom_from_file)
#define fftw_export_wisdom_to_file  FFTW_MANGLE(export_wisdom_to_file)

namespace P3M {
/***************************************************/
/* FORWARD DECLARATIONS OF INTERNAL FUNCTIONS */
/***************************************************/
static int
find_comm_groups(Communication &comm, int grid1[3], int grid2[3], int *node_list1,
		int *node_list2, int *group, int *pos, int *my_pos);

static int
calc_local_grid(int n_pos[3], int n_grid[3], int grid[3], p3m_float grid_off[3],
		int loc_grid[3], int start[3]);

static int
calc_send_block(int pos1[3], int grid1[3], int pos2[3], int grid2[3],
		int grid[3], p3m_float grid_off[3], int block[6]);

static void
pack_block_permute1(p3m_float *in, p3m_float *out, int start[3], int size[3],
		int dim[3], int element);

static void
pack_block_permute2(p3m_float *in, p3m_float *out, int start[3], int size[3],
		int dim[3], int element);

/***************************************************/
/* IMPLEMENTATION */
/***************************************************/
/** Debug function to print forward_plan structure.
 * \param pl fft/communication plan (see \ref forward_plan).
 */
void Parallel3DFFT::forward_plan::print() {
	printf("      forward_plan:\n");
	printf("        dir=%d, row_dir=%d, n_permute=%d, n_ffts=%d\n", dir,
			row_dir, n_permute, n_ffts);

	printf(
			"        local: old_grid=(%d,%d,%d), new_grid=(%d,%d,%d), start=(%d,%d,%d)\n",
			old_grid[0], old_grid[1], old_grid[2], new_grid[0], new_grid[1],
			new_grid[2], start[0], start[1], start[2]);
	printf("        new_size %d\n", new_size);

	printf("        g_size=%d group=(", g_size);
	for (int i = 0; i < g_size - 1; i++)
		printf("%d,", group[i]);
	printf("%d)\n", group[g_size - 1]);

	printf("        send=[");
	for (int i = 0; i < g_size; i++)
		printf("(%d,%d,%d)+(%d,%d,%d), ", send_block[6 * i + 0],
				send_block[6 * i + 1], send_block[6 * i + 2],
				send_block[6 * i + 3], send_block[6 * i + 4],
				send_block[6 * i + 5]);
	printf("]\n");
	printf("        recv=[");
	for (int i = 0; i < g_size; i++)
		printf("(%d,%d,%d)+(%d,%d,%d), ", recv_block[6 * i + 0],
				recv_block[6 * i + 1], recv_block[6 * i + 2],
				recv_block[6 * i + 3], recv_block[6 * i + 4],
				recv_block[6 * i + 5]);
	printf("]\n");
}

Parallel3DFFT::Parallel3DFFT(Communication &comm) : comm(comm) {
	for (int i = 0; i < 4; i++) {
		plan[i].group = new int[comm.size];
		plan[i].send_block = NULL;
		plan[i].send_size = NULL;
		plan[i].recv_block = NULL;
		plan[i].recv_size = NULL;
	}

	is_prepared = false;
	max_comm_size = 0;
	max_grid_size = 0;
	send_buf = NULL;
	recv_buf = NULL;
	data_buf = NULL;
}

Parallel3DFFT::~Parallel3DFFT() {
	for (int i = 0; i < 4; i++) {
		sdelete(plan[i].group);
		sfree(plan[i].send_block);
		sfree(plan[i].send_size);
		sfree(plan[i].recv_block);
		sfree(plan[i].recv_size);
	}
	sfree(send_buf);
	sfree(recv_buf);
	if (data_buf != NULL) {
		fftw_free(data_buf);
		data_buf = NULL;
	}
}

p3m_float* Parallel3DFFT::malloc_data() {
	return static_cast<p3m_float*>(fftw_malloc(
			max_grid_size * sizeof(p3m_float)));
}

void Parallel3DFFT::free_data(p3m_float* data) {
	if (data != NULL)
		fftw_free(data);
}

void Parallel3DFFT::prepare(p3m_int *local_grid_dim, p3m_int *local_grid_margin,
		p3m_int* global_grid_dim, p3m_float *global_grid_off,
		p3m_int *ks_pnum) {
	/* helpers */

	P3M_DEBUG(printf("    prepare() started...\n"));

	max_comm_size = 0;
	max_grid_size = 0;

	int *n_id[4]; /* linear node identity lists for the node grids. */
	int *n_pos[4]; /* positions of nodes in the node grids. */
	for (int i = 0; i < 4; i++) {
		n_id[i] = new int[comm.size];
		n_pos[i] = new int[3 * comm.size];
	}

	int n_grid[4][3]; /* The four node grids. */
	int my_pos[4][3]; /* The position of this_node in the node grids. */
	/* === node grids === */
	/* real space node grid (n_grid[0]) */
	for (int i = 0; i < 3; i++) {
		n_grid[0][i] = comm.node_grid[i];
		my_pos[0][i] = comm.node_pos[i];
	}

	for (int i = 0; i < comm.size; i++) {
		MPI_Cart_coords(comm.mpicomm, i, 3, &(n_pos[0][3 * i]));
		int lin_ind = get_linear_index(n_pos[0][3 * i], n_pos[0][3 * i + 1],
				n_pos[0][3 * i + 2], n_grid[0]);
		n_id[0][lin_ind] = i;
	}

	/* Calc 2D FFT node grids (n_grid[1 - 3]) */
	for (int i = static_cast<p3m_int>(sqrt(comm.size)); i >= 1; i--)
		if (comm.size % i == 0) {
			n_grid[1][0] = comm.size / i;
			n_grid[1][1] = i;
			n_grid[1][2] = 1;
			break;
		}

	/* Map 3D to 2D grid */
	plan[1].row_dir = -1;

	p3m_int *g3d = n_grid[0];
	p3m_int *g2d = n_grid[1];

	/* Using 1 node in z dir */
	if (g3d[2] == 1) {
		plan[1].row_dir = 2;
	} else {
		plan[1].row_dir = -1;
		if (g2d[0] % g3d[0] == 0) {
			if (g2d[1] % g3d[1] == 0)
				plan[1].row_dir = 2;
			else if (g2d[1] % g3d[2] == 0) {
				plan[1].row_dir = 1;
				g2d[2] = g2d[1];
				g2d[1] = 1;
			}
		} else if (g2d[0] % g3d[1] == 0) {
			if (g2d[1] % g3d[0] == 0) {
				plan[1].row_dir = 2;
				std::swap(g2d[0], g2d[1]);
			} else if (g2d[1] % g3d[2] == 0) {
				plan[1].row_dir = 0;
				g2d[2] = g2d[1];
				g2d[1] = g2d[0];
				g2d[0] = 1;
			}
		} else if (g2d[0] % g3d[2] == 0) {
			if (g2d[1] % g3d[0] == 0) {
				plan[1].row_dir = 1;
				g2d[2] = g2d[0];
				g2d[0] = g2d[1];
				g2d[1] = 1;
			} else if (g2d[1] % g3d[1] == 0) {
				plan[1].row_dir = 0;
				g2d[2] = g2d[0];
				g2d[0] = 1;
			}
		}
	}

	plan[0].n_permute = 0;
	for (int i = 1; i < 4; i++)
		plan[i].n_permute = (plan[1].row_dir + i) % 3;
	for (int i = 0; i < 3; i++) {
		n_grid[2][i] = n_grid[1][(i + 1) % 3];
		n_grid[3][i] = n_grid[1][(i + 2) % 3];
	}
	plan[2].row_dir = (plan[1].row_dir - 1) % 3;
	plan[3].row_dir = (plan[1].row_dir - 2) % 3;

	/* === communication groups === */
	/* copy local grid off real space charge assignment grid */
	for (int i = 0; i < 3; i++)
		plan[0].new_grid[i] = local_grid_dim[i];
	for (int i = 1; i < 4; i++) {
		plan[i].g_size = find_comm_groups(comm, n_grid[i - 1], n_grid[i],
				n_id[i - 1], n_id[i], plan[i].group, n_pos[i], my_pos[i]);
		if (plan[i].g_size == -1) {
			/* try permutation */
			std::swap(n_grid[i][(plan[i].row_dir + 1) % 3],
					n_grid[i][(plan[i].row_dir + 2) % 3]);
			plan[i].g_size = find_comm_groups(comm, n_grid[i - 1], n_grid[i],
					n_id[i - 1], n_id[i], plan[i].group, n_pos[i], my_pos[i]);
			assert(plan[i].g_size != -1);
		}

		plan[i].send_block = static_cast<p3m_int*>(realloc(plan[i].send_block,
				6 * plan[i].g_size * sizeof(p3m_int)));
		plan[i].send_size = static_cast<p3m_int*>(realloc(plan[i].send_size,
				1 * plan[i].g_size * sizeof(p3m_int)));
		plan[i].recv_block = static_cast<p3m_int*>(realloc(plan[i].recv_block,
				6 * plan[i].g_size * sizeof(p3m_int)));
		plan[i].recv_size = static_cast<p3m_int*>(realloc(plan[i].recv_size,
				1 * plan[i].g_size * sizeof(p3m_int)));

		plan[i].new_size = calc_local_grid(my_pos[i], n_grid[i],
				global_grid_dim, global_grid_off, plan[i].new_grid,
				plan[i].start);
		permute_ifield(plan[i].new_grid, 3, -(plan[i].n_permute));
		permute_ifield(plan[i].start, 3, -(plan[i].n_permute));
		plan[i].n_ffts = plan[i].new_grid[0] * plan[i].new_grid[1];

		/* === send/recv block specifications === */
		for (int j = 0; j < plan[i].g_size; j++) {
			/* send block: this_node to comm-group-node i (identity: node) */
			p3m_int node = plan[i].group[j];
			plan[i].send_size[j] = calc_send_block(my_pos[i - 1], n_grid[i - 1],
					&(n_pos[i][3 * node]), n_grid[i], global_grid_dim,
					global_grid_off, &(plan[i].send_block[6 * j]));
			permute_ifield(&(plan[i].send_block[6 * j]), 3,
					-(plan[i - 1].n_permute));
			permute_ifield(&(plan[i].send_block[6 * j + 3]), 3,
					-(plan[i - 1].n_permute));
			if (plan[i].send_size[j] > max_comm_size)
				max_comm_size = plan[i].send_size[j];
			/* First plan send blocks have to be adjusted, since the CA grid
			 may have an additional margin outside the actual domain of the
			 node */
			if (i == 1)
				for (int k = 0; k < 3; k++)
					plan[1].send_block[6 * j + k] += local_grid_margin[2 * k];
			/* recv block: this_node from comm-group-node i (identity: node) */
			plan[i].recv_size[j] = calc_send_block(my_pos[i], n_grid[i],
					&(n_pos[i - 1][3 * node]), n_grid[i - 1], global_grid_dim,
					global_grid_off, &(plan[i].recv_block[6 * j]));
			permute_ifield(&(plan[i].recv_block[6 * j]), 3,
					-(plan[i].n_permute));
			permute_ifield(&(plan[i].recv_block[6 * j + 3]), 3,
					-(plan[i].n_permute));
			if (plan[i].recv_size[j] > max_comm_size)
				max_comm_size = plan[i].recv_size[j];
		}

		for (int j = 0; j < 3; j++)
			plan[i].old_grid[j] = plan[i - 1].new_grid[j];
		if (i == 1) {
#ifdef P3M_INTERLACE
			plan[i].element = 2;
			for(int j=0; j<plan[i].g_size; j++) {
				plan[i].send_size[j] *= 2;
				plan[i].recv_size[j] *= 2;
			}
#else
			plan[i].element = 1;
#endif
		} else {
			plan[i].element = 2;
			for (int j = 0; j < plan[i].g_size; j++) {
				plan[i].send_size[j] *= 2;
				plan[i].recv_size[j] *= 2;
			}
		}
		/* DEBUG */
		P3M_DEBUG(plan[i].print());
	}

	/* Factor 2 for complex fields */
	max_comm_size *= 2;
#ifdef P3M_INTERLACE
	/* When using interlacing we need a complex charge grid which has double size */
	max_grid_size = 2*(local_grid_dim[0]*local_grid_dim[1]*local_grid_dim[2]);
#else
	max_grid_size = (local_grid_dim[0] * local_grid_dim[1] * local_grid_dim[2]);
#endif
	for (int i = 1; i < 4; i++)
		if (2 * plan[i].new_size > max_grid_size)
			max_grid_size = 2 * plan[i].new_size;

	P3M_DEBUG(
			printf("      max_comm_size = %d, max_grid_size = %d\n",
					max_comm_size, max_grid_size));

	/* === pack function === */
	for (int i = 1; i < 4; i++) {
		plan[i].pack_function = pack_block_permute2;
		P3M_DEBUG(printf("      forward plan[%d] permute 2 \n", i));
	}
	(*ks_pnum) = 6;
	if (plan[1].row_dir == 2) {
		plan[1].pack_function = pack_block;
		P3M_DEBUG(printf("      forward plan[%d] permute 0 \n", 1));
		(*ks_pnum) = 4;
	} else if (plan[1].row_dir == 1) {
		plan[1].pack_function = pack_block_permute1;
		P3M_DEBUG(printf("      forward plan[%d] permute 1 \n", 1));
		(*ks_pnum) = 5;
	}

	/* Factor 2 for complex numbers */
	send_buf = (p3m_float *) realloc(send_buf,
			max_comm_size * sizeof(p3m_float));
	recv_buf = (p3m_float *) realloc(recv_buf,
			max_comm_size * sizeof(p3m_float));

	if (data_buf != NULL)
		free_data(data_buf);
	data_buf = malloc_data();

	if (!data_buf || !recv_buf || !send_buf)
		throw std::logic_error("Could not allocate FFT data arrays");

	fftw_complex *c_data_buf = reinterpret_cast<fftw_complex*>(data_buf);

	/* FFTW WISDOM stuff. */
	/* @todo: Planning shouldn't write to file. */
	/* === FFT Routines (Using FFTW / RFFTW package)=== */
	for (int i = 1; i < 4; i++) {
		plan[i].dir = FFTW_FORWARD;
		/* FFT plan creation.
		 Attention: destroys contents of c_data/data and c_data_buf/data_buf. */
		int wisdom_status = FFTW_FAILURE;
		char wisdom_file_name[255];
		sprintf(wisdom_file_name, "fftw3_1d_wisdom_forw_n%d.file",
				plan[i].new_grid[2]);
		FILE* wisdom_file = fopen(wisdom_file_name, "r");
		if (wisdom_file != NULL) {
			wisdom_status = fftw_import_wisdom_from_file(wisdom_file);
			fclose(wisdom_file);
		}
		if (is_prepared)
			fftw_destroy_plan(plan[i].plan);
		//printf("plan[%d].n_ffts=%d\n",i,plan[i].n_ffts);
		plan[i].plan =
		fftw_plan_many_dft(1, &plan[i].new_grid[2], plan[i].n_ffts, c_data_buf,
				NULL, 1, plan[i].new_grid[2], c_data_buf, NULL, 1,
				plan[i].new_grid[2], plan[i].dir, FFTW_PATIENT);

		if (wisdom_status == FFTW_FAILURE) {
			FILE* wisdom_file = fopen(wisdom_file_name, "w");
			if (wisdom_file != NULL) {
				fftw_export_wisdom_to_file(wisdom_file);
				fclose(wisdom_file);
			}
		}
	}

	/* === The BACK Direction === */
	/* this is needed because slightly different functions are used */
	for (int i = 1; i < 4; i++) {
		back[i].dir = FFTW_BACKWARD;
		int wisdom_status = FFTW_FAILURE;
		char wisdom_file_name[255];
		sprintf(wisdom_file_name, "fftw3_1d_wisdom_back_n%d.file",
				plan[i].new_grid[2]);
		FILE* wisdom_file = fopen(wisdom_file_name, "r");
		if (wisdom_file != NULL) {
			wisdom_status = fftw_import_wisdom_from_file(wisdom_file);
			fclose(wisdom_file);
		}
		if (is_prepared)
			fftw_destroy_plan(back[i].plan);
		back[i].plan =
		fftw_plan_many_dft(1, &plan[i].new_grid[2], plan[i].n_ffts, c_data_buf,
				NULL, 1, plan[i].new_grid[2], c_data_buf, NULL, 1,
				plan[i].new_grid[2], back[i].dir, FFTW_PATIENT);
		if (wisdom_status == FFTW_FAILURE) {
			FILE* wisdom_file = fopen(wisdom_file_name, "w");
			if (wisdom_file != NULL) {
				fftw_export_wisdom_to_file(wisdom_file);
				fclose(wisdom_file);
			}
		}
		back[i].pack_function = pack_block_permute1;
		P3M_DEBUG(printf("      back plan[%d] permute 1 \n", i));
	}
	if (plan[1].row_dir == 2) {
		back[1].pack_function = pack_block;
		P3M_DEBUG(printf("      back plan[%d] permute 0 \n", 1));
	} else if (plan[1].row_dir == 1) {
		back[1].pack_function = pack_block_permute2;
		P3M_DEBUG(printf("      back plan[%d] permute 2 \n", 1));
	}
	is_prepared = true;

	for (int i = 0; i < 4; i++) {
		delete[] n_id[i];
		delete[] n_pos[i];
	}

	P3M_DEBUG(printf("    prepare() finished.\n"));
} /* end of PREPARE */

void Parallel3DFFT::forward(p3m_float *data) {
	p3m_int i;

	/* int m,n,o; */
	/* ===== first direction  ===== */
	P3M_DEBUG(printf("    %d: fft_perform_forward: dir 1:\n", comm.rank));

	fftw_complex *c_data = (fftw_complex *) data;
	fftw_complex *c_data_buf = (fftw_complex *) data_buf;

	/* communication to current dir row format (in is data) */
	forward_grid_comm(plan[1], data, data_buf);

	/*
	 printf("%d: start grid \n",comm.rank);
	 i=0;
	 for(m=0;m<8;m++) {
	 for(n=0;n<8;n++) {
	 for(o=0;o<8;o++) {
	 printf("%.3f ",data_buf[i++]);
	 }
	 printf("\n");
	 }
	 printf("\n");
	 }
	 */

	/* complexify the real data array (in is data_buf) */
#ifndef P3M_INTERLACE
	for (i = 0; i < plan[1].new_size; i++) {
		data[2 * i] = data_buf[i]; /* real value */
		data[(2 * i) + 1] = 0; /* complex value */
	}
#else
	for(i=0;i<(2*plan[1].new_size);i++)
	data[i] = data_buf[i]; /* real value */
#endif
	/* perform FFT (in/out is data)*/
	fftw_execute_dft(plan[1].plan, c_data, c_data);

	/* ===== second direction ===== */
	P3M_DEBUG_LOCAL(printf("    %d: fft_perform_forward: dir 2\n", comm.rank));
	/* communication to current dir row format (in is data) */
	forward_grid_comm(plan[2], data, data_buf);
	/* perform FFT (in/out is data_buf)*/
	fftw_execute_dft(plan[2].plan, c_data_buf, c_data_buf);
	/* ===== third direction  ===== */
	P3M_DEBUG_LOCAL(printf("    %d: fft_perform_forward: dir 3\n", comm.rank));
	/* communication to current dir row format (in is data_buf) */
	forward_grid_comm(plan[3], data_buf, data);
	/* perform FFT (in/out is data)*/
	fftw_execute_dft(plan[3].plan, c_data, c_data);
	//P3M_DEBUG_LOCAL(print_global_grid(comm, plan[3], data, 1, 0));

	// REMARK: Result has to be in data.
}

void Parallel3DFFT::backward(p3m_float *data) {
	p3m_int i;

	fftw_complex *c_data = (fftw_complex *) data;
	fftw_complex *c_data_buf = (fftw_complex *) data_buf;

	/* ===== third direction  ===== */
	P3M_DEBUG(printf("    %d: backward: dir 3\n", comm.rank));

	/* perform FFT (in is data) */
	fftw_execute_dft(back[3].plan, c_data, c_data);
	/* communicate (in is data)*/
	backward_grid_comm(plan[3], back[3], data, data_buf);

	/* ===== second direction ===== */
	P3M_DEBUG_LOCAL(printf("    %d: back: dir 2\n", comm.rank));
	/* perform FFT (in is data_buf) */
	fftw_execute_dft(back[2].plan, c_data_buf, c_data_buf);
	/* communicate (in is data_buf) */
	backward_grid_comm(plan[2], back[2], data_buf, data);

	/* ===== first direction  ===== */
	P3M_DEBUG_LOCAL(printf("    %d: backward: dir 1\n", comm.rank));
	/* perform FFT (in is data) */
	fftw_execute_dft(back[1].plan, c_data, c_data);
#ifndef P3M_INTERLACE
	/* throw away the (hopefully) empty complex component (in is data)*/
	for (i = 0; i < plan[1].new_size; i++)
		data_buf[i] = data[2 * i]; /* real value */
#ifdef ADDITIONAL_CHECKS
	for (i = 0; i < plan[1].new_size; i++)
		if (data[2 * i + 1] > 1e-5) {
			printf("    %d: Imaginary value is not zero (i=%d,data=%g)!!!\n",
					comm.rank, i, data[2 * i + 1]);
			if (i > 100)
				exit(-1);
		}
#endif
#else
	/* keep imaginary part */
	for (i=0; i<(2*plan[1].new_size); i++)
	data_buf[i] = data[i];
#endif
	/* communicate (in is data_buf) */
	backward_grid_comm(plan[1], back[1], data_buf, data);

	/* REMARK: Result has to be in data. */
}

/** communicate the grid data according to the given forward_plan.
 * \param plan communication plan (see \ref forward_plan).
 * \param in   input grid.
 * \param out  output grid.
 */
void Parallel3DFFT::forward_grid_comm(forward_plan plan,
        p3m_float *in, p3m_float *out) {
	for (int i = 0; i < plan.g_size; i++) {
		plan.pack_function(in, send_buf, &(plan.send_block[6 * i]),
				&(plan.send_block[6 * i + 3]), plan.old_grid, plan.element);

		if (plan.group[i] == comm.rank)
		    /* Self communication... */
		    std::swap(send_buf, recv_buf);
		else {
		    MPI_Sendrecv(send_buf, plan.send_size[i], P3M_MPI_FLOAT, plan.group[i],
		            REQ_FFT_FORW, recv_buf, plan.recv_size[i], P3M_MPI_FLOAT, plan.group[i],
		            REQ_FFT_FORW, comm.mpicomm, MPI_STATUS_IGNORE);
		}

		unpack_block(recv_buf, out, &(plan.recv_block[6 * i]),
		        &(plan.recv_block[6 * i + 3]), plan.new_grid, plan.element);
	}
}

/** communicate the grid data according to the given forward_plan/fft_bakc_plan.
 * \param plan_f communication plan (see \ref forward_plan).
 * \param plan_b additional backward plan (see \ref fft.backward_plan).
 * \param in     input grid.
 * \param out    output grid.
 */
void Parallel3DFFT::backward_grid_comm(forward_plan plan_f,
		backward_plan plan_b, p3m_float *in, p3m_float *out) {

	/* Back means: Use the send/recieve stuff from the forward plan but
	 replace the recieve blocks by the send blocks and vice
	 versa. Attention then also new_grid and old_grid are exchanged */

	for (int i = 0; i < plan_f.g_size; i++) {
		plan_b.pack_function(in, send_buf, &(plan_f.recv_block[6 * i]),
				&(plan_f.recv_block[6 * i + 3]), plan_f.new_grid,
				plan_f.element);

        if (plan_f.group[i] == comm.rank)
            /* Self communication... */
            std::swap(send_buf, recv_buf);
        else {
            MPI_Sendrecv(send_buf, plan_f.recv_size[i], P3M_MPI_FLOAT, plan_f.group[i],
                    REQ_FFT_BACK, recv_buf, plan_f.send_size[i], P3M_MPI_FLOAT, plan_f.group[i],
                    REQ_FFT_BACK, comm.mpicomm, MPI_STATUS_IGNORE);
        }

        unpack_block(recv_buf, out, &(plan_f.send_block[6 * i]),
                &(plan_f.send_block[6 * i + 3]), plan_f.old_grid, plan_f.element);
	}
}

/** This ugly function does the bookkepping which nodes have to
 *  communicate to each other, when you change the node grid.
 *  Changing the domain decomposition requieres communication. This
 *  function finds (hopefully) the best way to do this. As input it
 *  needs the two grids (grid1, grid2) and a linear list (node_list1)
 *  with the node identities for grid1. The linear list (node_list2)
 *  for the second grid is calculated. For the communication group of
 *  the calling node it calculates a list (group) with the node
 *  identities and the positions (pos1, pos2) of that nodes in grid1
 *  and grid2. The return value is the size of the communication
 *  group. It gives -1 if the two grids do not fit to each other
 *  (grid1 and grid2 have to be component wise multiples of each
 *  other. see e.g. \ref calc_2d_grid in \ref grid.c for how to do
 *  this.).
 *
 * \param grid1       The node grid you start with (Input).
 * \param grid2       The node grid you want to have (Input).
 * \param node_list1  Linear node index list for grid1 (Input).
 * \param node_list2  Linear node index list for grid2 (Output).
 * \param group       communication group (node identity list) for the calling node  (Output).
 * \param pos        positions of the nodes in in grid2 (Output).
 * \param my_pos      position of this_node in  grid2.
 * \return Size of the communication group (Output of course!).  */
p3m_int find_comm_groups(Communication &comm, p3m_int grid1[3], p3m_int grid2[3],
		p3m_int *node_list1, p3m_int *node_list2, p3m_int *group, p3m_int *pos,
		p3m_int *my_pos) {
	p3m_int i;
	/* communication group cell size on grid1 and grid2 */
	p3m_int s1[3], s2[3];
	/* The communication group cells build the same super grid on grid1 and grid2 */
	p3m_int ds[3];
	/* communication group size */
	p3m_int g_size = 1;
	/* comm. group cell index */
	p3m_int gi[3];
	/* position of a node in a grid */
	p3m_int p1[3], p2[3];
	/* node identity */
	p3m_int n;
	/* this_node position in the communication group. */
	p3m_int c_pos = -1;
	/* flag for group identification */
	p3m_int my_group = 0;

	P3M_DEBUG(
			printf(
					"     fft_find_comm_groups(): for grid1=(%d,%d,%d) and grids=(%d,%d,%d)\n",
					grid1[0], grid1[1], grid1[2], grid2[0], grid2[1],
					grid2[2]));

	/* calculate dimension of comm. group cells for both grids */
	if ((grid1[0] * grid1[1] * grid1[2]) != (grid2[0] * grid2[1] * grid2[2]))
		return -1; /* unlike number of nodes */
	for (i = 0; i < 3; i++) {
		s1[i] = grid1[i] / grid2[i];
		if (s1[i] == 0)
			s1[i] = 1;
		else if (grid1[i] != grid2[i] * s1[i])
			return -1; /* grids do not match!!! */

		s2[i] = grid2[i] / grid1[i];
		if (s2[i] == 0)
			s2[i] = 1;
		else if (grid2[i] != grid1[i] * s2[i])
			return -1; /* grids do not match!!! */

		ds[i] = grid2[i] / s2[i];
		g_size *= s2[i];
	}

	/* calc node_list2 */
	/* loop over all comm. group cells */
	for (gi[2] = 0; gi[2] < ds[2]; gi[2]++)
		for (gi[1] = 0; gi[1] < ds[1]; gi[1]++)
			for (gi[0] = 0; gi[0] < ds[0]; gi[0]++) {
				/* loop over all nodes in that comm. group cell */
				for (i = 0; i < g_size; i++) {
					p1[0] = (gi[0] * s1[0]) + (i % s1[0]);
					p1[1] = (gi[1] * s1[1]) + ((i / s1[0]) % s1[1]);
					p1[2] = (gi[2] * s1[2]) + (i / (s1[0] * s1[1]));

					p2[0] = (gi[0] * s2[0]) + (i % s2[0]);
					p2[1] = (gi[1] * s2[1]) + ((i / s2[0]) % s2[1]);
					p2[2] = (gi[2] * s2[2]) + (i / (s2[0] * s2[1]));

					n =
							node_list1[get_linear_index(p1[0], p1[1], p1[2],
									grid1)];
					node_list2[get_linear_index(p2[0], p2[1], p2[2], grid2)] =
							n;

					pos[3 * n + 0] = p2[0];
					pos[3 * n + 1] = p2[1];
					pos[3 * n + 2] = p2[2];
					if (my_group == 1)
						group[i] = n;
					if (n == comm.rank && my_group == 0) {
						my_group = 1;
						c_pos = i;
						my_pos[0] = p2[0];
						my_pos[1] = p2[1];
						my_pos[2] = p2[2];
						i = -1; /* restart the loop */
					}
				}
				my_group = 0;
			}

	/* permute comm. group according to the nodes position in the group */
	/* This is necessary to have matching node pairs during communication! */
	while (c_pos > 0) {
		n = group[g_size - 1];
		for (i = g_size - 1; i > 0; i--)
			group[i] = group[i - 1];
		group[0] = n;
		c_pos--;
	}
	return g_size;
}

/** Calculate the local fft grid.  Calculate the local grid (loc_grid)
 *  of a node at position (n_pos) in a node grid (n_grid) for a global
 *  grid of size (grid) and a grid offset (grid_off (in grid units))
 *  and store also the first point (start) of the local grid.
 *
 * \return size     number of grid points in local grid.
 * \param  n_pos    Position of the node in n_grid.
 * \param  n_grid   node grid.
 * \param  grid     global grid dimensions.
 * \param  grid_off global grid offset (see \ref p3m_struct).
 * \param  loc_grid local grid dimension (output).
 * \param  start    first point of local grid in global grid (output).
 */
p3m_int calc_local_grid(p3m_int n_pos[3], p3m_int n_grid[3], p3m_int grid[3],
		p3m_float grid_off[3], p3m_int loc_grid[3], p3m_int start[3]) {
	p3m_int last[3], size = 1;

	for (int i = 0; i < 3; i++) {
		start[i] = (p3m_int) ceil(
				(grid[i] / (p3m_float) n_grid[i]) * n_pos[i] - grid_off[i]);
		last[i] = (p3m_int) floor(
				(grid[i] / (p3m_float) n_grid[i]) * (n_pos[i] + 1)
						- grid_off[i]);
		/* correct round-off errors */
		if ((grid[i] / (p3m_float) n_grid[i]) * (n_pos[i] + 1) - grid_off[i]
				- last[i] < 1.0e-15)
			last[i]--;
		if (1.0 + (grid[i] / (p3m_float) n_grid[i]) * n_pos[i] - grid_off[i]
				- start[i] < 1.0e-15)
			start[i]--;
		loc_grid[i] = last[i] - start[i] + 1;
		size *= loc_grid[i];
	}
	return size;
}

/** Calculate a send (or recv.) block for grid communication during a
 *  decomposition change.  Calculate the send block specification
 *  (block = lower left corner and upper right corner) which a node at
 *  position (pos1) in the actual node grid (grid1) has to send to
 *  another node at position (pos2) in the desired node grid
 *  (grid2). The global grid, subject to communication, is specified
 *  via its size (grid) and its grid offset (grid_off (in grid
 *  units)).
 *
 *  For the calculation of a receive block you have to change the arguments in the following way: <br>
 *  pos1  - position of receiving node in the desired node grid. <br>
 *  grid1 - desired node grid. <br>
 *  pos2  - position of the node you intend to receive the data from in the actual node grid. <br>
 *  grid2 - actual node grid.  <br>
 2 *
 *  \return          size of the send block.
 *  \param  pos1     Position of send node in grid1.
 *  \param  node_grid1    node grid 1.
 *  \param  pos2     Position of recv node in grid2.
 *  \param  node_grid2    node grid 2.
 *  \param  grid     global grid dimensions.
 *  \param  grid_off global grid offset (see \ref p3m_struct).
 *  \param  block    send block specification.
 */
p3m_int calc_send_block(p3m_int pos1[3], p3m_int node_grid1[3], p3m_int pos2[3],
		p3m_int node_grid2[3], p3m_int grid[3], p3m_float grid_off[3],
		p3m_int block[6]) {
	p3m_int i, size = 1;
	p3m_int grid1[3], first1[3], last1[3];
	p3m_int grid2[3], first2[3], last2[3];

	calc_local_grid(pos1, node_grid1, grid, grid_off, grid1, first1);
	calc_local_grid(pos2, node_grid2, grid, grid_off, grid2, first2);

	for (i = 0; i < 3; i++) {
		last1[i] = first1[i] + grid1[i] - 1;
		last2[i] = first2[i] + grid2[i] - 1;
		block[i] = imax(first1[i], first2[i]) - first1[i];
		block[i + 3] = (imin(last1[i], last2[i]) - first1[i]) - block[i] + 1;
		size *= block[i + 3];
	}
	return size;
}

/** pack a block (size[3] starting at start[3]) of an input 3d-grid
 *  with dimension dim[3] into an output 3d-block with dimension size[3].
 *
 *    The block with dimensions (size[0], size[1], size[2]) is stored
 *    in 'row-major-order' or 'C-order', that means the first index is
 *    changing slowest when running through the linear array. The
 *    element (i0 (slow), i1 (mid), i2 (fast)) has the linear index
 *    li = i2 + size[2] * (i1 + (size[1]*i0))
 *
 *  \param in      pointer to input 3d-grid.
 *  \param out     pointer to output 3d-grid (block).
 *  \param start   start index of the block in the in-grid.
 *  \param size    size of the block (=dimension of the out-grid).
 *  \param dim     size of the in-grid.
 *  \param element size of a grid element (e.g. 1 for Real, 2 for Complex).
 */
void Parallel3DFFT::pack_block(p3m_float *in, p3m_float *out, p3m_int start[3],
		p3m_int size[3], p3m_int dim[3], p3m_int element) {
	/* mid and slow changing indices */
	p3m_int m, s;
	/* linear index of in grid, linear index of out grid */
	p3m_int li_in, li_out = 0;
	/* copy size */
	p3m_int copy_size;
	/* offsets for indizes in input grid */
	p3m_int m_in_offset, s_in_offset;
	/* offsets for indizes in output grid */
	p3m_int m_out_offset;

	copy_size = element * size[2] * sizeof(p3m_float);
	m_in_offset = element * dim[2];
	s_in_offset = element * (dim[2] * (dim[1] - size[1]));
	m_out_offset = element * size[2];
	li_in = element * (start[2] + dim[2] * (start[1] + dim[1] * start[0]));

	for (s = 0; s < size[0]; s++) {
		for (m = 0; m < size[1]; m++) {
			memcpy(&(out[li_out]), &(in[li_in]), copy_size);
			li_in += m_in_offset;
			li_out += m_out_offset;
		}
		li_in += s_in_offset;
	}
}

/** pack a block with dimensions (size[0] * size[1] * aize[2]) starting
 *  at start[3] of an input 3d-grid with dimension dim[3] into an
 *  output 3d-grid with dimensions (size[2] * size[0] * size[1]) with
 *  a simulatanous one-fold permutation of the indices.
 *
 * The permutation is defined as:
 * slow_in -> fast_out, mid_in ->slow_out, fast_in -> mid_out
 *
 * An element (i0_in , i1_in , i2_in ) is then
 * (i0_out = i1_in-start[1], i1_out = i2_in-start[2], i2_out = i0_in-start[0]) and
 * for the linear indices we have:                              <br>
 * li_in = i2_in + size[2] * (i1_in + (size[1]*i0_in))          <br>
 * li_out = i2_out + size[0] * (i1_out + (size[2]*i0_out))
 *
 * For index definition see \ref fft_pack_block.
 *
 *  \param in      pointer to input 3d-grid.
 *  \param out     pointer to output 3d-grid (block).
 *  \param start   start index of the block in the in-grid.
 *  \param size    size of the block (=dimension of the out-grid).
 *  \param dim     size of the in-grid.
 *  \param element size of a grid element (e.g. 1 for Real, 2 for Complex).
 */
static void pack_block_permute1(p3m_float *in, p3m_float *out, p3m_int start[3],
		p3m_int size[3], p3m_int dim[3], p3m_int element) {
	/* slow,mid and fast changing indices for input  grid */
	p3m_int s, m, f, e;
	/* linear index of in grid, linear index of out grid */
	p3m_int li_in, li_out = 0;
	/* offsets for indizes in input grid */
	p3m_int m_in_offset, s_in_offset;
	/* offset for mid changing indices of output grid */
	p3m_int m_out_offset;

	m_in_offset = element * (dim[2] - size[2]);
	s_in_offset = element * (dim[2] * (dim[1] - size[1]));
	m_out_offset = (element * size[0]) - element;
	li_in = element * (start[2] + dim[2] * (start[1] + dim[1] * start[0]));

	for (s = 0; s < size[0]; s++) { /* fast changing out */
		li_out = element * s;
		for (m = 0; m < size[1]; m++) { /* slow changing out */
			for (f = 0; f < size[2]; f++) { /* mid  changing out */
				for (e = 0; e < element; e++)
					out[li_out++] = in[li_in++];
				li_out += m_out_offset;
			}
			li_in += m_in_offset;
		}
		li_in += s_in_offset;
	}
}

/** pack a block with dimensions (size[0] * size[1] * aize[2]) starting
 *  at start[3] of an input 3d-grid with dimension dim[3] into an
 *  output 3d-grid with dimensions (size[2] * size[0] * size[1]), this
 *  is a simulatanous two-fold permutation of the indices.
 *
 * The permutation is defined as:
 * slow_in -> mid_out, mid_in ->fast_out, fast_in -> slow_out
 *
 * An element (i0_in , i1_in , i2_in ) is then
 * (i0_out = i2_in-start[2], i1_out = i0_in-start[0], i2_out = i1_in-start[1]) and
 * for the linear indices we have:                              <br>
 * li_in = i2_in + size[2] * (i1_in + (size[1]*i0_in))          <br>
 * li_out = i2_out + size[0] * (i1_out + (size[2]*i0_out))
 *
 * For index definition see \ref fft_pack_block.
 *
 *  \param in      pointer to input 3d-grid.
 *  \param out     pointer to output 3d-grid (block).
 *  \param start   start index of the block in the in-grid.
 *  \param size    size of the block (=dimension of the out-grid).
 *  \param dim     size of the in-grid.
 *  \param element size of a grid element (e.g. 1 for Real, 2 for Complex).
 */
static void pack_block_permute2(p3m_float *in, p3m_float *out, p3m_int start[3],
		p3m_int size[3], p3m_int dim[3], p3m_int element) {
	/* slow,mid and fast changing indices for input  grid */
	p3m_int s, m, f, e;
	/* linear index of in grid, linear index of out grid */
	p3m_int li_in, li_out = 0;
	/* offsets for indizes in input grid */
	p3m_int m_in_offset, s_in_offset;
	/* offset for slow changing index of output grid */
	p3m_int s_out_offset;
	/* start index for mid changing index of output grid */
	p3m_int m_out_start;

	m_in_offset = element * (dim[2] - size[2]);
	s_in_offset = element * (dim[2] * (dim[1] - size[1]));
	s_out_offset = (element * size[0] * size[1]) - element;
	li_in = element * (start[2] + dim[2] * (start[1] + dim[1] * start[0]));

	for (s = 0; s < size[0]; s++) { /* mid changing out */
		m_out_start = element * (s * size[1]);
		for (m = 0; m < size[1]; m++) { /* fast changing out */
			li_out = m_out_start + element * m;
			for (f = 0; f < size[2]; f++) { /* slow  changing out */
				for (e = 0; e < element; e++)
					out[li_out++] = in[li_in++];
				li_out += s_out_offset;
			}
			li_in += m_in_offset;
		}
		li_in += s_in_offset;
	}

}

//fixme: unpack_block and add_block are the same! Unify.
void Parallel3DFFT::unpack_block(p3m_float *in, p3m_float *out,
		p3m_int start[3], p3m_int size[3], p3m_int dim[3], p3m_int element) {

	/* copy size */
	p3m_int copy_size = element * (size[2] * sizeof(p3m_float));
    /* offsets for indices in output grid */
	p3m_int m_out_offset = element * dim[2];
	p3m_int s_out_offset = element * (dim[2] * (dim[1] - size[1]));
    /* offset for indices in input grid */
	p3m_int m_in_offset = element * size[2];

    /* linear index of in grid, linear index of out grid */
    p3m_int li_in = 0;
	p3m_int li_out = element * (start[2] + dim[2] * (start[1] + dim[1] * start[0]));

	for (int s = 0; s < size[0]; s++) {
		for (int m = 0; m < size[1]; m++) {
			memcpy(&(out[li_out]), &(in[li_in]), copy_size);
			li_in += m_in_offset;
			li_out += m_out_offset;
		}
		li_out += s_out_offset;
	}
}

void Parallel3DFFT::add_block(p3m_float *in, p3m_float *out,
        int start[3], int size[3], int dim[3]) {

    /* linear index of in grid, linear index of out grid */
    p3m_int li_in=0;
    p3m_int li_out = start[2] + ( dim[2]*( start[1] + (dim[1]*start[0]) ) );
    /* offsets for indices in output grid */
    p3m_int m_out_offset  = dim[2] - size[2];
    p3m_int s_out_offset  = (dim[2] * (dim[1] - size[1]));

    for(int s=0 ;s<size[0]; s++) {
        for(int m=0; m<size[1]; m++) {
            for(int f=0; f<size[2]; f++) {
                out[li_out++] += in[li_in++];
            }
            li_out += m_out_offset;
        }
        li_out += s_out_offset;
    }
}

/** Debug function to print global fft grid.
 Print a globaly distributed grid contained in data. Element size is element.
 * \param plan     fft/communication plan (see \ref forward_plan).
 * \param data     grid data.
 * \param element  element size.
 * \param num      element index to print.
 */

void Parallel3DFFT::print_global_grid(Communication &comm,
		Parallel3DFFT::forward_plan plan, p3m_float *data, p3m_int element,
		p3m_int num) {
	p3m_int st[3], en[3], si[3];
	for (int i = 0; i < 3; i++) {
		st[i] = plan.start[i];
		en[i] = plan.start[i] + plan.new_grid[i];
		si[i] = plan.new_grid[i];
	}

	int grid = plan.new_grid[2];
	MPI_Barrier(comm.mpicomm);

	if (comm.rank == 0)
		printf("    global grid: (%d of %d elements)\n", num + 1, element);
	MPI_Barrier(comm.mpicomm);

	for (int i = 0; i < comm.size; i++) {
		MPI_Barrier(comm.mpicomm);
		if (i == comm.rank)
			printf("    %d: range (%d,%d,%d)-(%d,%d,%d)\n", comm.rank, st[0],
					st[1], st[2], en[0], en[1], en[2]);
	}
	MPI_Barrier(comm.mpicomm);

	int divide = 0;
	int block1 = -1;
	for (int b = 1; divide == 0; b++)
		if (b * grid > 7) {
			block1 = b;
			divide = (p3m_int) ceil(static_cast<p3m_float>(grid) / block1);
		}

	for (int b = 0; b < divide; b++) {
		int start1 = b * block1;
		bool my;
		for (int i0 = grid - 1; i0 >= 0; i0--) {
			for (int i1 = start1; i1 < imin(start1 + block1, grid); i1++) {
				for (int i2 = 0; i2 < grid; i2++) {
					my = (i0 >= st[0] && i0 < en[0] && i1 >= st[1] && i1 < en[1]
							&& i2 >= st[2] && i2 < en[2]);

					MPI_Barrier(comm.mpicomm);
					if (my) {
						p3m_float tmp =
								data[num
										+ (element
												* ((i2 - st[2])
														+ si[2]
																* ((i1 - st[1])
																		+ si[1]
																				* (i0
																						- st[0]))))];
						if (fabs(tmp) > 1.0e-15) {
							if (tmp < 0)
								printf("    %1.2" P3M_LMOD_FLOAT "e", tmp);
							else
								printf("    %1.2" P3M_LMOD_FLOAT "e", tmp);
						} else {
							printf(" %1.2e", 0.0);
						}
					}
					MPI_Barrier(comm.mpicomm);
				}
				if (my)
					printf(" | ");
			}
			if (my)
				printf("\n");
		}
		if (my)
			printf("\n");
	}
}
}

