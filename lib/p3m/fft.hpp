/*
  Copyright (C) 2013 Olaf Lenz
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
/** \file fft.h
 *
 *  Routines, row decomposition, data structures and communication for the 3D-FFT. 
 *
 *  The 3D-FFT is split into 3 ond dimensional FFTs. The data is
 *  distributed in such a way, that for the actual direction of the
 *  FFT each node has a certain number of rows for which it performs a
 *  1D-FFT. After performing the FFT on theat direction the data is
 *  redistributed.
 *
 *  For simplicity at the moment I have implemented a full complex to
 *  complex FFT (even though a real to complex FFT would be
 *  sufficient)
 *
 */
#ifndef _P3M_FFT_HPP
#define _P3M_FFT_HPP

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "communication.hpp"
#include "fftw3.h"

/***************************************************/
/* DATA TYPES */
/***************************************************/

#ifndef FCS_USE_COMMON_FFTW
typedef fftw_plan fcs_fftw_plan;
#endif

/** Structure for performing a 1D FFT.  
 *
 *  This includes the information about the redistribution of the 3D
 *  FFT *grid before the actual FFT.  
*/
typedef struct {
  /** plan direction: 0 = Forward FFT, 1 = Backward FFT. */
  fcs_int dir;
  /** row direction of that FFT. */
  fcs_int row_dir;
  /** permutations from normal coordinate system. */
  fcs_int n_permute;
  /** number of 1D FFTs. */ 
  fcs_int n_ffts;
  /** plan for fft. */
  fcs_fftw_plan fftw_plan;
  /** function for fft. */
  void (*fft_function)(fcs_fftw_plan);

  /** size of local grid before communication. */
  fcs_int old_grid[3];
  /** size of local grid after communication, also used for actual FFT. */
  fcs_int new_grid[3];
  /** lower left point of local FFT grid in global FFT grid coordinates. */
  fcs_int start[3];
  /** size of new grid (number of grid points). */
  fcs_int new_size;

  /** number of nodes which have to communicate with each other. */ 
  fcs_int g_size;
  /** group of nodes which have to communicate with each other. */ 
  fcs_int *group;

  /** packing function for send blocks. */
  void (*pack_function)(double*, double*, int*, int*, int*, int);
  /** Send block specification. 6 integers for each node: start[3], size[3]. */ 
  fcs_int *send_block;
  /** Send block communication sizes. */ 
  fcs_int *send_size;
  /** Recv block specification. 6 integers for each node: start[3], size[3]. */ 
  fcs_int *recv_block;
  /** Recv block communication sizes. */ 
  fcs_int *recv_size;
  /** size of send block elements. */
  fcs_int element;
} ifcs_fft_forw_plan;

/** Additional information for backwards FFT.*/
typedef struct {
  /** plan direction. (e.g. fftw makro)*/
  fcs_int dir;
  /** plan for fft. */
  fcs_fftw_plan fftw_plan;
  /** function for fft. */
  void (*fft_function)(fcs_fftw_plan);

  /** packing function for send blocks. */
  void (*pack_function)(double*, double*, int*, int*, int*, int); 
} ifcs_fft_back_plan;

typedef struct {
  /** Information about the three one dimensional FFTs and how the nodes
   *  have to communicate in between.
   *
   * NOTE: FFT numbering starts with 1 for technical reasons (because we
   *       have 4 node grids, the index 0 is used for the real space
   *       charge assignment grid).  */
  ifcs_fft_forw_plan plan[4];
  /** Information for Back FFTs (see fft.plan). */
  ifcs_fft_back_plan back[4];

  /** Whether FFT is initialized or not. */
  fcs_int init_tag;

  /** Maximal size of the communication buffers. */
  fcs_int max_comm_size;

  /** Maximal local grid size. */
  fcs_int max_grid_size;

  /** send buffer. */
  fcs_float *send_buf;
  /** receive buffer. */
  fcs_float *recv_buf;
  /** Buffer for receive data. */
  fcs_float *data_buf;
} ifcs_fft_data_struct;

/***************************************************/
/* EXPORTED FUNCTIONS */
/***************************************************/
/** Initialize fft data structure. */
void 
ifcs_fft_init(ifcs_fft_data_struct *fft, 
	      ifcs_p3m_comm_struct *comm);

void 
ifcs_fft_destroy(ifcs_fft_data_struct *fft, fcs_float *data, 
		 fcs_float *ks_data);

/** Prepare the 3D-FFT for a given size.

 * \return Maximal size of local fft grid (needed for allocation of ca_grid).
 * \param data              Pointer Pointer to data array containing the rs grid.
 * \param ks_data           Pointer Pointer to data array containing the ks grid.
 * \param local_grid_dim    Pointer to local CA grid dimensions.
 * \param local_grid_margin Pointer to local CA grid margins.
 * \param global_grid_dim   Pointer to global CA grid dimensions.
 * \param global_grid_off   Pointer to global CA grid margins.
 * \param ks_pnum           Pointer to number of permutations in k-space.
 */
fcs_int 
ifcs_fft_prepare(ifcs_fft_data_struct *fft, ifcs_p3m_comm_struct *comm,
		 fcs_float **data, fcs_float **ks_data,  
		 fcs_int *local_grid_dim, fcs_int *local_grid_margin, 
		 fcs_int* global_grid_dim, fcs_float *global_grid_off, 
		 fcs_int *ks_pnum);

/** perform the forward 3D FFT.
    The assigned charges are in \a data. The result is also stored in \a data.
    \warning The content of \a data is overwritten.
    \param data Grid.
*/
void 
ifcs_fft_perform_forw(ifcs_fft_data_struct *fft, 
		     ifcs_p3m_comm_struct *comm,
		     fcs_float *data);

/** perform the backward 3D FFT.
    \warning The content of \a data is overwritten.
    \param data Grid.
*/
void 
ifcs_fft_perform_back(ifcs_fft_data_struct *fft, 
		      ifcs_p3m_comm_struct *comm,
		      fcs_float *data);

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
void 
ifcs_fft_pack_block(fcs_float *in, fcs_float *out, 
		    fcs_int start[3], fcs_int size[3], 
		    fcs_int dim[3], fcs_int element);

/** unpack a 3d-grid input block (size[3]) into an output 3d-grid
 *  with dimension dim[3] at start position start[3].
 *
 *  see also \ref fcs_fft_pack_block.
 *
 *  \param in      pointer to input 3d-grid.
 *  \param out     pointer to output 3d-grid (block).
 *  \param start   start index of the block in the in-grid.
 *  \param size    size of the block (=dimension of the out-grid).
 *  \param dim     size of the in-grid.
 *  \param element size of a grid element (e.g. 1 for Real, 2 for Complex).
 */
void 
ifcs_fft_unpack_block(fcs_float *in, fcs_float *out, 
		      fcs_int start[3], fcs_int size[3],
		      fcs_int dim[3], fcs_int element);

/* MPI tags for the fft communications: */
/** Tag for communication in forw_grid_comm() */
#define REQ_FFT_FORW   201
/** Tag for communication in back_grid_comm() */
#define REQ_FFT_BACK   202
/* Tag for wisdom file I/O */
#  define FFTW_FAILURE 0

#endif
