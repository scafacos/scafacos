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

#include <config.h>
#include "p3mconfig.hpp"
#include "utils.hpp"
#include "Communication.hpp"

namespace P3M {
/***************************************************/
/* DATA TYPES */
/***************************************************/

/* MPI tags for the fft communications: */
/** Tag for communication in forward_grid_comm() */
#define REQ_FFT_FORW   201
/** Tag for communication in backward_grid_comm() */
#define REQ_FFT_BACK   202
/* Tag for wisdom file I/O */
#define FFTW_FAILURE 0

/** Structure for performing a 1D FFT.
 *
 *  This includes the information about the redistribution of the 3D
 *  FFT-grid before the actual FFT.
 */
class Parallel3DFFT {

public:
	struct forward_plan {
		/** plan direction: 0 = Forward FFT, 1 = Backward FFT. */
		p3m_int dir;
		/** row direction of that FFT. */
		p3m_int row_dir;
		/** permutations from normal coordinate system. */
		p3m_int n_permute;
		/** number of 1D FFTs. */
		p3m_int n_ffts;
		/** plan for fft. */
		fftw_plan plan;

		/** size of local grid before communication. */
		p3m_int old_grid[3];
		/** size of local grid after communication, also used for actual FFT. */
		p3m_int new_grid[3];
		/** lower left point of local FFT grid in global FFT grid coordinates. */
		p3m_int start[3];
		/** size of new grid (number of grid points). */
		p3m_int new_size;

		/** number of nodes which have to communicate with each other. */
		p3m_int g_size;
		/** group of nodes which have to communicate with each other. */
		p3m_int *group;

		/** packing function for send blocks. */
		void (*pack_function)(fcs_float*, fcs_float*, int*, int*, int*, int);
		/** Send block specification. 6 integers for each node: start[3], size[3]. */
		p3m_int *send_block;
		/** Send block communication sizes. */
		p3m_int *send_size;
		/** Recv block specification. 6 integers for each node: start[3], size[3]. */
		p3m_int *recv_block;
		/** Recv block communication sizes. */
		p3m_int *recv_size;
		/** size of send block elements. */
		p3m_int element;

		void print();
	};

	/** Additional information for backwards FFT.*/
	struct backward_plan {
		/** plan direction. (e.g. fftw makro)*/
		p3m_int dir;
		/** plan for fft. */
		fftw_plan plan;

		/** packing function for send blocks. */
		void (*pack_function)(fcs_float*, fcs_float*, int*, int*, int*, int);
	};

	/***************************************************/
	/* DATA MEMBERS */
	/***************************************************/
public:

private:
  Communication& comm;

  /** Whether FFT is initialized or not. */
  bool is_prepared;
  
  /** Information about the three one dimensional FFTs and how the nodes
   *  have to communicate in between.
   *
   * NOTE: FFT numbering starts with 1 for technical reasons (because we
   *       have 4 node grids, the index 0 is used for the real space
   *       charge assignment grid).  */
  forward_plan plan[4];
  /** Information for Backward FFTs (see fft.plan). */
  backward_plan back[4];
  
  /** Maximal size of the communication buffers. */
  p3m_int max_comm_size;
  
  /** Maximal local grid size. */
  p3m_int max_grid_size;
  
  /** send buffer. */
  p3m_float *send_buf;
  /** receive buffer. */
  p3m_float *recv_buf;

public:
  /***************************************************/
  /* MEMBER FUNCTIONS */
  /***************************************************/
  Parallel3DFFT(Communication &comm);
  ~Parallel3DFFT();
  
  /** Prepare the 3D-FFT for a given size.
      
   * \return Maximal size of local fft grid (needed for allocation of ca_grid).
   * \param local_grid_dim    Pointer to local CA grid dimensions.
   * \param local_grid_margin Pointer to local CA grid margins.
   * \param global_grid_dim   Pointer to global CA grid dimensions.
   * \param global_grid_off   Pointer to global CA grid margins.
   * \param ks_pnum           Pointer to number of permutations in k-space.
   */
  void
  prepare(p3m_int *local_grid_dim, p3m_int *local_grid_margin,
          p3m_int* global_grid_dim, p3m_float *global_grid_off,
          p3m_int *ks_pnum);

    /** Perform the forward 3D FFT. buffer has to be of the same size as data
     * and will be used internally. */
	void forward(p3m_float *data, p3m_float* buffer);
    /** Perform the backward 3D FFT. buffer has to be of the same size as data
     * and will be used internally. */
	void backward(p3m_float *data, p3m_float* buffer);

	p3m_float *malloc_data();
	void free_data(p3m_float* data);

	int getKSSize() const;
	void getKSExtent(const p3m_int*& offset, const p3m_int*& size) const;

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
	static void
	pack_block(p3m_float *in, p3m_float *out, p3m_int start[3], p3m_int size[3],
			p3m_int dim[3], p3m_int element);

	/** unpack a 3d-grid input block (size[3]) into an output 3d-grid
	 *  with dimension dim[3] at start position start[3].
	 *
	 *  see also \ref p3m_fft_pack_block.
	 *
	 *  \param in      pointer to input 3d-grid.
	 *  \param out     pointer to output 3d-grid (block).
	 *  \param start   start index of the block in the in-grid.
	 *  \param size    size of the block (=dimension of the out-grid).
	 *  \param dim     size of the in-grid.
	 *  \param element size of a grid element (e.g. 1 for Real, 2 for Complex).
	 */
	static void
	unpack_block(p3m_float *in, p3m_float *out,
	        p3m_int start[3], p3m_int size[3], p3m_int dim[3],
	        p3m_int element);

	static void
	add_block(p3m_float *in, p3m_float *out,
	        int start[3], int size[3], int dim[3]);

private:
    p3m_float *_malloc_data();
	void forward_grid_comm(forward_plan plan, p3m_float *in, p3m_float *out);
	void backward_grid_comm(forward_plan plan_f, backward_plan plan_b,
			p3m_float *in, p3m_float *out);

	void print_global_grid(Communication &comm,
	        Parallel3DFFT::forward_plan plan, p3m_float *data, p3m_int element,
	        p3m_int num);
};
}
#endif
