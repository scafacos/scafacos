/*
 Copyright (C) 2014,2011 Olaf Lenz

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
#ifndef _P3M_COMMUNICATION_H
#define _P3M_COMMUNICATION_H

#include "p3mconfig.hpp"
#include <mpi.h>

namespace P3M {

struct Communication {
	Communication(MPI_Comm mpicomm);
	~Communication();

	void prepare(p3m_float box_l[3]);

	/* The MPI communicator to use (cartesian) */
	MPI_Comm mpicomm;
	/* The original MPI communicator to use (possibly not cartesian) */
	MPI_Comm mpicomm_orig;
	/* The size of the communicator */
	p3m_int size;
	/* The rank within the communicator */
	p3m_int rank;

	/** The number of nodes in each spatial dimension. */
	p3m_int node_grid[3];
	/** position of this node in the node grid */
	p3m_int node_pos[3];
	/** the six nearest neighbors of a node in the node grid. */
	p3m_int node_neighbors[6];

	/** Size of the local box. */
	p3m_float local_box_l[3];
	/** Left (bottom, front) corner of this nodes local box. */
	p3m_float my_left[3];
	/** Right (top, back) corner of this nodes local box. */
	p3m_float my_right[3];
};
}
#endif
