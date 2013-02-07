/*
 *    vmg - a versatile multigrid solver
 *    Copyright (C) 2012 Institute for Numerical Simulation, University of Bonn
 *
 *  vmg is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  vmg is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * @file   error_handler.cpp
 * @author Julian Iseringhausen <isering@ins.uni-bonn.de>
 * @date   Mon Nov 21 13:27:22 2011
 *
 * @brief  Convert MPI errors to C++ exceptions. Used to be
 *         able to call a debugger when MPI crashes
 *         internally.
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef DEBUG

#ifdef HAVE_MPI

#include <mpi.h>
#ifdef HAVE_MARMOT
#include <enhancempicalls.h>
#include <sourceinfompicalls.h>
#endif

#include "comm/mpi/error_handler.hpp"

void VMG::MPI::ConvertToException(MPI_Comm* comm, int* err, ...)
{
  char error_string[MPI_MAX_ERROR_STRING];
  int resultlen;
  MPI_Error_string(*err, error_string, &resultlen);
  throw Exception(error_string);
}

#endif /* DEBUG */

#endif /* HAVE_MPI */
