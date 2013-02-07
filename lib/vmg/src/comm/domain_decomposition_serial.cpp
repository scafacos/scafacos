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
 * @file   domain_decomposition_serial.cpp
 * @author Julian Iseringhausen <isering@ins.uni-bonn.de>
 * @date   Mon Jun 27 12:33:10 2011
 *
 * @brief  Domain decomposition which gives all grid points
 *         to one process
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "base/interface.hpp"
#include "comm/domain_decomposition_serial.hpp"
#include "grid/grid.hpp"
#include "grid/multigrid.hpp"

using namespace VMG;

void DomainDecompositionSerial::Compute(Comm* comm, const Interface* interface, std::vector<GlobalIndices>& global)
{
  global = interface->Global();
}
