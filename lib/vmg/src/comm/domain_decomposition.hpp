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
 * @file   domain_decomposition.hpp
 * @author Julian Iseringhausen <isering@ins.uni-bonn.de>
 * @date   Mon Jun 27 12:22:57 2011
 *
 * @brief  Computes a domain decomposition on a global
 *         multigrid
 *
 */

#ifndef DOMAIN_DECOMPOSITION_HPP_
#define DOMAIN_DECOMPOSITION_HPP_

#include <vector>

namespace VMG
{

class Comm;
class GlobalIndices;
class Interface;
class Multigrid;

class DomainDecomposition
{
public:
  DomainDecomposition() {}
  virtual ~DomainDecomposition() {}

  virtual void Compute(Comm* comm, const Interface* interface, std::vector<GlobalIndices>& global) = 0;
};

}

#endif /* DOMAIN_DECOMPOSITION_HPP_ */
