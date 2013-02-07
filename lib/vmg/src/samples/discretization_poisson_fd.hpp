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
 * @file   discretization_poisson_fd.hpp
 * @author Julian Iseringhausen <isering@ins.uni-bonn.de>
 * @date   Mon Apr 18 13:01:40 2011
 *
 * @brief  Finite difference discretization of the
 *         Poisson equation.
 *
 */

#ifndef DISCRETIZATION_POISSON_FD_HPP_
#define DISCRETIZATION_POISSON_FD_HPP_

#include "base/discretization.hpp"
#include "base/helper.hpp"

namespace VMG
{

class DiscretizationPoissonFD : public Discretization
{
public:
  DiscretizationPoissonFD(const int& order);

  vmg_float OperatorPrefactor(const Grid& grid) const
  {
    return 1.0 / Helper::pow_2(grid.Extent().MeshWidth().Max());
  }

  void ModifyRightHandSide();
};

}

#endif /* DISCRETIZATION_POISSON_FD_HPP_ */
