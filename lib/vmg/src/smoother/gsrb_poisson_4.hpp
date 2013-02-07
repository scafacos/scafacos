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
 * @file   gsrb_poisson_4.hpp
 * @author Julian Iseringhausen <isering@ins.uni-bonn.de>
 * @date   Fri May 11 18:37:36 2012
 *
 * @brief  Gauss-Seidel Red-Black method, specialized to
 *         the Poisson equation. Performance improved by
 *         explicit loop unrolling.
 *
 */

#ifndef GSRB_POISSON_4_HPP_
#define GSRB_POISSON_4_HPP_

#include "smoother/smoother.hpp"

namespace VMG
{

class Multigrid;

class GaussSeidelRBPoisson4 : public Smoother
{
public:
  GaussSeidelRBPoisson4(bool register_ = true) :
    Smoother(register_)
  {}

  virtual ~GaussSeidelRBPoisson4() {}

private:
  void Compute(Grid& sol, Grid& rhs);
};

}

#endif /* GSRB_POISSON_4_HPP_ */
