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
 * @file   discretization.hpp
 * @author Julian Iseringhausen <isering@ins.uni-bonn.de>
 * @date   Tue Apr  5 20:32:09 2011
 *
 * @brief  Base class for controlling the discretization of
 *         the continuous system of equations. The discretized
 *         operator must be specified as a stencil.
 *
 */

#ifndef DISCRETIZATION_HPP_
#define DISCRETIZATION_HPP_

#include <algorithm>
#include <string>

#include "base/object.hpp"
#include "base/stencil.hpp"
#include "grid/grid.hpp"
#include "grid/multigrid.hpp"
#include "mg.hpp"

namespace VMG
{

class Discretization : public Object
{
public:
  Discretization(bool register_ = true) :
    Object("DISCRETIZATION", register_),
    stencil(1.0),
    order(2)
  {}

  Discretization(const int& order, bool register_ = true) :
    Object("DISCRETIZATION", register_),
    stencil(1.0),
    order(order)
  {}

  Discretization(const Stencil& stencil_, const int& order, bool register_ = true) :
    Object("DISCRETIZATION", register_),
    stencil(stencil_),
    order(order)
  {}

  Discretization(std::string id) :
    Object(id),
    stencil(1.0),
    order(2)
  {}

  Discretization(std::string id, const int& order) :
    Object(id),
    stencil(1.0),
    order(order)
  {}

  Discretization(std::string id, const Stencil& stencil_, const int& order) :
    Object(id),
    stencil(stencil_),
    order(order)
  {}

  const Stencil& GetStencil() const {return stencil;} ///< Returns the stencil of the discretized operator.

  virtual vmg_float OperatorPrefactor(const Grid& grid) const = 0; ///< Returns the prefactor of the operator.

  virtual void ModifyRightHandSide() {}

  /**
   * This function gets called whenever boundary points at inner boundaries are needed.
   * Inner boundaries occur when using adaptive grid refinement.
   *
   * @param sol_fine Solution vector / fine level
   * @param rhs_fine Right handside vector / fine level
   * @param sol_coarse Solution vector / coarse level
   */
  void SetInnerBoundary(Grid& sol_fine, Grid& rhs_fine, Grid& sol_coarse) const
  {
    if (sol_fine.Global().BoundaryType() == LocallyRefined)
      SetInnerBoundaryCompute(sol_fine, rhs_fine, sol_coarse);
  }

private:
  virtual void SetInnerBoundaryCompute(Grid& sol_fine, Grid& rhs_fine, Grid& sol_coarse) const {}

protected:
  VMG::Stencil stencil;
  int order;
};

}

#endif /* DISCRETIZATION_HPP_ */
