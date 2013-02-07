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
 * @file   defs.hpp
 * @author Julian Iseringhausen <isering@ins.uni-bonn.de>
 * @date   Tue Apr  5 19:40:03 2011
 *
 * @brief  Very basic types get defined here.
 *
 *
 */

#ifndef DEFS_HPP_
#define DEFS_HPP_

#include <cassert>

namespace VMG
{

/**
 * Some commands may want to control the control flow of
 * the command list they belong to. This can be achieved
 * by returning one of these values.
 *
 */
enum Request {
  Continue,      ///< Continue execution normally
  StopCycleNow,  ///< Stop execution of command list immediately
  StopCycleLater ///< Stop execution of loop after execution of all commands.
};

/**
 * This enum specifies the available boundary conditions.
 *
 */
enum BC {
  Periodic,
  Dirichlet,
  Open
};

class Boundary
{
public:
  Boundary(const BC& boundary_x, const BC& boundary_y, const BC& boundary_z)
  {
    bc[0] = boundary_x;
    bc[1] = boundary_y;
    bc[2] = boundary_z;
  }

  const BC& operator[](const int& index) const {return bc[index];}

  const BC& X() const {return bc[0];}
  const BC& Y() const {return bc[1];}
  const BC& Z() const {return bc[2];}

private:
  BC bc[3];
};

/**
 * The boundaries at different grid levels may have to
 * be handled differently. This enum specifies the type
 * of the grid level.
 *
 */
enum BT {
  LocallyRefined,  ///< For adaptive grids. Level is above the finest global grid.
  GlobalMax,       ///< Finest global grid.
  GlobalCoarsened, ///< Coarse global grid.
  EmptyGrid,        ///< This grid does not contain any data.
  BTUndefined
};

}

#endif /* DEFS_HPP_ */
