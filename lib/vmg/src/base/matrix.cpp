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
 * @file   matrix.cpp
 * @author Julian Iseringhausen <isering@ins.uni-bonn.de>
 * @date   Mon Apr 18 12:20:22 2011
 *
 * @brief  Source file for the class VMG::Matrix, which represents a 3x3 matrix
 *         and defines some useful operations.
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "base/matrix.hpp"

using namespace VMG;

const Matrix abT(const Vector lhs, const Vector rhs)
{
  Matrix matrix;

  for (int i=0; i<3; i++)
    for (int j=0; j<3; j++)
      matrix(i,j) = lhs.vec()[i] * rhs.vec()[j];

  return matrix;
}
