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
 * @file   matrix.hpp
 * @author Julian Iseringhausen <isering@ins.uni-bonn.de>
 * @date   Mon Apr 18 12:20:22 2011
 *
 * @brief  Header file for the class VMG::Matrix, which represents a 3x3 matrix
 *         and defines some useful operations.
 *
 */

#ifndef MATRIX_HPP_
#define MATRIX_HPP_

#include "base/vector.hpp"

namespace VMG
{

class Matrix
{
public:
  Matrix()
  {
    for (int i=0; i<9; i++)
      mat[i] = 0.0;
  }

  vmg_float& operator() (int i, int j) {return mat[j+3*i];}
  const vmg_float& GetVal(int i, int j) const {return mat[j+3*i];}

  const Vector operator*(const Vector& rhs) const
  {
    Vector result;

    for (int i=0; i<3; i++) {
      result.X() += this->GetVal(0, i) * rhs.X();
      result.Y() += this->GetVal(1, i) * rhs.Y();
      result.Z() += this->GetVal(2, i) * rhs.Z();
    }

    return result;
  }

private:
  vmg_float mat[9];
};

const Matrix abT(const Vector lhs, const Vector rhs);

}

#endif /* MATRIX_HPP_ */
