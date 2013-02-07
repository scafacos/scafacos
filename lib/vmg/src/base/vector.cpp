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
 * @file   vector.cpp
 * @author Julian Iseringhausen <isering@ins.uni-bonn.de>
 * @date   Mon Apr 18 12:25:24 2011
 *
 * @brief  VMG::Vector
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "base/index.hpp"
#include "base/vector.hpp"

using namespace VMG;

Vector::Vector(const Index& index)
{
  i[0] = static_cast<vmg_float>(index.X());
  i[1] = static_cast<vmg_float>(index.Y());
  i[2] = static_cast<vmg_float>(index.Z());
}

Vector::Vector(const vmg_float& x, const vmg_float& y, const vmg_float& z)
{
  i[0]=x;
  i[1]=y;
  i[2]=z;
}

Vector::Vector(const vmg_float& val)
{
  std::fill(i, i+3, val);
}

Vector::Vector(const vmg_float* arr)
{
  std::memcpy(i, arr, 3*sizeof(vmg_float));
}

Vector::Vector()
{
  std::fill(i, i+3, 0.0);
}

std::ostream& VMG::operator<<(std::ostream& out, const Vector& vector)
{
  out << "{" << vector.X() << " " << vector.Y() << " " << vector.Z() << "}";

  return out;
}
