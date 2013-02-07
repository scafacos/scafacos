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
 * @file   index.cpp
 * @author Julian Iseringhausen <isering@ins.uni-bonn.de>
 * @date   Mon Apr 18 12:19:16 2011
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "base/index.hpp"
#include "base/stencil.hpp"
#include "base/vector.hpp"

using VMG::Index;
using VMG::Vector;

Index::Index(const Index& rhs)
{
  std::memcpy(this->i, rhs.i, 3*sizeof(int));
}

Index::Index(const Vector& vec)
{
  i[0] = static_cast<int>(vec.X());
  i[1] = static_cast<int>(vec.Y());
  i[2] = static_cast<int>(vec.Z());
}

Vector Index::operator+(const Vector& rhs) const
{
  return Vector(*this) += rhs;
}

Vector Index::operator-(const Vector& rhs) const
{
  return Vector(*this) -= rhs;
}

Vector Index::operator*(const Vector& rhs) const
{
  return Vector(*this) *= rhs;
}

Vector Index::operator/(const Vector& rhs) const
{
  return Vector(*this) /= rhs;
}

Vector Index::operator+(const vmg_float& rhs) const
{
  return Vector(*this) += rhs;
}

Vector Index::operator-(const vmg_float& rhs) const
{
  return Vector(*this) -= rhs;
}

Vector Index::operator*(const vmg_float& rhs) const
{
  return Vector(*this) *= rhs;
}

Vector Index::operator/(const vmg_float& rhs) const
{
  return Vector(*this) /= rhs;
}

std::ostream& VMG::operator<<(std::ostream& out, const Index& index)
{
  out << "{" << index.X() << " " << index.Y() << " " << index.Z() << "}";

  return out;
}
