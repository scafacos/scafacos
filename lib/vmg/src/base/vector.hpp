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
 * @file   vector.hpp
 * @author Julian Iseringhausen <isering@ins.uni-bonn.de>
 * @date   Mon Apr 18 12:26:13 2011
 *
 * @brief  VMG::Vector
 *
 */

#ifndef VECTOR_HPP_
#define VECTOR_HPP_

#include <algorithm>
#include <cmath>
#include <cstring>
#include <iostream>

namespace VMG
{

class Index;

class Vector
{
public:
  Vector(const Index& index);
  Vector(const vmg_float& x, const vmg_float& y, const vmg_float& z);
  Vector(const vmg_float& val);
  Vector(const vmg_float* arr);
  Vector();

  Vector(const Vector& other)
  {
    std::memcpy(this->i, other.i, 3*sizeof(vmg_float));
  }

  vmg_float& operator[](const int& index) {return i[index];}
  const vmg_float& operator[](const int& index) const {return i[index];}

  vmg_float& X() {return i[0];}
  vmg_float& Y() {return i[1];}
  vmg_float& Z() {return i[2];}

  const vmg_float& X() const {return i[0];}
  const vmg_float& Y() const {return i[1];}
  const vmg_float& Z() const {return i[2];}

  Vector Abs() const {return Vector(fabs(i[0]), fabs(i[1]), fabs(i[2]));}
  vmg_float Length() const {return sqrt(i[0]*i[0]+i[1]*i[1]+i[2]*i[2]);}
  vmg_float Max() const {return std::max(std::max(i[0],i[1]),i[2]);}
  vmg_float Min() const {return std::min(std::min(i[0],i[1]),i[2]);}
  vmg_float Product() const {return i[0]*i[1]*i[2];}

  Vector& Floor()
  {
    i[0] = std::floor(i[0]);
    i[1] = std::floor(i[1]);
    i[2] = std::floor(i[2]);
    return *this;
  }

  bool IsComponentwiseLess(const Vector& other) const
  {
    return this->i[0] < other.i[0]
        && this->i[1] < other.i[1]
        && this->i[2] < other.i[2];
  }

  bool IsComponentwiseLessOrEqual(const Vector& other) const
  {
    return this->i[0] <= other.i[0]
        && this->i[1] <= other.i[1]
        && this->i[2] <= other.i[2];
  }

  bool IsComponentwiseGreater(const Vector& other) const
  {
    return this->i[0] > other.i[0]
        && this->i[1] > other.i[1]
        && this->i[2] > other.i[2];
  }

  bool IsComponentwiseGreaterOrEqual(const Vector& other) const
  {
    return this->i[0] >= other.i[0]
        && this->i[1] >= other.i[1]
        && this->i[2] >= other.i[2];
  }

  bool IsInBounds(const Vector& begin, const Vector& end) const
  {
    return this->IsComponentwiseGreaterOrEqual(begin) &&
           this->IsComponentwiseLessOrEqual(end);
  }

  Vector MaxComponentwise(const Vector& rhs) const
  {
    return Vector(std::max(i[0], rhs.i[0]), std::max(i[1], rhs.i[1]), std::max(i[2], rhs.i[2]));
  }

  Vector MinComponentwise(const Vector& rhs) const
  {
    return Vector(std::min(i[0], rhs.i[0]), std::min(i[1], rhs.i[1]), std::min(i[2], rhs.i[2]));
  }

  vmg_float* vec() {return i;}
  const vmg_float* vec() const {return i;}

  Vector& operator=(const Vector& rhs)
  {
    i[0] = rhs.X();
    i[1] = rhs.Y();
    i[2] = rhs.Z();
    return *this;
  }

  Vector& operator=(const vmg_float& rhs)
  {
    i[0] = rhs;
    i[1] = rhs;
    i[2] = rhs;
    return *this;
  }

  bool operator==(const Vector& other) const
  {
    return (i[0]==other.X() && i[1]==other.Y() && i[2]==other.Z());
  }

  bool operator!=(const Vector& other) const
  {
    return !(*this == other);
  }

  Vector& operator+=(const Vector& rhs)
  {
    i[0] += rhs.X();
    i[1] += rhs.Y();
    i[2] += rhs.Z();
    return *this;
  }

  Vector& operator-=(const Vector& rhs)
  {
    i[0] -= rhs.X();
    i[1] -= rhs.Y();
    i[2] -= rhs.Z();
    return *this;
  }

  Vector& operator*=(const Vector& rhs)
  {
    i[0] *= rhs.X();
    i[1] *= rhs.Y();
    i[2] *= rhs.Z();
    return *this;
  }

  Vector& operator/=(const Vector& rhs)
  {
    i[0] /= rhs.X();
    i[1] /= rhs.Y();
    i[2] /= rhs.Z();
    return *this;
  }

  Vector& operator+=(const vmg_float& rhs)
  {
    i[0] += rhs;
    i[1] += rhs;
    i[2] += rhs;
    return *this;
  }

  Vector& operator-=(const vmg_float& rhs)
  {
    i[0] -= rhs;
    i[1] -= rhs;
    i[2] -= rhs;
    return *this;
  }

  Vector& operator*=(const vmg_float& rhs)
  {
    i[0] *= rhs;
    i[1] *= rhs;
    i[2] *= rhs;
    return *this;
  }

  Vector& operator/=(const vmg_float& rhs)
  {
    i[0] /= rhs;
    i[1] /= rhs;
    i[2] /= rhs;
    return *this;
  }

  Vector operator+(const Vector& rhs) const
  {
    return Vector(*this) += rhs;
  }

  Vector operator-(const Vector& rhs) const
  {
    return Vector(*this) -= rhs;
  }

  Vector operator*(const Vector& rhs) const
  {
    return Vector(*this) *= rhs;
  }

  Vector operator/(const Vector& rhs) const
  {
    return Vector(*this) /= rhs;
  }

  Vector operator+(const vmg_float& rhs) const
  {
    return Vector(*this) += rhs;
  }

  Vector operator-(const vmg_float& rhs) const
  {
    return Vector(*this) -= rhs;
  }

  Vector operator*(const vmg_float& rhs) const
  {
    return Vector(*this) *= rhs;
  }

  Vector operator/(const vmg_float& rhs) const
  {
    return Vector(*this) /= rhs;
  }

private:
  vmg_float i[3];
};

inline Vector operator+(const vmg_float& lhs, const Vector& rhs)
{
  return Vector(lhs) += rhs;
}

inline Vector operator-(const vmg_float& lhs, const Vector& rhs)
{
  return Vector(lhs) -= rhs;
}

inline Vector operator*(const vmg_float& lhs, const Vector& rhs)
{
  return Vector(lhs) *= rhs;
}

inline Vector operator/(const vmg_float& lhs, const Vector& rhs)
{
  return Vector(lhs) /= rhs;
}

inline vmg_float Norm_2(const Vector& vec)
{
  return sqrt(vec.X()*vec.X()+vec.Y()*vec.Y()+vec.Z()*vec.Z());
}

inline vmg_float InnerProduct(const Vector& vec1, const Vector& vec2)
{
  return vec1.X()*vec2.X()+vec1.Y()*vec2.Y()+vec1.Z()*vec2.Z();
}

std::ostream& operator<<(std::ostream& out, const Vector& vector);

}

#endif /* VECTOR_HPP_ */
