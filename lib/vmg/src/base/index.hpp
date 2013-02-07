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
 * @file   index.hpp
 * @author Julian Iseringhausen <isering@ins.uni-bonn.de>
 * @date   Mon Apr 18 12:19:49 2011
 *
 * @brief  Header file for the class VMG::Index.
 *
 */

#ifndef INDEX_HPP_
#define INDEX_HPP_

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iostream>

namespace VMG
{

class Vector;

class Index
{
public:
  Index(const Index& rhs);
  Index(const Vector& vec);
  Index(int x, int y, int z) {i[0]=x;i[1]=y;i[2]=z;}
  Index(int val) {i[0]=val;i[1]=val;i[2]=val;}
  Index(int* arr) {i[0]=arr[0];i[1]=arr[1];i[2]=arr[2];}
  Index() {i[0]=0;i[1]=0;i[2]=0;}

  int& operator[](const int& index) {return i[index];}
  const int& operator[](const int& index) const {return i[index];}

  int& X() {return i[0];}
  int& Y() {return i[1];}
  int& Z() {return i[2];}

  const int& X() const {return i[0];}
  const int& Y() const {return i[1];}
  const int& Z() const {return i[2];}

  Index Abs() {return Index(abs(i[0]), abs(i[1]), abs(i[2]));}
  int Max() const {return std::max(std::max(i[0],i[1]),i[2]);}
  int Min() const {return std::min(std::min(i[0],i[1]),i[2]);}
  int Product() const {return i[0]*i[1]*i[2];}
  int Sum() const {return i[0]+i[1]+i[2];}

  bool IsInBounds(const Index& begin, const Index& end) const
  {
    return i[0] >= begin[0] && i[0] < end[0] &&
           i[1] >= begin[1] && i[1] < end[1] &&
           i[2] >= begin[2] && i[2] < end[2];
  }

  bool IsComponentwiseLess(const Index& other) const
  {
    return this->i[0] < other.i[0]
        && this->i[1] < other.i[1]
        && this->i[2] < other.i[2];
  }

  bool IsComponentwiseLessOrEqual(const Index& other) const
  {
    return this->i[0] <= other.i[0]
        && this->i[1] <= other.i[1]
        && this->i[2] <= other.i[2];
  }

  bool IsComponentwiseGreater(const Index& other) const
  {
    return this->i[0] > other.i[0]
        && this->i[1] > other.i[1]
        && this->i[2] > other.i[2];
  }

  bool IsComponentwiseGreaterOrEqual(const Index& other) const
  {
    return this->i[0] >= other.i[0]
        && this->i[1] >= other.i[1]
        && this->i[2] >= other.i[2];
  }

  Index MaxComponentwise(const Index& rhs) const
  {
    return Index(std::max(i[0], rhs.i[0]), std::max(i[1], rhs.i[1]), std::max(i[2], rhs.i[2]));
  }

  Index MinComponentwise(const Index& rhs) const
  {
    return Index(std::min(i[0], rhs.i[0]), std::min(i[1], rhs.i[1]), std::min(i[2], rhs.i[2]));
  }

  Index Clamp(const Index& lower_bound, const Index& upper_bound) const
  {
    Index index(*this);
    for (int j=0; j<3; ++j) {
      index.i[j] = std::max(index.i[j], lower_bound[j]);
      index.i[j] = std::min(index.i[j], upper_bound[j]);
    }
    return index;
  }

  int* vec() {return i;}
  const int* vec() const {return i;}

  Index& operator=(const Index& rhs)
  {
    i[0] = rhs.X();
    i[1] = rhs.Y();
    i[2] = rhs.Z();
    return *this;
  }

  Index& operator=(const int& rhs)
  {
    i[0] = rhs;
    i[1] = rhs;
    i[2] = rhs;
    return *this;
  }

  bool operator==(const Index& other) const
  {
    return (i[0]==other.X() && i[1]==other.Y() && i[2]==other.Z());
  }

  bool operator<(const Index& other) const
  {
    for (int j=0; j<3; ++j) {
      if (this->i[j] < other.i[j]) return true;
      if (this->i[j] != other.i[j]) return false;
    }

    return false;
  }

  bool operator!=(const Index& other) const
  {
    return !(*this == other);
  }

  Index& operator+=(const Index& rhs)
  {
    i[0] += rhs.X();
    i[1] += rhs.Y();
    i[2] += rhs.Z();
    return *this;
  }

  Index& operator-=(const Index& rhs)
  {
    i[0] -= rhs.X();
    i[1] -= rhs.Y();
    i[2] -= rhs.Z();
    return *this;
  }

  Index& operator*=(const Index& rhs)
  {
    i[0] *= rhs.X();
    i[1] *= rhs.Y();
    i[2] *= rhs.Z();
    return *this;
  }

  Index& operator/=(const Index& rhs)
  {
    i[0] /= rhs.X();
    i[1] /= rhs.Y();
    i[2] /= rhs.Z();
    return *this;
  }

  Index& operator%=(const Index& rhs)
  {
    i[0] %= rhs.X();
    i[1] %= rhs.Y();
    i[2] %= rhs.Z();
    return *this;
  }

  Index& operator+=(const int& rhs)
  {
    i[0] += rhs;
    i[1] += rhs;
    i[2] += rhs;
    return *this;
  }

  Index& operator-=(const int& rhs)
  {
    i[0] -= rhs;
    i[1] -= rhs;
    i[2] -= rhs;
    return *this;
  }

  Index& operator*=(const int& rhs)
  {
    i[0] *= rhs;
    i[1] *= rhs;
    i[2] *= rhs;
    return *this;
  }

  Index& operator/=(const int& rhs)
  {
    i[0] /= rhs;
    i[1] /= rhs;
    i[2] /= rhs;
    return *this;
  }

  Index& operator%=(const int& rhs)
  {
    i[0] %= rhs;
    i[1] %= rhs;
    i[2] %= rhs;
    return *this;
  }

  Index operator+(const Index& rhs) const
  {
    return Index(*this) += rhs;
  }

  Index operator-(const Index& rhs) const
  {
    return Index(*this) -= rhs;
  }

  Index operator*(const Index& rhs) const
  {
    return Index(*this) *= rhs;
  }

  Index operator/(const Index& rhs) const
  {
    return Index(*this) /= rhs;
  }

  Index operator%(const Index& rhs) const
  {
    return Index(*this) %= rhs;
  }

  Index operator+(const int& rhs) const
  {
    return Index(*this) += rhs;
  }

  Index operator-(const int& rhs) const
  {
    return Index(*this) -= rhs;
  }

  Index operator*(const int& rhs) const
  {
    return Index(*this) *= rhs;
  }

  Index operator/(const int& rhs) const
  {
    return Index(*this) /= rhs;
  }

  Index operator%(const int& rhs) const
  {
    return Index(*this) %= rhs;
  }

  Vector operator+(const Vector& rhs) const;
  Vector operator-(const Vector& rhs) const;
  Vector operator*(const Vector& rhs) const;
  Vector operator/(const Vector& rhs) const;

  Vector operator+(const vmg_float& rhs) const;
  Vector operator-(const vmg_float& rhs) const;
  Vector operator*(const vmg_float& rhs) const;
  Vector operator/(const vmg_float& rhs) const;

private:
  int i[3];
};

inline Index operator+(const int& lhs, const Index& rhs)
{
  return Index(lhs) += rhs;
}

inline Index operator-(const int& lhs, const Index& rhs)
{
  return Index(lhs) -= rhs;
}

inline Index operator*(const int& lhs, const Index& rhs)
{
  return Index(lhs) *= rhs;
}

inline Index operator/(const int& lhs, const Index& rhs)
{
  return Index(lhs) /= rhs;
}

inline Index operator%(const int& lhs, const Index& rhs)
{
  return Index(lhs) %= rhs;
}

const Index Max(const Index& index1, const Index& index2);
std::ostream& operator<<(std::ostream& out, const Index& index);

}

#endif /* INDEX_HPP_ */
