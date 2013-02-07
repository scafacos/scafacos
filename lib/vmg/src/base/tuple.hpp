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
 * @file   tuple.hpp
 * @author Julian Iseringhausen <isering@ins.uni-bonn.de>
 * @date   Wed May 18 16:48:46 2011
 *
 * @brief  Templated class for storing tuples.
 *
 */

#ifndef TUPLE_HPP_
#define TUPLE_HPP_

namespace VMG
{

template<class T>
class Tuple
{
public:
  T& Left() {return left;}
  T& Right() {return right;}

  const T& Left() const {return left;}
  const T& Right() const {return right;}

private:
  T left, right;
};

template <class T>
class Tuple3
{
public:
  Tuple<T>& operator[](const int& index) {return tuples[index];}
  Tuple<T>& X() {return tuples[0];}
  Tuple<T>& Y() {return tuples[1];}
  Tuple<T>& Z() {return tuples[2];}

  const Tuple<T>& operator[](const int& index) const {return tuples[index];}
  const Tuple<T>& X() const {return tuples[0];}
  const Tuple<T>& Y() const {return tuples[1];}
  const Tuple<T>& Z() const {return tuples[2];}

private:
  Tuple<T> tuples[3];
};

}

#endif /* TUPLE_HPP_ */
