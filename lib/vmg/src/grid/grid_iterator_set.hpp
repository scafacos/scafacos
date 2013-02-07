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
 * @file   grid_iterator_set.hpp
 * @author Julian Iseringhausen <isering@ins.uni-bonn.de>
 * @date   Thu May 12 16:22:02 2011
 *
 * @brief  Classes that group grid iterators usefully.
 *
 */

#ifndef GRID_ITERATOR_SET_HPP_
#define GRID_ITERATOR_SET_HPP_

#include "base/index.hpp"
#include "base/object.hpp"
#include "grid/grid_iterator.hpp"

namespace VMG
{

class GridIteratorSet : public Object
{
public:
  GridIteratorSet()
  {}

  GridIteratorSet(const GridIteratorSet& rhs) :
    begin(rhs.begin)
  {}

  GridIteratorSet(const Index& begin, const Index& end)
  {
    SetBounds(begin, end);
  }

  virtual ~GridIteratorSet()
  {}

  const GridIterator& Begin() const {return begin;}
  static int End() {return end;}

  void SetBounds(const Index& begin_, const Index& end_)
  {
    if ((end_ - begin_).Product() > 0)
      begin = GridIterator(begin_, begin_, end_);
    else
      begin = GridIterator(-1, begin_, end_);
  }

private:
  GridIterator begin;
  static const int end = -1;
};

class GridIteratorSet3
{
public:
  GridIteratorSet3()
  {}

  GridIteratorSet3(const GridIteratorSet3& rhs)
  {
    for (int i=0; i<3; ++i)
      this->iterators[i] = rhs.iterators[i];
  }

  GridIteratorSet& operator[](const int& index) {return iterators[index];}
  const GridIteratorSet& operator[](const int& index) const {return iterators[index];}

  GridIteratorSet& X() {return iterators[0];}
  GridIteratorSet& Y() {return iterators[1];}
  GridIteratorSet& Z() {return iterators[2];}

  const GridIteratorSet& X() const {return iterators[0];}
  const GridIteratorSet& Y() const {return iterators[1];}
  const GridIteratorSet& Z() const {return iterators[2];}

private:
  GridIteratorSet iterators[3];
};

}

#endif /* GRID_ITERATOR_SET_HPP_ */
