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
 * @file   grid_iterator.hpp
 * @author Julian Iseringhausen <isering@ins.uni-bonn.de>
 * @date   Wed Apr 20 12:24:09 2011
 *
 * @brief  Iterator class to iterate over subgrids.
 *
 */

#ifndef GRID_ITERATOR_HPP_
#define GRID_ITERATOR_HPP_

#include <list>

#include "base/index.hpp"
#include "base/object.hpp"

namespace VMG
{

class GridIterator : public Object
{
public:
  GridIterator() :
    index(0),
    begin(0),
    end(0)
  {}

  GridIterator(const Index& index, const Index& begin, const Index& end) :
    index(index),
    begin(begin),
    end(end)
  {}

  GridIterator(const GridIterator& rhs) :
    index(rhs.index),
    begin(rhs.begin),
    end(rhs.end)
  {}

  virtual ~GridIterator()
  {}

  GridIterator& operator=(const GridIterator& rhs)
  {
    this->index = rhs.index;
    this->begin = rhs.begin;
    this->end = rhs.end;
    return *this;
  }

  Index& operator*();
  Index* operator->();

  bool operator==(const GridIterator& rhs) const;
  bool operator==(const Index& rhs) const;
  bool operator==(const int& rhs) const;

  bool operator!=(const GridIterator& rhs) const;
  bool operator!=(const Index& rhs) const;
  bool operator!=(const int& rhs) const;

  GridIterator& operator++();
  GridIterator operator++(int);

  GridIterator& operator--();
  GridIterator operator--(int);

  GridIterator& operator+=(const int& steps);

  const Index& GetIndex() const {return index;}
  const Index& GetBegin() const {return begin;}
  const Index& GetEnd() const {return end;}

private:
  Index index;
  Index begin, end;
};

inline Index& GridIterator::operator*()
{
  return this->index;
}

inline Index* GridIterator::operator->()
{
  return &*(static_cast<GridIterator>(*this));
}

inline bool GridIterator::operator==(const GridIterator& rhs) const
{
  return this->index == rhs.index;
}

inline bool GridIterator::operator==(const Index& rhs) const
{
  return this->index == rhs;
}

inline bool GridIterator::operator==(const int& rhs) const
{
  return this->index[0] == rhs &&
         this->index[1] == rhs &&
         this->index[2] == rhs;
}

inline bool GridIterator::operator!=(const GridIterator& rhs) const
{
  return this->index != rhs.index;
}

inline bool GridIterator::operator!=(const Index& rhs) const
{
  return this->index != rhs;
}

inline bool GridIterator::operator!=(const int& rhs) const
{
  return this->index[0] != rhs ||
         this->index[1] != rhs ||
         this->index[2] != rhs;
}

inline GridIterator& GridIterator::operator++()
{
  if (++index.Z() >= end.Z()) {
    index.Z() = begin.Z();
    if (++index.Y() >= end.Y()) {
      index.Y() = begin.Y();
      if (++index.X() >= end.X())
	index = -1;
    }
  }

  return *this;
}

inline GridIterator GridIterator::operator++(int)
{
  GridIterator result = *this;
  ++(*this);
  return result;
}

inline GridIterator& GridIterator::operator--()
{
  if (--index.Z() < begin.Z()) {
    index.Z() = end.Z() - 1;
    if (--index.Y() < begin.Y()) {
      index.Y() = end.Y() - 1;
      if (--index.X() < begin.X())
	index = -1;
    }
  }

  return *this;
}

inline GridIterator GridIterator::operator--(int)
{
  GridIterator result = *this;
  --(*this);
  return result;
}

inline GridIterator& GridIterator::operator+=(const int& steps)
{
  index.Z() += steps;
  if (index.Z() >= end.Z()) {
    index.Z() = begin.Z();
    index.Y() += steps;
    if (index.Y() >= end.Y()) {
      index.Y() = begin.Y();
      index.X() += steps;
      if (index.X() >= end.X())
	index = -1;
    }
  }

  return *this;
}

}

#endif /* GRID_ITERATOR_HPP_ */
