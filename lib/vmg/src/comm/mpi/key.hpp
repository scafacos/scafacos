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
 * @file   key.hpp
 * @author Julian Iseringhausen <isering@ins.uni-bonn.de>
 * @date   Mon Nov 21 13:27:22 2011
 *
 * @brief  Class to distinguish the different entities in maps.
 *
 */

#ifndef KEY_HPP_
#define KEY_HPP_

#ifndef HAVE_MPI
#error MPI is needed to compile this class
#endif /* HAVE_MPI */

#include <iostream>
#include <list>

#include "grid/grid.hpp"

namespace VMG
{

class Grid;
class Index;

namespace MPI
{

class KeyStorage
{
public:
  KeyStorage(const Grid& grid) :
    begin(grid.Global().LocalBegin()),
    end(grid.Global().LocalEnd()),
    size_local(grid.Local().SizeTotal()),
    size_global(grid.Global().GlobalSize()),
    level(grid.Level())
  {}

  KeyStorage(const Index& begin, const Index& end,
	     const Index& size_local, const Index& size_global,
	     const int& level) :
    begin(begin),
    end(end),
    size_local(size_local),
    size_global(size_global),
    level(level)
  {}

  KeyStorage(const KeyStorage& other) :
    begin(other.begin),
    end(other.end),
    size_local(other.size_local),
    size_global(other.size_global),
    level(other.level)
  {}

  ~KeyStorage() {}

  bool operator==(const KeyStorage& other) const
  {
    return this->begin == other.begin &&
           this->end == other.end &&
           this->size_local == other.size_local &&
           this->size_global == other.size_global &&
           this->level == other.level;
  }

  bool operator!=(const KeyStorage& other) const
  {
    return !(*this==other);
  }

  bool operator<(const KeyStorage& other) const
  {
    if (this->begin < other.begin) return true;
    if (this->begin != other.begin) return false;
    if (this->end < other.end) return true;
    if (this->end != other.end) return false;
    if (this->size_local < other.size_local) return true;
    if (this->size_local != other.size_local) return false;
    if (this->size_global < other.size_global) return true;
    if (this->size_global != other.size_global) return false;
    if (this->level < other.level) return true;
    return false;
  }

private:
  const Index begin, end, size_local, size_global;
  const int level;
};

class KeySorted {
public:
  KeySorted(const KeySorted& other)
  {
    std::list<KeyStorage>::const_iterator i;
    this->keys.clear();
    for (i=other.keys.begin(); i!=other.keys.end(); ++i)
      this->keys.push_back(*i);
  }

  KeySorted(const Grid& grid)
  {
    keys.push_back(KeyStorage(grid));
  }

  KeySorted(const Grid& grid_1, const Grid& grid_2)
  {
    keys.push_back(KeyStorage(grid_1));
    keys.push_back(KeyStorage(grid_2));
    keys.sort();
  }

  KeySorted(const Index& begin, const Index& end,
	    const Index& size_local, const Index& size_global,
	    const int& level)
  {
    keys.push_back(KeyStorage(begin, end, size_local, size_global, level));
  }

  ~KeySorted() {}

  bool operator<(const KeySorted& other) const
  {
    if (this->keys.size() < other.keys.size()) return true;
    if (this->keys.size() > other.keys.size()) return false;

    std::list<KeyStorage>::const_iterator i1, i2;
    for (i1=this->keys.begin(), i2=other.keys.begin();
	 i1!=this->keys.end() && i2!=other.keys.end();
	 ++i1, ++i2) {
      if (*i1 < *i2) return true;
      if (*i1 != *i2) return false;
    }

    return false;
  }

private:
  std::list<KeyStorage> keys;
};


  /*
   * direction: 0 - from multigrid to temporary grid
   *            1 - from temporary grid to multigrid
   *            for single grid datatypes always 0
   */
class KeyUnsorted {
public:
  KeyUnsorted(const KeyUnsorted& other)
  {
    std::list<KeyStorage>::const_iterator i;
    this->keys.clear();
    for (i=other.keys.begin(); i!=other.keys.end(); ++i)
      this->keys.push_back(*i);
    this->direction = other.direction;
  }

  KeyUnsorted(const Grid& grid, const int& direction)
  {
    keys.push_back(KeyStorage(grid));
    this->direction = direction;
  }

  KeyUnsorted(const Grid& grid_1, const Grid& grid_2, const int& direction)
  {
    keys.push_back(KeyStorage(grid_1));
    keys.push_back(KeyStorage(grid_2));
    this->direction = direction;
  }

  KeyUnsorted(const Index& begin, const Index& end,
	      const Index& size_local, const Index& size_global,
	      const int& level, const int& direction)
  {
    keys.push_back(KeyStorage(begin, end, size_local, size_global, level));
    this->direction = direction;
  }

  ~KeyUnsorted() {}

  bool operator<(const KeyUnsorted& other) const
  {
    if (this->keys.size() < other.keys.size()) return true;
    if (this->keys.size() > other.keys.size()) return false;

    std::list<KeyStorage>::const_iterator i1, i2;
    for (i1=this->keys.begin(), i2=other.keys.begin();
	 i1!=this->keys.end() && i2!=other.keys.end();
	 ++i1, ++i2) {
      if (*i1 < *i2) return true;
      if (*i1 != *i2) return false;
    }

    if (this->direction < other.direction) return true;

    return false;
  }

private:
  std::list<KeyStorage> keys;
  int direction;
};

}

}

#endif /* KEY_HPP_ */
