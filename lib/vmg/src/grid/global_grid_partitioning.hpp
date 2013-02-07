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
 * @file   global_grid_partitioning.hpp
 * @author Julian Iseringhausen <isering@ins.uni-bonn.de>
 * @date   Thu May 19 14:02:46 2011
 *
 * @brief  Class to store global grid partitioning.
 *
 */

#ifndef GLOBAL_GRID_PARTITIONING_HPP_
#define GLOBAL_GRID_PARTITIONING_HPP_

#include <map>

#include "base/index.hpp"

namespace VMG
{

class GlobalGridPartitioning
{
public:
  GlobalGridPartitioning(Index pos, Index procs, Index size, int points_min);

  const Index& Pos() const {return pos;}
  const Index& Procs() const {return procs;}
  const Index& Size() const {return size;}
  const int& PointsMin() const {return points_min;}

  const Index& GlobalBegin(const Index& pos) const;
  const Index& GlobalEnd(const Index& pos) const;

private:
  std::map<Index, Index> begin, end;
  Index pos, procs, size;
  int points_min;

};

}

#endif /* GLOBAL_GRID_PARTITIONING_HPP_ */
