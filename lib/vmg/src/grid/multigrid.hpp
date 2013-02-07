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
 * @file   multigrid.hpp
 * @author Julian Iseringhausen <isering@ins.uni-bonn.de>
 * @date   Mon Apr 18 12:54:51 2011
 *
 * @brief  VMG::Multigrid
 *
 */

#ifndef MULTIGRID_HPP_
#define MULTIGRID_HPP_

#include <vector>

#include "base/object.hpp"
#include "grid/grid_properties.hpp"

namespace VMG
{

class Comm;
class Grid;
class Index;
class Interface;

class Multigrid : public Object
{
public:
  Multigrid(Comm* comm, const Interface* interface);

  virtual ~Multigrid();

  Grid& operator()() {return *grids[levelIndex];}
  Grid& operator()(const int& level) {return *grids[levelMax-level];}

  const Grid& operator()() const {return *grids[levelIndex];}
  const Grid& operator()(const int& level) const {return *grids[levelMax-level];}

  const int& Level() const {return levelCurrent;}        ///< Current level
  int LevelIndex() const {return levelMax-levelCurrent;} ///< Index of current level

  const int& GlobalMaxLevel() const {return levelGlobalMax;}
  const int& MaxLevel() const {return levelMax;}             ///< Maximum level
  const int& MinLevel() const {return levelMin;}             ///< Minimum level
  const int& NumLevels() const {return numLevels;}           ///< Number of levels

  void SetLevel(int level);

  void ToCoarserLevel(); ///< Switch to next coarser level if possible
  void ToFinerLevel();   ///< Switch to next finer level if possible

  void ClearAll();             ///< Overwrites all grid points with zeros
  void ClearAllCoarseLevels(); ///< Overwrites all grid points on all levels except the finest one with zeros

  void SetCoarserDirichletValues();

protected:

  std::vector<Grid*> grids;

  int levelGlobalMax, levelMin, levelMax;
  int levelIndex, levelCurrent;
  int numLevels;
};

}

#endif /* MULTIGRID_HPP_ */
