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
 * @file   has_tempgrids.hpp
 * @author Julian Iseringhausen <isering@ins.uni-bonn.de>
 * @date   Tue Apr  5 21:00:56 2011
 *
 * @brief  Provides an arbitrary number of temporary grids,
 *         accessible by an integer key.
 *
 */

#ifndef HAS_TEMPGRIDS_HPP_
#define HAS_TEMPGRIDS_HPP_

#include <map>

#include "grid/tempgrid.hpp"
#include "mg.hpp"

namespace VMG
{

class HasTempGrids
{
public:
  virtual ~HasTempGrids()
  {
    for (std::map<int, TempGrid*>::iterator iter = grids_1.begin(); iter != grids_1.end(); ++iter)
      delete iter->second;

    for (std::map<const Grid*, TempGrid*>::iterator iter = grids_2.begin(); iter != grids_2.end(); ++iter)
      delete iter->second;
  }

protected:
  TempGrid* GetTempGrid(int key)
  {
    std::map<int, TempGrid*>::iterator iter = grids_1.find(key);

    if (iter == grids_1.end())
      iter = grids_1.insert(std::pair<int, TempGrid*>(key, new TempGrid())).first;

    return iter->second;
  }

  TempGrid* GetTempGrid(const Grid* key)
  {
    std::map<const Grid*, TempGrid*>::iterator iter = grids_2.find(key);

    if (iter == grids_2.end())
      iter = grids_2.insert(std::pair<const Grid*, TempGrid*>(key, new TempGrid())).first;

    return iter->second;
  }

private:
  std::map<int, TempGrid*> grids_1;
  std::map<const Grid*, TempGrid*> grids_2;
};

}

#endif /* HAS_TEMPGRIDS_HPP_ */
