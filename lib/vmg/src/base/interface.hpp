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
 * @file   interface.hpp
 * @author Julian Iseringhausen <isering@ins.uni-bonn.de>
 * @date   Mon Apr 18 12:56:05 2011
 *
 * @brief  VMG::Interface
 *
 */

#ifndef INTERFACE_HPP_
#define INTERFACE_HPP_

#include <vector>

#include "base/object.hpp"
#include "grid/grid_properties.hpp"

namespace VMG
{

class Grid;
class Multigrid;

class Interface : public Object
{
public:
  Interface(const Boundary& boundary, const int& levelMin, const int& levelMax,
	    const Vector& box_offset, const vmg_float box_size,
	    int coarseningSteps=0, vmg_float alpha=1.6,
            bool register_ = true) :
    Object("INTERFACE", register_),
    bc(boundary),
    levelMin(levelMin),
    levelMax(levelMax)
  {
    InitInterface(box_offset, box_size, coarseningSteps, alpha);
  }

  virtual ~Interface() {}

  virtual void ImportRightHandSide(Multigrid& grid) = 0;
  virtual void ExportSolution(Grid& grid) = 0;

  const std::vector<GlobalIndices>& Global() const {return global;}
  const std::vector<SpatialExtent>& Extent() const {return extent;}

  const GlobalIndices& Global(const int& level) const
  {
    return global[levelMax - level];
  }

  const SpatialExtent& Extent(const int& level) const
  {
    return extent[levelMax - level];
  }

  const Boundary& BoundaryConditions() const {return bc;}

  const int& MinLevel() const {return levelMin;}
  const int& MaxLevel() const {return levelMax;}

private:
  void InitInterface(const Vector& box_offset, const vmg_float& box_end,
		     const int& coarseningSteps, const vmg_float& alpha);

  std::vector<GlobalIndices> global;
  std::vector<SpatialExtent> extent;

  Boundary bc;

  int levelMin, levelMax;
};

}

#endif /* INTERFACE_HPP_ */
