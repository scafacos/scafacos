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
 * @file   interface.cpp
 * @author Julian Iseringhausen <isering@ins.uni-bonn.de>
 * @date   Mon Apr 18 12:55:48 2011
 *
 * @brief  VMG::Interface
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <algorithm>
#include <cmath>

#include "base/helper.hpp"
#include "base/interface.hpp"

using namespace VMG;

void Interface::InitInterface(const Vector& box_offset, const vmg_float& box_size,
			      const int& coarseningSteps, const vmg_float& alpha)
{
  int i;
  Index num_cells = Helper::intpow(2, levelMax);
  Index size_factor;

  const Vector box_center = box_offset + 0.5 * box_size;

  for (i=0; i<coarseningSteps; ++i) {

    global.push_back(GlobalIndices());
    extent.push_back(SpatialExtent());

    for (int j=0; j<3; ++j)
      size_factor[j] = (bc[j] == Periodic ? 1 : Helper::intpow(2, static_cast<int>(log(pow(alpha, i+1)) / log(2.0) + 1.0)));

    extent.back().Size() = box_size * static_cast<Vector>(size_factor);
    extent.back().Begin() = box_center - 0.5 * extent.back().Size();
    extent.back().End() = extent.back().Begin() + extent.back().Size();
    extent.back().MeshWidth() = pow(2.0, i-levelMax);

    num_cells = Helper::intpow(2,levelMax-i) * size_factor;

    global.back().GlobalSize() = num_cells + 1;
    global.back().LocalSize() = num_cells + 1;
    global.back().LocalBegin() = 0;
    global.back().LocalEnd() = num_cells + 1;

    global.back().FinestAbsSize() = Helper::intpow(2,i) * num_cells + 1;
    global.back().FinestAbsBegin() = ((global.back().FinestAbsSize()-1) * (1-size_factor)) / (2*size_factor);
    global.back().FinestAbsEnd() = global.back().FinestAbsBegin() + global.back().FinestAbsSize();

    if (i==0)
      global.back().GlobalFinerBegin() = (num_cells - Helper::intpow(2, levelMax-i))/2;
    else
      global.back().GlobalFinerBegin() = (global.back().FinestAbsSize() - (++global.rbegin())->FinestAbsSize()) / Helper::intpow(2,i+1);

    global.back().GlobalFinerEnd() = global.back().GlobalSize() - global.back().GlobalFinerBegin();
    global.back().GlobalFinerSize() = global.back().GlobalFinerEnd() - global.back().GlobalFinerBegin();

    global.back().LocalFinerBegin() = global.back().GlobalFinerBegin();
    global.back().LocalFinerEnd() = global.back().GlobalFinerEnd();
    global.back().LocalFinerSize() = global.back().GlobalFinerSize();
  }

  while (global.size() == 0 || global.back().GlobalSize().Min() > Helper::intpow(2, levelMin)+1) {

    if (global.size() > 0)
      num_cells /= 2;

    global.push_back(GlobalIndices());
    extent.push_back(SpatialExtent());

    if (global.size() == 1) {
      extent.back().Size() = box_size;
      extent.back().Begin() = box_offset;
      extent.back().End() = box_offset + box_size;
      extent.back().MeshWidth() = box_size / static_cast<Vector>(num_cells);
    }else {
      extent.back().Size() = (++extent.rbegin())->Size();
      extent.back().Begin() = (++extent.rbegin())->Begin();
      extent.back().End() = (++extent.rbegin())->End();
      extent.back().MeshWidth() = 2.0 * (++extent.rbegin())->MeshWidth();
    }

    global.back().GlobalSize() = num_cells + 1;
    global.back().LocalSize() = num_cells + 1;
    global.back().LocalBegin() = 0;
    global.back().LocalEnd() = num_cells + 1;

    if (global.size() == 1) {
      global.back().FinestAbsBegin() = 0;
      global.back().FinestAbsEnd() = global.back().GlobalSize();
      global.back().FinestAbsSize() = global.back().GlobalSize();
    }else {
      global.back().FinestAbsBegin() = (++global.rbegin())->FinestAbsBegin();
      global.back().FinestAbsEnd() = (++global.rbegin())->FinestAbsEnd();
      global.back().FinestAbsSize() = (++global.rbegin())->FinestAbsSize();
    }

    global.back().GlobalFinerBegin() = 0;
    global.back().GlobalFinerEnd() = global.back().GlobalSize();
    global.back().GlobalFinerSize() = global.back().GlobalSize();

    global.back().LocalFinerBegin() = global.back().GlobalFinerBegin();
    global.back().LocalFinerEnd() = global.back().GlobalFinerEnd();
    global.back().LocalFinerSize() = global.back().GlobalFinerSize();

  }

  for (i=0; i<3; ++i)
    if (bc[i] == Periodic)
      for (unsigned int j=0; j<global.size(); ++j) {
	global[j].GlobalSize()[i] -= 1;
	global[j].FinestAbsSize()[i] -= Helper::intpow(2, j);
	global[j].FinestAbsEnd()[i] -= Helper::intpow(2, j);
	global[j].LocalSize()[i] -= 1;
	global[j].LocalEnd()[i] -= 1;
	global[j].GlobalFinerSize()[i] -= 1;
	global[j].GlobalFinerEnd()[i] -= 1;
	global[j].LocalFinerSize()[i] -= 1;
	global[j].LocalFinerEnd()[i] -= 1;
      }

  levelMin = levelMax - global.size() + 1;

  global.back().BoundaryType() = GlobalCoarsened;

  for (i=global.size()-2; i>=0; --i) {
    if (global[i].FinestAbsSize().Product() >= global[i+1].FinestAbsSize().Product()) {
      global[i].BoundaryType() = GlobalCoarsened;
    }else {
      global[i].BoundaryType() = LocallyRefined;
      global[i+1].BoundaryType() = GlobalMax;
      break;
    }
  }

  for (; i>=0; --i)
    global[i].BoundaryType() = LocallyRefined;

  if (global.front().BoundaryType() != LocallyRefined &&
      global.front().BoundaryType() != GlobalMax)
    global.front().BoundaryType() = GlobalMax;

}
