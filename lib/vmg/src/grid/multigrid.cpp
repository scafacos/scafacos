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
 * @file   multigrid.cpp
 * @author Julian Iseringhausen <isering@ins.uni-bonn.de>
 * @date   Mon Apr 18 12:54:37 2011
 *
 * @brief  VMG::Multigrid
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>

#include "base/helper.hpp"
#include "base/index.hpp"
#include "base/interface.hpp"
#include "comm/comm.hpp"
#include "comm/domain_decomposition.hpp"
#include "grid/grid.hpp"
#include "grid/grid_properties.hpp"
#include "grid/multigrid.hpp"

using namespace VMG;

Multigrid::Multigrid(Comm* comm, const Interface* interface)
{
  Index points, remainder;
  LocalIndices local_l;
  std::vector<GlobalIndices> global_separated;

  const std::vector<GlobalIndices>& global = interface->Global();
  const std::vector<SpatialExtent>& extent = interface->Extent();

  comm->GetDomainDecomposition().Compute(comm, interface, global_separated);

  for (unsigned int i=0; i<global.size(); ++i) {

    /*
     * Check if this level is the finest global level
     */
    if (global[i].BoundaryType() == GlobalMax)
      levelGlobalMax = interface->MaxLevel() - i;


    if (global_separated[i].LocalSize().Product() > 0) {

      /*
       * Initialize some properties with zero
       */
      local_l.HaloBegin1() = 0;
      local_l.HaloEnd1() = 0;
      local_l.HaloBegin2() = 0;
      local_l.HaloEnd2() = 0;
      local_l.BoundaryBegin1() = 0;
      local_l.BoundaryEnd1() = 0;
      local_l.BoundaryBegin2() = 0;
      local_l.BoundaryEnd2() = 0;

      /*
       * Set boundary dependant properties
       */
      for (int j=0; j<3; ++j) {

	if (comm->BoundaryConditions()[j] == Dirichlet ||
	    comm->BoundaryConditions()[j] == Open) {

	  local_l.SizeTotal()[j] = global_separated[i].LocalSize()[j] +
	    (global_separated[i].BoundaryType() == LocallyRefined ? 2 : 0);

	  /*
	   * We have a boundary at the ends of the process grids
	   * and halo grid points otherwise
	   */
	  if (global_separated[i].LocalBegin()[j] == 0) {
	    local_l.BoundaryBegin1()[j] = 0;
	    local_l.BoundaryEnd1()[j] = 1;
	  }else {
	    ++local_l.SizeTotal()[j];
	    local_l.HaloBegin1()[j] = 0;
	    local_l.HaloEnd1()[j] = 1;
	  }

	  if (global_separated[i].LocalEnd()[j] == global_separated[i].GlobalSize()[j]) {
	    local_l.BoundaryBegin2()[j] = local_l.SizeTotal()[j] - 1;
	    local_l.BoundaryEnd2()[j] = local_l.SizeTotal()[j];
	  }else {
	    ++local_l.SizeTotal()[j];
	    local_l.HaloBegin2()[j] = local_l.SizeTotal()[j] - 1;
	    local_l.HaloEnd2()[j] = local_l.SizeTotal()[j];
	  }

	}else if (comm->BoundaryConditions()[j] == Periodic) {

	  local_l.SizeTotal()[j] = global_separated[i].LocalSize()[j] + 2;

	  /*
	   * No boundary
	   */
	  local_l.BoundaryBegin1()[j] = local_l.BoundaryEnd1()[j] = 0;
	  local_l.BoundaryBegin2()[j] = local_l.BoundaryEnd2()[j] = 0;

	  /*
	   * Halo grid points on all processes
	   */
	  local_l.HaloBegin1()[j] = 0;
	  local_l.HaloEnd1()[j] = 1;
	  local_l.HaloBegin2()[j] = local_l.SizeTotal()[j] - 1;
	  local_l.HaloEnd2()[j] = local_l.SizeTotal()[j];

	}

      }

      local_l.Begin() = 1;
      local_l.End() = local_l.SizeTotal() - 1;

      local_l.Size() = local_l.End() - local_l.Begin();
      local_l.HaloSize1() = local_l.HaloEnd1() - local_l.HaloBegin1();
      local_l.HaloSize2() = local_l.HaloEnd2() - local_l.HaloBegin2();
      local_l.BoundarySize1() = local_l.BoundaryEnd1() - local_l.BoundaryBegin1();
      local_l.BoundarySize2() = local_l.BoundaryEnd2() - local_l.BoundaryBegin2();

      local_l.FinerSize() = global_separated[i].LocalFinerSize() - local_l.BoundarySize1() - local_l.BoundarySize2();
      local_l.FinerBegin() = global_separated[i].LocalFinerBegin() - global_separated[i].LocalBegin() + local_l.Begin();
      local_l.FinerEnd() = local_l.FinerBegin() + local_l.FinerSize();


    }else {

      local_l.Begin() = 0;
      local_l.End() = 0;
      local_l.Size() = 0;
      local_l.SizeTotal() = 0;
      local_l.HaloBegin1() = 0;
      local_l.HaloEnd1() = 0;
      local_l.HaloSize1() = 0;
      local_l.HaloBegin2() = 0;
      local_l.HaloEnd2() = 0;
      local_l.HaloSize2() = 0;
      local_l.BoundaryBegin1() = 0;
      local_l.BoundaryEnd1() = 0;
      local_l.BoundarySize1() = 0;
      local_l.BoundaryBegin2() = 0;
      local_l.BoundaryEnd2() = 0;
      local_l.BoundarySize2() = 0;
      local_l.FinerBegin() = 0;
      local_l.FinerEnd() = 0;
      local_l.FinerSize() = 0;

    }

    grids.push_back(new Grid(global_separated[i], local_l, extent[i], interface->MaxLevel()-i, this));

  }

  numLevels = grids.size();

  levelMax = interface->MaxLevel();
  levelMin = levelMax - numLevels + 1;

  levelIndex = 0;
  levelCurrent = levelMax;

}

void Multigrid::SetLevel(int level)
{
  assert(level >= levelMin && level <= levelMax);
  levelCurrent = level;
  levelIndex = levelMax - level;
}

void Multigrid::ToCoarserLevel()
{
  assert(levelCurrent-1 >= levelMin);
  --levelCurrent;
  ++levelIndex;
}

void Multigrid::ToFinerLevel()
{
  assert(levelCurrent+1 <= levelMax);
  ++levelCurrent;
  --levelIndex;
}

void Multigrid::ClearAll()
{
  for (int i=MinLevel(); i<=MaxLevel(); ++i)
    (*this)(i).Clear();
}

// TODO: diff that in case that breaks QP
void Multigrid::ClearAllCoarseLevels()
{
  for (int i=MinLevel(); i<MaxLevel(); ++i)
    (*this)(i).Clear();
}


void Multigrid::SetCoarserDirichletValues()
{
  Index i_c, i_f;

  for (int i=GlobalMaxLevel(); i>MinLevel(); --i) {

    Grid& grid_f = (*this)(i);
    Grid& grid_c = (*this)(i-1);

    i_c.X() = grid_c.Local().BoundaryBegin1().X();
    i_f.X() = grid_f.Local().BoundaryBegin1().X();
    for (i_c.Y()=grid_c.Local().Begin().Y(); i_c.Y()<grid_c.Local().End().Y(); ++i_c.Y())
      for (i_c.Z()=grid_c.Local().Begin().Z(); i_c.Z()<grid_c.Local().End().Z(); ++i_c.Z())
	grid_c(i_c) = grid_f.GetVal(i_f.X(), 2*i_c.Y(), 2*i_c.Z());

    i_c.X() = grid_c.Local().BoundaryBegin2().X();
    i_f.X() = grid_f.Local().BoundaryBegin2().X();
    for (i_c.Y()=grid_c.Local().Begin().Y(); i_c.Y()<grid_c.Local().End().Y(); ++i_c.Y())
      for (i_c.Z()=grid_c.Local().Begin().Z(); i_c.Z()<grid_c.Local().End().Z(); ++i_c.Z())
	grid_c(i_c) = grid_f.GetVal(i_f.X(), 2*i_c.Y(), 2*i_c.Z());

    i_c.Y() = grid_c.Local().BoundaryBegin1().Y();
    i_f.Y() = grid_f.Local().BoundaryBegin1().Y();
    for (i_c.X()=0; i_c.X()<grid_c.Local().SizeTotal().X(); ++i_c.X())
      for (i_c.Z()=grid_c.Local().Begin().Z(); i_c.Z()<grid_c.Local().End().Z(); ++i_c.Z())
	grid_c(i_c) = grid_f.GetVal(2*i_c.X(), i_f.Y(), 2*i_c.Z());

    i_c.Y() = grid_c.Local().BoundaryBegin2().Y();
    i_f.Y() = grid_f.Local().BoundaryBegin2().Y();
    for (i_c.X()=0; i_c.X()<grid_c.Local().SizeTotal().X(); ++i_c.X())
      for (i_c.Z()=grid_c.Local().Begin().Z(); i_c.Z()<grid_c.Local().End().Z(); ++i_c.Z())
	grid_c(i_c) = grid_f.GetVal(2*i_c.X(), i_f.Y(), 2*i_c.Z());

    i_c.Z() = grid_c.Local().BoundaryBegin1().Z();
    i_f.Z() = grid_f.Local().BoundaryBegin1().Z();
    for (i_c.X()=0; i_c.X()<grid_c.Local().SizeTotal().X(); ++i_c.X())
      for (i_c.Y()=0; i_c.Y()<grid_c.Local().SizeTotal().Y(); ++i_c.Y())
	grid_c(i_c) = grid_f.GetVal(2*i_c.X(), 2*i_c.Y(), i_f.Z());

    i_c.Z() = grid_c.Local().BoundaryBegin2().Z();
    i_f.Z() = grid_f.Local().BoundaryBegin2().Z();
    for (i_c.X()=0; i_c.X()<grid_c.Local().SizeTotal().X(); ++i_c.X())
      for (i_c.Y()=0; i_c.Y()<grid_c.Local().SizeTotal().Y(); ++i_c.Y())
	grid_c(i_c) = grid_f.GetVal(2*i_c.X(), 2*i_c.Y(), i_f.Z());

  }
}

Multigrid::~Multigrid()
{
  for (unsigned int i=0; i<grids.size(); ++i)
    delete grids[i];
}
