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
 * @file   tempgrid.cpp
 * @author Julian Iseringhausen <isering@ins.uni-bonn.de>
 * @date   Mon Apr 18 12:55:05 2011
 *
 * @brief  VMG::TempGrid
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "base/discretization.hpp"
#include "base/interface.hpp"
#include "base/stencil.hpp"
#include "comm/comm.hpp"
#include "grid/grid_index_translations.hpp"
#include "grid/tempgrid.hpp"
#include "mg.hpp"

using namespace VMG;

void TempGrid::SetProperties(const Grid& rhs)
{
  local = rhs.Local();
  global = rhs.Global();
  extent = rhs.Extent();
  iterators.SetSubgrids(rhs.Local());
  level = rhs.Level();
  Allocate();
}

void TempGrid::SetProperties(const GlobalIndices& global_, const LocalIndices& local_, const SpatialExtent& extent_)
{
  local = local_;
  global = global_;
  extent = extent_;
  iterators.SetSubgrids(local_);
  Allocate();
}

void TempGrid::SetProperties(const Index& size, const Index& halo_size,
                   const Vector& spatial_begin, const Vector& spatial_end)
{
  global.LocalBegin() = 0;
  global.LocalEnd() = size;
  global.LocalSize() = size;
  global.GlobalFinerBegin() = 0;
  global.GlobalFinerEnd() = 0;
  global.GlobalFinerSize() = 0;
  global.LocalFinerBegin() = 0;
  global.LocalFinerEnd() = 0;
  global.LocalFinerSize() = 0;
  global.FinestAbsBegin() = 0;
  global.FinestAbsEnd() = size;
  global.FinestAbsSize() = size;
  global.GlobalSize() = size;
  global.BoundaryType() = BTUndefined;

  local.Begin() = halo_size;
  local.End() = this->local.Begin() + size;
  local.Size() = size;
  local.SizeTotal() = size + 2 * halo_size;
  local.HaloBegin1() = 0;
  local.HaloEnd1() = halo_size;
  local.HaloSize1() = halo_size;
  local.HaloBegin2() = this->local.End();
  local.HaloEnd2() = this->local.HaloBegin2() = halo_size;
  local.HaloSize2() = halo_size;
  local.BoundaryBegin1() = 0;
  local.BoundaryEnd1() = 0;
  local.BoundarySize1() = 0;
  local.BoundaryBegin2() = 0;
  local.BoundaryEnd2() = 0;
  local.BoundarySize2() = 0;
  local.FinerBegin() = 0;
  local.FinerEnd() = 0;
  local.FinerSize() = 0;

  extent.Begin() = spatial_begin;
  extent.End() = spatial_end;
  extent.Size() = spatial_end - spatial_begin;
  extent.MeshWidth() = this->extent.Size() / static_cast<Vector>(size-1);

  Allocate();
}

void TempGrid::SetPropertiesToFiner(const Grid& grid, const Boundary& boundary)
{
  assert(grid.Father() != NULL);
  assert(grid.Level() < grid.Father()->MaxLevel());

  const Grid& grid_finer = (*grid.Father())(grid.Level()+1);
  const Index off = GridIndexTranslations::EndOffset(boundary);

  level = grid.Level() + 1;

  global.GlobalFinerBegin() = 0;
  global.GlobalFinerEnd() = 0;
  global.GlobalFinerSize() = 0;

  global.LocalFinerBegin() = 0;
  global.LocalFinerEnd() = 0;
  global.LocalFinerSize() = 0;

  global.FinestAbsBegin() = grid_finer.Global().FinestAbsBegin();
  global.FinestAbsEnd() = grid_finer.Global().FinestAbsEnd();
  global.FinestAbsSize() = grid_finer.Global().FinestAbsSize();

  global.GlobalSize() = 2 * (grid.Global().GlobalFinerSize() - off) + off;
  global.BoundaryType() = grid_finer.Global().BoundaryType();

  global.LocalBegin() = 2 * grid.Global().LocalBegin().Clamp(grid.Global().GlobalFinerBegin(), grid.Global().GlobalFinerEnd());
  global.LocalEnd() = 2 * grid.Global().LocalEnd().Clamp(grid.Global().GlobalFinerBegin(), grid.Global().GlobalFinerEnd()) - off;

  global.LocalSize() = global.LocalEnd() - global.LocalBegin();
  if (global.LocalSize().Product() == 0) {
    global.LocalBegin() = 0;
    global.LocalEnd() = 0;
    global.LocalSize() = 0;
    global.BoundaryType() = EmptyGrid;
  }

  local.FinerBegin() = 0;
  local.FinerEnd() = 0;
  local.FinerSize() = 0;

  local.Begin() = 0;
  local.End() = global.LocalSize();
  local.Size() = global.LocalSize();
  local.SizeTotal() = global.LocalSize();

  local.HaloBegin1() = 0;
  local.HaloEnd1() = 0;
  local.HaloBegin2() = 0;
  local.HaloEnd2() = 0;

  local.BoundaryBegin1() = 0;
  local.BoundaryEnd1() = 0;
  local.BoundaryBegin2() = 0;
  local.BoundaryEnd2() = 0;

  for (int i=0; i<3; ++i) {

    if (grid.Local().HaloSize1()[i] > 0) {
      local.Begin()[i] += grid.Local().HaloSize1()[i];
      local.End()[i] += grid.Local().HaloSize1()[i];
      local.HaloEnd1()[i] = grid.Local().HaloSize1()[i];
    }

    if (grid.Local().BoundarySize1()[i] > 0) {
      local.Begin()[i] += grid.Local().BoundarySize1()[i];
      local.BoundaryEnd1()[i] = grid.Local().BoundarySize1()[i];
    }

    if (grid.Local().HaloSize2()[i] > 0) {
      local.HaloBegin2()[i] = local.End()[i];
      local.HaloEnd2()[i] = local.End()[i] + grid.Local().HaloSize2()[i];
    }

    if (grid.Local().BoundarySize2()[i] > 0) {
      local.End()[i] -= grid.Local().BoundarySize2()[i];
      local.BoundaryBegin2()[i] = local.End()[i];
      local.BoundaryEnd2()[i] = local.End()[i]+grid.Local().BoundarySize2()[i];
    }

  }

  local.HaloSize1() = local.HaloEnd1() - local.HaloBegin1();
  local.HaloSize2() = local.HaloEnd2() - local.HaloBegin2();
  local.BoundarySize1() = local.BoundaryEnd1() - local.BoundaryBegin1();
  local.BoundarySize2() = local.BoundaryEnd2() - local.BoundaryBegin2();
  local.Size() = local.End() - local.Begin();
  local.SizeTotal() = local.Size() +
    local.HaloSize1() + local.HaloSize2() +
    local.BoundarySize1() + local.BoundarySize2();

  extent.Size() = grid.Extent().Size();
  extent.Begin() = grid.Extent().Begin();
  extent.End() = grid.Extent().End();
  extent.MeshWidth() = 0.5 * grid.Extent().MeshWidth();

  iterators.SetSubgrids(local);

  Allocate();
}

void TempGrid::SetPropertiesToCoarser(const Grid& grid, const Boundary& boundary)
{
  assert(grid.Father() != NULL);
  assert(grid.Level() > grid.Father()->MinLevel());

  const Grid& grid_coarser = (*grid.Father())(grid.Level()-1);
  const Index off = GridIndexTranslations::EndOffset(boundary);

  level = grid.Level() - 1;

  global.GlobalFinerBegin() = grid_coarser.Global().GlobalFinerBegin();
  global.GlobalFinerEnd() = grid_coarser.Global().GlobalFinerEnd();
  global.GlobalFinerSize() = grid_coarser.Global().GlobalFinerSize();

  global.FinestAbsBegin() = grid.Global().FinestAbsBegin();
  global.FinestAbsEnd() = grid.Global().FinestAbsEnd();
  global.FinestAbsSize() = grid.Global().FinestAbsSize();

  global.GlobalSize() = (grid.Global().GlobalSize() - off)/2 + off;
  global.BoundaryType() = grid_coarser.Global().BoundaryType();

  global.LocalBegin() = grid.Global().LocalBegin();
  global.LocalEnd() = grid.Global().LocalEnd();
  GridIndexTranslations::FineToCoarse(global.LocalBegin(), global.LocalEnd());

  global.LocalBegin() += grid_coarser.Global().GlobalFinerBegin();
  global.LocalEnd() += grid_coarser.Global().GlobalFinerBegin();

  global.LocalSize() = global.LocalEnd() - global.LocalBegin();

  if (global.LocalSize().Product() == 0) {
    global.LocalBegin() = 0;
    global.LocalEnd() = 0;
    global.LocalSize() = 0;
    global.BoundaryType() = EmptyGrid;
  }

  global.LocalFinerBegin() = global.LocalBegin();
  global.LocalFinerEnd() = global.LocalEnd();
  global.LocalFinerSize() = global.LocalSize();

  local.SizeTotal() = global.LocalSize();
  local.Size() = global.LocalSize();
  local.Begin() = 0;
  local.End() = global.LocalSize();

  for (int i=0; i<3; ++i) {

    if (grid.Local().HaloSize1()[i] > 0) {
      local.SizeTotal()[i] += grid.Local().HaloSize1()[i];
      local.Begin()[i] += grid.Local().HaloSize1()[i];
      local.End()[i] += grid.Local().HaloSize1()[i];
      local.HaloBegin1()[i] = 0;
      local.HaloEnd1()[i] = grid.Local().HaloSize1()[i];
    }else {
      local.HaloBegin1()[i] = 0;
      local.HaloEnd1()[i] = 0;
    }

    if (grid.Local().BoundarySize1()[i]> 0) {
      local.Size()[i] -= grid.Local().BoundarySize1()[i];
      local.Begin()[i] += grid.Local().BoundarySize1()[i];
      local.BoundaryBegin1()[i] = 0;
      local.BoundaryEnd1()[i] = grid.Local().BoundarySize1()[i];
    }else {
      local.BoundaryBegin1()[i] = 0;
      local.BoundaryEnd1()[i] = 0;
    }

    if (grid.Local().HaloSize2()[i] > 0) {
      local.SizeTotal()[i] += grid.Local().HaloSize2()[i];
      local.HaloBegin2()[i] = local.End()[i];
      local.HaloEnd2()[i] = local.End()[i] + grid.Local().HaloSize2()[i];
    }else {
      local.HaloBegin2()[i] = 0;
      local.HaloEnd2()[i] = 0;
    }

    if (grid.Local().BoundarySize2()[i] > 0) {
      local.Size()[i] -= grid.Local().BoundarySize2()[i];
      local.End()[i] -= grid.Local().BoundarySize2()[i];
      local.BoundaryBegin2()[i] = local.End()[i];
      local.BoundaryEnd2()[i] = local.End()[i] + grid.Local().BoundarySize2()[i];
    }else {
      local.BoundaryBegin2()[i] = 0;
      local.BoundaryEnd2()[i] = 0;
    }

  }

  local.HaloSize1() = local.HaloEnd1() - local.HaloBegin1();
  local.HaloSize2() = local.HaloEnd2() - local.HaloBegin2();
  local.BoundarySize1() = local.BoundaryEnd1() - local.BoundaryBegin1();
  local.BoundarySize2() = local.BoundaryEnd2() - local.BoundaryBegin2();

  local.FinerBegin() = local.Begin();
  local.FinerEnd() = local.End();
  local.FinerSize() = local.Size();

  Extent().Size() = grid.Extent().Size();
  Extent().Begin() = grid.Extent().Begin();
  Extent().End() = grid.Extent().End();
  Extent().MeshWidth() = 2.0 * grid.Extent().MeshWidth();

  iterators.SetSubgrids(local);
  Allocate();
}

void TempGrid::ImportFromResidual(Grid& sol, Grid& rhs)
{
  Grid::iterator iter;

  const vmg_float prefactor = MG::GetDiscretization()->OperatorPrefactor(sol);
  const Stencil& A = MG::GetDiscretization()->GetStencil();

  this->Clear();

  MG::GetComm()->CommToGhosts(sol);

  for (iter=Iterators().Local().Begin(); iter!=Iterators().Local().End(); ++iter)
    (*this)(*iter) = rhs.GetVal(*iter) - prefactor * A.Apply(sol, *iter);
}

void TempGrid::Allocate()
{
  const int size = local.SizeTotal().Product();

  if (size > size_max) {
    size_max = size;
    delete [] grid;
    grid = new vmg_float[size];
  }
}

TempGrid::TempGrid() :
  size_max(0)
{
}

TempGrid::TempGrid(const Grid& rhs) :
  size_max(0)
{
  SetProperties(rhs);
}

TempGrid::TempGrid(const GlobalIndices& global, const LocalIndices& local, const SpatialExtent& extent) :
  size_max(0)
{
  SetProperties(global, local, extent);
}

TempGrid::TempGrid(const Index& size, const Index& halo_size,
                   const Vector& spatial_begin, const Vector& spatial_end) :
  size_max(0)
{
  SetProperties(size, halo_size, spatial_begin, spatial_end);
}

TempGrid::~TempGrid()
{
}
