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
 * @file   grid.cpp
 * @author Julian Iseringhausen <isering@ins.uni-bonn.de>
 * @date   Mon Apr 18 12:53:27 2011
 *
 * @brief  VMG::Grid
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <cstdio>
#include <cmath>
#include <limits>

#include "base/helper.hpp"
#include "base/stencil.hpp"
#include "comm/comm.hpp"
#include "grid/grid.hpp"
#include "grid/tempgrid.hpp"
#include "mg.hpp"

using namespace VMG;

void Grid::InitGrid()
{
  grid = new vmg_float[local.SizeTotal().Product()];
}

Grid::~Grid()
{
  delete [] grid;
}

Grid& Grid::operator=(const Grid& rhs)
{
  if (this != &rhs) {

    global = rhs.Global();
    local = rhs.Local();
    extent = rhs.Extent();
    iterators = rhs.Iterators();
    father = rhs.Father();

    index_translations = rhs.Indexing();
    level = rhs.Level();

    delete [] grid;
    grid = new vmg_float[local.SizeTotal().Product()];

    SetGrid(rhs);

  }

  return *this;
}

void Grid::Clear()
{
  for (int i=0; i<local.SizeTotal().Product(); ++i)
    grid[i] = 0.0;
}

void Grid::ClearInner()
{
  for (Grid::iterator iter=Iterators().Local().Begin(); iter!=Iterators().Local().End(); ++iter)
    (*this)(*iter) = 0.0;
}

void Grid::ClearHalo()
{
  Grid::iterator iter;

  for (int i=0; i<3; ++i) {

    for (iter = Iterators().Halo1()[i].Begin(); iter != Iterators().Halo1()[i].End(); ++iter)
      (*this)(*iter) = 0.0;

    for (iter = Iterators().Halo2()[i].Begin(); iter != Iterators().Halo2()[i].End(); ++iter)
      (*this)(*iter) = 0.0;

  }

}

void Grid::ClearBoundary()
{
  Grid::iterator iter;

  for (int i=0; i<3; ++i) {

    for (iter = Iterators().Boundary1()[i].Begin(); iter != Iterators().Boundary1()[i].End(); ++iter)
      (*this)(*iter) = 0.0;

    for (iter = Iterators().Boundary2()[i].Begin(); iter != Iterators().Boundary2()[i].End(); ++iter)
      (*this)(*iter) = 0.0;

  }

}

void Grid::SetAverageToZero()
{
  Grid::iterator iter;
  vmg_float avg = 0.0;

  for (iter = Iterators().Local().Begin(); iter != Iterators().Local().End(); ++iter)
    avg += GetVal(*iter);

  avg = MG::GetComm()->GlobalSum(avg);
  avg /= Global().GlobalSize().Product();

#ifdef DEBUG_OUTPUT
  MG::GetComm()->PrintStringOnce("Global constraint enforcement: %e", avg);
#endif

  if (std::abs(avg) > std::numeric_limits<vmg_float>::epsilon())
    for (iter = Iterators().Local().Begin(); iter != Iterators().Local().End(); ++iter)
      (*this)(*iter) -= avg;
}

void Grid::ForceDiscreteCompatibilityCondition()
{
  Grid::iterator iter;
  vmg_float val = 0.0;

  for (iter = Iterators().Local().Begin(); iter != Iterators().Local().End(); ++iter)
    val += GetVal(*iter);

  val = MG::GetComm()->GlobalSum(val) / Global().GlobalSize().Product();

  if (std::abs(val) > std::numeric_limits<vmg_float>::epsilon()) {

#ifdef DEBUG_OUTPUT
    MG::GetComm()->PrintStringOnce("WARNING: Right hand side does not satisfy the compatibility condition.");
#endif

    for (iter = Iterators().Local().Begin(); iter != Iterators().Local().End(); ++iter)
      (*this)(*iter) -= val;

#ifdef DEBUG_OUTPUT
    val = 0.0;
    for (iter = Iterators().Local().Begin(); iter != Iterators().Local().End(); ++iter)
      val += GetVal(*iter);
    val = MG::GetComm()->GlobalSumRoot(val);
    MG::GetComm()->PrintStringOnce("Sum of grid charges after forcing the discrete compatibility condition: %e", val);
#endif
  }

}

void Grid::SetGrid(const Grid& rhs)
{
#ifdef DEBUG
  IsCompatible(rhs);
#endif

  std::memcpy(grid, rhs.grid, local.SizeTotal().Product()*sizeof(vmg_float));
}

void Grid::SetBoundary(const Grid& rhs)
{
#ifdef DEBUG
  IsCompatible(rhs);
#endif
  Grid::iterator iter;

  for (int i=0; i<3; ++i) {

    for (iter = Iterators().Boundary1()[i].Begin(); iter != Iterators().Boundary1()[i].End(); ++iter)
      (*this)(*iter) = rhs.GetVal(*iter);

    for (iter = Iterators().Boundary2()[i].Begin(); iter != Iterators().Boundary2()[i].End(); ++iter)
      (*this)(*iter) = rhs.GetVal(*iter);

  }
}

void Grid::AddGrid(const Grid& rhs)
{
#ifdef DEBUG_MATRIX_CHECKS
  IsCompatible(rhs);
#endif

  for (Grid::iterator iter = Iterators().CompleteGrid().Begin(); iter != Iterators().CompleteGrid().End(); ++iter)
    (*this)(*iter) += rhs.GetVal(*iter);
}

void Grid::SubtractGrid(const Grid& rhs)
{
#ifdef DEBUG_MATRIX_CHECKS
  IsCompatible(rhs);
#endif

  for (Grid::iterator iter = Iterators().CompleteGrid().Begin(); iter != Iterators().CompleteGrid().End(); ++iter)
    (*this)(*iter) -= rhs.GetVal(*iter);
}

void Grid::MultiplyScalar(const vmg_float& scalar)
{
  for (Grid::iterator iter = Iterators().CompleteGrid().Begin(); iter != Iterators().CompleteGrid().End(); ++iter)
    (*this)(*iter) *= scalar;
}

void Grid::ApplyStencil(const Stencil& stencil)
{
  Grid::iterator grid_iter;
  Stencil::iterator stencil_iter;
  TempGrid *temp = MG::GetTempGrid();

  temp->SetProperties(*this);
  temp->SetGrid(*this);
  MG::GetComm()->CommToGhosts(*temp);

  for (grid_iter = Iterators().Local().Begin(); grid_iter != Iterators().Local().End(); ++grid_iter) {

	(*this)(*grid_iter) = stencil.GetDiag() * temp->GetVal(*grid_iter);

	for (stencil_iter=stencil.begin(); stencil_iter!=stencil.end(); ++stencil_iter)
	  (*this)(*grid_iter) += stencil_iter->Val() * temp->GetVal(*grid_iter + stencil_iter->Disp());
      }
}

bool Grid::IsCompatible(const Grid& rhs) const
{
  bool eq = true;

  eq &= Helper::IsEq(Local().Begin(), rhs.Local().Begin(), "Local().Begin");
  eq &= Helper::IsEq(Local().End(), rhs.Local().End(), "Local().End");
  eq &= Helper::IsEq(Local().Size(), rhs.Local().Size(), "Local().Size");
  eq &= Helper::IsEq(Local().SizeTotal(), rhs.Local().SizeTotal(), "Local().SizeTotal");

  eq &= Helper::IsEq(Global().LocalBegin(), rhs.Global().LocalBegin(), "Global().LocalBegin");
  eq &= Helper::IsEq(Global().LocalEnd(), rhs.Global().LocalEnd(), "Global().LocalEnd");
  eq &= Helper::IsEq(Global().GlobalSize(), rhs.Global().GlobalSize(), "Global().GlobalSize");
  eq &= Helper::IsEq(Global().LocalSize(), rhs.Global().LocalSize(), "Global().LocalSize");

  eq &= Helper::IsEq(Local().HaloBegin1(), rhs.Local().HaloBegin1(), "Local().HaloBegin1");
  eq &= Helper::IsEq(Local().HaloBegin2(), rhs.Local().HaloBegin2(), "Local().HaloBegin2");
  eq &= Helper::IsEq(Local().HaloEnd1(), rhs.Local().HaloEnd1(), "Local().HaloEnd1");
  eq &= Helper::IsEq(Local().HaloEnd2(), rhs.Local().HaloEnd2(), "Local().HaloEnd2");

  eq &= Helper::IsEq(Extent().MeshWidth(), rhs.Extent().MeshWidth(), "MeshWidth");

  return eq;
}

bool Grid::IsConsistent() const
{
  bool consistent = true;

  for (Grid::iterator iter=Iterators().Local().Begin(); iter!=Iterators().Local().End(); ++iter)
    consistent &= Helper::CheckNumber(GetVal(*iter));

  return consistent;
}
