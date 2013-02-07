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
 * @file   comm.hpp
 * @author Julian Iseringhausen <isering@ins.uni-bonn.de>
 * @date   Wed Jun 16 13:21:06 2010
 *
 * @brief  Base class for communication.
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef HAVE_BOOST_FILESYSTEM
#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;
#endif

#include "base/discretization.hpp"
#include "base/helper.hpp"
#include "base/stencil.hpp"
#include "comm/comm.hpp"
#include "comm/domain_decomposition.hpp"
#include "grid/grid.hpp"
#include "mg.hpp"

using namespace VMG;

vmg_float Comm::ComputeResidualNorm(Multigrid& sol_mg, Multigrid& rhs_mg)
{
#ifdef DEBUG_MATRIX_CHECKS
  sol().IsCompatible(rhs());
  sol().IsConsistent();
  rhs().IsConsistent();
#endif

  Stencil::iterator stencil_iter;
  vmg_float norm = 0.0;

  const vmg_float prefactor = MG::GetDiscretization()->OperatorPrefactor(sol_mg());
  const Stencil& A = MG::GetDiscretization()->GetStencil();

  this->CommToGhosts(sol_mg());

  if (sol_mg().Global().BoundaryType() == LocallyRefined)
    MG::GetDiscretization()->SetInnerBoundary(sol_mg(), rhs_mg(), sol_mg(sol_mg.Level()-1));

  const Grid& sol = sol_mg();
  const Grid& rhs = rhs_mg();

  for (int i=rhs.Local().Begin().X(); i<rhs.Local().End().X(); ++i)
    for (int j=rhs.Local().Begin().Y(); j<rhs.Local().End().Y(); ++j)
      for (int k=rhs.Local().Begin().Z(); k<rhs.Local().End().Z(); ++k) {
	vmg_float val = rhs.GetVal(i,j,k) - prefactor * A.GetDiag() * sol.GetVal(i,j,k);
	for (stencil_iter=A.begin(); stencil_iter!=A.end(); ++stencil_iter)
	  val -= prefactor * stencil_iter->Val() * sol.GetVal(i + stencil_iter->Disp().X(),
							      j + stencil_iter->Disp().Y(),
							      k + stencil_iter->Disp().Z());
	norm += val*val;
      }

  norm = GlobalSum(norm);
  norm = std::sqrt(sol.Extent().MeshWidth().Product() * norm);

  return norm;
}

Grid& Comm::GetParticleGrid()
{
  if (particle_grid != NULL)
    return *particle_grid;

  const Multigrid& multigrid = *MG::GetRhs();
  const Grid& grid = multigrid(multigrid.MaxLevel());
  LocalIndices local = grid.Local();

  const int& near_field_cells = MG::GetFactory().GetObjectStorageVal<int>("PARTICLE_NEAR_FIELD_CELLS");

  local.BoundaryBegin1() = 0;
  local.BoundaryEnd1() = 0;
  local.BoundarySize1() = 0;

  local.BoundaryBegin2() = 0;
  local.BoundaryEnd2() = 0;
  local.BoundarySize2() = 0;

  local.FinerBegin() = 0;
  local.FinerEnd() = 0;
  local.FinerSize() = 0;

  local.End() -= local.Begin();
  local.Begin() = 0;

  // Set grid size of intermediate temporary grid
  for (int i=0; i<3; ++i) {

    if (local.HaloSize1()[i] > 0) {
      local.HaloBegin1()[i] = 0;
      local.HaloEnd1()[i] = near_field_cells;
      local.HaloSize1()[i] = near_field_cells;
      local.Begin()[i] = near_field_cells;
      local.End()[i] = local.Begin()[i] + local.Size()[i];
    }

    if (local.HaloSize2()[i] > 0) {
      local.HaloBegin2()[i] = local.End()[i];
      local.HaloEnd2()[i] = local.HaloBegin2()[i] + near_field_cells;
      local.HaloSize2()[i] = near_field_cells;
    }

  }

  local.SizeTotal() = local.Size() + local.HaloSize1() + local.HaloSize2();

  particle_grid = new Grid(grid.Global(), local, grid.Extent());

  return *particle_grid;
}

Comm::~Comm()
{
  delete dd;
  delete particle_grid;
}

const std::string& Comm::OutputPath()
{
  if (!output_directory_is_created) {
    output_path_str = CreateOutputDirectory();
    output_directory_is_created = true;
  }

  return output_path_str;
}
