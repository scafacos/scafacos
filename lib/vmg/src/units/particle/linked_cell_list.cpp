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
 * @file   linked_cell_list.cpp
 * @author Julian Iseringhausen <isering@ins.uni-bonn.de>
 * @date   Mon Nov 21 13:27:22 2011
 *
 * @brief  A linked cell list implementation.
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "comm/comm.hpp"
#include "units/particle/linked_cell_list.hpp"
#include "mg.hpp"

using namespace VMG;

Particle::LinkedCellList::LinkedCellList(std::list<Particle>& particles,
					 const int& near_field_cells, const Grid& grid)
{
  std::list<Particle>::iterator iter;
  Index index_global, index_local;
  LocalIndices local = grid.Local();

  local.BoundaryBegin1() = 0;
  local.BoundaryEnd1() = 0;
  local.BoundarySize1() = 0;

  local.BoundaryBegin2() = 0;
  local.BoundaryEnd2() = 0;
  local.BoundarySize2() = 0;

  local.FinerBegin() = 0;
  local.FinerEnd() = 0;
  local.FinerSize() = 0;

  for (int i=0; i<3; ++i) {

    if (local.HaloSize1()[i] > 0) {
      local.HaloBegin1()[i] = 0;
      local.HaloEnd1()[i] = near_field_cells+1;
      local.HaloSize1()[i] = near_field_cells+1;
      local.Begin()[i] = near_field_cells+1;
      local.End()[i] = local.Begin()[i] + local.Size()[i];
      local.SizeTotal()[i] = local.Size()[i] +
	                     local.HaloEnd1()[i] - local.HaloBegin1()[i] +
	                     local.HaloEnd2()[i] - local.HaloBegin2()[i];
    }else {
      local.Begin()[i] = 0;
      local.End()[i] = local.Size()[i];
    }

    if (local.HaloSize2()[i] > 0) {
      local.HaloBegin2()[i] = local.End()[i];
      local.HaloEnd2()[i] = local.HaloBegin2()[i] + near_field_cells+1;
      local.HaloSize2()[i] = near_field_cells+1;
      local.SizeTotal()[i] = local.Size()[i] +
	                     local.HaloEnd1()[i] - local.HaloBegin1()[i] +
	                     local.HaloEnd2()[i] - local.HaloBegin2()[i];
    }
  }

  SetGridSize(grid.Global(), local, grid.Extent());

  for (iter=particles.begin(); iter!=particles.end(); ++iter)
    AddParticle(&(*iter));
}

Particle::LinkedCellList::~LinkedCellList()
{
  ClearHalo();
}

void Particle::LinkedCellList::AddParticle(Particle* p)
{
  const Index global_index = (p->Pos() - Extent().Begin()) / Extent().MeshWidth();
  const Index local_index = global_index - Global().LocalBegin() + Local().Begin();

  assert(local_index.IsInBounds(Local().Begin(), Local().End()));
  (*this)(local_index).push_back(p);
}

void Particle::LinkedCellList::AddParticleToHalo(const vmg_float* x, const vmg_float& q)
{
  const Index global_index = ((Vector(x) - Extent().Begin()) / Extent().MeshWidth()).Floor();
  const Index local_index = global_index - Global().LocalBegin() + Local().Begin();

  assert(local_index.IsInBounds(0, Local().SizeTotal()));
  assert((local_index[0] >= Local().HaloBegin1()[0] && local_index[0] < Local().HaloEnd1()[0]) ||
	 (local_index[1] >= Local().HaloBegin1()[1] && local_index[1] < Local().HaloEnd1()[1]) ||
	 (local_index[2] >= Local().HaloBegin1()[2] && local_index[2] < Local().HaloEnd1()[2]) ||
	 (local_index[0] >= Local().HaloBegin2()[0] && local_index[0] < Local().HaloEnd2()[0]) ||
	 (local_index[1] >= Local().HaloBegin2()[1] && local_index[1] < Local().HaloEnd2()[1]) ||
	 (local_index[2] >= Local().HaloBegin2()[2] && local_index[2] < Local().HaloEnd2()[2]));

  (*this)(local_index).push_back(new Particle(x, q));
}

void Particle::LinkedCellList::ClearHalo()
{
  Grid::iterator g_iter;
  std::list<Particle*>::iterator p_iter;

  for (int i=0; i<3; ++i) {

    for (g_iter = Iterators().Halo1()[i].Begin(); g_iter != Iterators().Halo1()[i].End(); ++g_iter)
      for (p_iter = (*this)(*g_iter).begin(); p_iter!=(*this)(*g_iter).end(); ++p_iter)
	delete *p_iter;

    for (g_iter = Iterators().Halo2()[i].Begin(); g_iter != Iterators().Halo2()[i].End(); ++g_iter)
      for (p_iter = (*this)(*g_iter).begin(); p_iter!=(*this)(*g_iter).end(); ++p_iter)
	delete *p_iter;

  }
}
