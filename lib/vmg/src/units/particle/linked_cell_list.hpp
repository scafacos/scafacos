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
 * @file   linked_cell_list.hpp
 * @author Julian Iseringhausen <isering@ins.uni-bonn.de>
 * @date   Mon Nov 21 13:27:22 2011
 *
 * @brief  A linked cell list implementation.
 *
 */

#ifndef LINKED_CELL_LIST_HPP_
#define LINKED_CELL_LIST_HPP_

#include <list>
#include <vector>

#include "base/index.hpp"
#include "grid/is_grid.hpp"
#include "units/particle/particle.hpp"

namespace VMG
{

namespace Particle
{

class LinkedCellList : public IsGrid< std::list<Particle*> >
{
public:
  typedef std::list<Particle*>::iterator iterator;
  typedef std::list<Particle*>::const_iterator const_iterator;
  typedef std::list<Particle*>::reverse_iterator reverse_iterator;
  typedef std::list<Particle*>::const_reverse_iterator const_reverse_iterator;

  LinkedCellList(std::list<Particle>& particles,
		 const int& near_field_cells, const Grid& grid);

  ~LinkedCellList();

  void AddParticle(Particle* particle);
  void AddParticleToHalo(const vmg_float* x, const vmg_float& q);

  void ClearHalo();

  const Index& NearFieldCells() const {return near_field_cells_;}
private:
  Index near_field_cells_;
};

}

}

#endif /* LINKED_CELL_LIST_HPP_ */
