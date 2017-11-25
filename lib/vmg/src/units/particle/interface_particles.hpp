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
 * @file   interface_particles.hpp
 * @author Julian Iseringhausen <isering@ins.uni-bonn.de>
 * @date   Tue Apr 12 17:40:36 2011
 *
 * @brief  Interface for computing forces and energies in
 *         particle systems.
 *
 */

#ifndef INTERFACE_PARTICLES_HPP
#define INTERFACE_PARTICLES_HPP

#include <list>

#include "base/defs.hpp"
#include "base/interface.hpp"
#include "units/particle/bspline.hpp"
#include "units/particle/particle.hpp"

namespace VMG
{

class Grid;
class Multigrid;

class InterfaceParticles : public Interface
{
public:
  InterfaceParticles(const Boundary& boundary, const int& levelMin, const int& levelMax,
		     const Vector& box_offset, const vmg_float& box_size,
		     const int& near_field_cells,
		     const int& coarsening_steps, const vmg_float& alpha) :
    Interface(boundary, levelMin, levelMax, box_offset, box_size, coarsening_steps, alpha),
    spl(near_field_cells, Extent(MaxLevel()).MeshWidth().Max())
  {}

  void ImportRightHandSide(Multigrid& multigrid);
  void ExportSolution(Grid& grid);

protected:
  Particle::BSpline spl;

private:
  std::list<Particle::Particle> particles;
};

}

#endif /* INTERFACE_PARTICLES_HPP */
