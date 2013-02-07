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
 * @file   interface_particles.cpp
 * @author Julian Iseringhausen <isering@ins.uni-bonn.de>
 * @date   Mon Apr 18 12:56:48 2011
 *
 * @brief  VMG::InterfaceParticles
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef HAVE_MPI
#include <mpi.h>
#ifdef HAVE_MARMOT
#include <enhancempicalls.h>
#include <sourceinfompicalls.h>
#endif
#endif

#include <algorithm>
#include <cmath>
#include <cstring>

#include "base/helper.hpp"
#include "base/index.hpp"
#include "base/math.hpp"
#include "base/vector.hpp"
#include "comm/comm.hpp"
#include "grid/grid.hpp"
#include "grid/multigrid.hpp"
#include "grid/tempgrid.hpp"
#include "units/particle/comm_mpi_particle.hpp"
#include "units/particle/interface_particles.hpp"
#include "units/particle/interpolation.hpp"
#include "units/particle/linked_cell_list.hpp"
#include "mg.hpp"

using namespace VMG;

void InterfaceParticles::ImportRightHandSide(Multigrid& multigrid)
{
  Index index_global, index_local, index;
  Vector pos_rel, pos_abs, grid_val;

  Factory& factory = MG::GetFactory();
  Particle::CommMPI& comm = *dynamic_cast<Particle::CommMPI*>(MG::GetComm());

  const int& near_field_cells = factory.GetObjectStorageVal<int>("PARTICLE_NEAR_FIELD_CELLS");

  Grid& grid = multigrid(multigrid.MaxLevel());
  Grid& particle_grid = comm.GetParticleGrid();

  particle_grid.Clear();

  assert(particle_grid.Global().LocalSize().IsComponentwiseGreater(near_field_cells));

  /*
   * Distribute particles to their processes
   */
  particles.clear();
  comm.CommParticles(grid, particles);

  /*
   * Charge assignment on the grid
   */
  std::list<Particle::Particle>::iterator iter;

#ifdef DEBUG_OUTPUT
  vmg_float particle_charges = 0.0;
  for (iter=particles.begin(); iter!=particles.end(); ++iter)
    particle_charges += iter->Charge();
  particle_charges = MG::GetComm()->GlobalSumRoot(particle_charges);
  comm.PrintStringOnce("Particle list charge sum: %e", particle_charges);
  comm.PrintString("Local number of particles: %d", particles.size());
#endif

  for (iter=particles.begin(); iter!=particles.end(); ++iter)
    spl.SetSpline(particle_grid, *iter);

  // Communicate charges over halo
  comm.CommFromGhosts(particle_grid);

  // Assign charge values to the right hand side
  for (int i=0; i<grid.Local().Size().X(); ++i)
    for (int j=0; j<grid.Local().Size().Y(); ++j)
      for (int k=0; k<grid.Local().Size().Z(); ++k)
	grid(grid.Local().Begin().X() + i,
	     grid.Local().Begin().Y() + j,
	     grid.Local().Begin().Z() + k) = 4.0 * Math::pi *
	  particle_grid.GetVal(particle_grid.Local().Begin().X() + i,
			       particle_grid.Local().Begin().Y() + j,
			       particle_grid.Local().Begin().Z() + k);

#ifdef DEBUG_OUTPUT
  Grid::iterator grid_iter;
  vmg_float charge_sum = 0.0;
  for (grid_iter=grid.Iterators().Local().Begin(); grid_iter!=grid.Iterators().Local().End(); ++grid_iter)
    charge_sum += grid.GetVal(*grid_iter);
  charge_sum = MG::GetComm()->GlobalSum(charge_sum);
  comm.PrintStringOnce("Grid charge sum: %e", charge_sum);
#endif
}

void InterfaceParticles::ExportSolution(Grid& grid)
{
  Index i;

#ifdef DEBUG_OUTPUT
  vmg_float e = 0.0;
  vmg_float e_long = 0.0;
  vmg_float e_self = 0.0;
  vmg_float e_short_peak = 0.0;
  vmg_float e_short_spline = 0.0;
#endif

  Factory& factory = MG::GetFactory();
  Particle::CommMPI& comm = *dynamic_cast<Particle::CommMPI*>(MG::GetComm());

  /*
   * Get parameters and arrays
   */
  const vmg_int& near_field_cells = factory.GetObjectStorageVal<int>("PARTICLE_NEAR_FIELD_CELLS");
  const vmg_int& interpolation_degree = factory.GetObjectStorageVal<int>("PARTICLE_INTERPOLATION_DEGREE");

  Particle::Interpolation ip(interpolation_degree);

  const vmg_float r_cut = near_field_cells * grid.Extent().MeshWidth().Max();

  /*
   * Copy potential values to a grid with sufficiently large halo size.
   * This may be optimized in future.
   * The parameters of this grid have been set in the import step.
   */
  Grid& particle_grid = comm.GetParticleGrid();

  for (i.X()=0; i.X()<grid.Local().Size().X(); ++i.X())
    for (i.Y()=0; i.Y()<grid.Local().Size().Y(); ++i.Y())
      for (i.Z()=0; i.Z()<grid.Local().Size().Z(); ++i.Z())
        particle_grid(i + particle_grid.Local().Begin()) = grid.GetVal(i + grid.Local().Begin());

  comm.CommToGhosts(particle_grid);

  /*
   * Compute potentials
   */
  Particle::LinkedCellList lc(particles, near_field_cells, grid);
  Particle::LinkedCellList::iterator p1, p2;
  Grid::iterator iter;

  comm.CommLCListToGhosts(lc);

  for (int i=lc.Local().Begin().X(); i<lc.Local().End().X(); ++i)
    for (int j=lc.Local().Begin().Y(); j<lc.Local().End().Y(); ++j)
      for (int k=lc.Local().Begin().Z(); k<lc.Local().End().Z(); ++k) {

	if (lc(i,j,k).size() > 0)
	  ip.ComputeCoefficients(particle_grid, Index(i,j,k) - lc.Local().Begin() + particle_grid.Local().Begin());

	for (p1=lc(i,j,k).begin(); p1!=lc(i,j,k).end(); ++p1) {

	  // Interpolate long-range part of potential and electric field
	  ip.Evaluate(**p1);

	  // Subtract self-induced potential
	  (*p1)->Pot() -= (*p1)->Charge() * spl.GetAntiDerivativeAtZero();

#ifdef DEBUG_OUTPUT
	  e_long += 0.5 * (*p1)->Charge() * ip.EvaluatePotentialLR(**p1);
	  e_self += 0.5 * (*p1)->Charge() * (*p1)->Charge() * spl.GetAntiDerivativeAtZero();
#endif

	  for (int dx=-1*near_field_cells; dx<=near_field_cells; ++dx)
	    for (int dy=-1*near_field_cells; dy<=near_field_cells; ++dy)
	      for (int dz=-1*near_field_cells; dz<=near_field_cells; ++dz) {

		for (p2=lc(i+dx,j+dy,k+dz).begin(); p2!=lc(i+dx,j+dy,k+dz).end(); ++p2)

		  if (*p1 != *p2) {

		    const Vector dir = (*p1)->Pos() - (*p2)->Pos();
		    const vmg_float length = dir.Length();

		    if (length < r_cut) {

		      (*p1)->Pot() += (*p2)->Charge() / length * (1.0 + spl.EvaluatePotential(length));
		      (*p1)->Field() += (*p2)->Charge() * dir * spl.EvaluateField(length);

#ifdef DEBUG_OUTPUT
		      e_short_peak += 0.5 * (*p1)->Charge() * (*p2)->Charge() / length;
		      e_short_spline += 0.5 * (*p1)->Charge() * (*p2)->Charge() / length * spl.EvaluatePotential(length);
#endif
		    }
		  }
	      }
	}
      }

  /* Remove average force term */
  Vector average_force = 0.0;
  for (std::list<Particle::Particle>::const_iterator iter=particles.begin(); iter!=particles.end(); ++iter)
    average_force += iter->Charge() * iter->Field();
  const vmg_int& npl = MG::GetFactory().GetObjectStorageVal<vmg_int>("PARTICLE_NUM_LOCAL");
  const vmg_int num_particles_global = comm.GlobalSum(npl);
  average_force /= num_particles_global;
  comm.GlobalSumArray(average_force.vec(), 3);
  for (std::list<Particle::Particle>::iterator iter=particles.begin(); iter!=particles.end(); ++iter)
    iter->Field() -= average_force / iter->Charge();

  comm.CommParticlesBack(particles);

#ifdef DEBUG_OUTPUT
  vmg_float* q = factory.GetObjectStorageArray<vmg_float>("PARTICLE_CHARGE_ARRAY");
  const vmg_int& num_particles_local = factory.GetObjectStorageVal<vmg_int>("PARTICLE_NUM_LOCAL");
  const vmg_float* p = factory.GetObjectStorageArray<vmg_float>("PARTICLE_POTENTIAL_ARRAY");
  const vmg_float* f = factory.GetObjectStorageArray<vmg_float>("PARTICLE_FIELD_ARRAY");


  e_long = comm.GlobalSumRoot(e_long);
  e_short_peak = comm.GlobalSumRoot(e_short_peak);
  e_short_spline = comm.GlobalSumRoot(e_short_spline);
  e_self = comm.GlobalSumRoot(e_self);

  for (int j=0; j<num_particles_local; ++j)
    e += 0.5 * p[j] * q[j];
  e = comm.GlobalSumRoot(e);

  comm.PrintStringOnce("E_long:         %e", e_long);
  comm.PrintStringOnce("E_short_peak:   %e", e_short_peak);
  comm.PrintStringOnce("E_short_spline: %e", e_short_spline);
  comm.PrintStringOnce("E_self:         %e", e_self);
  comm.PrintStringOnce("E_total:        %e", e);
  comm.PrintStringOnce("E_total*:       %e", e_long + e_short_peak + e_short_spline - e_self);

#endif /* DEBUG_OUTPUT */

}
