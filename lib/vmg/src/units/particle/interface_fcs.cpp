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
 * @file   interface_fcs.cpp
 * @author Julian Iseringhausen <isering@ins.uni-bonn.de>
 * @date   Mon Apr 18 12:56:20 2011
 *
 * @brief  VMG::InterfaceFCS
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifndef HAVE_MPI
#error MPI is needed to use the Scafacos interface.
#endif

#include <mpi.h>
#ifdef HAVE_MARMOT
#include <enhancempicalls.h>
#include <sourceinfompicalls.h>
#endif

#include "base/object.hpp"
#include "base/timer.hpp"
#include "comm/domain_decomposition_mpi.hpp"
#ifdef DEBUG
#include "comm/mpi/error_handler.hpp"
#endif
#include "cycles/cycle_cs_periodic.hpp"
#include "cycles/cycle_fas_dirichlet.hpp"
#include "level/level_operator_cs.hpp"
#include "level/level_operator_fas.hpp"
#include "level/stencils.hpp"
#include "samples/discretization_poisson_fd.hpp"
#include "samples/discretization_poisson_fv.hpp"
#include "smoother/gsrb_poisson_2.hpp"
#include "smoother/gsrb_poisson_4.hpp"
#include "solver/givens.hpp"
#include "solver/solver_regular.hpp"
#include "solver/solver_singular.hpp"
#include "units/particle/comm_mpi_particle.hpp"
#include "units/particle/interface_fcs.h"
#include "units/particle/interface_particles.hpp"


using namespace VMG;

namespace VMGBackupSettings
{
  static vmg_int level = -1;
  static vmg_int periodic[3] = {-1, -1, -1};
  static vmg_int max_iter = -1;
  static vmg_int smoothing_steps = -1;
  static vmg_int cycle_type = -1;
  static vmg_float precision = -1;
  static vmg_float box_offset[3];
  static vmg_float box_size = -1.0;
  static vmg_int near_field_cells = -1;
  static vmg_int interpolation_degree = -1;
  static vmg_int discretization_order = -1;
  static MPI_Comm mpi_comm;
}

static void VMG_fcs_init(vmg_int level, vmg_int* periodic,vmg_int max_iter,
			 vmg_int smoothing_steps, vmg_int cycle_type, vmg_float precision,
			 vmg_float* box_offset, vmg_float box_size,
			 vmg_int near_field_cells, vmg_int interpolation_degree,
                         vmg_int discretization_order, MPI_Comm mpi_comm)
{
  VMGBackupSettings::level = level;
  std::memcpy(VMGBackupSettings::periodic, periodic, 3*sizeof(vmg_int));
  VMGBackupSettings::max_iter = max_iter;
  VMGBackupSettings::smoothing_steps = smoothing_steps;
  VMGBackupSettings::cycle_type = cycle_type;
  VMGBackupSettings::precision = precision;
  std::memcpy(VMGBackupSettings::box_offset, box_offset, 3*sizeof(vmg_float));
  VMGBackupSettings::box_size = box_size;
  VMGBackupSettings::near_field_cells = near_field_cells;
  VMGBackupSettings::interpolation_degree = interpolation_degree;
  VMGBackupSettings::discretization_order = discretization_order;
  VMGBackupSettings::mpi_comm = mpi_comm;

#ifdef DEBUG
  MPI_Errhandler mpiErrorHandler;
  MPI_Comm_create_errhandler(VMG::MPI::ConvertToException, &mpiErrorHandler);
  MPI_Comm_set_errhandler(mpi_comm, mpiErrorHandler);
#endif

  const Boundary boundary(periodic[0] ? Periodic : Open,
			  periodic[1] ? Periodic : Open,
			  periodic[2] ? Periodic : Open);

  const bool singular = periodic[0] * periodic[1] * periodic[2];

  /*
   * Choose multigrid components
   */
  if (singular) {

    new Particle::CommMPI(boundary, new DomainDecompositionMPI(), mpi_comm);
    new DiscretizationPoissonFD(discretization_order);
    new InterfaceParticles(boundary, 2, level, Vector(box_offset), box_size, near_field_cells, 0, 1.0);
    new LevelOperatorCS(Stencils::RestrictionFullWeight, Stencils::InterpolationTrilinear);
    new Givens<SolverSingular>();
    new CycleCSPeriodic(cycle_type);

  }else {

    new Particle::CommMPI(boundary, new DomainDecompositionMPI(), mpi_comm);
    new DiscretizationPoissonFV(discretization_order);
    new InterfaceParticles(boundary, 2, level, Vector(box_offset), box_size, near_field_cells, 2, 1.6);
    new LevelOperatorFAS(Stencils::RestrictionFullWeight, Stencils::Injection, Stencils::InterpolationTrilinear);
    new Givens<SolverRegular>();
    new CycleFASDirichlet(cycle_type);

  }

  /*
   * Use Gauss-Seidel Red-Black ordering
   */
  if (discretization_order == 2)
    new GaussSeidelRBPoisson2();
  else
    new GaussSeidelRBPoisson4();

  /*
   * Register required parameters
   */
  new ObjectStorage<int>("PRESMOOTHSTEPS", smoothing_steps);
  new ObjectStorage<int>("POSTSMOOTHSTEPS", smoothing_steps);
  new ObjectStorage<vmg_float>("PRECISION", precision);
  new ObjectStorage<int>("MAX_ITERATION", max_iter);
  new ObjectStorage<int>("PARTICLE_NEAR_FIELD_CELLS", near_field_cells);
  new ObjectStorage<int>("PARTICLE_INTERPOLATION_DEGREE", interpolation_degree);

  /*
   * Post init
   */
  MG::PostInit();

  /*
   * Check whether the library is correctly initialized now.
   */
  MG::IsInitialized();
}

void VMG_fcs_setup(vmg_int level, vmg_int* periodic, vmg_int max_iter,
		   vmg_int smoothing_steps, vmg_int cycle_type, vmg_float precision,
		   vmg_float* box_offset, vmg_float box_size,
		   vmg_int near_field_cells, vmg_int interpolation_degree,
                   vmg_int discretization_order, MPI_Comm mpi_comm)
{
  if (VMGBackupSettings::level != level ||
      VMGBackupSettings::periodic[0] != periodic[0] ||
      VMGBackupSettings::periodic[1] != periodic[1] ||
      VMGBackupSettings::periodic[2] != periodic[2] ||
      VMGBackupSettings::max_iter != max_iter ||
      VMGBackupSettings::smoothing_steps != smoothing_steps ||
      VMGBackupSettings::cycle_type != cycle_type ||
      VMGBackupSettings::precision != precision ||
      VMGBackupSettings::box_offset[0] != box_offset[0] ||
      VMGBackupSettings::box_offset[1] != box_offset[1] ||
      VMGBackupSettings::box_offset[2] != box_offset[2] ||
      VMGBackupSettings::box_size != box_size ||
      VMGBackupSettings::near_field_cells != near_field_cells ||
      VMGBackupSettings::interpolation_degree != interpolation_degree ||
      VMGBackupSettings::discretization_order != discretization_order ||
      VMGBackupSettings::mpi_comm != mpi_comm) {

    VMG_fcs_destroy();
    VMG_fcs_init(level, periodic, max_iter,
		 smoothing_steps, cycle_type, precision,
		 box_offset, box_size, near_field_cells,
                 interpolation_degree, discretization_order,
		 mpi_comm);

  }
}

int VMG_fcs_check()
{
  const int& near_field_cells = MG::GetFactory().GetObjectStorageVal<int>("PARTICLE_NEAR_FIELD_CELLS");
  const Multigrid& multigrid = *MG::GetRhs();
  const Grid& grid = multigrid(multigrid.MaxLevel());

  vmg_int error_code = 0;

  if (!grid.Global().LocalSize().IsComponentwiseGreater(near_field_cells))
    error_code = 1;

  return MG::GetComm()->GlobalMax(error_code);
}

void VMG_fcs_run(vmg_float* x, vmg_float* q, vmg_float* p, vmg_float* f, vmg_int num_particles_local)
{
  /*
   * Register parameters for later use.
   */
  new ObjectStorage<vmg_float*>("PARTICLE_POS_ARRAY", x);
  new ObjectStorage<vmg_float*>("PARTICLE_CHARGE_ARRAY", q);
  new ObjectStorage<vmg_float*>("PARTICLE_POTENTIAL_ARRAY", p);
  new ObjectStorage<vmg_float*>("PARTICLE_FIELD_ARRAY", f);
  new ObjectStorage<vmg_int>("PARTICLE_NUM_LOCAL", num_particles_local);

  /*
   * Start the multigrid solver
   */
  MG::Solve();

#ifdef DEBUG_MEASURE_TIME
  Timer::Print();
#endif
}

void VMG_fcs_print_timer()
{
  Timer::Print();
}

void VMG_fcs_destroy(void)
{
  /*
   * Delete all data.
   */
  MG::Destroy();
}
