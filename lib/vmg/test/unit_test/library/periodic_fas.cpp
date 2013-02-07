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

/*
 * periodic_cs.cpp
 *
 *  Created on: Oct 20, 2010
 *      Author: Julian Iseringhausen
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include "base/factory.hpp"
#include "base/math.hpp"
#include "comm/comm_serial.hpp"
#include "cycles/cycle_fas_periodic.hpp"
#include "level/level_operator_fas.hpp"
#include "samples/discretization_poisson_fd.hpp"
#include "smoother/gsrb.hpp"
#include "solver/givens.hpp"
#include "solver/solver_singular.hpp"
#include "mg.hpp"

#include "interface_sinus.hpp"

using namespace VMG;

const vmg_float sine_factor = static_cast<vmg_float>(2.0 * Math::pi);

struct LibraryPeriodicFASFixture
{
  LibraryPeriodicFASFixture()
  {
    const Boundary boundary(Periodic, Periodic, Periodic);

    new CommSerial(boundary);
    new VMGInterfaces::InterfaceSinus(sine_factor, boundary, 2, 6, 0.0, 1.0);
    new DiscretizationPoissonFD(2);
    new LevelOperatorFAS(Stencils::RestrictionFullWeight, Stencils::Injection, Stencils::InterpolationTrilinear);
    new GaussSeidelRB();
    new Givens<SolverSingular>();
    new CycleFASPeriodic(2);

    new ObjectStorage<int>("PRESMOOTHSTEPS", 3);
    new ObjectStorage<int>("POSTSMOOTHSTEPS", 3);
    new ObjectStorage<int>("MAX_ITERATION", 7);
    new ObjectStorage<vmg_float>("PRECISION", 1.0e-10);

    MG::PostInit();

    MG::IsInitialized();
  }

  ~LibraryPeriodicFASFixture()
  {
    MG::Destroy();
  }
};

BOOST_FIXTURE_TEST_CASE(LibraryPeriodicFASTest, LibraryPeriodicFASFixture)
{
  MG::Solve();

  double res_init = MG::GetFactory().Get("INITIAL_RESIDUAL")->Cast< ObjectStorage<double> >()->Val();
  double res = MG::GetComm()->ComputeResidualNorm(*MG::GetSol(), *MG::GetRhs());

  BOOST_CHECK_SMALL(res/res_init, 1e-10);
}
