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
 * dirichlet_cs_mpi.cpp
 *
 *  Created on: 25.01.2011
 *      Author: Julian Iseringhausen
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#ifdef HAVE_MPI
#include <mpi.h>
#ifdef HAVE_MARMOT
#include <enhancempicalls.h>
#include <sourceinfompicalls.h>
#endif

#include "base/factory.hpp"
#include "base/math.hpp"
#include "base/vector.hpp"
#include "comm/comm_mpi.hpp"
#include "comm/domain_decomposition_mpi.hpp"
#include "cycles/cycle_cs_dirichlet.hpp"
#include "level/level_operator_cs.hpp"
#include "level/level_operator.hpp"
#include "samples/discretization_poisson_fd.hpp"
#include "smoother/gsrb.hpp"
#include "solver/givens.hpp"
#include "solver/solver_regular.hpp"
#include "mg.hpp"

#include "interface_sinus.hpp"

using namespace VMG;

const vmg_float sine_factor = static_cast<vmg_float>(2.0 * Math::pi);

struct LibraryDirichletCSMPIFixture
{
  LibraryDirichletCSMPIFixture()
  {
    const Boundary boundary(Dirichlet, Dirichlet, Dirichlet);

    new CommMPI(boundary, new DomainDecompositionMPI());
    new VMGInterfaces::InterfaceSinus(sine_factor, boundary, 2, 6, 0.0, 1.0);
    new DiscretizationPoissonFD(2);
    new LevelOperatorCS(Stencils::RestrictionFullWeight, Stencils::InterpolationTrilinear);
    new GaussSeidelRB();
    new Givens<SolverRegular>();
    new CycleCSDirichlet(2);

    new ObjectStorage<int>("PRESMOOTHSTEPS", 3);
    new ObjectStorage<int>("POSTSMOOTHSTEPS", 3);
    new ObjectStorage<int>("MAX_ITERATION", 7);
    new ObjectStorage<vmg_float>("PRECISION", 1.0e-10);

    MG::PostInit();

    MG::IsInitialized();
  }

  ~LibraryDirichletCSMPIFixture()
  {
    MG::Destroy();
  }
};

BOOST_FIXTURE_TEST_CASE(LibraryDirichletCSMPITest, LibraryDirichletCSMPIFixture)
{
  MG::Solve();

  double res_init = MG::GetFactory().Get("INITIAL_RESIDUAL")->Cast< ObjectStorage<double> >()->Val();
  double res = MG::GetComm()->ComputeResidualNorm(*MG::GetSol(), *MG::GetRhs());

  BOOST_CHECK_SMALL(res/res_init, 1e-10);
}

#endif /* HAVE_MPI */
