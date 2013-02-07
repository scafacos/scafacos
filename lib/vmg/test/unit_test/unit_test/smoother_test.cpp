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
 * smoother_test.cpp
 *
 *  Created on: 20.09.2010
 *      Author: Julian Iseringhausen
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include "base/interface.hpp"
#include "base/math.hpp"
#include "comm/comm_serial.hpp"
#include "samples/discretization_poisson_fd.hpp"
#include "smoother/gs.hpp"
#include "smoother/gsrb.hpp"
#include "mg.hpp"

#include "interface_sinus.hpp"

using namespace VMG;

struct SmootherFixture
{
  SmootherFixture()
  {
    const Boundary boundary(Dirichlet, Dirichlet, Dirichlet);

    new CommSerial(boundary);
    new DiscretizationPoissonFD(2);
    new VMGInterfaces::InterfaceSinus(2.0*Math::pi, boundary, 4, 4, 0.0, 1.0);

    MG::PostInit();

    MG::GetInterface()->ImportRightHandSide(*MG::GetRhs());

    gs = new GaussSeidel(false);
    gsrb = new GaussSeidelRB(false);
  }

  ~SmootherFixture()
  {
    MG::Destroy();

    delete gs;
    delete gsrb;
  }

  Smoother* gs;
  Smoother* gsrb;
};

BOOST_FIXTURE_TEST_CASE(SmootherGSLexTest, SmootherFixture)
{
  double norm;
  Multigrid& sol = *MG::GetSol();
  Multigrid& rhs = *MG::GetRhs();
  Comm& comm = *MG::GetComm();

  sol.ClearAll();

  for (int i=0; i<20; ++i) {

    gs->Run(sol, rhs, 50);
    norm = comm.ComputeResidualNorm(sol, rhs);

    if (norm < 1e-10)
      break;
  }

  BOOST_CHECK_SMALL(norm, 1.0e-10);
}

BOOST_FIXTURE_TEST_CASE(SmootherGSRBTest, SmootherFixture)
{
  double norm;

  MG::GetSol()->ClearAll();

  for (int i=0; i<20; ++i) {
    gsrb->Run(*MG::GetSol(), *MG::GetRhs(), 50);
    norm = MG::GetComm()->ComputeResidualNorm(*MG::GetSol(), *MG::GetRhs());
    if (norm < 1e-10)
      break;
  }

  BOOST_CHECK_SMALL(norm, 1e-10);
}
