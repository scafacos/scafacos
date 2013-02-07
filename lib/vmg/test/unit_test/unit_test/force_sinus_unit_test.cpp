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
 * force_sinus.cpp
 *
 *  Created on: 19.06.2012
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

#include "grid/temp_grid.hpp"
#include "units/particle/interpolation.hpp"
#include "units/particle/particle.hpp"

#include "interface_sinus.hpp"

using namespace VMG;

const vmg_float sine_factor = static_cast<vmg_float>(2.0 * Math::pi);

struct LibraryForceSinusFixture
{
  LibraryForceSinusFixture()
  {
    grid = new TempGrid(129, 0, 0.0, 1.0);
  }

  ~LibraryForceSinusFixture()
  {
    delete grid;
  }

  TempGrid* grid;
};

BOOST_FIXTURE_TEST_CASE(LibraryForceSinusTest, LibraryForceSinusFixture)
{
  Particle::Interpolation ip(5);
  Particle::Particle p;
  Vector pos;

  MG::Solve();

  Index i;
  for (i.X()=grid->Local().Begin().X(); i.X()<grid->Local().End().X(); ++i.X())
      for (i.Y()=grid->Local().Begin().Y(); i.Y()<grid->Local().End().Y(); ++i.Y())
        for (i.Z()=grid->Local().Begin().Z(); i.Z()<grid->Local().End().Z(); ++i.Z()) {
          (*grid)(i) = std::sin(sine_factor * (grid->Extent().Begin().X() + (i.X() - grid->Local().Begin().X() + grid->Global().LocalBegin().X()) * grid->Extent().MeshWidth().X())) *
                       std::sin(sine_factor * (grid->Extent().Begin().Y() + (i.Y() - grid->Local().Begin().Y() + grid->Global().LocalBegin().Y()) * grid->Extent().MeshWidth().Y())) *
                       std::sin(sine_factor * (grid->Extent().Begin().Z() + (i.Z() - grid->Local().Begin().Z() + grid->Global().LocalBegin().Z()) * grid->Extent().MeshWidth().Z()));
        }

  for (pos.X()=0.25; pos.X()<=0.75; pos.X()+=0.01)
    for (pos.Y()=0.25; pos.Y()<=0.75; pos.Y()+=0.01)
      for (pos.Z()=0.25; pos.Z()<=0.75; pos.Z()+=0.01) {

	p.Pos() = pos;

	const Index index = (p.Pos()-grid->.Extent().Begin())/grid->.Extent().MeshWidth()-sol.Global().LocalBegin()+sol.Local().Begin();

	ip.ComputeCoefficients(*grid, index);
	ip.Evaluate(p);

	Vector ref = Vector(std::cos(sine_factor*p.Pos().X())*std::sin(sine_factor*p.Pos().Y())*std::sin(sine_factor*p.Pos().Z()),
			    std::sin(sine_factor*p.Pos().X())*std::cos(sine_factor*p.Pos().Y())*std::sin(sine_factor*p.Pos().Z()),
			    std::sin(sine_factor*p.Pos().X())*std::sin(sine_factor*p.Pos().Y())*std::cos(sine_factor*p.Pos().Z()));
	ref *= -1.0 * sine_factor;

	for (int i=0; i<3; ++i)
	  BOOST_CHECK_SMALL(p.Field()[i] - ref[i], 0.01);

      }

  pos.X() = 0.25;
  pos.Y() = 0.25;
  for (pos.Z() = 0.25; pos.Z()<=0.75; pos.Z() += 0.001) {
	p.Pos() = pos;

	const Index index = (p.Pos()-sol.Extent().Begin())/sol.Extent().MeshWidth()-sol.Global().LocalBegin()+sol.Local().Begin();

	ip.ComputeCoefficients(sol, index);
	ip.Evaluate(p);

	Vector ref = Vector(std::cos(sine_factor*p.Pos().X())*std::sin(sine_factor*p.Pos().Y())*std::sin(sine_factor*p.Pos().Z()),
			    std::sin(sine_factor*p.Pos().X())*std::cos(sine_factor*p.Pos().Y())*std::sin(sine_factor*p.Pos().Z()),
			    std::sin(sine_factor*p.Pos().X())*std::sin(sine_factor*p.Pos().Y())*std::cos(sine_factor*p.Pos().Z()));
	ref *= -1.0 * sine_factor;
        std::cout << pos.Z() << " " << std::abs((ref.Z() - p.Field().Z()) / ref.Z()) << std::endl;
  }
}

#endif /* HAVE_MPI */
