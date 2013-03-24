/*
 * interpolate_test.cpp
 *
 *  Created on: 17.04.2012
 *      Author: Julian Iseringhausen
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include "base/index.hpp"
#include "base/math.hpp"
#include "base/vector.hpp"
#include "comm/comm_serial.hpp"
#include "grid/grid.hpp"
#include "grid/multigrid.hpp"
#include "grid/tempgrid.hpp"
#include "samples/discretization_poisson_fd.hpp"
#include "units/particle/interpolation.hpp"
#include "units/particle/particle.hpp"
#include "mg.hpp"

#include "interface_sinus.hpp"

using namespace VMG;

const vmg_float sine_factor = 2.0 * Math::pi;

static inline vmg_float f(const Vector& pos)
{
  return std::sin(sine_factor*pos.X())
       * std::sin(sine_factor*pos.Y())
       * std::sin(sine_factor*pos.Z());
}

static inline Vector grad_f(const Vector& pos)
{
  return -1.0 * sine_factor * Vector(std::cos(sine_factor*pos.X()) * std::sin(sine_factor*pos.Y()) * std::sin(sine_factor*pos.Z()),
                                     std::sin(sine_factor*pos.X()) * std::cos(sine_factor*pos.Y()) * std::sin(sine_factor*pos.Z()),
                                     std::sin(sine_factor*pos.X()) * std::sin(sine_factor*pos.Y()) * std::cos(sine_factor*pos.Z()));
}

struct InterpolateFixture
{
  InterpolateFixture()
  {
    grid = new TempGrid(64, 5, 0.0, 1.0);

    Index i;
    for (i.X() = 0; i.X()<grid->Local().SizeTotal().X(); ++i.X()) {
      for (i.Y() = 0; i.Y()<grid->Local().SizeTotal().Y(); ++i.Y()) {
        for (i.Z() = 0; i.Z()<grid->Local().SizeTotal().Z(); ++i.Z()) {
          const Vector pos = grid->Extent().Begin() + (i - grid->Local().HaloSize1()) * grid->Extent().MeshWidth();
          (*grid)(i) = f(pos);
        }
      }
    }
  }

  ~InterpolateFixture()
  {
    delete grid;
  }

  Grid* grid;
};

BOOST_FIXTURE_TEST_SUITE(InterpolateSuite, InterpolateFixture)

BOOST_AUTO_TEST_CASE(InterpolateTest)
{
  Index i;
  Particle::Interpolation ip(5);
  Particle::Particle p;

  for (i.X() = grid->Local().Begin().X(); i.X()<grid->Local().End().X(); ++i.X())
    for (i.Y() = grid->Local().Begin().Y(); i.Y()<grid->Local().End().Y(); ++i.Y())
      for (i.Z() = grid->Local().Begin().Z(); i.Z()<grid->Local().End().Z(); ++i.Z()) {

        ip.ComputeCoefficients(*grid, i);

        const Vector pos_begin = grid->Extent().Begin()
          + (i - grid->Local().Begin() + grid->Global().LocalBegin()) * grid->Extent().MeshWidth();
        const Vector pos_end = pos_begin + grid->Extent().MeshWidth();

        p.Pos() = pos_begin;

        ip.Evaluate(p);
        BOOST_CHECK_SMALL(p.Pot() - grid->GetVal(i), vmg_float(1.0e-12));

        p.Pos() += 0.5 * grid->Extent().MeshWidth();
        ip.Evaluate(p);
        BOOST_CHECK_SMALL(p.Pot() - f(p.Pos()), vmg_float(1.0e-7));

      }

}

/*
BOOST_AUTO_TEST_CASE(InterpolateFieldTest)
{
  Index i;
  Index i_old = -1;
  Particle::Interpolation ip(5);
  Particle::Particle p;

  vmg_float max_error = 0.0;

  for(p.Pos().X()=grid->Extent().Begin().X(); p.Pos().X()<grid->Extent().End().X(); p.Pos().X() += 0.01 * (grid->Extent().End().X() - grid->Extent().Begin().X()))
    for(p.Pos().Y()=grid->Extent().Begin().Y(); p.Pos().Y()<grid->Extent().End().Y(); p.Pos().Y() += 0.01 * (grid->Extent().End().Y() - grid->Extent().Begin().Y()))
      for(p.Pos().Z()=grid->Extent().Begin().Z(); p.Pos().Z()<grid->Extent().End().Z(); p.Pos().Z() += 0.01 * (grid->Extent().End().Z() - grid->Extent().Begin().Z())) {

        i = (p.Pos()-grid->Extent().Begin())/grid->Extent().MeshWidth()-grid->Global().LocalBegin()+grid->Local().Begin();

        if (i != i_old) {
          ip.ComputeCoefficients(*grid, i);
          i_old = i;
        }

        ip.Evaluate(p);

        for (int j=0; j<3; ++j)
          max_error = std::max(max_error, p.Field()[j] - grad_f(p.Pos())[j]);

      }

  std::printf("Max error: %e\n", max_error);

  p.Pos().X() = 0.25;
  p.Pos().Y() = 0.25;
  for(p.Pos().Z()=0.25; p.Pos().Z()<0.75; p.Pos().Z() += 0.001) {

    i = (p.Pos()-grid->Extent().Begin())/grid->Extent().MeshWidth()-grid->Global().LocalBegin()+grid->Local().Begin();

    if (i != i_old) {
      ip.ComputeCoefficients(*grid, i);
      i_old = i;
    }

    ip.Evaluate(p);
  }
}
*/

BOOST_AUTO_TEST_SUITE_END()
