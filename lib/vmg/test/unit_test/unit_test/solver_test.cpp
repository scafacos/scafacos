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
 * solver_test.cpp
 *
 *  Created on: 21.07.2010
 *      Author: Julian Iseringhausen
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#ifdef HAVE_LAPACK
#include "solver/dgesv.hpp"
#include "solver/dsysv.hpp"
#endif
#include "solver/givens.hpp"
#include "solver/solver_regular.hpp"

using namespace VMG;

#ifdef HAVE_LAPACK
class DGESVTest : public DGESV<SolverRegular>
{
public:
  DGESVTest() :
    DGESV<SolverRegular>(5)
  {
    this->Mat(0, 0) = 1.0;
    this->Mat(0, 1) = -4.0;
    this->Mat(0, 2) = 0.0;
    this->Mat(0, 3) = -5.0;
    this->Mat(0, 4) = 0.0;

    this->Mat(1, 0) = 0.0;
    this->Mat(1, 1) = 2.0;
    this->Mat(1, 2) = 2.0;
    this->Mat(1, 3) = 5.0;
    this->Mat(1, 4) = 4.0;

    this->Mat(2, 0) = 4.0;
    this->Mat(2, 1) = -16.0;
    this->Mat(2, 2) = 0.0;
    this->Mat(2, 3) = -25.0;
    this->Mat(2, 4) = 0.0;

    this->Mat(3, 0) = -1.0;
    this->Mat(3, 1) = -2.0;
    this->Mat(3, 2) = -4.0;
    this->Mat(3, 3) = -5.0;
    this->Mat(3, 4) = -8.0;

    this->Mat(4, 0) = 0.0;
    this->Mat(4, 1) = 0.0;
    this->Mat(4, 2) = 2.0;
    this->Mat(4, 3) = 0.0;
    this->Mat(4, 4) = 5.0;

    this->Rhs(0) = 12.5;
    this->Rhs(1) = 8.0;
    this->Rhs(2) = -10.0;
    this->Rhs(3) = 3.5;
    this->Rhs(4) = -16.0;

    this->Compute();

    for (int i=0; i<5; ++i)
      solution_array[i] = this->Sol(i);
  }

  double solution_array[5];
};

class DSYSVTest : public DSYSV<SolverRegular>
{
public:
  DSYSVTest() :
    DSYSV<SolverRegular>(10)
  {
    this->Mat(0, 0) = 4.0;
    this->Mat(0, 1) = -1.0;
    this->Mat(0, 2) = -1.0;
    this->Mat(0, 3) = -1.0;
    this->Mat(0, 4) = 0.0;
    this->Mat(0, 5) = 0.0;
    this->Mat(0, 6) = -1.0;
    this->Mat(0, 7) = 0.0;
    this->Mat(0, 8) = 0.0;
    this->Mat(0, 9) = 1.0;

    this->Mat(1, 0) = -1.0;
    this->Mat(1, 1) = 4.0;
    this->Mat(1, 2) = -1.0;
    this->Mat(1, 3) = 0.0;
    this->Mat(1, 4) = -1.0;
    this->Mat(1, 5) = 0.0;
    this->Mat(1, 6) = 0.0;
    this->Mat(1, 7) = -1.0;
    this->Mat(1, 8) = 0.0;
    this->Mat(1, 9) = 1.0;

    this->Mat(2, 0) = -1.0;
    this->Mat(2, 1) = -1.0;
    this->Mat(2, 2) = 4.0;
    this->Mat(2, 3) = 0.0;
    this->Mat(2, 4) = 0.0;
    this->Mat(2, 5) = -1.0;
    this->Mat(2, 6) = 0.0;
    this->Mat(2, 7) = 0.0;
    this->Mat(2, 8) = -1.0;
    this->Mat(2, 9) = 1.0;

    this->Mat(3, 0) = -1.0;
    this->Mat(3, 1) = 0.0;
    this->Mat(3, 2) = 0.0;
    this->Mat(3, 3) = 4.0;
    this->Mat(3, 4) = -1.0;
    this->Mat(3, 5) = -1.0;
    this->Mat(3, 6) = -1.0;
    this->Mat(3, 7) = 0.0;
    this->Mat(3, 8) = 0.0;
    this->Mat(3, 9) = 1.0;

    this->Mat(4, 0) = 0.0;
    this->Mat(4, 1) = -1.0;
    this->Mat(4, 2) = 0.0;
    this->Mat(4, 3) = -1.0;
    this->Mat(4, 4) = 4.0;
    this->Mat(4, 5) = -1.0;
    this->Mat(4, 6) = 0.0;
    this->Mat(4, 7) = -1.0;
    this->Mat(4, 8) = 0.0;
    this->Mat(4, 9) = 1.0;

    this->Mat(5, 0) = 0.0;
    this->Mat(5, 1) = 0.0;
    this->Mat(5, 2) = -1.0;
    this->Mat(5, 3) = -1.0;
    this->Mat(5, 4) = -1.0;
    this->Mat(5, 5) = 4.0;
    this->Mat(5, 6) = 0.0;
    this->Mat(5, 7) = 0.0;
    this->Mat(5, 8) = -1.0;
    this->Mat(5, 9) = 1.0;

    this->Mat(6, 0) = -1.0;
    this->Mat(6, 1) = 0.0;
    this->Mat(6, 2) = 0.0;
    this->Mat(6, 3) = -1.0;
    this->Mat(6, 4) = 0.0;
    this->Mat(6, 5) = 0.0;
    this->Mat(6, 6) = 4.0;
    this->Mat(6, 7) = -1.0;
    this->Mat(6, 8) = -1.0;
    this->Mat(6, 9) = 1.0;

    this->Mat(7, 0) = 0.0;
    this->Mat(7, 1) = -1.0;
    this->Mat(7, 2) = 0.0;
    this->Mat(7, 3) = 0.0;
    this->Mat(7, 4) = -1.0;
    this->Mat(7, 5) = 0.0;
    this->Mat(7, 6) = -1.0;
    this->Mat(7, 7) = 4.0;
    this->Mat(7, 8) = -1.0;
    this->Mat(7, 9) = 1.0;

    this->Mat(8, 0) = 0.0;
    this->Mat(8, 1) = 0.0;
    this->Mat(8, 2) = -1.0;
    this->Mat(8, 3) = 0.0;
    this->Mat(8, 4) = 0.0;
    this->Mat(8, 5) = -1.0;
    this->Mat(8, 6) = -1.0;
    this->Mat(8, 7) = -1.0;
    this->Mat(8, 8) = 4.0;
    this->Mat(8, 9) = 1.0;

    this->Mat(9, 0) = 1.0;
    this->Mat(9, 1) = 1.0;
    this->Mat(9, 2) = 1.0;
    this->Mat(9, 3) = 1.0;
    this->Mat(9, 4) = 1.0;
    this->Mat(9, 5) = 1.0;
    this->Mat(9, 6) = 1.0;
    this->Mat(9, 7) = 1.0;
    this->Mat(9, 8) = 1.0;
    this->Mat(9, 9) = 0.0;

    this->Rhs(0) = -0.125;
    this->Rhs(1) = -0.125;
    this->Rhs(2) = -0.125;
    this->Rhs(3) = -0.125;
    this->Rhs(4) = 1.0;
    this->Rhs(5) = -0.125;
    this->Rhs(6) = -0.125;
    this->Rhs(7) = -0.125;
    this->Rhs(8) = -0.125;
    this->Rhs(9) = 0.0;

    this->Compute();

    for (int i=0; i<10; ++i)
      solution_array[i] = this->Sol(i);
  }

  double solution_array[10];
};
#endif /* HAVE_LAPACK */

class GivensTest : public Givens<SolverRegular>
{
public:
  GivensTest() :
    Givens<SolverRegular>(5)
  {
    this->Mat(0, 0) = 1.0;
    this->Mat(0, 1) = -4.0;
    this->Mat(0, 2) = 0.0;
    this->Mat(0, 3) = -5.0;
    this->Mat(0, 4) = 0.0;

    this->Mat(1, 0) = 0.0;
    this->Mat(1, 1) = 2.0;
    this->Mat(1, 2) = 2.0;
    this->Mat(1, 3) = 5.0;
    this->Mat(1, 4) = 4.0;

    this->Mat(2, 0) = 4.0;
    this->Mat(2, 1) = -16.0;
    this->Mat(2, 2) = 0.0;
    this->Mat(2, 3) = -25.0;
    this->Mat(2, 4) = 0.0;

    this->Mat(3, 0) = -1.0;
    this->Mat(3, 1) = -2.0;
    this->Mat(3, 2) = -4.0;
    this->Mat(3, 3) = -5.0;
    this->Mat(3, 4) = -8.0;

    this->Mat(4, 0) = 0.0;
    this->Mat(4, 1) = 0.0;
    this->Mat(4, 2) = 2.0;
    this->Mat(4, 3) = 0.0;
    this->Mat(4, 4) = 5.0;

    this->Rhs(0) = 12.5;
    this->Rhs(1) = 8.0;
    this->Rhs(2) = -10.0;
    this->Rhs(3) = 3.5;
    this->Rhs(4) = -16.0;

    this->Compute();

    for (int i=0; i<5; ++i)
      solution_array[i] = this->Sol(i);
  }

  double solution_array[5];
};

BOOST_AUTO_TEST_SUITE(SolverTestSuite)

BOOST_AUTO_TEST_CASE(SolverDGESVTest)
{
#ifdef HAVE_LAPACK
  DGESVTest* dgesv = new DGESVTest();
  BOOST_CHECK_CLOSE(  8.5, dgesv->solution_array[0], 1.0e-12);
  BOOST_CHECK_CLOSE(-16.0, dgesv->solution_array[1], 1.0e-12);
  BOOST_CHECK_CLOSE(-18.0, dgesv->solution_array[2], 1.0e-12);
  BOOST_CHECK_CLOSE( 12.0, dgesv->solution_array[3], 1.0e-12);
  BOOST_CHECK_CLOSE(  4.0, dgesv->solution_array[4], 1.0e-12);
#endif /* HAVE_LAPACK */
}

BOOST_AUTO_TEST_CASE(SolverDSYSVTest)
{
#ifdef HAVE_LAPACK
  DSYSVTest* dsysv = new DSYSVTest();
  BOOST_CHECK_CLOSE(-0.0625, dsysv->solution_array[0], 1.0e-12);
  BOOST_CHECK_SMALL(         dsysv->solution_array[1], 1.0e-12);
  BOOST_CHECK_CLOSE(-0.0625, dsysv->solution_array[2], 1.0e-12);
  BOOST_CHECK_SMALL(         dsysv->solution_array[3], 1.0e-12);
  BOOST_CHECK_CLOSE(   0.25, dsysv->solution_array[4], 1.0e-12);
  BOOST_CHECK_SMALL(         dsysv->solution_array[5], 1.0e-12);
  BOOST_CHECK_CLOSE(-0.0625, dsysv->solution_array[6], 1.0e-12);
  BOOST_CHECK_SMALL(         dsysv->solution_array[7], 1.0e-12);
  BOOST_CHECK_CLOSE(-0.0625, dsysv->solution_array[8], 1.0e-12);
#endif
}

BOOST_AUTO_TEST_CASE(SolverGivensTest)
{
  GivensTest* givens = new GivensTest();
  BOOST_CHECK_CLOSE(  8.5, givens->solution_array[0], 1.0e-12);
  BOOST_CHECK_CLOSE(-16.0, givens->solution_array[1], 1.0e-12);
  BOOST_CHECK_CLOSE(-18.0, givens->solution_array[2], 1.0e-12);
  BOOST_CHECK_CLOSE( 12.0, givens->solution_array[3], 1.0e-12);
  BOOST_CHECK_CLOSE(  4.0, givens->solution_array[4], 1.0e-12);
}

BOOST_AUTO_TEST_SUITE_END()
