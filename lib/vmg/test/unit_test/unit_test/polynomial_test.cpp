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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "base/polynomial.hpp"

using namespace VMG;

BOOST_AUTO_TEST_SUITE(PolynomialTest)

BOOST_AUTO_TEST_CASE(PolynomialTest1)
{
  Polynomial p(2, 0.0, 0.0, 1.0);
  BOOST_REQUIRE_CLOSE(p(2.0), 4.0, 0.0000001);
}

BOOST_AUTO_TEST_CASE(PolynomialTest2)
{
  Polynomial p(0, 4.0);
  BOOST_REQUIRE_CLOSE(p(0.0), 4.0, 1.0e-16);
  BOOST_REQUIRE_CLOSE(p(1.0), 4.0, 1.0e-16);
  BOOST_REQUIRE_CLOSE(p(1000.0), 4.0, 1.0e-16);
}

BOOST_AUTO_TEST_CASE(PolynomialTest3)
{
  Polynomial p(4, -5.0, 7.0, -3.0, 4.5, -4.0);
  BOOST_REQUIRE_CLOSE(p(1.5), -6.3125, 1.0e-16);
  BOOST_REQUIRE_CLOSE(p(-7.5), -14780.9375, 1.0e-16);
}

BOOST_AUTO_TEST_CASE(PolynomialTestCopyConstructor)
{
  Polynomial p1(3, 2.5, 33.2, 7.0, 3.0);
  Polynomial p2(p1);
  BOOST_REQUIRE_CLOSE(p1(0.0), p2(0.0), 1.0e-16);
  BOOST_REQUIRE_CLOSE(p1(32.3), p2(32.3), 1.0e-16);
}

BOOST_AUTO_TEST_CASE(PolynomialTestAssignment)
{
  Polynomial p1(3, 2.1, 5.3, -2.3, -4.3);
  Polynomial p2;
  p2 = p1;
  BOOST_REQUIRE_CLOSE(p1(0.0), p2(0.0), 1.0e-16);
  BOOST_REQUIRE_CLOSE(p1(32.3), p2(32.3), 1.0e-16);
}

BOOST_AUTO_TEST_SUITE_END()
