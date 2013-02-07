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

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#ifndef BOOST_TEST_NO_MAIN
#define BOOST_TEST_NO_MAIN
#endif

#ifdef HAVE_BOOST_UNIT_TEST_FRAMEWORK_LIB
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <boost/test/framework.hpp>
#include <boost/test/detail/unit_test_parameters.hpp>
#include <boost/test/impl/unit_test_main.ipp>

#ifdef BOOST_TEST_ALTERNATIVE_INIT_API
bool init_unit_test()
#else
::boost::unit_test::test_suite* init_unit_test_suite(int argc, char* argv[])
#endif
{
  boost::unit_test::framework::master_test_suite().p_name.value = "VMG Test Suite";
#ifdef BOOST_TEST_ALTERNATIVE_INIT_API
  return true;
#else
  return NULL;
#endif
}

int BOOST_TEST_CALL_DECL
main(int argc, char* argv[])
{
#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
#endif

  // prototype for user's unit test init function
#ifdef BOOST_TEST_ALTERNATIVE_INIT_API
  boost::unit_test::init_unit_test_func init_func = &init_unit_test;
#else
  boost::unit_test::init_unit_test_func init_func = &init_unit_test_suite;
#endif

  int rval = ::boost::unit_test::unit_test_main(init_func, argc, argv);

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return rval;
}
