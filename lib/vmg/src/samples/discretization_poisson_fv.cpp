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
 * @file   discretization_poisson_fv.cpp
 * @author Julian Iseringhausen <isering@ins.uni-bonn.de>
 * @date   Mon Apr 18 13:03:47 2011
 *
 * @brief  Finite volume discretization for the Poisson
 *         equation. Absolutely equivalent to the finite
 *         difference discretization unless you use
 *         hierarchically coarsened grids.
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "samples/discretization_poisson_fv.hpp"

using namespace VMG;

void DiscretizationPoissonFV::SetInnerBoundaryCompute(Grid& sol_f, Grid& rhs_f, Grid& sol_c) const
{
  Index i_c, i_f;

  const Vector h2_inv = 0.5 / sol_f.Extent().MeshWidth();

  const Index b1_c = sol_c.Local().FinerBegin();
  const Index b2_c = sol_c.Local().FinerEnd() - 1;
  const Index b1_f = 0;
  const Index b2_f = rhs_f.Local().SizeTotal() - 1;

  const Index begin_f = sol_f.Local().Begin();
  const Index end_f = sol_f.Local().End();
  const Index begin_c = sol_c.Local().FinerBegin();

  const vmg_float c_1_3 = 1.0 / 3.0;
  const Vector c_2_3_sp = 2.0 / 3.0 * sol_f.Extent().MeshWidth();
  const vmg_float c_4_3 = 4.0 / 3.0;

  //
  // X-direction
  //
  for (i_f.Y()=begin_f.Y(), i_c.Y()=begin_c.Y(); i_f.Y()<end_f.Y(); i_f.Y()+=2, ++i_c.Y())
    for (i_f.Z()=begin_f.Z(), i_c.Z()=begin_c.Z(); i_f.Z()<end_f.Z(); i_f.Z()+=2, ++i_c.Z()) {
      rhs_f(b1_f.X(),i_f.Y(),i_f.Z()) = (sol_c(b1_c.X()-1,i_c.Y(),i_c.Z()) - sol_f(b1_f.X()+1,i_f.Y(),i_f.Z())) * h2_inv.X();
      rhs_f(b2_f.X(),i_f.Y(),i_f.Z()) = (sol_c(b2_c.X()+1,i_c.Y(),i_c.Z()) - sol_f(b2_f.X()-1,i_f.Y(),i_f.Z())) * h2_inv.X();
    }

  for (i_f.Y()=begin_f.Y()+1; i_f.Y()<end_f.Y(); i_f.Y()+=2)
    for (i_f.Z()=begin_f.Z(); i_f.Z()<end_f.Z(); i_f.Z()+=2) {
      rhs_f(b1_f.X(),i_f.Y(),i_f.Z()) = 0.5 * (rhs_f(b1_f.X(),i_f.Y()-1,i_f.Z()) + rhs_f(b1_f.X(),i_f.Y()+1,i_f.Z()));
      rhs_f(b2_f.X(),i_f.Y(),i_f.Z()) = 0.5 * (rhs_f(b2_f.X(),i_f.Y()-1,i_f.Z()) + rhs_f(b2_f.X(),i_f.Y()+1,i_f.Z()));
    }

  for (i_f.Y()=begin_f.Y(); i_f.Y()<end_f.Y(); i_f.Y()+=2)
    for (i_f.Z()=begin_f.Z()+1; i_f.Z()<end_f.Z(); i_f.Z()+=2) {
      rhs_f(b1_f.X(),i_f.Y(),i_f.Z()) = 0.5 * (rhs_f(b1_f.X(),i_f.Y(),i_f.Z()-1) + rhs_f(b1_f.X(),i_f.Y(),i_f.Z()+1));
      rhs_f(b2_f.X(),i_f.Y(),i_f.Z()) = 0.5 * (rhs_f(b2_f.X(),i_f.Y(),i_f.Z()-1) + rhs_f(b2_f.X(),i_f.Y(),i_f.Z()+1));
    }

  for (i_f.Y()=begin_f.Y()+1; i_f.Y()<end_f.Y(); i_f.Y()+=2)
    for (i_f.Z()=begin_f.Z()+1; i_f.Z()<end_f.Z(); i_f.Z()+=2) {

      rhs_f(b1_f.X(),i_f.Y(),i_f.Z()) = 0.25 * (rhs_f(b1_f.X(),i_f.Y()-1,i_f.Z()-1) +
						rhs_f(b1_f.X(),i_f.Y()+1,i_f.Z()-1) +
						rhs_f(b1_f.X(),i_f.Y()-1,i_f.Z()+1) +
						rhs_f(b1_f.X(),i_f.Y()+1,i_f.Z()+1));

      rhs_f(b2_f.X(),i_f.Y(),i_f.Z()) = 0.25 * (rhs_f(b2_f.X(),i_f.Y()-1,i_f.Z()-1) +
						rhs_f(b2_f.X(),i_f.Y()+1,i_f.Z()-1) +
						rhs_f(b2_f.X(),i_f.Y()-1,i_f.Z()+1) +
						rhs_f(b2_f.X(),i_f.Y()+1,i_f.Z()+1));
    }

  for (i_f.Y()=begin_f.Y(); i_f.Y()<end_f.Y(); ++i_f.Y())
    for (i_f.Z()=begin_f.Z(); i_f.Z()<end_f.Z(); ++i_f.Z()) {

      rhs_f(b1_f.X(),i_f.Y(),i_f.Z()) = sol_f(b1_f.X(),i_f.Y(),i_f.Z()) = c_2_3_sp.X() * rhs_f(b1_f.X(),i_f.Y(),i_f.Z()) +
	c_4_3 * sol_f(b1_f.X()+1,i_f.Y(),i_f.Z()) -
	c_1_3 * sol_f(b1_f.X()+2,i_f.Y(),i_f.Z());

      rhs_f(b2_f.X(),i_f.Y(),i_f.Z()) = sol_f(b2_f.X(),i_f.Y(),i_f.Z()) = c_2_3_sp.X() * rhs_f(b2_f.X(),i_f.Y(),i_f.Z()) +
	c_4_3 * sol_f(b2_f.X()-1,i_f.Y(),i_f.Z()) -
	c_1_3 * sol_f(b2_f.X()-2,i_f.Y(),i_f.Z());
    }

  //
  // Y-direction
  //

  for (i_f.X()=begin_f.X(), i_c.X()=begin_c.X(); i_f.X()<end_f.X(); i_f.X()+=2, ++i_c.X())
    for (i_f.Z()=begin_f.Z(), i_c.Z()=begin_c.Z(); i_f.Z()<end_f.Z(); i_f.Z()+=2, ++i_c.Z()) {
      rhs_f(i_f.X(),b1_f.Y(),i_f.Z()) = (sol_c(i_c.X(),b1_c.Y()-1,i_c.Z()) - sol_f(i_f.X(),b1_f.Y()+1,i_f.Z())) * h2_inv.Y();
      rhs_f(i_f.X(),b2_f.Y(),i_f.Z()) = (sol_c(i_c.X(),b2_c.Y()+1,i_c.Z()) - sol_f(i_f.X(),b2_f.Y()-1,i_f.Z())) * h2_inv.Y();
    }

  for (i_f.X()=begin_f.X()+1; i_f.X()<end_f.X(); i_f.X()+=2)
    for (i_f.Z()=begin_f.Z(); i_f.Z()<end_f.Z(); i_f.Z()+=2) {
      rhs_f(i_f.X(),b1_f.Y(),i_f.Z()) = 0.5 * (rhs_f(i_f.X()-1,b1_f.Y(),i_f.Z()) + rhs_f(i_f.X()+1,b1_f.Y(),i_f.Z()));
      rhs_f(i_f.X(),b2_f.Y(),i_f.Z()) = 0.5 * (rhs_f(i_f.X()-1,b2_f.Y(),i_f.Z()) + rhs_f(i_f.X()+1,b2_f.Z(),i_f.Z()));
    }

  for (i_f.X()=begin_f.X(); i_f.X()<end_f.X(); i_f.X()+=2)
    for (i_f.Z()=begin_f.Z()+1; i_f.Z()<end_f.Z(); i_f.Z()+=2) {
      rhs_f(i_f.X(),b1_f.Y(),i_f.Z()) = 0.5 * (rhs_f(i_f.X(),b1_f.Y(),i_f.Z()-1) + rhs_f(i_f.X(),b1_f.Y(),i_f.Z()+1));
      rhs_f(i_f.X(),b2_f.Y(),i_f.Z()) = 0.5 * (rhs_f(i_f.X(),b2_f.Y(),i_f.Z()-1) + rhs_f(i_f.X(),b2_f.Y(),i_f.Z()+1));
    }

  for (i_f.X()=begin_f.X()+1; i_f.X()<end_f.X(); i_f.X()+=2)
    for (i_f.Z()=begin_f.Z()+1; i_f.Z()<end_f.Z(); i_f.Z()+=2) {

      rhs_f(i_f.X(),b1_f.Y(),i_f.Z()) = 0.25 * (rhs_f(i_f.X()-1,b1_f.Y(),i_f.Z()-1) +
						rhs_f(i_f.X()+1,b1_f.Y(),i_f.Z()-1) +
						rhs_f(i_f.X()-1,b1_f.Y(),i_f.Z()+1) +
						rhs_f(i_f.X()+1,b1_f.Y(),i_f.Z()+1));

      rhs_f(i_f.X(),b2_f.Y(),i_f.Z()) = 0.25 * (rhs_f(i_f.X()-1,b2_f.Y(),i_f.Z()-1) +
						rhs_f(i_f.X()+1,b2_f.Y(),i_f.Z()-1) +
						rhs_f(i_f.X()-1,b2_f.Y(),i_f.Z()+1) +
						rhs_f(i_f.X()+1,b2_f.Y(),i_f.Z()+1));
    }

  for (i_f.X()=begin_f.X(); i_f.X()<end_f.X(); ++i_f.X())
    for (i_f.Z()=begin_f.Z(); i_f.Z()<end_f.Z(); ++i_f.Z()) {

      rhs_f(i_f.X(),b1_f.Y(),i_f.Z()) = sol_f(i_f.X(),b1_f.Y(),i_f.Z()) = c_2_3_sp.Y() * rhs_f(i_f.X(),b1_f.Y(),i_f.Z()) +
	c_4_3 * sol_f(i_f.X(),b1_f.Y()+1,i_f.Z()) -
	c_1_3 * sol_f(i_f.X(),b1_f.Y()+2,i_f.Z());

      rhs_f(i_f.X(),b2_f.Y(),i_f.Z()) = sol_f(i_f.X(),b2_f.Y(),i_f.Z()) = c_2_3_sp.Y() * rhs_f(i_f.X(),b2_f.Y(),i_f.Z()) +
	c_4_3 * sol_f(i_f.X(),b2_f.Y()-1,i_f.Z()) -
	c_1_3 * sol_f(i_f.X(),b2_f.Y()-2,i_f.Z());
    }

  //
  // Z-direction
  //

  for (i_f.X()=begin_f.X(), i_c.X()=begin_c.X(); i_f.X()<end_f.X(); i_f.X()+=2, ++i_c.X())
    for (i_f.Y()=begin_f.Y(), i_c.Y()=begin_c.Y(); i_f.Y()<end_f.Y(); i_f.Y()+=2, ++i_c.Y()) {
      rhs_f(i_f.X(),i_f.Y(),b1_f.Z()) = (sol_c(i_c.X(),i_c.Y(),b1_c.Z()-1) - sol_f(i_f.X(),i_f.Y(),b1_f.Z()+1)) * h2_inv.Z();
      rhs_f(i_f.X(),i_f.Y(),b2_f.Z()) = (sol_c(i_c.X(),i_c.Y(),b2_c.Z()+1) - sol_f(i_f.X(),i_f.Y(),b2_f.Z()-1)) * h2_inv.Z();
    }

  for (i_f.X()=begin_f.X()+1; i_f.X()<end_f.X(); i_f.X()+=2)
    for (i_f.Y()=begin_f.Y(); i_f.Y()<end_f.Y(); i_f.Y()+=2) {
      rhs_f(i_f.X(),i_f.Y(),b1_f.Z()) = 0.5 * (rhs_f(i_f.X()-1,i_f.Y(),b1_f.Z()) + rhs_f(i_f.X()+1,i_f.Y(),b1_f.Z()));
      rhs_f(i_f.X(),i_f.Y(),b2_f.Z()) = 0.5 * (rhs_f(i_f.X()-1,i_f.Y(),b2_f.Z()) + rhs_f(i_f.X()+1,i_f.Y(),b2_f.Z()));
    }

  for (i_f.X()=begin_f.X(); i_f.X()<end_f.X(); i_f.X()+=2)
    for (i_f.Y()=begin_f.Y()+1; i_f.Y()<end_f.Y(); i_f.Y()+=2) {
      rhs_f(i_f.X(),i_f.Y(),b1_f.Z()) = 0.5 * (rhs_f(i_f.X(),i_f.Y()-1,b1_f.Z()) + rhs_f(i_f.X(),i_f.Y()+1,b1_f.Z()));
      rhs_f(i_f.X(),i_f.Y(),b2_f.Z()) = 0.5 * (rhs_f(i_f.X(),i_f.Y()-1,b2_f.Z()) + rhs_f(i_f.X(),i_f.Y()+1,b2_f.Z()));
    }

  for (i_f.X()=begin_f.X()+1; i_f.X()<end_f.X(); i_f.X()+=2)
    for (i_f.Y()=begin_f.Y()+1; i_f.Y()<end_f.Y(); i_f.Y()+=2) {
      rhs_f(i_f.X(),i_f.Y(),b1_f.Z()) = 0.25 * (rhs_f(i_f.X()-1,i_f.Y()-1,b1_f.Z()) +
						rhs_f(i_f.X()+1,i_f.Y()-1,b1_f.Z()) +
						rhs_f(i_f.X()-1,i_f.Y()+1,b1_f.Z()) +
						rhs_f(i_f.X()+1,i_f.Y()+1,b1_f.Z()));

      rhs_f(i_f.X(),i_f.Y(),b2_f.Z()) = 0.25 * (rhs_f(i_f.X()-1,i_f.Y()-1,b2_f.Z()) +
						rhs_f(i_f.X()+1,i_f.Y()-1,b2_f.Z()) +
						rhs_f(i_f.X()-1,i_f.Y()+1,b2_f.Z()) +
						rhs_f(i_f.X()+1,i_f.Y()+1,b2_f.Z()));
    }

  for (i_f.X()=begin_f.X(); i_f.X()<end_f.X(); ++i_f.X())
    for (i_f.Y()=begin_f.Y(); i_f.Y()<end_f.Y(); ++i_f.Y()) {

      rhs_f(i_f.X(),i_f.Y(),b1_f.Z()) = sol_f(i_f.X(),i_f.Y(),b1_f.Z()) = c_2_3_sp.Z() * rhs_f(i_f.X(),i_f.Y(),b1_f.Z()) +
	c_4_3 * sol_f(i_f.X(),i_f.Y(),b1_f.Z()+1) -
	c_1_3 * sol_f(i_f.X(),i_f.Y(),b1_f.Z()+2);

      rhs_f(i_f.X(),i_f.Y(),b2_f.Z()) = sol_f(i_f.X(),i_f.Y(),b2_f.Z()) = c_2_3_sp.Z() * rhs_f(i_f.X(),i_f.Y(),b2_f.Z()) +
	c_4_3 * sol_f(i_f.X(),i_f.Y(),b2_f.Z()-1) -
	c_1_3 * sol_f(i_f.X(),i_f.Y(),b2_f.Z()-2);
    }

#ifdef DEBUG_MATRIX_CHECKS
  rhs_f.IsConsistent();
  sol_f.IsConsistent();
#endif
}
