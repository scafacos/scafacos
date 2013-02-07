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
 * @file   bspline.cpp
 * @author Julian Iseringhausen <isering@ins.uni-bonn.de>
 * @date   Mon Nov 21 13:27:22 2011
 *
 * @brief  B-Splines for molecular dynamics.
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifndef BSPLINE_DEGREE
#error BSPLINE_DEGREE not defined.
#endif

#if BSPLINE_DEGREE < 3 || BSPLINE_DEGREE > 6
#error Please choose a B-Spline degree between 3 and 6
#endif

#include "base/helper.hpp"
#include "base/math.hpp"
#include "units/particle/bspline.hpp"

#define POW(x,y) Helper::pow(x,y)

using namespace VMG;

Particle::BSpline::BSpline(const int& near_field_cells, const vmg_float& h) :
  spline_nom((BSPLINE_DEGREE+1)/2),
  spline_denom((BSPLINE_DEGREE+1)/2),
  potential_nom((BSPLINE_DEGREE+1)/2+1),
  potential_denom((BSPLINE_DEGREE+1)/2+1),
  field_nom((BSPLINE_DEGREE+1)/2),
  field_denom((BSPLINE_DEGREE+1)/2),
  intervals((BSPLINE_DEGREE+1)/2),
  R(near_field_cells*h),
  near_field_cells(near_field_cells)
{
  for (unsigned int i=0; i<intervals.size(); ++i)
    intervals[i] = R * ( -1.0 + 2.0 / static_cast<vmg_float>(BSPLINE_DEGREE) * (i + BSPLINE_DEGREE/2 + 1));

#if BSPLINE_DEGREE == 3

  spline_nom[0] = Polynomial(2, 81.0*POW(R,2), 0.0, -243.0);
  spline_nom[1] = Polynomial(2, 243.0*POW(R,2), -486.0*R, 243.0);

  spline_denom[0] = Polynomial(0, 16.0 * Math::pi * POW(R,5));
  spline_denom[1] = Polynomial(0, 32.0 * Math::pi * POW(R,5));

  potential_nom[0] = Polynomial(5,
				0.0,
				-195.0*POW(R,4),
				0.0,
				270.0*POW(R,2),
				0.0,
				-243.0);

  potential_nom[1] = Polynomial(5,
				2.0     * POW(R,5),
				-405.0  * POW(R,4),
				0.0,
				810.0   * POW(R,2),
				-810.0  * R,
				243.0);

  potential_nom[2] = Polynomial(0, -1.0);

  potential_denom[0] = Polynomial(0, 80.0*POW(R,5));
  potential_denom[1] = Polynomial(0, 160.0*POW(R,5));
  potential_denom[2] = Polynomial(0, 1.0);

  field_nom[0] = Polynomial(5,
			    20.0 * POW(R,5),
			    0.0,
			    0.0,
			    -135.0 * POW(R,2),
			    0.0,
			    243.0);

  field_nom[1] = Polynomial(5,
			    81.0 * POW(R,5),
			    0.0,
			    0.0,
			    -810.0 * POW(R,2),
			    1215.0 * R,
			    -486.0);

  field_denom[0] = Polynomial(3,
			      0.0,
			      0.0,
			      0.0,
			      20.0 * POW(R,5));

  field_denom[1] = Polynomial(3,
			      0.0,
			      0.0,
			      0.0,
			      80.0 * POW(R,5));

  antid = 39.0 / (16.0 * R);

#elif BSPLINE_DEGREE == 6

  spline_nom[0] = Polynomial(5,
                             297.0     * POW(R,5),
                             0.0,
                             -2430.0   * POW(R,3),
                             0.0,
                             10935.0   * R,
                             -10935.0);

  spline_nom[1] = Polynomial(5,
                             459.0    * POW(R,5),
                             2025.0   * POW(R,4),
                             -17010.0 * POW(R,3),
                             36450.0  * POW(R,2),
                             -32805.0 * R,
                             10935.0);

  spline_nom[2] = Polynomial(5,
                             2187.0   * POW(R,5),
                             -10935.0 * POW(R,4),
                             21870.0  * POW(R,3),
                             -21870.0 * POW(R,2),
                             10935.0  * R,
                             -2187.0);

  spline_denom[0] = Polynomial(0, 20.0 * Math::pi * POW(R,8));
  spline_denom[1] = Polynomial(0, 40.0 * Math::pi * POW(R,8));
  spline_denom[2] = Polynomial(0, 40.0 * Math::pi * POW(R,8));

  potential_nom[0] = Polynomial(8,
                                0.0,
                                -956.0   * POW(R,7),
                                0.0,
                                2772.0   * POW(R,5),
                                0.0,
                                -6804.0  * POW(R,3),
                                0.0,
                                14580.0  * R,
                                -10935.0);

  potential_nom[1] = Polynomial(8,
                                -5.0      * POW(R,8),
                                -5676.0   * POW(R,7),
                                0.0,
                                12852.0   * POW(R,5),
                                28350.0   * POW(R,4),
                                -142884.0 * POW(R,3),
                                204120.0  * POW(R,2),
                                -131220.0 * R,
                                32805.0);

  potential_nom[2] = Polynomial(8,
                                169.0    * POW(R,8),
                                -2916.0  * POW(R,7),
                                0.0,
                                20412.0  * POW(R,5),
                                -51030.0 * POW(R,4),
                                61236.0  * POW(R,3),
                                -40824.0 * POW(R,2),
                                14580.0  * R,
                                -2187.0);

  potential_nom[3]  = Polynomial(0, -1.0);

  potential_denom[0] = Polynomial(0, 280.0 * POW(R,8));
  potential_denom[1] = Polynomial(0, 1680.0 * POW(R,8));
  potential_denom[2] = Polynomial(0, 560.0 * POW(R,8));
  potential_denom[3] = Polynomial(0, 0.0);

  field_nom[0] = Polynomial(8,
			    280.0    * POW(R,8),
			    0.0,
			    0.0,
			    -5544.0  * POW(R,5),
			    0.0,
			    27216.0  * POW(R,3),
			    0.0,
			    -87480.0 * R,
			    76545.0);

  field_nom[1] = Polynomial(8,
			    1675.0 * POW(R,8),
			    0.0,
			    0.0,
			    -25704.0 * POW(R,5),
			    -85050.0 * POW(R,4),
			    571536.0 * POW(R,3),
			    -1020600.0 * POW(R,2),
			    787320.0 * R,
			    -229635.0);

  field_nom[2] = Polynomial(8,
			    729.0 * POW(R,8),
			    0.0,
			    0.0,
			    -40824.0 * POW(R,5),
			    153090.0 * POW(R,4),
			    -244944.0 * POW(R,3),
			    204120.0 * POW(R,2),
			    -87480.0 * R,
			    15309.0);

  field_denom[0] = Polynomial(3,
			      0.0,
			      0.0,
			      0.0,
			      280.0 * POW(R,8));
  field_denom[1] = Polynomial(3,
			      0.0,
			      0.0,
			      0.0,
			      1680.0 * POW(R,8));

  field_denom[2] = Polynomial(3,
			      0.0,
			      0.0,
			      0.0,
			      560.0 * POW(R,8));

  antid = 239.0 / (70.0 * R);

#else
#error B-Spline degree not supported. Choose 3 or 6.
#endif /* BSPLINE_DEGREE */
}
