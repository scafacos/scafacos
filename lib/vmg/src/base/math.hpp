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
 *  @file   math.hpp
 *
 *  @author Julian Iseringhausen
 *
 *  @brief Mathematical constants
 *
 */

#ifndef MATH_HPP_
#define MATH_HPP_

namespace VMG
{

namespace Math
{

  const double e              = 2.7182818284590452354;	/* e */
  const double log2e          = 1.4426950408889634074;	/* log_2 e */
  const double log10e         = 0.43429448190325182765;	/* log_10 e */
  const double ln2            = 0.69314718055994530942;	/* log_e 2 */
  const double ln10           = 2.30258509299404568402;	/* log_e 10 */
  const double pi             = 3.14159265358979323846;	/* pi */
  const double pi_2           = 1.57079632679489661923;	/* pi/2 */
  const double pi_4           = 0.78539816339744830962;	/* pi/4 */
  const double pi_inv         = 0.31830988618379067154;	/* 1/pi */
  const double pi_2_inv       = 0.63661977236758134308;	/* 2/pi */
  const double sqrtpi_2_inv   = 1.12837916709551257390;	/* 2/sqrt(pi) */
  const double sqrt2          = 1.41421356237309504880;	/* sqrt(2) */
  const double sqrt2_inv      = 0.70710678118654752440;	/* 1/sqrt(2) */

}

}

#endif /* MATH_HPP_ */
