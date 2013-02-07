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

#ifndef INTERPOLATION_HPP_
#define INTERPOLATION_HPP_

#include <vector>

#include "base/index.hpp"
#include "base/vector.hpp"

namespace VMG
{

class Grid;

namespace Particle
{

class Particle;

class Interpolation
{
public:
  Interpolation(const int& degree);
  ~Interpolation();

  void ComputeCoefficients(const Grid& grid, const Index& index);
  void Evaluate(Particle& p);

  vmg_float EvaluatePotentialLR(const Particle& p);

private:
  vmg_float& _access_coeff(const Index& index)
  {
    return coeff[index.Z() + deg_1 * (index.Y() + deg_1 * index.X())];
  }

  vmg_float& _access_coeff(const int& i, const int& j, const int& k)
  {
    return coeff[k + deg_1 * (j + deg_1 * i)];
  }

  void _compute_coefficients_1d(const Index& index, const unsigned int& direction);

  vmg_float* coeff;
  vmg_float* coeff_buffer;
  int deg, deg_1;
  Vector pos_begin;
  Vector h;
  std::vector<Vector> buffer;
  std::vector< std::vector<Vector> > buffer_diff;
};

}

}

#endif /* INTERPOLATION_HPP_ */
