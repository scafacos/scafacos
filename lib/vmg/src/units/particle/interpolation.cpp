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

#include <vector>

#include "base/helper.hpp"
#include "base/index.hpp"
#include "base/vector.hpp"
#include "grid/grid.hpp"
#include "units/particle/interpolation.hpp"
#include "units/particle/particle.hpp"

#include "mg.hpp"
#include "comm/comm.hpp"

using namespace VMG;

Particle::Interpolation::Interpolation(const int& degree) :
  deg(degree),
  deg_1(degree+1),
  buffer(degree+1),
  buffer_diff(degree)
{
  coeff = new vmg_float[Helper::intpow(deg_1, 3)];
  coeff_buffer = new vmg_float[deg_1];
  for (int i=0; i<degree; ++i)
    buffer_diff[i].resize(i+1);
}

Particle::Interpolation::~Interpolation()
{
  delete [] coeff;
  delete [] coeff_buffer;
}

void Particle::Interpolation::ComputeCoefficients(const Grid& grid, const Index& index)
{
  Index i;

  const Index begin = index - deg/2;

  h = grid.Extent().MeshWidth();

  for (i[0]=0; i[0]<deg_1; ++i[0])
    for (i[1]=0; i[1]<deg_1; ++i[1])
      for (i[2]=0; i[2]<deg_1; ++i[2])
        _access_coeff(i) = grid.GetVal(i+begin);

  pos_begin = grid.Extent().Begin()
    + (begin - grid.Local().Begin() + grid.Global().LocalBegin()) * grid.Extent().MeshWidth();

  // compute coefficients x-direction
  for (i=0; i[1]<deg_1; ++i[1])
    for (i[2]=0; i[2]<deg_1; ++i[2])
      _compute_coefficients_1d(i, 0);

  // compute coefficients y-direction
  for (i=0; i[0]<deg_1; ++i[0])
    for (i[2]=0; i[2]<deg_1; ++i[2])
      _compute_coefficients_1d(i, 1);

  // compute coefficients z-direction
  for (i=0; i[0]<deg_1; ++i[0])
    for (i[1]=0; i[1]<deg_1; ++i[1])
      _compute_coefficients_1d(i, 2);
}

void Particle::Interpolation::_compute_coefficients_1d(const Index& index, const unsigned int& direction)
{
  vmg_float power = 1.0;
  unsigned long faculty = 1;
  int c;
  Index i;

  for (i=index, c=0; c<deg_1; ++i[direction], ++c)
    coeff_buffer[c] = _access_coeff(i);

  i=index;
  ++i[direction];
  for (c=1; c<deg_1; ++i[direction], ++c) {
    for (int j=0; j<deg_1-c; ++j)
      coeff_buffer[j] = coeff_buffer[j+1] - coeff_buffer[j];
    faculty *= c;
    power *= h[direction];
    _access_coeff(i) = coeff_buffer[0] / (faculty*power);
  }
}

void Particle::Interpolation::Evaluate(Particle& p)
{
  const Vector& pos = p.Pos();
  vmg_float& pot = p.Pot();
  Vector& field = p.Field();

  pot = 0.0;
  field = 0.0;

  Vector offset = pos - pos_begin;
  buffer[0] = 1.0;
  for (int i=0; i<deg; ++i) {
    buffer[i+1] = buffer[i] * offset;
    for (int j=0; j<i; ++j)
      buffer_diff[i][j] = buffer_diff[i-1][j] * offset;
    buffer_diff[i][i] = buffer[i];
    offset -= h;
  }

  for (int i=1; i<deg; ++i)
    for (int j=1; j<=i; ++j)
      buffer_diff[i][0] += buffer_diff[i][j];

  for (int i=0; i<deg_1; ++i)
    for (int j=0; j<deg_1; ++j)
      for (int k=0; k<deg_1; ++k)
 	pot += _access_coeff(i,j,k) * buffer[i][0] * buffer[j][1] * buffer[k][2];

  for (int i=0; i<deg_1; ++i)
    for (int j=0; j<deg_1; ++j)
      for (int k=0; k<deg; ++k) {
	field[0] -= _access_coeff(k+1, i, j) * buffer[i][1] * buffer[j][2] * buffer_diff[k][0][0];
	field[1] -= _access_coeff(i, k+1, j) * buffer[i][0] * buffer[j][2] * buffer_diff[k][0][1];
	field[2] -= _access_coeff(i, j, k+1) * buffer[i][0] * buffer[j][1] * buffer_diff[k][0][2];
      }
}

vmg_float Particle::Interpolation::EvaluatePotentialLR(const Particle& p)
{
  vmg_float result = 0.0;
  Vector prod, offset;
  Index i;

  const Vector& pos = p.Pos();

  prod[0] = 1.0;
  offset[0] = pos[0] - pos_begin[0];
  for (i[0]=0; i[0]<deg_1; ++i[0]) {
    prod[1] = prod[0];
    offset[1] = pos[1] - pos_begin[1];
    for (i[1]=0; i[1]<deg_1; ++i[1]) {
      prod[2] = prod[1];
      offset[2] = pos[2] - pos_begin[2];
      for (i[2]=0; i[2]<deg_1; ++i[2]) {
        result += _access_coeff(i) * prod[2];
        prod[2] *= offset[2];
        offset[2] -= h[2];
      }
      prod[1] *= offset[1];
      offset[1] -= h[1];
    }
    prod[0] *= offset[0];
    offset[0] -= h[0];
  }

  return result;
}
