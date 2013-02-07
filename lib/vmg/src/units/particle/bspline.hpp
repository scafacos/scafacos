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
 * @file   bspline.hpp
 * @author Julian Iseringhausen <isering@ins.uni-bonn.de>
 * @date   Mon Nov 21 13:27:22 2011
 *
 * @brief  B-Splines for molecular dynamics.
 *
 */

#ifndef BSPLINE_HPP_
#define BSPLINE_HPP_

#include "base/helper.hpp"
#include "base/index.hpp"
#include "base/polynomial.hpp"
#include "base/vector.hpp"
#include "grid/grid.hpp"
#include "units/particle/particle.hpp"

namespace VMG
{

namespace Particle
{

class BSpline
{
public:
  BSpline(const int& near_field_cells, const vmg_float& h);

  vmg_float EvaluateSpline(const vmg_float& val) const
  {
    for (unsigned int i=0; i<intervals.size(); ++i)
      if (val < intervals[i])
	return spline_nom[i](val) / spline_denom[i](val);
    return 0.0;
  }

  void SetSpline(Grid& grid, const Particle& p) const
  {
    assert(p.Pos().X() >= grid.Extent().Begin().X() && p.Pos().X() < grid.Extent().End().X());
    assert(p.Pos().Y() >= grid.Extent().Begin().Y() && p.Pos().Y() < grid.Extent().End().Y());
    assert(p.Pos().Z() >= grid.Extent().Begin().Z() && p.Pos().Z() < grid.Extent().End().Z());

    vmg_float* vals = new vmg_float[Helper::intpow(2*near_field_cells+1,3)];

    vmg_float temp_val;
    vmg_float int_val = 0.0;
    int c = 0;

    const int index_global_x = (p.Pos().X() - grid.Extent().Begin().X()) / grid.Extent().MeshWidth().X();
    const int index_global_y = (p.Pos().Y() - grid.Extent().Begin().Y()) / grid.Extent().MeshWidth().Y();
    const int index_global_z = (p.Pos().Z() - grid.Extent().Begin().Z()) / grid.Extent().MeshWidth().Z();

    assert(index_global_x >= grid.Global().LocalBegin().X() && index_global_x < grid.Global().LocalEnd().X());
    assert(index_global_y >= grid.Global().LocalBegin().Y() && index_global_y < grid.Global().LocalEnd().Y());
    assert(index_global_z >= grid.Global().LocalBegin().Z() && index_global_z < grid.Global().LocalEnd().Z());

    const int index_local_x = index_global_x - grid.Global().LocalBegin().X() + grid.Local().Begin().X();
    const int index_local_y = index_global_y - grid.Global().LocalBegin().Y() + grid.Local().Begin().Y();
    const int index_local_z = index_global_z - grid.Global().LocalBegin().Z() + grid.Local().Begin().Z();

    assert(index_local_x >= grid.Local().Begin().X() && index_local_x < grid.Local().End().X());
    assert(index_local_y >= grid.Local().Begin().Y() && index_local_y < grid.Local().End().Y());
    assert(index_local_z >= grid.Local().Begin().Z() && index_local_z < grid.Local().End().Z());

    const vmg_float pos_beg_x = p.Pos().X() - grid.Extent().Begin().X() - grid.Extent().MeshWidth().X() * (index_global_x - near_field_cells);
    const vmg_float pos_beg_y = p.Pos().Y() - grid.Extent().Begin().Y() - grid.Extent().MeshWidth().Y() * (index_global_y - near_field_cells);
    const vmg_float pos_beg_z = p.Pos().Z() - grid.Extent().Begin().Z() - grid.Extent().MeshWidth().Z() * (index_global_z - near_field_cells);

    const vmg_float& h_x = grid.Extent().MeshWidth().X();
    const vmg_float& h_y = grid.Extent().MeshWidth().Y();
    const vmg_float& h_z = grid.Extent().MeshWidth().Z();

    // Iterate over all grid points which lie in the support of the interpolating B-Spline
    vmg_float dir_x = pos_beg_x;
    for (int i=-1*near_field_cells; i<=near_field_cells; ++i) {
      vmg_float dir_y = pos_beg_y;
      for (int j=-1*near_field_cells; j<=near_field_cells; ++j) {
	vmg_float dir_z = pos_beg_z;
	for (int k=-1*near_field_cells; k<=near_field_cells; ++k) {

	  // Compute distance from grid point to particle
	  temp_val = EvaluateSpline(std::sqrt(dir_x*dir_x+dir_y*dir_y+dir_z*dir_z));
	  vals[c++] = temp_val * p.Charge();
	  int_val += temp_val;

	  dir_z -= h_z;
	}
	dir_y -= h_y;
      }
      dir_x -= h_x;
    }

    // Reciprocal value of the numerically integrated spline
    int_val = 1.0 / (int_val * h_x * h_y * h_z);

    c = 0;
    for (int i=-1*near_field_cells; i<=near_field_cells; ++i)
      for (int j=-1*near_field_cells; j<=near_field_cells; ++j)
	for (int k=-1*near_field_cells; k<=near_field_cells; ++k)
	  grid(index_local_x + i,
	       index_local_y + j,
	       index_local_z + k) += vals[c++] * int_val;

    delete [] vals;
  }

  vmg_float EvaluatePotential(const vmg_float& val) const
  {
    for (unsigned int i=0; i<intervals.size(); ++i)
      if (val < intervals[i])
	return potential_nom[i](val) / potential_denom[i](val);
    return potential_nom.back()(val) / potential_denom.back()(val);
  }

  vmg_float EvaluateField(const vmg_float& val) const
  {
    for (unsigned int i=0; i<intervals.size(); ++i)
      if (val < intervals[i])
	return field_nom[i](val) / field_denom[i](val);
    return 0.0;
  }

  const vmg_float& GetAntiDerivativeAtZero() const
  {
    return antid;
  }

private:
  std::vector<Polynomial> spline_nom, spline_denom;
  std::vector<Polynomial> potential_nom, potential_denom;
  std::vector<Polynomial> field_nom, field_denom;
  vmg_float antid;
  std::vector<vmg_float> intervals;

  const vmg_float R;
  const int near_field_cells;
};

}

}

#endif /* BSPLINE_HPP_ */
