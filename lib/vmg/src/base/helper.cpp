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
 * @file   helper.cpp
 * @author Julian Iseringhausen <isering@ins.uni-bonn.de>
 * @date   Tue Apr  5 21:03:47 2011
 *
 * @brief  Provides various helper functions.
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <cstring>
#include <string>

#include "base/helper.hpp"
#include "base/index.hpp"
#include "base/vector.hpp"
#include "grid/grid.hpp"

using namespace VMG;

char* Helper::GetCharArray(const std::string& str)
{
  int size = str.size();

  char* rval = new char[size+1];

  strncpy(rval, str.c_str(), size);

  rval[size] = '\0';

  return rval;
}

std::string Helper::ReplaceWhitespaces(const char* buffer, const char* replace)
{
  size_t pos;

  std::string str(buffer);

  while ((pos = str.find(' ')) != std::string::npos)
    str.replace(pos, 1, replace);

  return str;
}

vmg_float Helper::InterpolateTrilinear(const Vector& point, const Grid& grid)
{
  vmg_float interpolate_vals[4], grid_vals[8];

  const Index index_global = (point - grid.Extent().Begin()) / grid.Extent().MeshWidth();
  const Index index_local = index_global - grid.Global().LocalBegin() + grid.Local().Begin();
  const Vector coord = (point - grid.Extent().Begin() - index_global * grid.Extent().MeshWidth()) / grid.Extent().MeshWidth();

  grid_vals[0] = grid.GetVal(index_local.X()  , index_local.Y()  , index_local.Z()  );
  grid_vals[1] = grid.GetVal(index_local.X()+1, index_local.Y()  , index_local.Z()  );
  grid_vals[2] = grid.GetVal(index_local.X()  , index_local.Y()+1, index_local.Z()  );
  grid_vals[3] = grid.GetVal(index_local.X()+1, index_local.Y()+1, index_local.Z()  );
  grid_vals[4] = grid.GetVal(index_local.X()  , index_local.Y()  , index_local.Z()+1);
  grid_vals[5] = grid.GetVal(index_local.X()+1, index_local.Y()  , index_local.Z()+1);
  grid_vals[6] = grid.GetVal(index_local.X()  , index_local.Y()+1, index_local.Z()+1);
  grid_vals[7] = grid.GetVal(index_local.X()+1, index_local.Y()+1, index_local.Z()+1);

  for (int i=0; i<4; ++i)
    interpolate_vals[i] = (1.0 - coord.X()) * grid_vals[2*i] + coord.X() * grid_vals[2*i+1];

  for (int i=0; i<2; ++i)
    interpolate_vals[i] = (1.0 - coord.Y()) * interpolate_vals[2*i] + coord.Y() * interpolate_vals[2*i+1];

  return (1.0 - coord.Z()) * interpolate_vals[0] + coord.Z() * interpolate_vals[1];
}
