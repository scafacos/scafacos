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
 * @file   stencil.hpp
 * @author Julian Iseringhausen <isering@ins.uni-bonn.de>
 * @date   Mon Apr 18 12:25:05 2011
 *
 * @brief  VMG::Stencil
 *
 */

#ifndef STENCIL_HPP_
#define STENCIL_HPP_

#include <cstddef>
#include <vector>

#include "grid/grid.hpp"
#include "base/index.hpp"

namespace VMG
{

class Displacement
{
public:
  Displacement() :
    disp(0),
    val(0.0)
  {}

  Displacement(const Index& disp_, const vmg_float& val_) :
    disp(disp_),
    val(val_)
  {}

  const Index& Disp() const {return disp;}
  const vmg_float& Val() const {return val;}

private:
  Index disp;
  vmg_float val;
};

class Stencil
{
public:
  typedef std::vector<Displacement>::const_iterator iterator;

  Stencil(const vmg_float& diag_) :
    diag(diag_)
  {}

  Stencil(const Stencil& stencil_)
  {
    this->diag = stencil_.GetDiag();

    for (Stencil::iterator iter=stencil_.begin(); iter!=stencil_.end(); iter++)
      this->push_back(*iter);
  }

  const vmg_float& GetDiag() const {return diag;}
  void SetDiag(const vmg_float& diag_) {diag = diag_;}

  const Displacement& operator[](const int& index) const {return disp[index];}

  void push_back(const int& x, const int& y, const int& z, const vmg_float& val)
  {
    disp.push_back(Displacement(Index(x,y,z), val));
  }

  void push_back(const Displacement& displacement)
  {
    disp.push_back(displacement);
  }

  iterator begin() const
  {
    return disp.begin();
  }

  iterator end() const
  {
    return disp.end();
  }

  size_t size() const
  {
    return disp.size();
  }

  void clear()
  {
    disp.clear();
  }

  vmg_float Apply(const Grid& grid, const Index& index) const
  {
    vmg_float result = diag * grid.GetVal(index);
    for (Stencil::iterator iter=disp.begin(); iter!=disp.end(); ++iter)
      result += iter->Val() * grid.GetVal(index.X() + iter->Disp().X(),
					  index.Y() + iter->Disp().Y(),
					  index.Z() + iter->Disp().Z());
    return result;
  }

  void Apply(Grid& grid) const;

private:
  std::vector<Displacement> disp;
  vmg_float diag;
};

}

#endif /* STENCIL_HPP_ */
