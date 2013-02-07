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
 * @file   level_operator.hpp
 * @author Julian Iseringhausen <isering@ins.uni-bonn.de>
 * @date   Mon Apr 18 12:59:28 2011
 *
 * @brief  VMG::LevelOperator
 *
 */

#ifndef LEVEL_OPERATOR_HPP_
#define LEVEL_OPERATOR_HPP_

#include "base/object.hpp"
#include "base/stencil.hpp"
#include "level/stencils.hpp"

namespace VMG
{

class Grid;
class Multigrid;

class LevelOperator : public Object
{
public:
  LevelOperator(const Stencil& operatorRestrict_, const Stencil& operatorProlongate_,
                bool register_ = true) :
    Object("LEVEL_OPERATOR", register_),
    operatorRestrict(operatorRestrict_),
    operatorProlongate(operatorProlongate_)
  {}

  virtual void Restrict(Multigrid& sol, Multigrid& rhs) = 0;
  virtual void Prolongate(Multigrid& sol, Multigrid& rhs) = 0;

protected:
  Stencil& OperatorRestrict() {return operatorRestrict;}
  Stencil& OperatorProlongate() {return operatorProlongate;}

private:
  Stencil operatorRestrict, operatorProlongate;
};

}

#endif
