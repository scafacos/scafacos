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
 * @file   level_operator_fas.hpp
 * @author Julian Iseringhausen <isering@ins.uni-bonn.de>
 * @date   Mon Apr 18 13:00:41 2011
 *
 * @brief  VMG::LevelOperatorFAS
 *
 */

#ifndef LEVEL_OPERATOR_FAS_HPP_
#define LEVEL_OPERATOR_FAS_HPP_

#include "base/has_tempgrids.hpp"
#include "level/level_operator.hpp"

namespace VMG
{

  class LevelOperatorFAS : public LevelOperator, public VMG::HasTempGrids
{
public:
  LevelOperatorFAS(const Stencil& restrictDefect_, const Stencil& restrictSol_,
                   const Stencil& prolongate_, bool register_ = true) :
    LevelOperator(restrictDefect_, prolongate_, register_),
    operatorSol(restrictSol_)
  {}

  void Restrict(Multigrid& sol, Multigrid& rhs);
  void Prolongate(Multigrid& sol, Multigrid& rhs);

private:
  const Stencil& OperatorSol() const {return operatorSol;}

  Stencil operatorSol;
};

}

#endif /* LEVEL_OPERATOR_FAS_HPP_ */
