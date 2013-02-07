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
 * @file   level_operator_cs.hpp
 * @author Julian Iseringhausen <isering@ins.uni-bonn.de>
 * @date   Mon Apr 18 13:00:01 2011
 *
 * @brief  VMG::LevelOperatorCS
 *
 */

#ifndef LEVEL_OPERATOR_CS_HPP_
#define LEVEL_OPERATOR_CS_HPP_

#include "level/level_operator.hpp"

namespace VMG
{

class LevelOperatorCS : public LevelOperator
{
public:
  LevelOperatorCS(const Stencil& operatorRestrict_, const Stencil& operatorProlongate_,
                  bool register_ = true) :
    LevelOperator(operatorRestrict_, operatorProlongate_, register_)
  {}

  void Restrict(Multigrid& sol, Multigrid& rhs);
  void Prolongate(Multigrid& sol, Multigrid& rhs);
};

}

#endif /* LEVEL_OPERATOR_CS_HPP_ */
