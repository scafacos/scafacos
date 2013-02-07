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

/*
 * interface_sinus.hpp
 *
 *  Created on: 31.03.2011
 *      Author: Julian Iseringhausen
 */

#ifndef INTERFACE_SINUS_HPP_
#define INTERFACE_SINUS_HPP_

#include "base/interface.hpp"

namespace VMG
{
class MGGrid;
class MGMultigrid;
}

namespace VMGInterfaces
{

class InterfaceSinus : public VMG::Interface
{
public:
  InterfaceSinus(vmg_float sine_factor,
                 VMG::Boundary boundary, int levelMin, int levelMax,
		 vmg_float box_begin, vmg_float box_end,
		 int coarseningSteps=0, double alpha=1.6) :
    VMG::Interface(boundary, levelMin, levelMax,
		   box_begin, box_end, coarseningSteps, alpha),
    sine_factor(sine_factor)
  {}

  virtual ~InterfaceSinus() {}

  void ImportRightHandSide(VMG::Multigrid& multigrid);
  void ExportSolution(VMG::Grid& grid);

private:
  vmg_float sine_factor;
};

}

#endif /* INTERFACE_SINUS_HPP_ */
