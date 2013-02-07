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
 * @file   solver_regular.hpp
 * @author Julian Iseringhausen <isering@ins.uni-bonn.de>
 * @date   Mon Apr 18 13:11:47 2011
 *
 * @brief  VMG::SolverRegular
 *
 */

#ifndef SOLVER_REGULAR_HPP_
#define SOLVER_REGULAR_HPP_

#include "solver.hpp"

namespace VMG
{

class Multigrid;

class SolverRegular : public Solver
{
public:
  SolverRegular(bool register_ = true) :
    Solver(register_)
  {}

  SolverRegular(int size, bool register_ = true) :
    Solver(size, register_)
  {}

private:
  void AssembleMatrix(const Grid& rhs); ///< Assembles all matrices and vectors.
  void ExportSol(Grid& sol, Grid& rhs); ///< Exports the solution back to a given mesh.
};

}

#endif /* SOLVER_REGULAR_HPP_ */
