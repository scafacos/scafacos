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
 * @file   solver_singular.hpp
 * @author Julian Iseringhausen <isering@ins.uni-bonn.de>
 * @date   Mon Apr 18 13:12:18 2011
 *
 * @brief  VMG::SolverSingular
 *
 */

#ifndef SOLVER_SINGULAR_HPP_
#define SOLVER_SINGULAR_HPP_

#include "solver/solver.hpp"

namespace VMG
{

class Multigrid;

class SolverSingular : public Solver
{
public:
  SolverSingular(bool register_ = true) :
    Solver(register_)
  {}

  SolverSingular(int size, bool register_ = true) :
    Solver(size, register_)
  {}

private:
  void AssembleMatrix(const Grid& rhs); ///< Assembles all matrices and vectors.
  void ExportSol(Grid& sol, Grid& rhs); ///< Exports the solution back to a given mesh.
};

}

#endif /* SOLVER_SINGULAR_HPP_ */
