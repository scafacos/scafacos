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
 * @file   solver.hpp
 * @author Julian Iseringhausen <isering@ins.uni-bonn.de>
 * @date   Mon Apr 18 13:11:14 2011
 *
 * @brief  VMG::Solver
 *
 */

#ifndef SOLVER_HPP_
#define SOLVER_HPP_

#include <vector>

#include "base/discretization.hpp"
#include "base/object.hpp"
#include "base/stencil.hpp"
#include "comm/comm.hpp"
#include "grid/grid.hpp"
#include "mg.hpp"

namespace VMG
{

class Solver : public Object
{
public:
  Solver(bool register_ = true) :
    Object("SOLVER", register_),
    size(0)
  {}

  Solver(int size, bool register_ = true) :
    Object("SOLVER", register_),
    size(size)
  {
    this->Realloc(size);
  }

  virtual ~Solver() {}

  void Run(Grid& sol, Grid& rhs);

  void Realloc(int n);
  void Realloc(Grid& x);

  vmg_float& Mat(int i, int j)
  {
    return A[j+size*i];
  }

  vmg_float& Rhs(int i)
  {
    return b[i];
  }

  vmg_float& Sol(int i)
  {
    return x[i];
  }

  const int& Size() const {return size;}

protected:
  virtual void Compute() = 0; ///< Solves the system of equations

private:
  virtual void AssembleMatrix(const Grid& rhs) = 0; ///< Assembles all matrices and vectors.
  virtual void ExportSol(Grid& sol, Grid& rhs) = 0; ///< Exports the solution back to a given mesh.

  std::vector<vmg_float> A, b, x;
  int size;
};

}

#endif
