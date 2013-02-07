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
 * @file   com_print_residual_norm.cpp
 * @author Julian Iseringhausen <isering@ins.uni-bonn.de>
 * @date   Mon Apr 18 12:47:22 2011
 *
 * @brief  Writes the current discrete L2-norm of the residual
 *         to standard output.
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <sstream>

#include "base/command.hpp"
#include "comm/comm.hpp"
#include "grid/grid.hpp"
#include "grid/multigrid.hpp"
#include "level/level_operator.hpp"
#include "mg.hpp"

using namespace VMG;

class VMGCommandPrintResidualNorm : public Command
{
public:
  Request Run(Command::argument_vector arguments)
  {
    MPE_EVENT_BEGIN()

    Multigrid* sol = MG::GetSol();
    Multigrid* rhs = MG::GetRhs();
    Comm* comm = MG::GetComm();

    if ((*sol)(sol->MaxLevel()).IsActive()) {
      vmg_float residual = comm->ComputeResidualNorm(*sol, *rhs);
      comm->PrintStringOnce("Residual: %e", residual);
    }

    MPE_EVENT_END()

    return Continue;
  }

  static const char* Name() {return "PrintResidualNorm";}
  static int Arguments() {return 0;}
};

CREATE_INITIALIZER(VMGCommandPrintResidualNorm)
