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
 * @file   com_initialize_residual_norm.cpp
 * @author Julian Iseringhausen <isering@ins.uni-bonn.de>
 * @date   Mon Apr 18 12:40:12 2011
 *
 * @brief  Computes the residual in the discrete L2-norm
 *         and hands the value over to the factory.
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "base/command.hpp"
#include "base/object.hpp"
#include "comm/comm.hpp"
#include "grid/multigrid.hpp"
#include "mg.hpp"

using namespace VMG;

class VMGCommandInitializeResidualNorm : public Command
{
public:
  Request Run(Command::argument_vector arguments)
  {
    MPE_EVENT_BEGIN()

    vmg_float residual = MG::GetComm()->ComputeResidualNorm(*MG::GetSol(), *MG::GetRhs());
    new ObjectStorage<vmg_float>(arguments[0], residual);

#ifdef DEBUG_OUTPUT
    MG::GetComm()->PrintStringOnce("Initial residual: %e", residual);
#endif /* DEBUG_OUTPUT */

    MPE_EVENT_END()

    return Continue;
  }

  static const char* Name() {return "InitializeResidualNorm";}
  static int Arguments() {return 1;}
};

CREATE_INITIALIZER(VMGCommandInitializeResidualNorm)
