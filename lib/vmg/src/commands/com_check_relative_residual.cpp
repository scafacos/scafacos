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
 * @file   com_check_relative_residual.cpp
 * @author Julian Iseringhausen <isering@ins.uni-bonn.de>
 * @date   Mon Apr 18 12:30:18 2011
 *
 * @brief  Checks if the relative residual is small enough
 *         to stop execution.
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <cmath>

#include "base/command.hpp"
#include "base/factory.hpp"
#include "base/object.hpp"
#include "comm/comm.hpp"
#include "grid/multigrid.hpp"
#include "mg.hpp"

using namespace VMG;

class VMGCommandCheckRelativeResidual : public Command
{
public:
  Request Run(Command::argument_vector arguments)
  {
    MPE_EVENT_BEGIN()

    VMG::Factory& factory = MG::GetFactory();

    const vmg_float& res = factory.GetObjectStorageVal<vmg_float>(arguments[0]);
    const vmg_float& init_res = factory.GetObjectStorageVal<vmg_float>(arguments[1]);
    const vmg_float& precision = factory.GetObjectStorageVal<vmg_float>("PRECISION");
    const vmg_float rel_res = std::fabs(res / init_res);

#ifdef DEBUG_OUTPUT
    MG::GetComm()->PrintStringOnce("Relative residual: %e", rel_res);
#endif /* DEBUG_OUTPUT */

    MPE_EVENT_END()

    if (rel_res < precision)
      return StopCycleLater;
    else
      return Continue;
  }

  static const char* Name() {return "CheckRelativeResidual";}
  static int Arguments() {return 2;}
};

CREATE_INITIALIZER(VMGCommandCheckRelativeResidual)

