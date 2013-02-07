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
 * @file   com_check_residual.cpp
 * @author Julian Iseringhausen <isering@ins.uni-bonn.de>
 * @date   Mon Apr 18 12:31:37 2011
 *
 * @brief  Checks if the residual is small enough to stop
 *         execution.
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

class VMGCommandCheckResidual : public Command
{
public:
  Request Run(Command::argument_vector arguments)
  {
    MPE_EVENT_BEGIN()

    const vmg_float& res = MG::GetFactory().GetObjectStorageVal<vmg_float>(arguments[0]);
    const vmg_float& precision = MG::GetFactory().Get("PRECISION")->Cast< ObjectStorage<vmg_float> >()->Val();

    MPE_EVENT_END()

    if (std::fabs(res) < precision)
      return StopCycleLater;
    else
      return Continue;
  }

  static const char* Name() {return "CheckResidual";}
  static int Arguments() {return 1;}
};

CREATE_INITIALIZER(VMGCommandCheckResidual)

