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
 * @file   cycle_cs_periodic.cpp
 * @author Julian Iseringhausen <isering@ins.uni-bonn.de>
 * @date   Wed Jun 27 14:48:42 2012
 *
 * @brief  Generates a Correction Scheme cycle for
 *         periodic boundary conditions.
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "base/command_list.hpp"
#include "base/interface.hpp"
#include "cycles/cycle_cs_periodic.hpp"
#include "mg.hpp"

using namespace VMG;

void CycleCSPeriodic::Generate()
{
  CommandList* init = new CommandList();
  CommandList* loop = new CommandList();
  CommandList* finalize = new CommandList();

  const int min_level = MG::GetInterface()->MinLevel();
  const int max_level = MG::GetInterface()->MaxLevel();

  init->AddCommand("ClearGrid", "RHS");
  init->AddCommand("ClearGrid", "SOL");
  init->AddCommand("ImportRightHandSide");
  init->AddCommand("ForceDiscreteCompatibility");
  init->AddCommand("CheckConsistency", "RHS");
  init->AddCommand("InitializeIterationCounter");
  init->AddCommand("InitializeResidualNorm", "INITIAL_RESIDUAL");

  loop->AddCommand("ClearCoarseLevels", "SOL");
  loop->AddCommand("ClearCoarseLevels", "RHS");

  AddCycle(*loop, max_level-min_level+1, gamma);

  loop->AddCommand("ComputeResidualNorm", "RESIDUAL");
  loop->AddCommand("CheckResidual", "RESIDUAL");
  loop->AddCommand("CheckRelativeResidual", "RESIDUAL:INITIAL_RESIDUAL");
  loop->AddCommand("CheckIterationCounter");

  finalize->AddCommand("SetAverageToZero", "SOL");
  finalize->AddCommand("ExportSolution");

  init->Register("COMMANDLIST_INIT");
  loop->Register("COMMANDLIST_LOOP");
  finalize->Register("COMMANDLIST_FINALIZE");

}
