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
 * @file   mg.cpp
 * @author Julian Iseringhausen <isering@ins.uni-bonn.de>
 * @date   Sat Jun 12 20:36:24 2010
 *
 * @brief  A multigrid solver
 *
 * This file contains the implementation of the main class for
 * a multigrid solver.
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef DEBUG_MEASURE_TIME
#ifdef HAVE_MPI
#include <mpi.h>
#ifdef HAVE_MARMOT
#include <enhancempicalls.h>
#include <sourceinfompicalls.h>
#endif
#endif
#endif

#include <cmath>
#include <cstdio>
#include <cfloat>
#include <ctime>
#include <string>
#include <fstream>
#include <sstream>

#include "base/command_list.hpp"
#include "base/discretization.hpp"
#include "base/factory.hpp"
#include "base/interface.hpp"
#include "base/timer.hpp"
#include "comm/comm.hpp"
#include "cycles/cycle.hpp"
#include "grid/grid.hpp"
#include "grid/tempgrid.hpp"
#include "level/level_operator.hpp"
#include "smoother/smoother.hpp"
#include "solver/solver.hpp"
#include "mg.hpp"

using namespace VMG;

#define REGISTER_COMMAND(a) extern void Initialize##a();Initialize##a();

VMG::CommandFactory MG::command_factory;

static void VMGRegisterBuiltinCommands()
{
  REGISTER_COMMAND(VMGCommandCheckConsistency);
  REGISTER_COMMAND(VMGCommandCheckIterationCounter);
  REGISTER_COMMAND(VMGCommandCheckRelativeResidual);
  REGISTER_COMMAND(VMGCommandCheckResidual);
  REGISTER_COMMAND(VMGCommandClearCoarseLevels);
  REGISTER_COMMAND(VMGCommandClearGrid);
  REGISTER_COMMAND(VMGCommandComputeResidualNorm);
  REGISTER_COMMAND(VMGCommandCopyBoundary);
  REGISTER_COMMAND(VMGCommandExecuteCycle);
  REGISTER_COMMAND(VMGCommandExecuteCycleLoop);
  REGISTER_COMMAND(VMGCommandExecuteFullCycle);
  REGISTER_COMMAND(VMGCommandExecuteFullCycleLoop);
  REGISTER_COMMAND(VMGCommandExportSolution);
  REGISTER_COMMAND(VMGCommandForceDiscreteCompatibility);
  REGISTER_COMMAND(VMGCommandImportRightHandSide);
  REGISTER_COMMAND(VMGCommandInterpolateFMG);
  REGISTER_COMMAND(VMGCommandInitializeIterationCounter);
  REGISTER_COMMAND(VMGCommandInitializeResidualNorm);
  REGISTER_COMMAND(VMGCommandNOP);
  REGISTER_COMMAND(VMGCommandPrintAllSettings);
  REGISTER_COMMAND(VMGCommandPrintDefect);
  REGISTER_COMMAND(VMGCommandPrintGridStructure);
  REGISTER_COMMAND(VMGCommandPrintGrid);
  REGISTER_COMMAND(VMGCommandPrintResidualNorm);
  REGISTER_COMMAND(VMGCommandPrintRunningTime);
  REGISTER_COMMAND(VMGCommandProlongate);
  REGISTER_COMMAND(VMGCommandRestrict);
  REGISTER_COMMAND(VMGCommandSetAverageToZero);
  REGISTER_COMMAND(VMGCommandSetCoarserDirichletValues);
  REGISTER_COMMAND(VMGCommandSetLevel);
  REGISTER_COMMAND(VMGCommandSmooth);
  REGISTER_COMMAND(VMGCommandSolve);
}

MG::MG()
{
  state = 0;
  VMGRegisterBuiltinCommands();
}

MG::~MG()
{
  MG::Destroy();
}

// Brings Multigrid back to starting state.
void MG::Destroy()
{
  MG::Instance()->factories.clear();
  MG::Instance()->state = 0;
  Timer::Clear();
}

/*
 * Post init communication class
 */
void MG::PostInit()
{
  if (GetFactory().TestObject("COMM") && GetFactory().TestObject("INTERFACE")) {

    Multigrid* sol = new Multigrid(GetComm(), GetInterface());
    sol->Register("SOL");

    Multigrid* rhs = new Multigrid(GetComm(), GetInterface());
    rhs->Register("RHS");

    TempGrid* temp = new TempGrid();
    temp->Register("TEMPGRID");

    new ObjectStorage<int>("GLOBAL_MAXLEVEL", sol->GlobalMaxLevel());
    new ObjectStorage<int>("MINLEVEL", sol->MinLevel());
    new ObjectStorage<int>("MAXLEVEL", sol->MaxLevel());

    if (GetFactory().TestObject("CYCLE"))
      GetCycle()->Generate();

    if (GetFactory().TestObject("COMM"))
      GetComm()->PostInit(*GetSol(), *GetRhs());

  }
}

/**
 * Solves a given system with a multigrid method
 *
 */
void MG::Solve()
{
#ifdef DEBUG_MEASURE_TIME
#ifdef HAVE_MPI
  GetComm()->Barrier();
#endif
  Timer::Start("CompleteRunningTime");
#endif

  CommandList* cl_init = MG::GetFactory().Get("COMMANDLIST_INIT")->Cast<CommandList>();
  CommandList* cl_loop = MG::GetFactory().Get("COMMANDLIST_LOOP")->Cast<CommandList>();
  CommandList* cl_finalize = MG::GetFactory().Get("COMMANDLIST_FINALIZE")->Cast<CommandList>();

  cl_init->ExecuteList();

  while (cl_loop->ExecuteList() == Continue);

  cl_finalize->ExecuteList();

#ifdef DEBUG_MEASURE_TIME
#ifdef HAVE_MPI
  GetComm()->Barrier();
#endif
  Timer::Stop("CompleteRunningTime");
#ifdef DEBUG_MEASURE_TIME_OUTPUT
#ifdef HAVE_MPI
  Timer::PrintGlobal();
#else
  Timer::Print();
#endif
#endif
#endif
}

void MG::SetState(const int& state_)
{
  MG::Instance()->state = state_;
}

VMG::Factory& MG::GetFactory()
{
  std::map<int, VMG::Factory>::iterator iter = MG::Instance()->factories.find(MG::Instance()->state);

  if (iter == MG::Instance()->factories.end())
    iter = MG::Instance()->factories.insert(std::make_pair(MG::Instance()->state, Factory())).first;

  assert(iter != MG::Instance()->factories.end());

  return iter->second;
}

VMG::CommandFactory& MG::GetCommands()
{
  return MG::command_factory;
}

Comm* MG::GetComm()
{
  return MG::GetFactory().Get("COMM")->Cast<VMG::Comm>();
}

Cycle* MG::GetCycle()
{
  return MG::GetFactory().Get("CYCLE")->Cast<VMG::Cycle>();
}

Discretization* MG::GetDiscretization()
{
  return MG::GetFactory().Get("DISCRETIZATION")->Cast<VMG::Discretization>();
}

LevelOperator* MG::GetLevelOperator()
{
  return MG::GetFactory().Get("LEVEL_OPERATOR")->Cast<VMG::LevelOperator>();
}

Multigrid* MG::GetRhs()
{
  return MG::GetFactory().Get("RHS")->Cast<VMG::Multigrid>();
}

Multigrid* MG::GetSol()
{
  return MG::GetFactory().Get("SOL")->Cast<VMG::Multigrid>();
}

VMG::Grid& MG::GetRhsMaxLevel()
{
  return (*MG::GetRhs())(MG::GetRhs()->MaxLevel());
}

VMG::Grid& MG::GetSolMaxLevel()
{
  return (*MG::GetSol())(MG::GetSol()->MaxLevel());
}

Smoother* MG::GetSmoother()
{
  return MG::GetFactory().Get("SMOOTHER")->Cast<VMG::Smoother>();
}

Solver* MG::GetSolver()
{
  return MG::GetFactory().Get("SOLVER")->Cast<VMG::Solver>();
}

TempGrid* MG::GetTempGrid()
{
  return MG::GetFactory().Get("TEMPGRID")->Cast<VMG::TempGrid>();
}

Interface* MG::GetInterface()
{
  return MG::GetFactory().Get("INTERFACE")->Cast<VMG::Interface>();
}

bool MG::IsInitialized()
{
  const Factory& f = GetFactory();

  bool init = true;

  init &= f.TestObject("COMM");
  init &= f.TestObject("LEVEL_OPERATOR");
  init &= f.TestObject("RHS");
  init &= f.TestObject("SOL");
  init &= f.TestObject("SOLVER");
  init &= f.TestObject("SMOOTHER");
  init &= f.TestObject("DISCRETIZATION");
  init &= f.TestObject("MAX_ITERATION");
  init &= f.TestObject("PRECISION");
  init &= f.TestObject("PRESMOOTHSTEPS");
  init &= f.TestObject("POSTSMOOTHSTEPS");
  init &= f.TestObject("COMMANDLIST_INIT");
  init &= f.TestObject("COMMANDLIST_LOOP");
  init &= f.TestObject("COMMANDLIST_FINALIZE");
  init &= f.TestObject("MINLEVEL");
  init &= f.TestObject("MAXLEVEL");
  init &= f.TestObject("GLOBAL_MAXLEVEL");
  init &= f.TestObject("INTERFACE");

  return init;
}
