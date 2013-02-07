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
 * @file   mg.hpp
 * @author Julian Iseringhausen <isering@ins.uni-bonn.de>
 * @date   Mon Apr 18 13:01:13 2011
 *
 * @brief  Header file for the VMG main class.
 *
 */

#ifndef MG_HPP_
#define MG_HPP_

#include <map>

#include "base/command_factory.hpp"
#include "base/factory.hpp"

namespace VMG
{

class Comm;
class Cycle;
class Discretization;
class Grid;
class Interface;
class LevelOperator;
class Multigrid;
class Smoother;
class Solver;
class Stencil;
class TempGrid;

class MG
{
public:
  static MG* Instance()
  {
    static MG mgs;
    return &mgs;
  }

  static VMG::Comm* GetComm();
  static VMG::Cycle* GetCycle();
  static VMG::Discretization* GetDiscretization();
  static VMG::Interface* GetInterface();
  static VMG::LevelOperator* GetLevelOperator();
  static VMG::Multigrid* GetRhs();
  static VMG::Multigrid* GetSol();
  static VMG::Grid& GetRhsMaxLevel();
  static VMG::Grid& GetSolMaxLevel();
  static VMG::Smoother* GetSmoother();
  static VMG::Solver* GetSolver();
  static VMG::TempGrid* GetTempGrid();

  static VMG::Factory& GetFactory();
  static VMG::CommandFactory& GetCommands();

  static void PostInit();
  static void Solve();
  static void Destroy();
  static bool IsInitialized();

  static void SetState(const int& key);

private:
  MG();
  virtual ~MG();

  void RegisterLibraryCommands();

  vmg_float ComputeVectorNorm(const Multigrid& vec);

  static VMG::CommandFactory command_factory;
  std::map<int, VMG::Factory> factories;
  int state;
};

}

#endif
