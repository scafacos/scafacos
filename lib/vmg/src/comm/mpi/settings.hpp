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

#ifndef SETTINGS_HPP_
#define SETTINGS_HPP_

#ifndef HAVE_MPI
#error You need MPI in order to compile VMG::MPI::Settings
#endif

#include <map>
#include <set>

#include "comm/mpi/datatypes_global.hpp"
#include "comm/mpi/datatypes_local.hpp"
#include "comm/mpi/key.hpp"

namespace VMG
{

class TempGrid;

namespace MPI
{

class DatatypesGlobal;

class Settings
{
public:
  Settings();
  ~Settings();

  void ComputeSettings(Multigrid& sol, Multigrid& rhs, MPI_Comm& comm);

  Grid& FinerGrid(const Grid& grid);
  Grid& CoarserGrid(const Grid& grid);

  MPI_Comm CommunicatorGlobal(const Grid& grid) const;
  MPI_Comm CommunicatorLocal(const Grid& grid) const;

  MPI_Datatype& Datatype(const Index& begin, const Index& end,
			 const Index& size_local, const Index& size_global,
			 const int& level);
  VMG::MPI::DatatypesGlobal& DatatypesGlobal(const Grid& grid_old, const Grid& grid_new, const int& direction);
  VMG::MPI::DatatypesLocal& DatatypesLocal(const Grid& grid);

private:
  Index GlobalDims(MPI_Comm comm, Index pos);
  void AddDatatypeGlobal(const Grid& grid_old, const Grid& grid_new, const int& direction);

  void CreateGlobalCommunicator(MPI_Comm& comm_global, const Grid* grid_1, const Grid* grid_2=NULL, const Grid* grid_3=NULL);
  void CreateLocalCommunicator(MPI_Comm& comm_global, const Grid& grid);

  std::map<int, MPI_Comm> communicators_global;
  std::map<KeyUnsorted, MPI_Comm> communicators_local;
  std::set<MPI_Comm> communicators_local_unique;
  std::map<const Grid*, Grid*> finer_grids, coarser_grids;
  std::map<KeyUnsorted, MPI_Datatype> datatypes;
  std::map<KeyUnsorted, VMG::MPI::DatatypesGlobal> datatypes_global;
  std::map<KeyUnsorted, VMG::MPI::DatatypesLocal> datatypes_local;
};

}

}

#endif /* SETTINGS_HPP_ */
