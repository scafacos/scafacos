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
 * @file   comm_mpi.hpp
 * @author Julian Iseringhausen <isering@ins.uni-bonn.de>
 * @date   Wed Jun 16 13:21:06 2010
 *
 * @brief  Class for MPI-based communication.
 *
 */

#ifndef COMM_MPI_HPP_
#define COMM_MPI_HPP_

#ifndef HAVE_MPI
#error You need MPI in order to compile CommMPI
#endif

#include <cstdarg>
#include <cstdio>
#include <list>
#include <map>
#include <vector>

#include "base/has_tempgrids.hpp"
#include "comm/comm.hpp"
#include "comm/mpi/settings.hpp"
#include "comm/mpi/has_request_vec.hpp"

namespace VMG
{

class DomainDecomposition;
class GridIteratorSet;
class TempGrid;

class CommMPI : public Comm, public HasTempGrids, public HasRequestVec
{
public:
  CommMPI(const Boundary& boundary, DomainDecomposition* domain_dec, const MPI_Comm& mpi_comm, bool register_ = true) :
    Comm(boundary, domain_dec, register_)
#ifdef VMG_ONE_SIDED
    ,win_created(false)
#endif
  {
    InitCommMPI(mpi_comm);
  }

  CommMPI(const Boundary& boundary, DomainDecomposition* domain_dec, bool register_ = true) :
    Comm(boundary, domain_dec, register_)
#ifdef VMG_ONE_SIDED
    ,win_created(false)
#endif
  {
    InitCommMPI(MPI_COMM_WORLD);
  }

  virtual ~CommMPI();

  Grid& GetCoarserGrid(Multigrid& multigrid);
  Grid& GetFinerGrid(Multigrid& multigrid);

  void CommFromGhosts(Grid& grid);
  void CommToGhosts(Grid& grid);
  void CommSubgrid(Grid& grid_old, Grid& grid_new, const int& direction);
  void CommAddSubgrid(Grid& grid_old, Grid& grid_new, const int& direction);

  void CommToGhostsAsyncStart(Grid& grid);
  void CommToGhostsAsyncFinish(Grid& grid);
  void CommFromGhostsAsyncStart(Grid& grid);
  void CommFromGhostsAsyncFinish(Grid& grid);

  void Barrier();

  vmg_float GlobalSum(vmg_float value);
  vmg_float GlobalSumRoot(vmg_float value);
  void GlobalSumArray(vmg_float* array, const vmg_int& size);
  vmg_float GlobalMax(vmg_float value);
  vmg_float GlobalMaxRoot(vmg_float value);
  void GlobalMaxArray(vmg_float* array, const vmg_int& size);
  void GlobalBroadcast(vmg_float& value);
  void GlobalGather(vmg_float& value, vmg_float* array);

  vmg_int GlobalSum(vmg_int value);
  vmg_int GlobalSumRoot(vmg_int value);
  void GlobalSumArray(vmg_int* array, const vmg_int& size);
  vmg_int GlobalMax(vmg_int value);
  vmg_int GlobalMaxRoot(vmg_int value);
  void GlobalMaxArray(vmg_int* array, const vmg_int& size);
  void GlobalBroadcast(vmg_int& value);
  void GlobalGather(vmg_int& value, vmg_int* array);

  void GlobalBroadcast(char* str);

  vmg_float LevelSum(const Grid& grid, vmg_float value);
  vmg_float LevelSumRoot(const Grid& grid, vmg_float value);
  void LevelSumArray(const Grid& grid,  vmg_float* array, const vmg_int& size);

  vmg_int LevelSum(const Grid& grid, vmg_int value);
  vmg_int LevelSumRoot(const Grid& grid, vmg_int value);
  void LevelSumArray(const Grid& grid,  vmg_int* array, const vmg_int& size);

  void PrintString(const char* format, ...);
  void PrintStringOnce(const char* format, ...);
  void PrintXML(const std::string& filename, const std::string& xml_data);
  void PrintXMLAll(const std::string& filename, const std::string& xml_data);
  void PrintAllSettings();
  void PrintGrid(Grid& grid, const char* information);
  void PrintDefect(Grid& sol, Grid& rhs, const char* information);

  virtual int GlobalRank() const;
  virtual int GlobalSize() const;
  virtual Index GlobalPos() const;
  virtual Index GlobalProcs() const;

  virtual int Rank(const Grid& grid) const;
  virtual int Size(const Grid& grid) const;
  virtual Index Pos(const Grid& grid) const;
  virtual Index Procs(const Grid& grid) const;

  void PostInit(Multigrid& sol, Multigrid& rhs)
  {
    settings.ComputeSettings(sol, rhs, comm_global);
  }

  virtual void DebugPrintError(const Grid& sol, const char* information) {}
  virtual void DebugPrintErrorNorm(Grid& sol) {}
  virtual void DebugPrintGridStructure(Multigrid& multigrid) {}

private:
  void InitCommMPI(const MPI_Comm& comm);

  void CreateOutputFiles(const Grid& grid, const std::stringstream& serial_data, const char* information,
			 const Index& begin_global, const Index& end_global,
			 const Index& begin_local, const Index& end_local,
			 const int& output_count);

  void CreateParallelOutputFile(const Grid& grid, MPI_Comm& comm,
				const int& output_count, const char* information,
				const Index& begin_global, const Index& end_global,
				const Index& begin_local, const Index& end_local);

  MPI_File CreateSerialOutputFile(const Grid& grid, MPI_Comm& comm,
				  const int& output_count, const char* information,
				  const Index& begin_global, const Index& end_global,
				  const Index& begin_local, const Index& end_local);

  void FinalizeSerialOutputFile(MPI_File& file);

  std::string CreateOutputDirectory();

protected:
  void IsendAll(Grid& grid, std::vector<VMG::MPI::Datatype>& types, const MPI_Comm& comm, const int& tag_start);
  void IrecvAll(Grid& grid, std::vector<VMG::MPI::Datatype>& types, const MPI_Comm& comm, const int& tag_start);

  void IsendAllBuffered(const Grid& grid, std::vector<VMG::MPI::Datatype>& types, const MPI_Comm& comm, const int& tag_start);
  void IrecvAllBuffered(std::vector<VMG::MPI::Datatype>& types, const MPI_Comm& comm, const int& tag_start);

  void ReplaceBufferAll(Grid& grid, const std::vector<VMG::MPI::Datatype>& types);
  void AddBufferAll(Grid& grid, const std::vector<VMG::MPI::Datatype>& types);

  void PrintGridInformation(const Grid& grid, char* filename, const std::string& name);

  VMG::MPI::Settings settings;

  MPI_Comm comm_global;
  MPI_Info info;

#ifdef VMG_ONE_SIDED
  bool win_created;
  MPI_Win win;
#endif

};

}

#endif /* COMM_MPI_HPP_ */
