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
 * @file   comm.hpp
 * @author Julian Iseringhausen <isering@ins.uni-bonn.de>
 * @date   Wed Jun 16 13:21:06 2010
 *
 * @brief  Base class for communication.
 *
 */

#ifndef COMM_HPP_
#define COMM_HPP_

#include <list>

#include "base/defs.hpp"
#include "base/index.hpp"
#include "base/object.hpp"
#include "thirdparty/pugixml/pugixml.hpp"

namespace VMG
{

class DomainDecomposition;
class Grid;
class Multigrid;
class Stencil;

class Comm : public Object
{
public:
  Comm(const Boundary& boundary, DomainDecomposition* domain_dec, bool register_ = true) :
    Object("COMM", register_),
    bound(boundary),
    dd(domain_dec),
    particle_grid(NULL),
    output_directory_is_created(false),
    output_count(0)
  {}

  virtual ~Comm();

  virtual Grid& GetCoarserGrid(Multigrid& multigrid) = 0;
  virtual Grid& GetFinerGrid(Multigrid& multigrid) = 0;
  Grid& GetParticleGrid();

  virtual void CommFromGhosts(Grid& grid) = 0;
  virtual void CommToGhosts(Grid& grid) = 0;
  virtual void CommSubgrid(Grid& grid_old, Grid& grid_new, const int& direction) = 0;
  virtual void CommAddSubgrid(Grid& grid_old, Grid& grid_new, const int& direction) = 0;

  virtual void CommToGhostsAsyncStart(Grid& grid) = 0;
  virtual void CommToGhostsAsyncFinish(Grid& grid) = 0;
  virtual void CommFromGhostsAsyncStart(Grid& grid) = 0;
  virtual void CommFromGhostsAsyncFinish(Grid& grid) = 0;

  virtual void Barrier() {}

  virtual vmg_float GlobalSum(vmg_float value) {return value;}
  virtual vmg_float GlobalSumRoot(vmg_float value) {return value;}
  virtual void GlobalSumArray(vmg_float* array, const vmg_int& size) {}
  virtual vmg_float GlobalMax(vmg_float value) {return value;}
  virtual vmg_float GlobalMaxRoot(vmg_float value) {return value;}
  virtual void GlobalMaxArray(vmg_float* array, const vmg_int& size) {}
  virtual void GlobalBroadcast(vmg_float& value) {}
  virtual void GlobalGather(vmg_float& value, vmg_float* array) {array[0] = value;}

  virtual vmg_int GlobalSum(vmg_int value) {return value;}
  virtual vmg_int GlobalSumRoot(vmg_int value) {return value;}
  virtual void GlobalSumArray(vmg_int* array, const vmg_int& size) {}
  virtual vmg_int GlobalMax(vmg_int value) {return value;}
  virtual vmg_int GlobalMaxRoot(vmg_int value) {return value;}
  virtual void GlobalMaxArray(vmg_int* array, const vmg_int& size) {}
  virtual void GlobalBroadcast(vmg_int& value) {}
  virtual void GlobalGather(vmg_int& value, vmg_int* array) {array[0] = value;}

  virtual void GlobalBroadcast(char* str) {}

  virtual vmg_float LevelSum(const Grid& grid, vmg_float value) {return value;}
  virtual vmg_float LevelSumRoot(const Grid& grid, vmg_float value) {return value;}
  virtual void LevelSumArray(const Grid& grid,  vmg_float* array, const vmg_int& size) {}

  virtual vmg_int LevelSum(const Grid& grid, vmg_int value) {return value;}
  virtual vmg_int LevelSumRoot(const Grid& grid, vmg_int value) {return value;}
  virtual void LevelSumArray(const Grid& grid, vmg_int* array, const vmg_int& size) {}

  virtual void Print(const OutputLevel level, const char* format, ...) = 0;
  virtual void PrintOnce(const OutputLevel level, const char* format, ...) = 0;

  virtual void PrintXML(const std::string& filename, const std::string& xml_data) = 0;
  virtual void PrintXMLAll(const std::string& filename, const std::string& xml_data) = 0;
  virtual void PrintAllSettings() = 0;
  virtual void PrintGrid(Grid& grid, const char* information) = 0;
  virtual void PrintDefect(Grid& sol, Grid& rhs, const char* information) = 0;

  virtual void DebugPrintError(const Grid& sol, const char* information) {}
  virtual void DebugPrintErrorNorm(Grid& sol) {}
  virtual void DebugPrintGridStructure(Multigrid& multigrid) {}

  virtual int GlobalRank() const {return 0;}
  virtual int GlobalSize() const {return 1;}
  virtual Index GlobalPos() const {return Index(0);}
  virtual Index GlobalProcs() const {return Index(1);}

  virtual int Rank(const Grid& grid) const {return 0;}
  virtual int Size(const Grid& grid) const {return 1;}
  virtual Index Pos(const Grid& grid) const {return Index(0);}
  virtual Index Procs(const Grid& grid) const {return Index(1);}

  virtual void PostInit(Multigrid& sol, Multigrid& rhs) = 0;

  const Boundary& BoundaryConditions() const {return bound;}

  DomainDecomposition& GetDomainDecomposition() {return *dd;}

  vmg_float ComputeResidualNorm(Multigrid& sol, Multigrid& rhs);

protected:
  const std::string& OutputPath();
  int OutputCount() {return ++output_count;}

  Boundary bound;
  DomainDecomposition* dd;

private:
  virtual std::string CreateOutputDirectory() = 0;

  Grid* particle_grid;
  std::string output_path_str;
  bool output_directory_is_created;
  int output_count;
};

}

#endif /* COMM_HPP_ */
