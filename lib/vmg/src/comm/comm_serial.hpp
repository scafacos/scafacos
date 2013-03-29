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
 * @file   comm_serial.hpp
 * @author Julian Iseringhausen <isering@ins.uni-bonn.de>
 * @date   Mon Apr 18 12:28:29 2011
 *
 * @brief  VMG::CommSerial
 *
 */

#ifndef COMM_SERIAL_HPP_
#define COMM_SERIAL_HPP_

#include <fstream>
#include <map>

#include "comm/comm.hpp"
#include "comm/domain_decomposition_serial.hpp"

namespace VMG
{

class DomainDecomposition;
class Grid;
class Multigrid;

class CommSerial : public Comm
{
public:
  CommSerial(const Boundary& boundary, bool register_ = true) :
    Comm(boundary, new DomainDecompositionSerial(), register_)
  {
    InitCommSerial();
  }

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

  void Print(const OutputLevel level, const char* format, ...);
  void PrintOnce(const OutputLevel level, const char* format, ...);
  void PrintXML(const std::string& filename, const std::string& xml_data);
  void PrintXMLAll(const std::string& filename, const std::string& xml_data);
  void PrintAllSettings();
  void PrintGrid(Grid& grid, const char* information);
  void PrintDefect(Grid& sol, Grid& rhs, const char* information);

  void DebugPrintGridStructure(Multigrid& multigrid);

  void PostInit(Multigrid& sol, Multigrid& rhs) {}

private:
  void InitCommSerial();

  void OpenFileAndPrintHeader(std::ofstream& out, const Grid& mesh, const char* information);
  void PrintGridStructureLevel(Grid& grid, std::ofstream& out);
  std::string CreateOutputDirectory();

  int error_norm_count, residual_count;

  std::map< const Grid*, Grid*> finer, coarser;

};

}

#endif /* COMM_SERIAL_HPP_ */
