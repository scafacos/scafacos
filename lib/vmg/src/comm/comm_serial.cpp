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
 * @file   comm_serial.cpp
 * @author Julian Iseringhausen <isering@ins.uni-bonn.de>
 * @date   Mon Apr 18 12:28:12 2011
 *
 * @brief  VMG::CommSerial
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef HAVE_BOOST_FILESYSTEM
#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;
#endif

#include <cstdarg>
#include <cstdio>
#include <sstream>

#include "base/discretization.hpp"
#include "base/helper.hpp"
#include "base/index.hpp"
#include "base/stencil.hpp"
#include "base/vector.hpp"
#include "comm/comm_serial.hpp"
#include "grid/multigrid.hpp"
#include "grid/tempgrid.hpp"
#include "mg.hpp"

using namespace VMG;

static char print_buffer[512];

void CommSerial::InitCommSerial()
{
  error_norm_count = 0;
  residual_count = 0;
}

void CommSerial::CommSubgrid(Grid& grid_old, Grid& grid_new, const int& direction)
{
  Grid::iterator iter;

  if (&grid_old == &grid_new)
    return;

  for (iter = grid_old.Iterators().CompleteGrid().Begin(); iter != grid_old.Iterators().CompleteGrid().End(); ++iter)
    grid_new(*iter) = grid_old.GetVal(*iter);
}

void CommSerial::CommAddSubgrid(Grid& grid_old, Grid& grid_new, const int& direction)
{
  Grid::iterator iter;

  if (&grid_old == &grid_new)
    return;

  for (iter = grid_old.Iterators().CompleteGrid().Begin(); iter != grid_old.Iterators().CompleteGrid().End(); ++iter)
    grid_new(*iter) += grid_old.GetVal(*iter);
}

Grid& CommSerial::GetFinerGrid(Multigrid& multigrid)
{
  return multigrid(multigrid.Level()+1);
}

Grid& CommSerial::GetCoarserGrid(Multigrid& multigrid)
{
  return multigrid(multigrid.Level()-1);
}

void CommSerial::CommFromGhosts(Grid& grid)
{
  Grid::iterator iter;

  const GridIteratorSuite& iterators = grid.Iterators();
  const Index& size = grid.Local().Size();

  if (BoundaryConditions().X() == Periodic) {

    for (iter = iterators.Halo1().X().Begin(); iter != iterators.Halo1().X().End(); ++iter)
      grid((*iter).X()+size.X(), (*iter).Y(), (*iter).Z()) += grid.GetVal(*iter);

    for (iter = iterators.Halo2().X().Begin(); iter != iterators.Halo2().X().End(); ++iter)
      grid((*iter).X()-size.X(), (*iter).Y(), (*iter).Z()) += grid.GetVal(*iter);
  }

  if (BoundaryConditions().Y() == Periodic) {

    for (iter = iterators.Halo1().Y().Begin(); iter != iterators.Halo1().Y().End(); ++iter)
      grid((*iter).X(), (*iter).Y()+size.Y(), (*iter).Z()) += grid.GetVal(*iter);

    for (iter = iterators.Halo2().Y().Begin(); iter != iterators.Halo2().Y().End(); ++iter)
      grid((*iter).X(), (*iter).Y()-size.Y(), (*iter).Z()) += grid.GetVal(*iter);

  }

  if (BoundaryConditions().Z() == Periodic) {

    for (iter = iterators.Halo1().Z().Begin(); iter != iterators.Halo1().Z().End(); ++iter)
      grid((*iter).X(), (*iter).Y(), (*iter).Z()+size.Z()) += grid.GetVal(*iter);

    for (iter = iterators.Halo2().Z().Begin(); iter != iterators.Halo2().Z().End(); ++iter)
      grid((*iter).X(), (*iter).Y(), (*iter).Z()-size.Z()) += grid.GetVal(*iter);

  }
}

void CommSerial::CommToGhosts(Grid& grid)
{
  Grid::iterator iter;

  const GridIteratorSuite& iterators = grid.Iterators();
  const Index& size = grid.Local().Size();

  if (BoundaryConditions().Z() == Periodic) {

    for (iter = iterators.Halo1().Z().Begin(); iter != iterators.Halo1().Z().End(); ++iter)
      grid(*iter) = grid.GetVal((*iter).X(), (*iter).Y(), (*iter).Z()+size.Z());

    for (iter = iterators.Halo2().Z().Begin(); iter != iterators.Halo2().Z().End(); ++iter)
      grid(*iter) = grid.GetVal((*iter).X(), (*iter).Y(), (*iter).Z()-size.Z());

  }

  if (BoundaryConditions().Y() == Periodic) {

    for (iter = iterators.Halo1().Y().Begin(); iter != iterators.Halo1().Y().End(); ++iter)
      grid(*iter) = grid.GetVal((*iter).X(), (*iter).Y()+size.Y(), (*iter).Z());

    for (iter = iterators.Halo2().Y().Begin(); iter != iterators.Halo2().Y().End(); ++iter)
      grid(*iter) = grid.GetVal((*iter).X(), (*iter).Y()-size.Y(), (*iter).Z());

  }

  if (BoundaryConditions().X() == Periodic) {

    for (iter = iterators.Halo1().X().Begin(); iter != iterators.Halo1().X().End(); ++iter)
      grid(*iter) = grid.GetVal((*iter).X()+size.X(), (*iter).Y(), (*iter).Z());

    for (iter = iterators.Halo2().X().Begin(); iter != iterators.Halo2().X().End(); ++iter)
      grid(*iter) = grid.GetVal((*iter).X()-size.X(), (*iter).Y(), (*iter).Z());

  }
}

void CommSerial::CommToGhostsAsyncStart(Grid& grid)
{
  CommToGhosts(grid);
}

void CommSerial::CommToGhostsAsyncFinish(Grid& grid)
{
}

void CommSerial::CommFromGhostsAsyncStart(Grid& grid)
{
  CommFromGhosts(grid);
}

void CommSerial::CommFromGhostsAsyncFinish(Grid& grid)
{
}

void CommSerial::Print(const OutputLevel level, const char* format, ...)
{
  bool print = (level == Output);

#ifdef OUTPUT_INFO
  print |= (level == Info);
#endif
#ifdef OUTPUT_DEBUG
  print |= (level == Debug);
#endif
#ifdef OUTPUT_TIMING
  print |= (level == Timing);
#endif

  if (print) {
    va_list args;
    va_start(args, format);
    vsprintf(print_buffer, format, args);
    printf("VMG: %s\n", print_buffer);
    va_end(args);
  }
}

void CommSerial::PrintOnce(const OutputLevel level, const char* format, ...)
{
  bool print = (level == Output);

#ifdef OUTPUT_INFO
  print |= (level == Info);
#endif
#ifdef OUTPUT_DEBUG
  print |= (level == Debug);
#endif
#ifdef OUTPUT_TIMING
  print |= (level == Timing);
#endif

  if (print) {
    va_list args;
    va_start(args, format);
    vsprintf(print_buffer, format, args);
    printf("VMG: %s\n", print_buffer);
    va_end(args);
  }
}

void CommSerial::PrintXML(const std::string& filename, const std::string& xml_data)
{
  pugi::xml_document doc;
  doc.load(xml_data.c_str());

  pugi::xml_node node_global = doc.prepend_child("Global");
  pugi::xml_node node_num_processes = node_global.append_child("NumProcesses");
  pugi::xml_node node_num_processes_data = node_num_processes.append_child(pugi::node_pcdata);
  node_num_processes_data.set_value(Helper::ToString(GlobalProcs().Product()).c_str());

  doc.save_file(filename.c_str());
}

void CommSerial::PrintXMLAll(const std::string& filename, const std::string& xml_data)
{
  PrintXML(filename, xml_data);
}

void CommSerial::PrintGrid(Grid& grid, const char* information)
{
  Index i;
  std::stringstream out;
  std::ofstream out_file;

  OpenFileAndPrintHeader(out_file, grid, information);

  for (i.Z()=grid.Local().Begin().Z(); i.Z()<grid.Local().End().Z(); ++i.Z())
    for (i.Y()=grid.Local().Begin().Y(); i.Y()<grid.Local().End().Y(); ++i.Y())
      for (i.X()=grid.Local().Begin().X(); i.X()<grid.Local().End().X(); ++i.X())
	out << std::scientific << grid.GetVal(i) << std::endl;

  out_file << out.str();
  out_file.close();
}

void CommSerial::PrintDefect(Grid& sol, Grid& rhs, const char* information)
{
  Grid::iterator iter;
  std::ofstream out_file;
  std::stringstream out;

  TempGrid *temp = MG::GetTempGrid();
  temp->SetProperties(sol);
  temp->ImportFromResidual(sol, rhs);

  OpenFileAndPrintHeader(out_file, *temp, information);

  for (iter=temp->Iterators().Local().Begin(); iter!=temp->Iterators().Local().End(); ++iter)
    out << std::scientific << temp->GetVal(*iter) << std::endl;

  out_file << out.str();
  out_file.close();
}

void CommSerial::OpenFileAndPrintHeader(std::ofstream& out, const Grid& grid, const char* information)
{
  char path_str[129];
  int count = OutputCount();

  sprintf(path_str, "%s%04d.dat", OutputPath().c_str(), count);

  out.open(path_str, std::ios::trunc);

  out << "# vtk DataFile Version 2.0" << std::endl
      << count << ": " << information << std::endl
      << "ASCII" << std::endl
      << "DATASET STRUCTURED_POINTS" << std::endl
      << grid.Local().Size().X() << " "
      << grid.Local().Size().Y() << " "
      << grid.Local().Size().Z() << std::endl
      << "ORIGIN 0 0 0" << std::endl
      << "SPACING " << grid.Extent().MeshWidth().X() << " "
      << grid.Extent().MeshWidth() << " "
      << grid.Extent().MeshWidth() << std::endl
      << "POINT_DATA " << grid.Local().Size().Product() << std::endl
      << "SCALARS residual double 1" << std::endl
      << "LOOKUP_TABLE default" << std::endl;
}

void CommSerial::PrintAllSettings()
{
  Multigrid* mg = MG::GetFactory().Get("SOL")->Cast<Multigrid>();
  std::stringstream buf;
  std::ofstream out;

  buf << OutputPath() << "settings.txt";

  out.open(buf.str().c_str(), std::ios::trunc);

  for (int i=mg->MinLevel(); i<=mg->MaxLevel(); ++i) {

    out << "###########################################################" << std::endl
	<< "LEVEL:                " << i << std::endl
	<< "GLOBAL:" << std::endl
	<< "  LOCAL_BEGIN:        " << (*mg)(i).Global().LocalBegin() << std::endl
	<< "  LOCAL_END:          " << (*mg)(i).Global().LocalEnd() << std::endl
	<< "  LOCAL_SIZE:         " << (*mg)(i).Global().LocalSize() << std::endl
	<< "  GLOBAL_FINER_BEGIN: " << (*mg)(i).Global().GlobalFinerBegin() << std::endl
	<< "  GLOBAL_FINER_END:   " << (*mg)(i).Global().GlobalFinerEnd() << std::endl
	<< "  GLOBAL_FINER_SIZE:  " << (*mg)(i).Global().GlobalFinerSize() << std::endl
	<< "  LOCAL_FINER_BEGIN:  " << (*mg)(i).Global().LocalFinerBegin() << std::endl
	<< "  LOCAL_FINER_END:    " << (*mg)(i).Global().LocalFinerEnd() << std::endl
	<< "  LOCAL_FINER_SIZE:   " << (*mg)(i).Global().LocalFinerSize() << std::endl
	<< "  FINEST_ABS_BEGIN:   " << (*mg)(i).Global().FinestAbsBegin() << std::endl
	<< "  FINEST_ABS_END:     " << (*mg)(i).Global().FinestAbsEnd() << std::endl
	<< "  FINEST_ABS_SIZE:    " << (*mg)(i).Global().FinestAbsSize() << std::endl
	<< "  GLOBAL_SIZE:        " << (*mg)(i).Global().GlobalSize() << std::endl
	<< "  BOUNDARY_TYPE:      " << (*mg)(i).Global().BoundaryType() << std::endl
	<< "LOCAL:" << std::endl
	<< "  BEGIN:              " << (*mg)(i).Local().Begin() << std::endl
	<< "  END:                " << (*mg)(i).Local().End() << std::endl
	<< "  SIZE:               " << (*mg)(i).Local().Size() << std::endl
	<< "  SIZE_TOTAL:         " << (*mg)(i).Local().SizeTotal() << std::endl
	<< "  HALO_BEGIN_1:       " << (*mg)(i).Local().HaloBegin1() << std::endl
	<< "  HALO_END_1:         " << (*mg)(i).Local().HaloEnd1() << std::endl
	<< "  HALO_SIZE_1:        " << (*mg)(i).Local().HaloSize1() << std::endl
	<< "  HALO_BEGIN_2:       " << (*mg)(i).Local().HaloBegin2() << std::endl
	<< "  HALO_END_2:         " << (*mg)(i).Local().HaloEnd2() << std::endl
	<< "  HALO_SIZE_2:        " << (*mg)(i).Local().HaloSize2() << std::endl
	<< "  BOUNDARY_BEGIN_1:   " << (*mg)(i).Local().BoundaryBegin1() << std::endl
	<< "  BOUNDARY_END_1:     " << (*mg)(i).Local().BoundaryEnd1() << std::endl
	<< "  BOUNDARY_SIZE_1:    " << (*mg)(i).Local().BoundarySize1() << std::endl
	<< "  BOUNDARY_BEGIN_2:   " << (*mg)(i).Local().BoundaryBegin2() << std::endl
	<< "  BOUNDARY_END_2:     " << (*mg)(i).Local().BoundaryEnd2() << std::endl
	<< "  BOUNDARY_SIZE_2:    " << (*mg)(i).Local().BoundarySize2() << std::endl
	<< "  FINER_BEGIN:        " << (*mg)(i).Local().FinerBegin() << std::endl
	<< "  FINER_END:          " << (*mg)(i).Local().FinerEnd() << std::endl
	<< "  FINER_SIZE:         " << (*mg)(i).Local().FinerSize() << std::endl
	<< "EXTENT:" << std::endl
	<< "  BEGIN:              " << (*mg)(i).Extent().Begin() << std::endl
	<< "  END:                " << (*mg)(i).Extent().End() << std::endl
	<< "  SIZE:               " << (*mg)(i).Extent().Size() << std::endl
	<< "  MESH_WIDTH:         " << (*mg)(i).Extent().MeshWidth() << std::endl
	<< "###########################################################" << std::endl;

  }

  assert(out.good());
  out.close();
}

inline int GetIndex(const Grid& grid, int i, int j, int k)
{
  if (grid.Global().BoundaryType() == LocallyRefined)
    return k + grid.Local().Size().Z() * (j + grid.Local().Size().Y() * i);
  else
    return k + grid.Local().SizeTotal().Z() * (j + grid.Local().SizeTotal().Y() * i);
}

void CommSerial::PrintGridStructureLevel(Grid& grid, std::ofstream& out)
{
  const Vector& sp = grid.Extent().MeshWidth();
  int numLines;

  if (grid.Global().BoundaryType() == LocallyRefined)
    numLines = grid.Local().Size().X() * grid.Local().Size().Y() +
      grid.Local().Size().Y() * grid.Local().Size().Z() +
      grid.Local().Size().X() * grid.Local().Size().Z();
  else
    numLines = grid.Local().SizeTotal().X() * grid.Local().SizeTotal().Y() +
      grid.Local().SizeTotal().Y() * grid.Local().SizeTotal().Z() +
      grid.Local().SizeTotal().X() * grid.Local().SizeTotal().Z();

  out << "    <Piece";

  if (grid.Global().BoundaryType() == LocallyRefined) {
    out << " NumberOfPoints=\"" << grid.Local().Size().Product()  << "\""
	<< " NumberOfVerts=\"" << grid.Local().Size().Product()  << "\"";
  }else {
    out << " NumberOfPoints=\"" << grid.Local().SizeTotal().Product()  << "\""
	<< " NumberOfVerts=\"" << grid.Local().SizeTotal().Product()  << "\"";
  }

  out << " NumberOfLines=\"" << numLines << "\""
      << " NumberOfStrips=\"0\""
      << " NumberOfPolys=\"0\"" << ">" << std::endl
      << "      <Points>" << std::endl
      << "        <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">" << std::endl
      << "          ";

  if (grid.Global().BoundaryType() == LocallyRefined) {
    for (int i=0; i<grid.Local().Size().X(); i++)
      for (int j=0; j<grid.Local().Size().Y(); j++)
	for (int k=0; k<grid.Local().Size().Z(); k++)
	  out << grid.Extent().Begin().X() + i * sp.X() << " "
	      << grid.Extent().Begin().Y() + j * sp.Y() << " "
	      << grid.Extent().Begin().Z() + k * sp.Z() << " ";
  }else {
    for (int i=0; i<grid.Local().SizeTotal().X(); i++)
      for (int j=0; j<grid.Local().SizeTotal().Y(); j++)
	for (int k=0; k<grid.Local().SizeTotal().Z(); k++)
	  out << grid.Extent().Begin().X() + i * sp.X() << " "
	      << grid.Extent().Begin().Y() + j * sp.Y() << " "
	      << grid.Extent().Begin().Z() + k * sp.Z() << " ";
  }

  out << std::endl
      << "        </DataArray>" << std::endl
      << "      </Points>" << std::endl
      << "      <Verts>" << std::endl
      << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">" << std::endl
      << "          ";

  if (grid.Global().BoundaryType() == LocallyRefined) {
    for (int i=0; i<grid.Local().Size().Product(); i++)
      out << i << " ";
  }else {
    for (int i=0; i<grid.Local().SizeTotal().Product(); i++)
      out << i << " ";
  }

  out << std::endl
      << "        </DataArray>" << std::endl
      << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">" << std::endl
      << "          ";

  if (grid.Global().BoundaryType() == LocallyRefined) {
    for (int i=1; i<=grid.Local().Size().Product(); i++)
      out  << i << " ";
  }else {
    for (int i=1; i<=grid.Local().SizeTotal().Product(); i++)
      out  << i << " ";
  }

  out << std::endl
      << "        </DataArray>" << std::endl
      << "      </Verts>" << std::endl
      << "      <Lines>" << std::endl
      << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">" << std::endl
      << "          ";

  if (grid.Global().BoundaryType() == LocallyRefined) {
    for (int i=0; i<grid.Local().Size().X(); i++)
      for (int j=0; j<grid.Local().Size().Y(); j++)
	out << GetIndex(grid,i,j,0) << " "
	    << GetIndex(grid,i,j,grid.Local().Size().Z()-1) << " ";

    for (int j=0; j<grid.Local().Size().Y(); j++)
      for (int k=0; k<grid.Local().Size().Z(); k++)
	out << GetIndex(grid,0,j,k) << " "
	    << GetIndex(grid,grid.Local().Size().X()-1,j,k) << " ";

    for (int i=0; i<grid.Local().Size().X(); i++)
      for (int k=0; k<grid.Local().Size().Z(); k++)
	out << GetIndex(grid,i,0,k) << " "
	    << GetIndex(grid,i,grid.Local().Size().Y()-1,k) << " ";
  }else {
    for (int i=0; i<grid.Local().SizeTotal().X(); i++)
      for (int j=0; j<grid.Local().SizeTotal().Y(); j++)
	out << GetIndex(grid,i,j,0) << " "
	    << GetIndex(grid,i,j,grid.Local().SizeTotal().Z()-1) << " ";

    for (int j=0; j<grid.Local().SizeTotal().Y(); j++)
      for (int k=0; k<grid.Local().SizeTotal().Z(); k++)
	out << GetIndex(grid,0,j,k) << " "
	    << GetIndex(grid,grid.Local().SizeTotal().X()-1,j,k) << " ";

    for (int i=0; i<grid.Local().SizeTotal().X(); i++)
      for (int k=0; k<grid.Local().SizeTotal().Z(); k++)
	out << GetIndex(grid,i,0,k) << " "
	    << GetIndex(grid,i,grid.Local().SizeTotal().Y()-1,k) << " ";
  }

  out << std::endl
      << "        </DataArray>" << std::endl
      << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">" << std::endl
      << "          ";

  for (int i=1; i<=numLines; i++)
    out << 2*i << " ";

  out << std::endl
      << "        </DataArray>" << std::endl
      << "      </Lines>" << std::endl
      << "    </Piece>" << std::endl;
}

void CommSerial::DebugPrintGridStructure(Multigrid& grid)
{
  std::ofstream out;
  char path_str[129];

  sprintf(path_str, "%sgrid.vtp", OutputPath().c_str());

  out.open(path_str, std::ios::trunc);

  if (!out.good()) {
    Print(Info, "File %s not accessible.", path_str);
    return;
  }

  out << "<?xml version=\"1.0\"?>" << std::endl
      << "<VTKFile type=\"PolyData\">" << std::endl
      << "  <PolyData>" << std::endl;

  for (int i=grid.MinLevel(); i<=grid.MaxLevel(); i++)
    PrintGridStructureLevel(grid(i), out);

  out << "  </PolyData>" << std::endl
      << "</VTKFile>" << std::endl;

  out.close();
}

std::string CommSerial::CreateOutputDirectory()
{
#ifdef HAVE_BOOST_FILESYSTEM
  std::string path, unique_path;
  std::stringstream unique_suffix;
  int suffix_counter = 0;
  char buffer[129];
  time_t rawtime;
  struct tm *timeinfo;

  time(&rawtime);
  timeinfo = localtime(&rawtime);
  strftime(buffer, 128, "./output/%Y_%m_%d_%H_%M_%S/", timeinfo);
  path = buffer;

  if (!fs::exists("output"))
    fs::create_directory("output");

  unique_path = path;

  while (fs::exists(unique_path.c_str())) {

    unique_suffix.str("");
    unique_suffix << "_" << suffix_counter++ << "/";

    unique_path = path;
    unique_path.replace(unique_path.size()-1, 1, unique_suffix.str());

  }

  fs::create_directory(unique_path.c_str());

  return unique_path;

#else

  return "./";

#endif
}
