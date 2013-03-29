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
 * @file   comm_mpi.cpp
 * @author Julian Iseringhausen <isering@ins.uni-bonn.de>
 * @date   Wed Jun 16 13:21:06 2010
 *
 * @brief  Class for MPI-based communication.
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef HAVE_MPI

#include <mpi.h>
#ifdef HAVE_MARMOT
#include <enhancempicalls.h>
#include <sourceinfompicalls.h>
#endif

#ifdef HAVE_BOOST_FILESYSTEM
#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;
#endif

#ifdef HAVE_VTK
#include <vtkAbstractArray.h>
#include <vtkImageData.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkXMLImageDataWriter.h>
#endif

#include <cstring>
#include <sstream>

#include "base/helper.hpp"
#include "base/tuple.hpp"
#include "comm/comm_mpi.hpp"
#include "comm/mpi/datatypes_local.hpp"
#include "grid/grid.hpp"
#include "grid/multigrid.hpp"
#include "grid/tempgrid.hpp"
#include "mg.hpp"
#include "base/timer.hpp"

static char print_buffer[512];

using namespace VMG;

void CommMPI::IsendAll(Grid& grid, std::vector<VMG::MPI::Datatype>& types, const MPI_Comm& comm, const int& tag_start)
{
  for (std::vector<VMG::MPI::Datatype>::const_iterator iter=types.begin(); iter!=types.end(); ++iter)
    iter->Isend(grid, tag_start, comm, Request());
}

void CommMPI::IrecvAll(Grid& grid, std::vector<VMG::MPI::Datatype>& types, const MPI_Comm& comm, const int& tag_start)
{
  for (std::vector<VMG::MPI::Datatype>::const_iterator iter=types.begin(); iter!=types.end(); ++iter)
    iter->Irecv(grid, tag_start, comm, Request());
}

void CommMPI::IsendAllBuffered(const Grid& grid, std::vector<VMG::MPI::Datatype>& types, const MPI_Comm& comm, const int& tag_start)
{
  for (std::vector<VMG::MPI::Datatype>::iterator iter=types.begin(); iter!=types.end(); ++iter)
    iter->IsendBuffered(grid, tag_start, comm, Request());
}

void CommMPI::IrecvAllBuffered(std::vector<VMG::MPI::Datatype>& types, const MPI_Comm& comm, const int& tag_start)
{
  for (std::vector<VMG::MPI::Datatype>::iterator iter=types.begin(); iter!=types.end(); ++iter)
    iter->IrecvBuffered(tag_start, comm, Request());
}

void CommMPI::ReplaceBufferAll(Grid& grid, const std::vector<VMG::MPI::Datatype>& types)
{
  for (std::vector<VMG::MPI::Datatype>::const_iterator iter=types.begin(); iter!= types.end(); ++iter)
    iter->GridReplace(grid);
}

void CommMPI::AddBufferAll(Grid& grid, const std::vector<VMG::MPI::Datatype>& types)
{
  for (std::vector<VMG::MPI::Datatype>::const_iterator iter=types.begin(); iter!= types.end(); ++iter)
    iter->GridSum(grid);
}

void CommMPI::CommSubgrid(Grid& grid_old, Grid& grid_new, const int& direction)
{
  MPI_Comm comm = settings.CommunicatorGlobal(grid_old);
  if (comm != MPI_COMM_NULL) {
    VMG::MPI::DatatypesGlobal& datatypes = settings.DatatypesGlobal(grid_old, grid_new, direction);
    IrecvAllBuffered(datatypes.Receive(), comm, 0411);
    IsendAllBuffered(grid_old, datatypes.Send(), comm, 0411);
    WaitAll();
    ReplaceBufferAll(grid_new, datatypes.Receive());
  }
}

void CommMPI::CommAddSubgrid(Grid& grid_old, Grid& grid_new, const int& direction)
{
  MPI_Comm comm = settings.CommunicatorGlobal(grid_old);
  if (comm != MPI_COMM_NULL) {
    VMG::MPI::DatatypesGlobal& datatypes = settings.DatatypesGlobal(grid_old, grid_new, direction);
    IrecvAllBuffered(datatypes.Receive(), comm, 1806);
    IsendAllBuffered(grid_old, datatypes.Send(), comm, 1806);
    WaitAll();
    AddBufferAll(grid_new, datatypes.Receive());
  }
}

void CommMPI::CommToGhosts(Grid& grid)
{
  MPI_Comm comm = settings.CommunicatorLocal(grid);
  if (comm != MPI_COMM_NULL) {
    VMG::MPI::DatatypesLocal& types = settings.DatatypesLocal(grid);
    IrecvAllBuffered(types.Halo(), comm, 2310);
    IsendAllBuffered(grid, types.NB(), comm, 2310);
    WaitAll();
    ReplaceBufferAll(grid, types.Halo());
  }
}

void CommMPI::CommToGhostsAsyncStart(Grid& grid)
{
  MPI_Comm comm = settings.CommunicatorLocal(grid);
  if (comm != MPI_COMM_NULL) {
    VMG::MPI::DatatypesLocal& types = settings.DatatypesLocal(grid);
    IrecvAllBuffered(types.Halo(), comm, 2412);
    IsendAllBuffered(grid, types.NB(), comm, 2412);
    TestAll();
  }
}

void CommMPI::CommToGhostsAsyncFinish(Grid& grid)
{
  WaitAll();
  ReplaceBufferAll(grid, settings.DatatypesLocal(grid).Halo());
}

void CommMPI::CommFromGhosts(Grid& grid)
{
  MPI_Comm comm = settings.CommunicatorLocal(grid);
  if (comm != MPI_COMM_NULL) {
    VMG::MPI::DatatypesLocal& types = settings.DatatypesLocal(grid);
    IrecvAllBuffered(types.NB(), comm, 1337);
    IsendAllBuffered(grid, types.Halo(), comm, 1337);
    WaitAll();
    AddBufferAll(grid, types.NB());
  }
}

void CommMPI::CommFromGhostsAsyncStart(Grid& grid)
{
 MPI_Comm comm = settings.CommunicatorLocal(grid);
 if (comm != MPI_COMM_NULL) {
   VMG::MPI::DatatypesLocal& types = settings.DatatypesLocal(grid);
   IrecvAllBuffered(types.NB(), comm, 0xc0ffee);
   IsendAllBuffered(grid, types.Halo(), comm, 0xc0ffee);
   TestAll();
 }
}

void CommMPI::CommFromGhostsAsyncFinish(Grid& grid)
{
  WaitAll();
  AddBufferAll(grid, settings.DatatypesLocal(grid).NB());
}

void CommMPI::Barrier()
{
  MPI_Barrier(comm_global);
}

vmg_float CommMPI::GlobalSum(vmg_float value)
{
  vmg_float result = 0;
  MPI_Allreduce(&value, &result, 1, MPI_DOUBLE, MPI_SUM, comm_global);
  return result;
}

vmg_float CommMPI::GlobalSumRoot(vmg_float value)
{
  vmg_float result = 0;
  MPI_Reduce(&value, &result, 1, MPI_DOUBLE, MPI_SUM, 0, comm_global);
  return result;
}

void CommMPI::GlobalSumArray(vmg_float* array, const vmg_int& size)
{
  MPI_Allreduce(MPI_IN_PLACE, array, size, MPI_DOUBLE, MPI_SUM, comm_global);
}

void CommMPI::GlobalBroadcast(vmg_float& value)
{
  MPI_Bcast(&value, 1, MPI_DOUBLE, 0, comm_global);
}

void CommMPI::GlobalGather(vmg_float& value, vmg_float* array)
{
  MPI_Gather(&value, 1, MPI_DOUBLE, array, 1, MPI_DOUBLE, 0, comm_global);
}

vmg_int CommMPI::GlobalSum(vmg_int value)
{
  vmg_int result = 0;
  MPI_Allreduce(&value, &result, 1, MPI_INT, MPI_SUM, comm_global);
  return result;
}

vmg_int CommMPI::GlobalSumRoot(vmg_int value)
{
  vmg_int result = 0;
  MPI_Reduce(&value, &result, 1, MPI_INT, MPI_SUM, 0, comm_global);
  return result;
}

void CommMPI::GlobalSumArray(vmg_int* array, const vmg_int& size)
{
  MPI_Allreduce(MPI_IN_PLACE, array, size, MPI_INT, MPI_SUM, comm_global);
}

vmg_float CommMPI::GlobalMax(vmg_float value)
{
  vmg_float result = 0.0;
  MPI_Allreduce(&value, &result, 1, MPI_DOUBLE, MPI_MAX, comm_global);
  return result;
}

vmg_float CommMPI::GlobalMaxRoot(vmg_float value)
{
  vmg_float result = 0.0;
  MPI_Reduce(&value, &result, 1, MPI_DOUBLE, MPI_MAX, 0, comm_global);
  return result;
}

void CommMPI::GlobalMaxArray(vmg_float* array, const vmg_int& size)
{
  MPI_Allreduce(MPI_IN_PLACE, array, size, MPI_DOUBLE, MPI_MAX, comm_global);
}

vmg_int CommMPI::GlobalMax(vmg_int value)
{
  vmg_int result = 0;
  MPI_Allreduce(&value, &result, 1, MPI_INT, MPI_MAX, comm_global);
  return result;
}

vmg_int CommMPI::GlobalMaxRoot(vmg_int value)
{
  vmg_int result = 0;
  MPI_Reduce(&value, &result, 1, MPI_INT, MPI_MAX, 0, comm_global);
  return result;
}

void CommMPI::GlobalMaxArray(vmg_int* array, const vmg_int& size)
{
  MPI_Allreduce(MPI_IN_PLACE, array, size, MPI_INT, MPI_MAX, comm_global);
}

void CommMPI::GlobalBroadcast(vmg_int& value)
{
  MPI_Bcast(&value, 1, MPI_INT, 0, comm_global);
}

void CommMPI::GlobalGather(vmg_int& value, vmg_int* array)
{
  MPI_Gather(&value, 1, MPI_INT, array, 1, MPI_INT, 0, comm_global);
}

void CommMPI::GlobalBroadcast(char* str)
{
  int size = std::strlen(str) + 1;
  MPI_Bcast(&size, 1, MPI_INT, 0, comm_global);
  MPI_Bcast(str, size, MPI_CHAR, 0, comm_global);
}

vmg_float CommMPI::LevelSum(const Grid& grid, vmg_float value)
{
  vmg_float result = 0.0;
  MPI_Comm comm = settings.CommunicatorLocal(grid);
  assert(comm != MPI_COMM_NULL);
  MPI_Allreduce(&value, &result, 1, MPI_DOUBLE, MPI_SUM, comm);
  return result;
}

vmg_float CommMPI::LevelSumRoot(const Grid& grid, vmg_float value)
{
  vmg_float result = 0.0;
  MPI_Comm comm = settings.CommunicatorLocal(grid);
  assert(comm != MPI_COMM_NULL);
  MPI_Reduce(&value, &result, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
  return result;
}

void CommMPI::LevelSumArray(const Grid& grid, vmg_float* array, const vmg_int& size)
{
  MPI_Comm comm = settings.CommunicatorLocal(grid);
  assert(comm != MPI_COMM_NULL);
  MPI_Allreduce(MPI_IN_PLACE, array, size, MPI_DOUBLE, MPI_SUM, comm);
}

vmg_int CommMPI::LevelSum(const Grid& grid, vmg_int value)
{
  vmg_int result = 0;
  MPI_Comm comm = settings.CommunicatorLocal(grid);
  assert(comm != MPI_COMM_NULL);
  MPI_Allreduce(&value, &result, 1, MPI_INT, MPI_SUM, comm);
  return result;
}

vmg_int CommMPI::LevelSumRoot(const Grid& grid, vmg_int value)
{
  vmg_int result = 0;
  MPI_Comm comm = settings.CommunicatorLocal(grid);
  assert(comm != MPI_COMM_NULL);
  MPI_Reduce(&value, &result, 1, MPI_INT, MPI_SUM, 0, comm);
  return result;
}

void CommMPI::LevelSumArray(const Grid& grid, vmg_int* array, const vmg_int& size)
{
  MPI_Comm comm = settings.CommunicatorLocal(grid);
  assert(comm != MPI_COMM_NULL);
  MPI_Allreduce(MPI_IN_PLACE, array, size, MPI_INT, MPI_SUM, comm);
}

void CommMPI::Print(const OutputLevel level, const char* format, ...)
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
    printf("VMG: Rank %d: %s\n", GlobalRank(), print_buffer);
    va_end(args);
  }
}

void CommMPI::PrintOnce(const OutputLevel level, const char* format, ...)
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

  if (GlobalRank() == 0 && print) {
    va_list args;
    va_start(args, format);
    vsprintf(print_buffer, format, args);
    printf("VMG: Rank %d: %s\n", GlobalRank(), print_buffer);
    va_end(args);
  }
}

void CommMPI::PrintXML(const std::string& filename, const std::string& xml_data)
{
  MPI_File file;
  std::stringstream path, xml_header;

  pugi::xml_document doc;
  pugi::xml_node node_data = doc.append_child("Global").append_child("NumProcesses").append_child(pugi::node_pcdata);
  node_data.set_value(Helper::ToString(GlobalProcs().Product()).c_str());
  doc.save(xml_header);

  path << OutputPath() << filename;

  char* filename_array = Helper::GetCharArray(path.str());
  char* xml_header_array = Helper::GetCharArray(xml_header.str());
  char* str_array = Helper::GetCharArray(xml_data);

  MPI_File_open(MPI_COMM_SELF, filename_array, MPI_MODE_WRONLY|MPI_MODE_CREATE, MPI_INFO_NULL, &file);
  MPI_File_set_size(file, 0);
  MPI_File_write(file, xml_header_array, xml_header.str().size(), MPI_CHAR, MPI_STATUS_IGNORE);
  MPI_File_write(file, str_array, xml_data.size(), MPI_CHAR, MPI_STATUS_IGNORE);
  MPI_File_close(&file);

  delete [] filename_array;
  delete [] xml_header_array;
  delete [] str_array;
}

void CommMPI::PrintXMLAll(const std::string& filename, const std::string& xml_data)
{
  MPI_File file;
  std::stringstream path;

  path << OutputPath() << filename;

  char* filename_array = Helper::GetCharArray(path.str());
  char* str_array = Helper::GetCharArray(xml_data);

  MPI_File_open(comm_global, filename_array, MPI_MODE_WRONLY|MPI_MODE_CREATE, MPI_INFO_NULL, &file);
  MPI_File_set_size(file, 0);

  if (GlobalRank() == 0) {
    std::stringstream xml_header;
    pugi::xml_document doc;
    pugi::xml_node node_data = doc.append_child("Global").append_child("NumProcesses").append_child(pugi::node_pcdata);
    node_data.set_value(Helper::ToString(GlobalProcs().Product()).c_str());
    doc.save(xml_header);

    char* xml_header_array = Helper::GetCharArray(xml_header.str());

    MPI_File_write_shared(file, xml_header_array, xml_header.str().size(), MPI_CHAR, MPI_STATUS_IGNORE);

    delete [] xml_header_array;
  }

  MPI_File_write_ordered(file, str_array, xml_data.size(), MPI_CHAR, MPI_STATUS_IGNORE);
  MPI_File_close(&file);

  delete [] filename_array;
  delete [] str_array;
}

void CommMPI::PrintGridInformation(const Grid& grid, char* filename, const std::string& name)
{
  std::stringstream buf;
  MPI_File file;
  int rank, size;
  int size_local, size_local_max;

  MPI_Comm comm = settings.CommunicatorGlobal(grid);
  MPI_Comm comm_local = settings.CommunicatorLocal(grid);

  if (comm_local != MPI_COMM_NULL)
    MPI_Comm_size(comm_local, &size_local);
  else
    size_local = 0;

  if (comm != MPI_COMM_NULL) {

    MPI_Reduce(&size_local, &size_local_max, 1, MPI_INT, MPI_MAX, 0, comm);

    MPI_File_open(comm, filename, MPI_MODE_WRONLY|MPI_MODE_CREATE|MPI_MODE_APPEND, MPI_INFO_NULL, &file);

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    if (rank == 0) {

      buf << "###########################################################" << std::endl
	  << "GLOBAL INFORMATION:" << std::endl
	  << "  NAME:                 " << name << std::endl
	  << "  LEVEL:                " << grid.Level() << std::endl
	  << "  COMM_SIZE_GLOBAL:     " << size << std::endl
	  << "  COMM_SIZE_LOCAL:      " << size_local_max << std::endl
	  << "  GLOBAL:" << std::endl
	  << "    GLOBAL_FINER_BEGIN: " << grid.Global().GlobalFinerBegin() << std::endl
	  << "    GLOBAL_FINER_END:   " << grid.Global().GlobalFinerEnd() << std::endl
	  << "    GLOBAL_FINER_SIZE:  " << grid.Global().GlobalFinerSize() << std::endl
	  << "    FINEST_ABS_BEGIN:   " << grid.Global().FinestAbsBegin() << std::endl
	  << "    FINEST_ABS_END:     " << grid.Global().FinestAbsEnd() << std::endl
	  << "    FINEST_ABS_SIZE:    " << grid.Global().FinestAbsSize() << std::endl
	  << "    GLOBAL_SIZE:        " << grid.Global().GlobalSize() << std::endl
	  << "  EXTENT:" << std::endl
	  << "    BEGIN:              " << grid.Extent().Begin() << std::endl
	  << "    END:                " << grid.Extent().End() << std::endl
	  << "    SIZE:               " << grid.Extent().Size() << std::endl
	  << "    MESH_WIDTH:         " << grid.Extent().MeshWidth() << std::endl
	  << std::endl
	  << "LOCAL INFORMATION:" << std::endl;
    }

    buf << "RANK " << rank << ":" << std::endl
	<< "  GLOBAL:" << std::endl
	<< "    LOCAL_BEGIN:          " << grid.Global().LocalBegin() << std::endl
	<< "    LOCAL_END:            " << grid.Global().LocalEnd() << std::endl
	<< "    LOCAL_SIZE:           " << grid.Global().LocalSize() << std::endl
	<< "    LOCAL_FINER_BEGIN:    " << grid.Global().LocalFinerBegin() << std::endl
	<< "    LOCAL_FINER_END:      " << grid.Global().LocalFinerEnd() << std::endl
	<< "    LOCAL_FINER_SIZE:     " << grid.Global().LocalFinerSize() << std::endl
	<< "    BOUNDARY_TYPE:        " << grid.Global().BoundaryType() << std::endl
	<< "  LOCAL:" << std::endl
	<< "    BEGIN:                " << grid.Local().Begin() << std::endl
	<< "    END:                  " << grid.Local().End() << std::endl
	<< "    SIZE:                 " << grid.Local().Size() << std::endl
	<< "    SIZE_TOTAL:           " << grid.Local().SizeTotal() << std::endl
	<< "    HALO_BEGIN_1:         " << grid.Local().HaloBegin1() << std::endl
	<< "    HALO_END_1:           " << grid.Local().HaloEnd1() << std::endl
	<< "    HALO_SIZE_1:          " << grid.Local().HaloSize1() << std::endl
	<< "    HALO_BEGIN_2:         " << grid.Local().HaloBegin2() << std::endl
	<< "    HALO_END_2:           " << grid.Local().HaloEnd2() << std::endl
	<< "    HALO_SIZE_2:          " << grid.Local().HaloSize2() << std::endl
	<< "    BOUNDARY_BEGIN_1:     " << grid.Local().BoundaryBegin1() << std::endl
	<< "    BOUNDARY_END_1:       " << grid.Local().BoundaryEnd1() << std::endl
	<< "    BOUNDARY_SIZE_1:      " << grid.Local().BoundarySize1() << std::endl
	<< "    BOUNDARY_BEGIN_2:     " << grid.Local().BoundaryBegin2() << std::endl
	<< "    BOUNDARY_END_2:       " << grid.Local().BoundaryEnd2() << std::endl
	<< "    BOUNDARY_SIZE_2:      " << grid.Local().BoundarySize2() << std::endl
	<< "    FINER_BEGIN:          " << grid.Local().FinerBegin() << std::endl
	<< "    FINER_END:            " << grid.Local().FinerEnd() << std::endl
	<< "    FINER_SIZE:           " << grid.Local().FinerSize() << std::endl;

    if (rank == size-1)
      buf << "###########################################################" << std::endl;

    char* char_buf = Helper::GetCharArray(buf.str());
    MPI_File_write_ordered(file, char_buf, buf.str().size(), MPI_CHAR, MPI_STATUS_IGNORE);
    delete [] char_buf;

    MPI_File_close(&file);

  }
}

void CommMPI::PrintAllSettings()
{
  std::stringstream buf;
  MPI_File file;

  const Multigrid& mg = *MG::GetSol();

  buf << OutputPath() << "settings.txt";
  char *filename = Helper::GetCharArray(buf.str());

  MPI_File_open(comm_global, filename, MPI_MODE_WRONLY|MPI_MODE_CREATE, MPI_INFO_NULL, &file);
  MPI_File_set_size(file, 0);
  MPI_File_close(&file);

  for (int i=mg.MinLevel(); i<=mg.MaxLevel(); ++i)
    PrintGridInformation(mg(i), filename, "MULTIGRID");

  for (int i=mg.MinLevel()+1; i<=mg.MaxLevel(); ++i)
    PrintGridInformation(settings.CoarserGrid(mg(i)), filename, "COARSER_GRID");

  for (int i=mg.MinLevel(); i<mg.MaxLevel(); ++i)
    PrintGridInformation(settings.FinerGrid(mg(i)), filename, "FINER_GRID");

  delete [] filename;

}

void CommMPI::PrintGrid(Grid& grid, const char* information)
{
  int output_count = OutputCount();

#ifdef HAVE_VTK

  if (settings.CommunicatorLocal(grid) != MPI_COMM_NULL) {

    Index end, end_global;

    for (int i=0; i<3; ++i) {
      end[i] = grid.Local().End()[i];
      end_global[i] = grid.Global().LocalEnd()[i];
    }

    vtkSmartPointer<vtkImageData> image = vtkSmartPointer<vtkImageData>::New();
    image->SetExtent(grid.Global().LocalBegin().X(), end_global.X()-1,
		     grid.Global().LocalBegin().Y(), end_global.Y()-1,
		     grid.Global().LocalBegin().Z(), end_global.Z()-1);
    image->SetSpacing(grid.Extent().MeshWidth().vec());
    image->SetOrigin(grid.Extent().Begin().vec());
    image->SetScalarTypeToDouble();
    image->SetNumberOfScalarComponents(1);
    image->AllocateScalars();
    image->GetPointData()->GetScalars()->SetName(information);

    Index i;
    for (i.X()=grid.Local().Begin().X(); i.X()<end.X(); ++i.X())
      for (i.Y()=grid.Local().Begin().Y(); i.Y()<end.Y(); ++i.Y())
	for (i.Z()=grid.Local().Begin().Z(); i.Z()<end.Z(); ++i.Z())
	  image->SetScalarComponentFromDouble(i.X() - grid.Local().Begin().X() + grid.Global().LocalBegin().X(),
					      i.Y() - grid.Local().Begin().Y() + grid.Global().LocalBegin().Y(),
					      i.Z() - grid.Local().Begin().Z() + grid.Global().LocalBegin().Z(),
					      0, grid.GetVal(i));

    image->Update();

    int rank, size;
    MPI_Comm_rank(comm_global, &rank);
    MPI_Comm_size(comm_global, &size);

    std::stringstream filename;
    filename << OutputPath() << std::setw(4) << std::setfill('0') << output_count << "_" << rank << ".vti";

    vtkSmartPointer<vtkXMLImageDataWriter> writer = vtkSmartPointer<vtkXMLImageDataWriter>::New();
    writer->SetFileName(filename.str().c_str());
    writer->SetCompressorTypeToNone();
    writer->SetDataModeToAscii();
    writer->SetInput(image);
    writer->Update();
    writer->Write();

  }

#else /* HAVE_VTK */
  Index i;
  std::stringstream buf;

  Index begin, end;
  Index begin_local, end_local, begin_global, end_global;

  CommToGhosts(grid);

  for (int i=0; i<3; ++i) {
    end[i] = grid.Local().End()[i] + (grid.Global().LocalEnd()[i] == grid.Global().GlobalSize()[i] ? 0 : grid.Local().HaloSize1()[i]);
    end_local[i] = grid.Global().LocalEnd()[i] - (grid.Global().LocalEnd()[i] == grid.Global().GlobalSize()[i] ? 1 : 0);
  }

  begin = grid.Local().Begin();
  begin_local = grid.Global().LocalBegin();
  begin_global = 0;
  end_global = grid.Global().GlobalSize()-1;

  for (i.Z()=begin.Z(); i.Z()<end.Z(); ++i.Z())
    for (i.Y()=begin.Y(); i.Y()<end.Y(); ++i.Y())
      for (i.X()=begin.X(); i.X()<end.X(); ++i.X())
	buf << std::scientific << grid.GetVal(i) << " ";

  CreateOutputFiles(grid, buf, information,
		    begin_global, end_global,
		    begin_local, end_local,
		    output_count);
#endif /* HAVE_VTK */
}

void CommMPI::PrintDefect(Grid& sol, Grid& rhs, const char* information)
{
  TempGrid *temp = MG::GetTempGrid();
  temp->SetProperties(sol);
  temp->ImportFromResidual(sol, rhs);
  PrintGrid(*temp, information);
}

void CommMPI::CreateOutputFiles(const Grid& grid, const std::stringstream& serial_data, const char* information,
				const Index& begin_global, const Index& end_global,
				const Index& begin_local, const Index& end_local,
				const int& output_count)
{
  MPI_Comm comm = settings.CommunicatorGlobal(grid);

  if (comm != MPI_COMM_NULL) {

    MPI_File file;
    std::string conv_information = Helper::ReplaceWhitespaces(information, "_");

    CreateParallelOutputFile(grid, comm, output_count, conv_information.c_str(),
			     begin_global, end_global, begin_local, end_local);

    file = CreateSerialOutputFile(grid, comm, output_count, conv_information.c_str(),
				  begin_global, end_global, begin_local, end_local);

    char *char_buf = Helper::GetCharArray(serial_data.str());
    MPI_File_write(file, char_buf, serial_data.str().size(), MPI_CHAR, MPI_STATUS_IGNORE);
    delete [] char_buf;

    FinalizeSerialOutputFile(file);

  }
}

void CommMPI::CreateParallelOutputFile(const Grid& grid, MPI_Comm& comm,
				       const int& output_count, const char* information,
				       const Index& begin_global, const Index& end_global,
				       const Index& begin_local, const Index& end_local)
{
  int rank;
  MPI_File file;
  char parallel_filename[513], serial_filename[513];
  std::stringstream buf;

  MPI_Comm_rank(comm, &rank);

  sprintf(parallel_filename, "%s%04d.pvti", OutputPath().c_str(), output_count);
  sprintf(serial_filename, "%04d_%d.vti", output_count, rank);

  MPI_File_open(comm, parallel_filename, MPI_MODE_WRONLY|MPI_MODE_CREATE, MPI_INFO_NULL, &file);
  MPI_File_set_size(file, 0);

  if (rank == 0) {

    buf << "<?xml version=\"1.0\"?>" << std::endl
	<< "<VTKFile type=\"PImageData\" version=\"0.1\" byte_order=\"LittleEndian\">" << std::endl
	<< "  <PImageData WholeExtent=\"";

    for (int i=0; i<3; ++i)
      buf << begin_global[i] << " " << end_global[i] << " ";

    buf << "\"" << std::endl
	<< "              GhostLevel=\"0\" Origin=\"0 0 0\" Spacing=\"";

    for (int i=0; i<3; ++i)
      buf << grid.Extent().MeshWidth()[i] << " ";

    buf << "\">" << std::endl
	<< "    <PPointData Scalars=\"" << information << "\">" << std::endl
	<< "      <PDataArray type=\"Float32\" Name=\"" << information << "\"/>" << std::endl
	<< "    </PPointData>" << std::endl;

    char* char_buf = Helper::GetCharArray(buf.str());

    MPI_File_write_shared(file, char_buf, buf.str().size(), MPI_CHAR, MPI_STATUS_IGNORE);

    delete [] char_buf;
  }

  buf.str("");

  if ((end_local-begin_local).Product() > 0) {
    buf << "    <Piece Extent=\"";

    for (int i=0; i<3; ++i)
      buf << begin_local[i] << " " << end_local[i]  << " ";

    buf << "\" Source=\"" << serial_filename << "\"/>" << std::endl;
  }

  char* char_buf = Helper::GetCharArray(buf.str());

  MPI_File_write_ordered(file, char_buf, buf.str().size(), MPI_CHAR, MPI_STATUS_IGNORE);

  delete [] char_buf;

  if (rank == 0) {

    buf.str("");

    buf << "  </PImageData>" << std::endl
	<< "</VTKFile>" << std::endl;

    char* char_buf = Helper::GetCharArray(buf.str());

    MPI_File_write_shared(file, char_buf, buf.str().size(), MPI_CHAR, MPI_STATUS_IGNORE);

    delete [] char_buf;
  }

  MPI_File_close(&file);
}

MPI_File CommMPI::CreateSerialOutputFile(const Grid& grid, MPI_Comm& comm,
					 const int& output_count, const char* information,
					 const Index& begin_global, const Index& end_global,
					 const Index& begin_local, const Index& end_local)
{
  char serial_filename[513];
  int rank;
  MPI_File file;
  std::stringstream buf;

  MPI_Comm_rank(comm_global, &rank);

  sprintf(serial_filename, "%s%04d_%d.vti", OutputPath().c_str(), output_count, rank);

  MPI_File_open(MPI_COMM_SELF, serial_filename, MPI_MODE_WRONLY|MPI_MODE_CREATE, MPI_INFO_NULL, &file);

  buf << "<?xml version=\"1.0\"?>" << std::endl
      << "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\">" << std::endl
      << "  <ImageData WholeExtent=\"";

  for (int i=0; i<3; ++i)
    buf << begin_global[i] << " " << end_global[i] << " ";

  buf << "\"" << std::endl
      << "             Origin=\"0 0 0\" Spacing=\"";

  for (int i=0; i<3; ++i)
    buf << grid.Extent().MeshWidth()[i] << " ";

  buf << "\">" << std::endl
      << "    <Piece Extent=\"";

  for (int i=0; i<3; ++i)
    buf << begin_local[i] << " " << end_local[i] << " ";

  buf << "\">" << std::endl
      << "      <PointData Scalars=\"" << information << "\">" << std::endl
      << "        <DataArray type=\"Float32\" Name=\"" << information << "\" format=\"ascii\">" << std::endl
      << "          ";

  char* char_buf = Helper::GetCharArray(buf.str());
  MPI_File_write(file, char_buf, buf.str().size(), MPI_CHAR, MPI_STATUS_IGNORE);
  delete [] char_buf;

  return file;
}

void CommMPI::FinalizeSerialOutputFile(MPI_File& file)
{
  std::stringstream buf;

  buf << std::endl
      << "        </DataArray>" << std::endl
      << "      </PointData>" << std::endl
      << "    </Piece>" << std::endl
      << "  </ImageData>" << std::endl
      << "</VTKFile>" << std::endl;

  char* char_buf = Helper::GetCharArray(buf.str());
  MPI_File_write(file, char_buf, buf.str().size(), MPI_CHAR, MPI_STATUS_IGNORE);
  delete [] char_buf;

  MPI_File_close(&file);
}

int CommMPI::GlobalRank() const
{
  int rank;
  MPI_Comm_rank(comm_global, &rank);
  return rank;
}

int CommMPI::GlobalSize() const
{
  int size;
  MPI_Comm_size(comm_global, &size);
  return size;
}

Index CommMPI::GlobalPos() const
{
  Index dims, periods, coords;
  MPI_Cart_get(comm_global, 3, dims.vec(), periods.vec(), coords.vec());
  return coords;
}

Index CommMPI::GlobalProcs() const
{
  Index dims, periods, coords;
  MPI_Cart_get(comm_global, 3, dims.vec(), periods.vec(), coords.vec());
  return dims;
}

int CommMPI::Rank(const Grid& grid) const
{
  int rank;
  MPI_Comm comm = settings.CommunicatorLocal(grid);
  assert(comm != MPI_COMM_NULL);
  MPI_Comm_rank(comm, &rank);
  return rank;
}

int CommMPI::Size(const Grid& grid) const
{
  int size;
  MPI_Comm comm = settings.CommunicatorLocal(grid);
  assert(comm != MPI_COMM_NULL);
  MPI_Comm_size(comm, &size);
  return size;
}

Index CommMPI::Pos(const Grid& grid) const
{
  Index dims, periods, coords;
  MPI_Comm comm = settings.CommunicatorLocal(grid);
  assert(comm != MPI_COMM_NULL);
  MPI_Cart_get(comm, 3, dims.vec(), periods.vec(), coords.vec());
  return coords;
}

Index CommMPI::Procs(const Grid& grid) const
{
  Index dims, periods, coords;
  MPI_Comm comm = settings.CommunicatorLocal(grid);
  assert(comm != MPI_COMM_NULL);
  MPI_Cart_get(comm, 3, dims.vec(), periods.vec(), coords.vec());
  return dims;
}

void CommMPI::InitCommMPI(const MPI_Comm& comm)
{
  int status, size, rank;
  int dims[3] = {0, 0, 0};
  int periods[3];

  for (int i=0; i<3; ++i)
    periods[i] = (BoundaryConditions()[i] == Periodic ? 1 : 0);

  MPI_Comm_size(comm, &size);
  MPI_Comm_rank(comm, &rank);

  MPI_Topo_test(comm, &status);

  if (status == MPI_CART) {

    comm_global = comm;

  }else {

    const int log2 = Helper::log_2(size);

    if (Helper::intpow(2, log2) == size) {
      for (int i=0; i<3; ++i)
        dims[i] = Helper::intpow(2, log2 / 3 + (log2%3 > i ? 1 : 0));
    }else {
      MPI_Dims_create(size, 3, dims);
    }

#ifdef OUTPUT_DEBUG
    if (rank == 0)
      std::printf("Process grid: %d %d %d\n", dims[0], dims[1], dims[2]);
#endif

    MPI_Cart_create(comm, 3, dims, periods, 1, &comm_global);

  }

  MPI_Info_create(&info);
  char key[] = "no_locks";
  char val[] = "true";
  MPI_Info_set(info, key, val);

}

CommMPI::~CommMPI()
{
  MPI_Comm_free(&comm_global);
#ifdef VMG_ONE_SIDED
  if (win_created)
    MPI_Win_free(&win);
#endif
  MPI_Info_free(&info);
}

Grid& CommMPI::GetCoarserGrid(Multigrid& multigrid)
{
  return settings.CoarserGrid(multigrid(multigrid.Level()));
}

Grid& CommMPI::GetFinerGrid(Multigrid& multigrid)
{
  return settings.FinerGrid(multigrid(multigrid.Level()));
}

std::string CommMPI::CreateOutputDirectory()
{
#ifdef HAVE_BOOST_FILESYSTEM
  std::string path, unique_path;
  std::stringstream unique_suffix;
  int suffix_counter = 0;
  char buffer[129];
  time_t rawtime;
  struct tm *timeinfo;
  int path_size;
  char* path_array;

  if (GlobalRank() == 0) {

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

    path_size = unique_path.size() + 1;
    path_array = Helper::GetCharArray(unique_path);

    MPI_Bcast(&path_size, 1, MPI_INT, 0, comm_global);

  }else {

    MPI_Bcast(&path_size, 1, MPI_INT, 0, comm_global);
    path_array = new char[path_size];

  }

  MPI_Bcast(path_array, path_size, MPI_CHAR, 0, comm_global);

  unique_path = path_array;

  delete [] path_array;

  return unique_path;

#else

  return "./";

#endif
}


#endif /* HAVE_MPI */
