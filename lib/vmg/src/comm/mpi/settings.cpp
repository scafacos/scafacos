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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef HAVE_MPI
#include <mpi.h>
#ifdef HAVE_MARMOT
#include <enhancempicalls.h>
#include <sourceinfompicalls.h>
#endif
#else
#error MPI is needed to compile CommMPISettings
#endif

#include <cassert>

#include "comm/comm.hpp"
#include "comm/mpi/settings.hpp"
#include "grid/multigrid.hpp"
#include "grid/tempgrid.hpp"
#include "mg.hpp"

using namespace VMG;

Grid& VMG::MPI::Settings::FinerGrid(const Grid& grid)
{
  assert(finer_grids.find(&grid) != finer_grids.end());
  return *finer_grids.find(&grid)->second;
}

Grid& VMG::MPI::Settings::CoarserGrid(const Grid& grid)
{
  assert(coarser_grids.find(&grid) != coarser_grids.end());
  return *coarser_grids.find(&grid)->second;
}

MPI_Comm VMG::MPI::Settings::CommunicatorGlobal(const Grid& grid) const
{
  assert(communicators_global.find(grid.Level()) != communicators_global.end());
  return communicators_global.find(grid.Level())->second;
}

MPI_Comm VMG::MPI::Settings::CommunicatorLocal(const Grid& grid) const
{
  KeyUnsorted k(grid, 0);
  assert(communicators_local.find(k) != communicators_local.end());
  return communicators_local.find(k)->second;
}

VMG::MPI::DatatypesGlobal& VMG::MPI::Settings::DatatypesGlobal(const Grid& grid_old, const Grid& grid_new, const int& direction)
{
  KeyUnsorted k(grid_old, grid_new, direction);
  assert(datatypes_global.find(k) != datatypes_global.end());
  return datatypes_global.find(k)->second;
}

VMG::MPI::DatatypesLocal& VMG::MPI::Settings::DatatypesLocal(const Grid& grid)
{
  KeyUnsorted k(grid, 0);
  assert(datatypes_local.find(k) != datatypes_local.end());
  return datatypes_local.find(k)->second;
}

void VMG::MPI::Settings::ComputeSettings(Multigrid& sol, Multigrid& rhs, MPI_Comm& comm_global)
{

  std::map<const Grid*, Grid*>::const_iterator grid_iter;

  Comm& comm = *MG::GetComm();

  /*
   * Create coarser grids
   */
  for (int i=sol.MaxLevel(); i>sol.MinLevel(); --i) {

    TempGrid* temp_grid = new TempGrid();
    temp_grid->SetPropertiesToCoarser(sol(i), comm.BoundaryConditions());

    if (temp_grid->Global().LocalBegin().IsComponentwiseGreaterOrEqual(sol(i-1).Global().LocalBegin()) &&
	temp_grid->Global().LocalBegin().IsComponentwiseLessOrEqual(sol(i-1).Global().LocalEnd()) &&
	temp_grid->Global().LocalEnd().IsComponentwiseGreaterOrEqual(sol(i-1).Global().LocalBegin()) &&
	temp_grid->Global().LocalEnd().IsComponentwiseLessOrEqual(sol(i-1).Global().LocalEnd())) {
      delete temp_grid;
      coarser_grids.insert(std::make_pair(&sol(i), &sol(i-1)));
    }else {
      coarser_grids.insert(std::make_pair(&sol(i), temp_grid));
    }

  }

  for (int i=rhs.MaxLevel(); i>rhs.MinLevel(); --i) {

    TempGrid* temp_grid = new TempGrid();
    temp_grid->SetPropertiesToCoarser(rhs(i), comm.BoundaryConditions());

    if (temp_grid->Global().LocalBegin().IsComponentwiseGreaterOrEqual(sol(i-1).Global().LocalBegin()) &&
	temp_grid->Global().LocalBegin().IsComponentwiseLessOrEqual(sol(i-1).Global().LocalEnd()) &&
	temp_grid->Global().LocalEnd().IsComponentwiseGreaterOrEqual(sol(i-1).Global().LocalBegin()) &&
	temp_grid->Global().LocalEnd().IsComponentwiseLessOrEqual(sol(i-1).Global().LocalEnd())) {
      delete temp_grid;
      coarser_grids.insert(std::make_pair(&rhs(i), &rhs(i-1)));
    }else {
      coarser_grids.insert(std::make_pair(&rhs(i), temp_grid));
    }

  }

  /*
   * Create finer grids
   */
  for (int i=sol.MinLevel(); i<sol.MaxLevel(); ++i) {

    TempGrid* temp_grid = new TempGrid();
    temp_grid->SetPropertiesToFiner(sol(i), comm.BoundaryConditions());

    if (temp_grid->Global().LocalBegin() == sol(i+1).Global().LocalBegin() &&
	temp_grid->Global().LocalEnd() == sol(i+1).Global().LocalEnd()) {
      delete temp_grid;
      finer_grids.insert(std::make_pair(&sol(i), &sol(i+1)));
    }else {
      finer_grids.insert(std::make_pair(&sol(i), temp_grid));
    }

  }

  for (int i=rhs.MinLevel(); i<rhs.MaxLevel(); ++i) {

    TempGrid* temp_grid = new TempGrid();
    temp_grid->SetPropertiesToFiner(rhs(i), comm.BoundaryConditions());

    if (temp_grid->Global().LocalBegin() == rhs(i+1).Global().LocalBegin() &&
	temp_grid->Global().LocalEnd() == rhs(i+1).Global().LocalEnd()) {
      delete temp_grid;
      finer_grids.insert(std::make_pair(&rhs(i), &rhs(i+1)));
    }else {
      finer_grids.insert(std::make_pair(&rhs(i), temp_grid));
    }

  }

  /*
   * Create global communicators
   */
  for (int i=sol.MinLevel()+1; i<sol.MaxLevel(); ++i)
    CreateGlobalCommunicator(comm_global, &sol(i), &CoarserGrid(sol(i+1)), &FinerGrid(sol(i-1)));

  CreateGlobalCommunicator(comm_global, &sol(sol.MinLevel()), &CoarserGrid(sol(sol.MinLevel()+1)));
  CreateGlobalCommunicator(comm_global, &sol(sol.MaxLevel()), &FinerGrid(sol(sol.MaxLevel()-1)));

  MPI_Comm my_comm_global;
  MPI_Comm_dup(comm_global, &my_comm_global);
  communicators_global.insert(std::make_pair(0, my_comm_global));

  /*
   * Create local communicators
   */
  for (int i=sol.MinLevel(); i<=sol.MaxLevel(); ++i)
    CreateLocalCommunicator(comm_global, sol(i));

  for (int i=sol.MinLevel(); i<sol.MaxLevel(); ++i)
    CreateLocalCommunicator(comm_global, FinerGrid(sol(i)));

  for (int i=sol.MaxLevel(); i>sol.MinLevel(); --i)
    CreateLocalCommunicator(comm_global, CoarserGrid(sol(i)));

  if (MG::GetFactory().TestObject("PARTICLE_NEAR_FIELD_CELLS"))
    CreateLocalCommunicator(comm_global, comm.GetParticleGrid());

  /*
   * Create single grid datatypes
   */
  for (int i=sol.MinLevel(); i<=sol.MaxLevel(); ++i)
    datatypes_local.insert(std::make_pair(KeyUnsorted(sol(i), 0), VMG::MPI::DatatypesLocal(sol(i), CommunicatorLocal(sol(i)), true)));

  for (grid_iter=finer_grids.begin(); grid_iter!=finer_grids.end(); ++grid_iter)
    datatypes_local.insert(std::make_pair(KeyUnsorted(*grid_iter->second, 0), VMG::MPI::DatatypesLocal(*grid_iter->second, CommunicatorLocal(*grid_iter->second), true)));

  for (grid_iter=coarser_grids.begin(); grid_iter!=coarser_grids.end(); ++grid_iter)
    datatypes_local.insert(std::make_pair(KeyUnsorted(*grid_iter->second, 0), VMG::MPI::DatatypesLocal(*grid_iter->second, CommunicatorLocal(*grid_iter->second), true)));

  if (MG::GetFactory().TestObject("PARTICLE_NEAR_FIELD_CELLS"))
    datatypes_local.insert(std::make_pair(KeyUnsorted(comm.GetParticleGrid(), 0), VMG::MPI::DatatypesLocal(comm.GetParticleGrid(), CommunicatorLocal(comm.GetParticleGrid()), true)));

  /*
   * Create two grid datatypes
   */
  for (int i=sol.MinLevel(); i<sol.MaxLevel(); ++i) {
    AddDatatypeGlobal(sol(i), CoarserGrid(sol(i+1)), 0);
    AddDatatypeGlobal(CoarserGrid(sol(i+1)), sol(i), 1);
  }

  for (int i=sol.MaxLevel(); i>sol.MinLevel(); --i) {
    AddDatatypeGlobal(sol(i), FinerGrid(sol(i-1)), 0);
    AddDatatypeGlobal(FinerGrid(sol(i-1)), sol(i), 1);
  }
}

void VMG::MPI::Settings::CreateGlobalCommunicator(MPI_Comm& comm_global, const Grid* grid_1, const Grid* grid_2, const Grid* grid_3)
{
  int rank;
  MPI_Comm comm_new;

  const bool in_communicator = (grid_1->Global().LocalSize().Product() > 0) ||
                               (grid_2 && grid_2->Global().LocalSize().Product() > 0) ||
                               (grid_3 && grid_3->Global().LocalSize().Product() > 0);

  MPI_Comm_rank(comm_global, &rank);

  if (in_communicator) {
    Index dims, periods, coords;
    MPI_Comm comm_temp;
    MPI_Cart_get(comm_global, 3, dims.vec(), periods.vec(), coords.vec());
    MPI_Comm_split(comm_global, 1, rank, &comm_temp);
    dims = GlobalDims(comm_temp, coords);
    MPI_Cart_create(comm_temp, 3, dims.vec(), periods.vec(), 0, &comm_new);
    MPI_Comm_free(&comm_temp);
  }else {
    MPI_Comm_split(comm_global, MPI_UNDEFINED, rank, &comm_new);
  }

  communicators_global.insert(std::make_pair(grid_1->Level(), comm_new));
}

void VMG::MPI::Settings::CreateLocalCommunicator(MPI_Comm& comm_global, const Grid& grid)
{
  int rank, comm_equal;
  MPI_Comm comm_new;
  std::set<MPI_Comm>::iterator iter;

  MPI_Comm_rank(comm_global, &rank);

  if (grid.Global().LocalSize().Product() > 0) {
    Index dims, periods, coords;
    MPI_Comm comm_temp;
    MPI_Cart_get(comm_global, 3, dims.vec(), periods.vec(), coords.vec());
    MPI_Comm_split(comm_global, 1, rank, &comm_temp);
    dims = GlobalDims(comm_temp, coords);
    MPI_Cart_create(comm_temp, 3, dims.vec(), periods.vec(), 0, &comm_new);
    MPI_Comm_free(&comm_temp);
  }else {
    MPI_Comm_split(comm_global, MPI_UNDEFINED, rank, &comm_new);
  }

  if (comm_new != MPI_COMM_NULL) {
    for (iter=communicators_local_unique.begin(); iter!=communicators_local_unique.end(); ++iter) {
      if (*iter != MPI_COMM_NULL) {
	MPI_Comm_compare(comm_new, *iter, &comm_equal);
	assert(comm_equal != MPI_SIMILAR);
	if (comm_equal == MPI_IDENT || comm_equal == MPI_CONGRUENT) {
	  MPI_Comm_free(&comm_new);
	  comm_new = *iter;
	  break;
	}
      }
    }
  }

  std::pair<std::set<MPI_Comm>::iterator, bool> insert_result = communicators_local_unique.insert(comm_new);
  communicators_local.insert(std::make_pair(KeyUnsorted(grid, 0), *insert_result.first));
}

void VMG::MPI::Settings::AddDatatypeGlobal(const Grid& grid_old, const Grid& grid_new, const int& direction)
{
  MPI_Comm comm = CommunicatorGlobal(grid_old);
  bool dt_is_new = true;

    // Insert into map
  std::pair< std::map<VMG::MPI::KeyUnsorted, VMG::MPI::DatatypesGlobal>::iterator, bool > insert_result =
    datatypes_global.insert(std::make_pair(VMG::MPI::KeyUnsorted(grid_old, grid_new, direction), VMG::MPI::DatatypesGlobal()));
  VMG::MPI::DatatypesGlobal& dt_global = insert_result.first->second;
  dt_is_new = insert_result.second;


  if (comm != MPI_COMM_NULL) {

    Index begin, end, offset_old, offset_new;
    int rank, size;

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    std::vector<int> buffer;
    buffer.resize(6*size);

    // Compute local offset for ghost cells
    for (int i=0; i<3; ++i) {
      offset_old[i] = (grid_old.Local().HaloSize1()[i] > 0 ? grid_old.Local().Begin()[i] : 0);
      offset_new[i] = (grid_new.Local().HaloSize1()[i] > 0 ? grid_new.Local().Begin()[i] : 0);
    }

    // Publish which grid part this process can offer
    if (&grid_old == &grid_new) {
      for (int i=0; i<6; ++i)
	buffer[6*rank+i] = 0;
    }else {
      for (int i=0; i<3; ++i) {
	buffer[6*rank+i] = grid_old.Global().LocalBegin()[i];
	buffer[6*rank+i+3] = grid_old.Global().LocalEnd()[i];
      }
    }

    MPI_Allgather(MPI_IN_PLACE, 6, MPI_INT, &buffer.front(), 6, MPI_INT, comm);

    if (dt_is_new) {

      // Decide who offers a useful grid part
      for (int i=0; i<size; ++i) {
	for (int j=0; j<3; ++j) {
	  begin[j] = buffer[6*i+j];
	  end[j] = buffer[6*i+j+3];
	}

	begin = begin.Clamp(grid_new.Global().LocalBegin(), grid_new.Global().LocalEnd());
	end = end.Clamp(grid_new.Global().LocalBegin(), grid_new.Global().LocalEnd());

	if ((end-begin).Product() > 0) {
	  // This process has a useful part
	  dt_global.Receive().push_back(VMG::MPI::Datatype(grid_new.Local().SizeTotal(),
							   end - begin,
							   begin - grid_new.Global().LocalBegin() + offset_new,
							   i, 0, 0, true));
	}
      }
    }

    // Publish which grid parts this process needs
    for (int i=0; i<3; ++i) {
      buffer[6*rank+i] = grid_new.Global().LocalBegin()[i];
      buffer[6*rank+i+3] = grid_new.Global().LocalEnd()[i];
    }

    MPI_Allgather(MPI_IN_PLACE, 6, MPI_INT, &buffer.front(), 6, MPI_INT, comm);

    if (dt_is_new) {

      // Decide who needs a part of my grid
      for (int i=0; i<size; ++i) {

	if ((i == rank) && (&grid_old == &grid_new))
	  continue;

	for (int j=0; j<3; ++j) {
	  begin[j] = buffer[6*i+j];
	  end[j] = buffer[6*i+j+3];
	}

	begin = begin.Clamp(grid_old.Global().LocalBegin(), grid_old.Global().LocalEnd());
	end = end.Clamp(grid_old.Global().LocalBegin(), grid_old.Global().LocalEnd());

	if ((end-begin).Product() > 0) {
	  // This process needs one of my parts
 	  dt_global.Send().push_back(VMG::MPI::Datatype(grid_old.Local().SizeTotal(),
							end - begin,
							begin - grid_old.Global().LocalBegin() + offset_old,
							i, 0, 0, true));
	}
      }
    }
  }
}

MPI_Datatype& VMG::MPI::Settings::Datatype(const Index& begin, const Index& end,
					   const Index& size_local, const Index& size_global,
					   const int& level)
{
  KeyUnsorted k(begin, end, size_local, size_global, level, 0);
  std::map<KeyUnsorted, MPI_Datatype>::iterator iter = datatypes.find(k);

  if (iter != datatypes.end())
    return iter->second;

  MPI_Datatype dt;
  Index sizes =  size_local;
  Index subsizes = end - begin;
  Index starts = begin;

  MPI_Type_create_subarray(3, sizes.vec(), subsizes.vec(), starts.vec(), MPI_ORDER_C, MPI_DOUBLE, &dt);
  MPI_Type_commit(&dt);

  return datatypes.insert(std::make_pair(k, dt)).first->second;
}

Index VMG::MPI::Settings::GlobalDims(MPI_Comm comm, Index pos)
{
  std::set<int> unique_set[3];
  Index dims;

  int size;
  MPI_Comm_size(comm, &size);

  int* coordinates = new int[3*size];

  MPI_Allgather(pos.vec(), 3, MPI_INT, coordinates, 3, MPI_INT, comm);

  for (int i=0; i<size; ++i)
    for (int j=0; j<3; ++j)
      unique_set[j].insert(coordinates[3*i+j]);

  for (int j=0; j<3; ++j)
    dims[j] = static_cast<int>(unique_set[j].size());

  delete [] coordinates;

  return dims;
}

VMG::MPI::Settings::Settings()
{
}

VMG::MPI::Settings::~Settings()
{
  std::map<int, MPI_Comm>::iterator iter_comm_global;
  for (iter_comm_global=communicators_global.begin(); iter_comm_global!=communicators_global.end(); ++iter_comm_global)
    if (iter_comm_global->second != MPI_COMM_NULL)
      MPI_Comm_free(&iter_comm_global->second);

  /*
   * We simply copied some communicators so we have to make sure that we free
   * each communicator exactly once
   */
  std::set<MPI_Comm>::iterator iter_comm_set;
  for (iter_comm_set=communicators_local_unique.begin(); iter_comm_set!=communicators_local_unique.end(); ++iter_comm_set)
    if (*iter_comm_set != MPI_COMM_NULL) {
      MPI_Comm comm_temp = *iter_comm_set;
      MPI_Comm_free(&comm_temp);
    }

  std::map<KeyUnsorted, MPI_Datatype>::iterator iter_dt;
  for (iter_dt=datatypes.begin(); iter_dt!=datatypes.end(); ++iter_dt)
    if (iter_dt->second != MPI_DATATYPE_NULL)
      MPI_Type_free(&iter_dt->second);

  std::map<const Grid*, Grid*>::iterator iter_grid;
  for (iter_grid=finer_grids.begin(); iter_grid!=finer_grids.end(); ++iter_grid)
    if (iter_grid->second->Father() == NULL)
      delete iter_grid->second;

  for (iter_grid=coarser_grids.begin(); iter_grid!=coarser_grids.end(); ++iter_grid)
    if (iter_grid->second->Father() == NULL)
      delete iter_grid->second;
}
