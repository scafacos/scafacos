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
 * @file   grid_properties.hpp
 * @author Julian Iseringhausen <isering@ins.uni-bonn.de>
 * @date   Mon Apr 18 12:54:05 2011
 *
 * @brief  VMG::GlobalIndices, VMG::LocalIndices and VMG::SpatialExtent
 *
 */

#ifndef GRID_PROPERTIES_HPP_
#define GRID_PROPERTIES_HPP_

#include "base/defs.hpp"
#include "base/index.hpp"
#include "base/vector.hpp"

namespace VMG
{

class GlobalIndices
{
public:
  GlobalIndices() :
    local_begin(0), local_end(0), local_size(0),
    global_finer_begin(0), global_finer_end(0), global_finer_size(0),
    local_finer_begin(0), local_finer_end(0), local_finer_size(0),
    finest_abs_begin(0), finest_abs_end(0), finest_abs_size(0),
    global_size(0),
    boundary(EmptyGrid)
  {}

  GlobalIndices(const GlobalIndices& other) :
    local_begin(other.local_begin), local_end(other.local_end), local_size(other.local_size),
    global_finer_begin(other.global_finer_begin), global_finer_end(other.global_finer_end), global_finer_size(other.global_finer_size),
    local_finer_begin(other.local_finer_begin), local_finer_end(other.local_finer_end), local_finer_size(other.local_finer_size),
    finest_abs_begin(other.finest_abs_begin), finest_abs_end(other.finest_abs_end), finest_abs_size(other.finest_abs_size),
    global_size(other.global_size),
    boundary(other.boundary)
  {}

  Index& LocalBegin() {return local_begin;}
  Index& LocalEnd() {return local_end;}
  Index& LocalSize() {return local_size;}

  const Index& LocalBegin() const {return local_begin;}
  const Index& LocalEnd() const {return local_end;}
  const Index& LocalSize() const {return local_size;}

  Index& GlobalFinerBegin() {return global_finer_begin;}
  Index& GlobalFinerEnd() {return global_finer_end;}
  Index& GlobalFinerSize() {return global_finer_size;}

  const Index& GlobalFinerBegin() const {return global_finer_begin;}
  const Index& GlobalFinerEnd() const {return global_finer_end;}
  const Index& GlobalFinerSize() const {return global_finer_size;}

  Index& LocalFinerBegin() {return local_finer_begin;}
  Index& LocalFinerEnd() {return local_finer_end;}
  Index& LocalFinerSize() {return local_finer_size;}

  const Index& LocalFinerBegin() const {return local_finer_begin;}
  const Index& LocalFinerEnd() const {return local_finer_end;}
  const Index& LocalFinerSize() const {return local_finer_size;}

  Index& FinestAbsBegin() {return finest_abs_begin;}
  Index& FinestAbsEnd() {return finest_abs_end;}
  Index& FinestAbsSize() {return finest_abs_size;}

  const Index& FinestAbsBegin() const {return finest_abs_begin;}
  const Index& FinestAbsEnd() const {return finest_abs_end;}
  const Index& FinestAbsSize() const {return finest_abs_size;}

  Index& GlobalSize() {return global_size;}
  const Index& GlobalSize() const {return global_size;}

  BT& BoundaryType() {return boundary;}
  const BT& BoundaryType() const {return boundary;}

private:
  Index local_begin, local_end, local_size;
  Index global_finer_begin, global_finer_end, global_finer_size;
  Index local_finer_begin, local_finer_end, local_finer_size;
  Index finest_abs_begin, finest_abs_end, finest_abs_size;
  Index global_size;
  BT boundary;
};

class LocalIndices
{
public:
  LocalIndices() :
    begin(0), end(0), size(0), size_total(0),
    halo_begin_1(0), halo_end_1(0), halo_size_1(0),
    halo_begin_2(0), halo_end_2(0), halo_size_2(0),
    boundary_begin_1(0), boundary_end_1(0), boundary_size_1(0),
    boundary_begin_2(0), boundary_end_2(0), boundary_size_2(0),
    finer_begin(0), finer_end(0), finer_size(0)
  {}

  LocalIndices(const LocalIndices& other) :
    begin(other.begin), end(other.end), size(other.size), size_total(other.size_total),
    halo_begin_1(other.halo_begin_1), halo_end_1(other.halo_end_1), halo_size_1(other.halo_size_1),
    halo_begin_2(other.halo_begin_2), halo_end_2(other.halo_end_2), halo_size_2(other.halo_size_2),
    boundary_begin_1(other.boundary_begin_1), boundary_end_1(other.boundary_end_1), boundary_size_1(other.boundary_size_1),
    boundary_begin_2(other.boundary_begin_2), boundary_end_2(other.boundary_end_2), boundary_size_2(other.boundary_size_2),
    finer_begin(other.finer_begin), finer_end(other.finer_end), finer_size(other.finer_size)
  {}

  Index& Begin() {return begin;}    	                 ///< Index of first local grid point
  Index& End() {return end;}                             ///< Index of first non-local grid point
  Index& Size() {return size;}                           ///< Local grid size excluding halo
  Index& SizeTotal() {return size_total;}                ///< Local grid size including halo and boundary

  const Index& Begin() const {return begin;}             ///< Index of first local grid point
  const Index& End() const {return end;}                 ///< Index of first non-local grid point
  const Index& Size() const {return size;}               ///< Local grid size excluding halo
  const Index& SizeTotal() const {return size_total;}    ///< Local grid size including halo and boundary

  Index& HaloBegin1() {return halo_begin_1;}             ///< Index of first halo point
  Index& HaloEnd1() {return halo_end_1;}                 ///< Index of first non-halo point
  Index& HaloSize1() {return halo_size_1;}               ///< Size of halo

  Index& HaloBegin2() {return halo_begin_2;}             ///< Index of first halo point
  Index& HaloEnd2() {return halo_end_2;}                 ///< Index of first non-halo point
  Index& HaloSize2() {return halo_size_2;}               ///< Size of halo

  const Index& HaloBegin1() const {return halo_begin_1;} ///< Index of first halo point
  const Index& HaloEnd1() const {return halo_end_1;}     ///< Index of first non-halo point
  const Index& HaloSize1() const {return halo_size_1;}   ///< Size of halo

  const Index& HaloBegin2() const {return halo_begin_2;} ///< Index of first halo point
  const Index& HaloEnd2() const {return halo_end_2;}     ///< Index of first non-halo point
  const Index& HaloSize2() const {return halo_size_2;}   ///< Size of halo

  Index& BoundaryBegin1() {return boundary_begin_1;}     ///< Index of first boundary point
  Index& BoundaryEnd1() {return boundary_end_1;}         ///< Index of first non-boundary point
  Index& BoundarySize1() {return boundary_size_1;}       ///< Size of boundary

  Index& BoundaryBegin2() {return boundary_begin_2;}     ///< Index of first boundary point
  Index& BoundaryEnd2() {return boundary_end_2;}         ///< Index of first non-boundary point
  Index& BoundarySize2() {return boundary_size_2;}       ///< Size of boundary

  const Index& BoundaryBegin1() const {return boundary_begin_1;} ///< Index of first boundary point
  const Index& BoundaryEnd1() const {return boundary_end_1;}     ///< Index of first non-boundary point
  const Index& BoundarySize1() const {return boundary_size_1;}   ///< Size of boundary

  const Index& BoundaryBegin2() const {return boundary_begin_2;} ///< Index of first boundary point
  const Index& BoundaryEnd2() const {return boundary_end_2;}     ///< Index of first non-boundary point
  const Index& BoundarySize2() const {return boundary_size_2;}   ///< Size of boundary

  Index& FinerBegin() {return finer_begin;}
  Index& FinerEnd() {return finer_end;}
  Index& FinerSize() {return finer_size;}

  const Index& FinerBegin() const {return finer_begin;}
  const Index& FinerEnd() const {return finer_end;}
  const Index& FinerSize() const {return finer_size;}

  bool HasFinerGrid() const {return static_cast<bool>(finer_size.Product());}

private:
  Index begin, end, size, size_total;
  Index halo_begin_1, halo_end_1, halo_size_1;
  Index halo_begin_2, halo_end_2, halo_size_2;
  Index boundary_begin_1, boundary_end_1, boundary_size_1;
  Index boundary_begin_2, boundary_end_2, boundary_size_2;
  Index finer_begin, finer_end, finer_size;
};

class SpatialExtent
{
public:
  SpatialExtent() :
    begin(0.0),
    end(0.0),
    size(0.0),
    mesh_width(0.0)
  {}

  SpatialExtent(const SpatialExtent& other) :
    begin(other.begin),
    end(other.end),
    size(other.size),
    mesh_width(other.mesh_width)
  {}

  Vector& Begin() {return begin;}
  const Vector& Begin() const {return begin;}

  Vector& End() {return end;}
  const Vector& End() const {return end;}

  Vector& Size() {return size;}
  const Vector& Size() const {return size;}

  Vector& MeshWidth() {return mesh_width;}
  const Vector& MeshWidth() const {return mesh_width;}

private:
  Vector begin, end, size;
  Vector mesh_width;
};

}

#endif /* GRID_PROPERTIES_HPP_ */
