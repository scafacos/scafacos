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
 * @file   grid_index_translations.hpp
 * @author Julian Iseringhausen <isering@ins.uni-bonn.de>
 * @date   Tue May 17 11:46:37 2011
 *
 * @brief  Class to convert different representations of grid
 *         indices.
 *
 */

#ifndef GRID_INDEX_TRANSLATIONS_HPP_
#define GRID_INDEX_TRANSLATIONS_HPP_

#include "base/defs.hpp"
#include "base/index.hpp"

namespace VMG
{

class Comm;
class Grid;

class GridIndexTranslations
{
public:
  GridIndexTranslations(const Grid* father_) :
    father(father_)
  {}

  Index GlobalToLocal(const Index& index_global) const;
  static Index GlobalToFiner(const Index& index_global);
  static Index GlobalToCoarser(const Index& index_global);

  static Index GlobalToFiner(const Index& index_global, const int& num_levels);
  static Index GlobalToCoarser(const Index& index_global, const int& num_levels);

  Index LocalToGlobal(const Index& index_local) const;
  Index LocalToFiner(const Index& index_local) const;
  Index LocalToCoarser(const Index& index_local) const;

  Index FinestGlobalToLocal(const Index& index_finest) const;
  Index FinestGlobalToGlobal(const Index& index_finest) const;

  static void FineToCoarse(Index& begin, Index& end);
  static void CoarseToFine(Index& begin, Index& end, const Index& size_global);

  static Index EndOffset(const Boundary& boundary)
  {
    Index offset;

    for (int i=0; i<3; ++i)
      offset[i] = (boundary[i] == Open || boundary[i] == Dirichlet ? 1 : 0);

    return offset;
  }

protected:
  const Grid* Father() const {return father;}

private:
  const Grid* father;
};

}

#endif /* GRID_INDEX_TRANSLATIONS_HPP_ */
