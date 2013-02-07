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
 * @file   grid_index_translations.cpp
 * @author Julian Iseringhausen <isering@ins.uni-bonn.de>
 * @date   Tue May 17 11:46:37 2011
 *
 * @brief  Class to convert different representations of grid
 *         indices.
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <cassert>

#include "base/helper.hpp"
#include "grid/grid_index_translations.hpp"
#include "grid/grid.hpp"
#include "grid/multigrid.hpp"

using namespace VMG;

Index GridIndexTranslations::GlobalToLocal(const Index& index_global) const
{
  const Index index_local = index_global - Father()->Global().LocalBegin() + Father()->Local().Begin();
  return index_local;
}

Index GridIndexTranslations::GlobalToFiner(const Index& index_global)
{
  const Index index_finer = 2 * index_global;
  return index_finer;
}

Index GridIndexTranslations::GlobalToCoarser(const Index& index_global)
{
  assert(index_global % 2 == 0);
  const Index index_coarser = index_global / 2;
  return index_coarser;
}

Index GridIndexTranslations::GlobalToFiner(const Index& index_global, const int& num_levels)
{
  const Index index_finer = Helper::intpow(2, num_levels) * index_global;
  return index_finer;
}

Index GridIndexTranslations::GlobalToCoarser(const Index& index_global, const int& num_levels)
{
  int quotient = Helper::intpow(2, num_levels);
  assert(index_global % quotient == 0);
  const Index index_coarser = index_global / quotient;
  return index_coarser;
}

Index GridIndexTranslations::LocalToGlobal(const Index& index_local) const
{
  const Index index_global = index_local - Father()->Local().Begin() + Father()->Global().LocalBegin();
  return index_global;
}

Index GridIndexTranslations::LocalToFiner(const Index& index_local) const
{
  assert(Father() != NULL);
  assert(Father()->Father() != NULL);

  const Multigrid& multigrid = *(Father()->Father());

  assert(Father()->Level() < multigrid.MaxLevel());

  const Index index_global_fine = GlobalToFiner(LocalToGlobal(index_local));
  const Index index_local_fine = multigrid(Father()->Level()+1).Indexing().GlobalToLocal(index_global_fine);

  return index_local_fine;
}

Index GridIndexTranslations::LocalToCoarser(const Index& index_local) const
{
  assert(Father() != NULL);
  assert(Father()->Father() != NULL);

  const Multigrid& multigrid = *(Father()->Father());

  assert(Father()->Level() > multigrid.MinLevel());

  const Index index_global_coarse = GlobalToCoarser(LocalToGlobal(index_local));
  const Index index_local_coarse = multigrid(Father()->Level()-1).Indexing().GlobalToLocal(index_global_coarse);

  return index_local_coarse;
}

Index GridIndexTranslations::FinestGlobalToLocal(const Index& index_finest) const
{
  const Index index_local = GlobalToLocal(FinestGlobalToGlobal(index_finest));
  return index_local;
}

Index GridIndexTranslations::FinestGlobalToGlobal(const Index& index_finest) const
{
  const int quotient = Helper::intpow(2, Father()->Father()->MaxLevel() - Father()->Level());
  assert(index_finest % quotient == 0);
  const Index index_global = index_finest / quotient;
  return index_global;
}

void GridIndexTranslations::FineToCoarse(Index& begin, Index& end)
{
  Index last_point = end - 1;

  for (int j=0; j<3; ++j) {

    if (begin[j] % 2 == 0)
      begin[j] /= 2;
    else
      begin[j] = (begin[j]+1) / 2;

    if (last_point[j] % 2 == 0)
      last_point[j] /= 2;
    else
      last_point[j] = (last_point[j]-1) / 2;

  }


  end = last_point + 1;
}

void GridIndexTranslations::CoarseToFine(Index& begin, Index& end, const Index& size_global)
{
  for (int i=0; i<3; ++i) {

    if (size_global[i] % 2 == 0) {

      begin[i] = 2*begin[i];
      end[i] = 2*end[i];

    }else {

      if (begin[i] > 0)
	begin[i] = 2*begin[i] - 1;
      end[i] = std::max(2*end[i] - 1, begin[i]);

    }
  }
}
