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
 * @file   cycle.hpp
 * @author Julian Iseringhausen <isering@ins.uni-bonn.de>
 * @date   Wed Jun 27 14:45:29 2012
 *
 * @brief  Base class for multigrid cycles.
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "base/command_list.hpp"
#include "cycles/cycle.hpp"

using namespace VMG;

void Cycle::AddCycle(CommandList& list, int numLevels, int gamma)
{
  for (int i=0; i<numLevels-1; i++) {
    list.AddCommand("Smooth", "PRESMOOTHSTEPS");
    list.AddCommand("Restrict");
  }

  list.AddCommand("Solve");

  for (int i=1; i<gamma; i++) {

    for (int j=0; j<numLevels-2; j++) {

      for (int k=0; k<=j; k++) {
        list.AddCommand("Prolongate");
        list.AddCommand("Smooth", "PRESMOOTHSTEPS");
      }

      for (int k=0; k<=j; k++) {
        list.AddCommand("Smooth", "PRESMOOTHSTEPS");
        list.AddCommand("Restrict");
      }

      list.AddCommand("Solve");

    }

    for (int j=numLevels-4; j>=0; j--) {

      for (int k=0; k<=j; k++) {
        list.AddCommand("Prolongate");
        list.AddCommand("Smooth", "POSTSMOOTHSTEPS");
      }

      for (int k=0; k<=j; k++) {
        list.AddCommand("Smooth", "PRESMOOTHSTEPS");
        list.AddCommand("Restrict");
      }

      list.AddCommand("Solve");
    }

  }

  for (int i=0; i<numLevels-1; i++) {
    list.AddCommand("Prolongate");
    list.AddCommand("Smooth", "POSTSMOOTHSTEPS");
  }
}

void Cycle::AddCycleDebug(CommandList& list, int numLevels, int gamma)
{
  for (int i=0; i<numLevels-1; i++) {
    list.AddCommand("PrintGrid", "RHS");
    list.AddCommand("Smooth", "PRESMOOTHSTEPS");
    list.AddCommand("PrintGrid", "SOL");
    list.AddCommand("PrintDefect");
    list.AddCommand("Restrict");
  }

  list.AddCommand("PrintGrid", "RHS");
  list.AddCommand("Solve");
  list.AddCommand("PrintGrid", "SOL");
  list.AddCommand("PrintDefect");

  for (int i=1; i<gamma; i++) {

    for (int j=0; j<numLevels-2; j++) {

      for (int k=0; k<=j; k++) {
        list.AddCommand("Prolongate");
        list.AddCommand("PrintGrid", "RHS");
        list.AddCommand("Smooth", "PRESMOOTHSTEPS");
        list.AddCommand("PrintGrid", "SOL");
        list.AddCommand("PrintDefect");
      }

      for (int k=0; k<=j; k++) {
        list.AddCommand("PrintGrid", "RHS");
        list.AddCommand("Smooth", "PRESMOOTHSTEPS");
        list.AddCommand("PrintGrid", "SOL");
        list.AddCommand("PrintDefect");
        list.AddCommand("Restrict");
      }

      list.AddCommand("PrintGrid", "RHS");
      list.AddCommand("Solve");
      list.AddCommand("PrintGrid", "SOL");
      list.AddCommand("PrintDefect");

    }

    for (int j=numLevels-4; j>=0; j--) {

      for (int k=0; k<=j; k++) {
        list.AddCommand("Prolongate");
        list.AddCommand("PrintGrid", "RHS");
        list.AddCommand("Smooth", "POSTSMOOTHSTEPS");
        list.AddCommand("PrintGrid", "SOL");
        list.AddCommand("PrintDefect");
      }

      for (int k=0; k<=j; k++) {
        list.AddCommand("PrintGrid", "RHS");
        list.AddCommand("Smooth", "PRESMOOTHSTEPS");
        list.AddCommand("PrintGrid", "SOL");
        list.AddCommand("PrintDefect");
        list.AddCommand("Restrict");
      }

      list.AddCommand("PrintGrid", "RHS");
      list.AddCommand("Solve");
      list.AddCommand("PrintGrid", "SOL");
      list.AddCommand("PrintDefect");
    }

  }

  for (int i=0; i<numLevels-1; i++) {
    list.AddCommand("Prolongate");
    list.AddCommand("PrintGrid", "RHS");
    list.AddCommand("Smooth", "POSTSMOOTHSTEPS");
    list.AddCommand("PrintGrid", "SOL");
    list.AddCommand("PrintDefect");
  }
}
