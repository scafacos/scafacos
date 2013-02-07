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
 * @file   object.cpp
 * @author Julian Iseringhausen <isering@ins.uni-bonn.de>
 * @date   Mon Apr 18 12:21:55 2011
 *
 * @brief  Source file for the class VMG::Object.
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "base/object.hpp"
#include "mg.hpp"

using namespace VMG;

void Object::Register(std::string name_)
{
  if (!registered) {
    name = name_;
    MG::GetFactory().Register(this);
    registered = true;
  }else
    assert(0 == "This object has already been registered.");
}

void Object::ObjectInit()
{
  MG::GetFactory().Register(this);
}
