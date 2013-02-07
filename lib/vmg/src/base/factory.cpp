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
 * @file   factory.cpp
 * @author Julian Iseringhausen <isering@ins.uni-bonn.de>
 * @date   Tue Apr  5 20:40:05 2011
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <cassert>
#include <cstdio>
#include <iostream>

#include "base/command.hpp"
#include "base/discretization.hpp"
#include "base/factory.hpp"
#include "base/object.hpp"
#include "comm/comm.hpp"
#include "grid/multigrid.hpp"
#include "grid/tempgrid.hpp"
#include "level/level_operator.hpp"
#include "smoother/smoother.hpp"
#include "solver/solver.hpp"
#include "mg.hpp"

using namespace VMG;

Factory::Factory()
{
}

Factory::~Factory()
{
  for (std::map<std::string, Object*>::iterator iter=object_map.begin(); iter!=object_map.end(); ++iter)
    delete iter->second;
}

void Factory::Register(Object* object)
{
  Delete(object->Name());
  object_map.insert(std::make_pair(object->Name(), object));
}

Object* Factory::Get(std::string id)
{
  std::map<std::string, Object*>::iterator iter = object_map.find(id);

  if (iter != object_map.end())
    return iter->second;

  MG::GetComm()->PrintStringOnce("Error: Object %s is not registered", id.c_str());
  assert(0);

  return NULL;
}

void Factory::Delete(std::string id)
{
  std::map<std::string, Object*>::iterator iter = object_map.find(id);

  if (iter != object_map.end()) {
    delete iter->second;
    object_map.erase(iter);
  }
}

void Factory::PrintAvailableObjects()
{
  MG::GetComm()->PrintString("Registered objects:");
  for (std::map<std::string, Object*>::iterator iter=object_map.begin(); iter!=object_map.end(); ++iter)
    MG::GetComm()->PrintString("%s", iter->second->Name().c_str());
}

bool Factory::TestObject(std::string id) const
{
  return object_map.find(id) != object_map.end();
}
