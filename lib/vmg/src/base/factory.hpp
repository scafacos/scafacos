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
 * @file   factory.hpp
 * @author Julian Iseringhausen <isering@ins.uni-bonn.de>
 * @date   Tue Apr  5 20:40:41 2011
 *
 * @brief  Factory class that holds commands and arbitrary objects.
 *
 *
 */

#ifndef FACTORY_HPP_
#define FACTORY_HPP_

#include <map>
#include <string>

#include "base/object.hpp"

namespace VMG
{

class MG;

class Factory
{
public:
  virtual ~Factory();

  void Register(Object* object); ///< Registers an object
  template <class T> T& RegisterObjectStorage(std::string id, const T& val);
  template <class T> T* RegisterObjectStorageArray(std::string id, const vmg_int& size);

  Object* Get(std::string id);   ///< Returns an object.
  template <class T> T& GetObjectStorageVal(std::string id);
  template <class T> T* GetObjectStorageArray(std::string id);

  void Delete(std::string id);   ///< Deletes an object

  void PrintAvailableObjects();  ///< Prints the name of all objects that have been registered to the factory.

  bool TestObject(std::string id) const; ///< Checks whether an object exists or not.

private:
  friend class MG;

  Factory();
  std::map<std::string, Object*> object_map;
};

template <class T>
T& Factory::RegisterObjectStorage(std::string id, const T& val)
{
  Object* object = new ObjectStorage<T>(id, val);
  return object->Cast< ObjectStorage<T> >()->Val();
}

template <class T>
T* Factory::RegisterObjectStorageArray(std::string id, const vmg_int& size)
{
  Object* object = new ObjectStorageArray<T>(id, size);
  return object->Cast< ObjectStorage<T*> >()->Val();
}

template <class T>
T& Factory::GetObjectStorageVal(std::string id)
{
  return Get(id)->Cast< ObjectStorage<T> >()->Val();
}

template <class T>
T* Factory::GetObjectStorageArray(std::string id)
{
  return Get(id)->Cast< ObjectStorage<T*> >()->Val();
}

}

#endif /* FACTORY_HPP_ */
