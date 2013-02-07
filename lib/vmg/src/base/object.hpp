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
 * @file   object.hpp
 * @author Julian Iseringhausen <isering@ins.uni-bonn.de>
 * @date   Mon Apr 18 12:22:29 2011
 *
 * @brief  Header file for the class VMG::Object.
 *
 */

#ifndef OBJECT_HPP_
#define OBJECT_HPP_

#include <cassert>
#include <string>

namespace VMG
{

class Object
{
public:
  Object() :
    registered(false)
  {
  }

  Object(std::string name_, bool register_ = true) :
    name(name_),
    registered(register_)
  {
    if (register_)
      ObjectInit();
  }

  Object(const Object& other) :
    name(other.name),
    registered(other.registered)
  {
  }

  virtual ~Object() {}

  void Register(std::string name_);

  template <class T>
  T* Cast()
  {
    T* casted = dynamic_cast<T*>(this);
    assert(casted != NULL);
    return casted;
  }

  std::string Name() {return name;}

private:
  void ObjectInit();

  std::string name;
  bool registered;
};

template <class T>
class ObjectStorage : public Object
{
public:
  ObjectStorage(const T& val) :
    val(val)
  {}

  ObjectStorage(std::string name, const T& val) :
    Object(name),
    val(val)
  {}

  T& Val() {return val;}

protected:
  T val;
};

template <class T>
class ObjectStorageArray : public ObjectStorage<T*>
{
public:
  ObjectStorageArray(const vmg_int& size) :
    ObjectStorage<T*>(new T[size])
  {}

  ObjectStorageArray(std::string name, const vmg_int& size) :
    ObjectStorage<T*>(name, new T[size])
  {}

  virtual ~ObjectStorageArray()
  {
    delete [] this->val;
  }
};

}

#endif /* OBJECT_HPP_ */
