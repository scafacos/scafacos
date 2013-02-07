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
 * @file   helper.hpp
 * @author Julian Iseringhausen <isering@ins.uni-bonn.de>
 * @date   Tue Apr  5 21:03:47 2011
 *
 * @brief  Provides various helper functions.
 *
 */

#ifndef HELPER_HPP_
#define HELPER_HPP_

#include <cassert>
#include <cstdio>
#include <limits>
#include <sstream>
#include <string>

namespace VMG
{

class Grid;
class Vector;

namespace Helper
{
  char* GetCharArray(const std::string& str);
  std::string ReplaceWhitespaces(const char* buffer, const char* replace);

  template <class T>
  std::string ToString(const T& val)
  {
    std::stringstream str;
    str << std::scientific << val;
    return str.str();
  }

  template <class T>
  T ToVal(const char* val_str)
  {
    T val;
    std::stringstream str;
    str << val_str;
    str >> val;

    return val;
  }

  template <class T>
  T ToValWithDefault(const char* val_str, const T& def)
  {
    T val;
    std::stringstream str(val_str);
    str >> val;

    if (str.fail() || str.bad() || !str.eof()) {
#ifdef DEBUG_VERBOSE
      std::printf("VMG::Helper::ToValWithDefault: Using default value.\n");
#endif
      val = def;
    }

    return val;
  }

  /**
   * Checks a number for validity, i.e. it is neither nan nor inf.
   */
  inline bool CheckNumber(const vmg_float& number)
  {
    bool valid = true;

    if (std::numeric_limits<vmg_float>::has_quiet_NaN) {
      valid &= number != std::numeric_limits<vmg_float>::quiet_NaN();
      assert(number != std::numeric_limits<vmg_float>::quiet_NaN());
    }

    if (std::numeric_limits<vmg_float>::has_signaling_NaN) {
      valid &= number != std::numeric_limits<vmg_float>::signaling_NaN();
      assert(number != std::numeric_limits<vmg_float>::signaling_NaN());
    }

    if (std::numeric_limits<vmg_float>::has_infinity) {
      valid &= number != std::numeric_limits<vmg_float>::infinity();
      valid &= number != -1 * std::numeric_limits<vmg_float>::infinity();
      assert(number != std::numeric_limits<vmg_float>::infinity());
      assert(number != -1 * std::numeric_limits<vmg_float>::infinity());
    }

    return valid;
  }

  inline int intpow(int base, unsigned int power)
  {
    int result = 1;
    while (power != 0) {
      if (power & 1)
	result *= base;
      base *= base;
      power >>= 1;
    }
    return result;
  }

  inline vmg_float pow(vmg_float base, unsigned int power)
  {
    vmg_float result = 1.0;
    while (power != 0) {
      if (power & 1)
	result *= base;
      base *= base;
      power >>= 1;
    }
    return result;
  }

  inline unsigned int fact(unsigned int number)
  {
    unsigned int result = 1;
    for (unsigned int i=2; i<=number; ++i)
      result *= i;
    return result;
  }

  inline vmg_float pow_2(vmg_float val)
  {
    return val*val;
  }

  inline vmg_float pow_3(vmg_float val)
  {
    return val*val*val;
  }

  inline int log_2(int val)
  {
    assert(val > 0);
    int log2 = 0;
    int x = 1;

    while (x < val) {
      x <<= 1;
      ++log2;
    }

    return log2;
  }


  /**
   * Tests two arbitrary objects for equality and prints
   * a warning if they differ.
   */
  template <class T>
  bool IsEq(const T& val1, const T& val2, const char name[])
  {
    bool rval = (val1 == val2);

#ifdef DEBUG_OUTPUT
    if (!rval)
      printf("WARNING: Values are not equal (%s)\n", name);
#endif /* DEBUG_OUTPUT */

    assert(rval);

    return rval;
  }

  vmg_float InterpolateTrilinear(const Vector& point, const Grid& grid);

}

}

#endif /* HELPER_HPP_ */
