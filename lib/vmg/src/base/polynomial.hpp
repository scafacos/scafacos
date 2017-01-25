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
 * @file   polynomial.hpp
 * @author Julian Iseringhausen <isering@ins.uni-bonn.de>
 * @date   Mon Nov 21 13:26:56 2011
 *
 * @brief  Class to represent and evaluate a polynomial.
 *
 */

#ifndef POLYNOMIAL_HPP_
#define POLYNOMIAL_HPP_

#include <cstdarg>
#include <vector>

namespace VMG
{

class Polynomial
{
public:
  Polynomial()
  {}

  Polynomial(const int degree,...) :
    coeff(degree+1)
  {
    va_list vl;
    va_start(vl, degree);
    for (int i=0; i<=degree; ++i)
      coeff[i] = va_arg(vl, vmg_float);
    va_end(vl);
  }

  Polynomial(const Polynomial& rhs) :
    coeff(rhs.coeff.size())
  {
    for (unsigned int i=0; i<coeff.size(); ++i)
      coeff[i] = rhs.coeff[i];
  }

  vmg_float operator()(const vmg_float& val) const
  {
    switch (coeff.size())
      {
      case 1:
	return coeff[0];
      case 2:
	return coeff[1] * val + coeff[0];
      case 3:
	return (coeff[2] * val + coeff[1]) * val + coeff[0];
      case 4:
	return ((coeff[3] * val + coeff[2]) * val + coeff[1]) * val + coeff[0];
      case 5:
	return (((coeff[4] * val + coeff[3]) * val + coeff[2]) * val + coeff[1]) * val + coeff[0];
      case 6:
	return ((((coeff[5] * val + coeff[4]) * val + coeff[3]) * val + coeff[2]) * val + coeff[1]) * val + coeff[0];
      case 7:
	return (((((coeff[6]*val+coeff[5])*val+coeff[4])*val+coeff[3])*val+coeff[2])*val+coeff[1])*val+coeff[0];
      case 8:
	return ((((((coeff[7]*val+coeff[6])*val+coeff[5])*val+coeff[4])*val+coeff[3])*val+coeff[2])*val+coeff[1])*val+coeff[0];
      default:
	vmg_float result = coeff.back();
	for (int i=coeff.size()-2; i>=0; --i)
	  result = coeff[i] + result * val;
	return result;
      }

  }

  Polynomial& operator=(const Polynomial& rhs)
  {
    coeff.resize(rhs.coeff.size());
    for (unsigned int i=0; i<coeff.size(); ++i)
      coeff[i] = rhs.coeff[i];

    return *this;
  }

private:
  std::vector<vmg_float> coeff;
};

}

#endif /* POLYNOMIAL_HPP_ */
