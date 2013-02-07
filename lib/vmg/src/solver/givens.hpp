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
 * @file   givens.hpp
 * @author Julian Iseringhausen <isering@ins.uni-bonn.de>
 * @date   Mon Apr 18 13:10:15 2011
 *
 * @brief  Solves system of linear equations using Givens rotations.
 *         Compare to Meister, Numerik lineare Gleichungssysteme.
 *
 */

#ifndef GIVENS_HPP_
#define GIVENS_HPP_

#include <cfloat>

#include "solver/solver.hpp"

namespace VMG
{

template<class T>
class Givens : public T
{
public:
  Givens(bool register_ = true) :
    T(register_)
  {}

  Givens(int size, bool register_ = true) :
    T(size, register_)
  {}

protected:
  void Compute();
};

template<class T>
void Givens<T>::Compute()
{
  int n = this->Size();
  vmg_float c,s,t;

  for (int i=0; i<n-1; i++)
    for (int j=i+1; j<n; j++)
      if (fabs(this->Mat(j,i)) > DBL_EPSILON) {

	t = 1.0 / sqrt(this->Mat(i,i)*this->Mat(i,i) + this->Mat(j,i)*this->Mat(j,i));
	s = t * this->Mat(j,i);
	c = t * this->Mat(i,i);

	for (int k=i; k<n; k++) {

	  t = c * this->Mat(i,k) + s * this->Mat(j,k);

	  if (k != i)
	    this->Mat(j,k) = c * this->Mat(j,k) - s * this->Mat(i,k);

	  this->Mat(i,k) = t;

	}

	t = c * this->Rhs(i) + s * this->Rhs(j);

	this->Rhs(j) = c * this->Rhs(j) - s * this->Rhs(i);

	this->Rhs(i) = t;

	this->Mat(j,i) = 0.0;

      }

  for (int i=n-1; i>=0; i--) {

    for (int j=i+1; j<n; j++)
      this->Rhs(i) -= this->Mat(i,j) * this->Sol(j);

    this->Sol(i) = this->Rhs(i) / this->Mat(i,i);

  }
}

}

#endif /* GIVENS_HPP_ */
