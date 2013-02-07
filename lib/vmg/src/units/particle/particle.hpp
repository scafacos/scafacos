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
 * @file   particle.hpp
 * @author Julian Iseringhausen <isering@ins.uni-bonn.de>
 * @date   Sat Sep 17 14:45:56 2011
 *
 * @brief  Class to represent a particle
 *
 */

#ifndef PARTICLE_HPP_
#define PARTICLE_HPP_

#include "base/vector.hpp"

namespace VMG
{

namespace Particle
{

class Particle
{
public:
  Particle() :
    x_(0.0),
    f_(0.0),
    q_(0.0),
    p_(0.0),
    rank_(-1),
    index_(-1)
  {}

  Particle(const vmg_float* x, const vmg_float& q) :
    x_(x),
    f_(0.0),
    q_(q),
    p_(0.0),
    rank_(-1),
    index_(-1)
  {}

  Particle(const vmg_float* x, const vmg_float& q, const vmg_float& p, const vmg_float* f, const int& rank, const vmg_int& index) :
    x_(x),
    f_(f),
    q_(q),
    p_(p),
    rank_(rank),
    index_(index)
  {}

  Particle(const Vector& x, const vmg_float& q, const vmg_float& p, const Vector& f, const int& rank, const vmg_int& index) :
    x_(x),
    f_(f),
    q_(q),
    p_(p),
    rank_(rank),
    index_(index)
  {}

  Particle(const Particle& other) :
    x_(other.x_),
    f_(other.f_),
    q_(other.q_),
    p_(other.p_),
    rank_(other.rank_),
    index_(other.index_)
  {}

  Vector& Pos() {return x_;}
  const Vector& Pos() const {return x_;}

  vmg_float& Charge() {return q_;}
  const vmg_float& Charge() const {return q_;}

  vmg_float& Pot() {return p_;}
  const vmg_float& Pot() const {return p_;}

  Vector& Field() {return f_;}
  const Vector& Field() const {return f_;}

  int& Rank() {return rank_;}
  const int& Rank() const {return rank_;}

  vmg_int& Index() {return index_;}
  const vmg_int& Index() const {return index_;}

  bool operator==(const Particle& rhs)
  {
    return (this->rank_ == rhs.rank_) && (this->index_ == rhs.index_);
  }

  bool operator!=(const Particle& rhs)
  {
    return (this->rank_ != rhs.rank_) || (this->index_ != rhs.index_);
  }

private:
  Vector x_, f_;
  vmg_float q_, p_;
  int rank_;
  vmg_int index_;
};

}

}

#endif /* PARTICLE_HPP_ */
