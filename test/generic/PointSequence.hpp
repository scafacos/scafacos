/*
  Copyright (C) 2012 Michael Hofmann
  
  This file is part of ScaFaCoS.
  
  ScaFaCoS is free software: you can redistribute it and/or modify
  it under the terms of the GNU Lesser Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  ScaFaCoS is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU Lesser Public License for more details.
  
  You should have received a copy of the GNU Lesser Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>. 
*/

#ifndef _POINTSEQUENCE_HPP
#define _POINTSEQUENCE_HPP


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef POINTSEQUENCEINT_T
typedef POINTSEQUENCEINT_T PointSequenceInt;
#else
typedef long PointSequenceInt;
#endif

static const PointSequenceInt primes[] = {
        2,    3,    5,    7,   11,   13,   17,   19,   23,   29,
};
static const PointSequenceInt nprimes = sizeof(primes) / sizeof(fcs_int);

template<int D, typename T>
class PointSequence {
  public:
    typedef PointSequenceInt int_t;
    PointSequence(int_t _n, int_t _step):n(_n),step(_step) { }
    PointSequence():n(0), step(0) { }
    void get_point(T *point);
  protected:
    int_t n, step;
};

template<int D, typename T>
class RandomPointSequence: public PointSequence<D, T> {
  public:
    typedef typename PointSequence<D, T>::int_t int_t;
    RandomPointSequence(int_t _n, int_t _step):PointSequence<D, T>(_n, _step) { }
    void get_next(T *point)
    {
      for (int_t i = 0; i < D; ++i)
        point[i] = (T) random() / (T) RAND_MAX;
      ++PointSequence<D, T>::step;
    }
};

template<int D, typename T, int H>
class HHPointSequence: public PointSequence<D, T> {
  public:
    typedef typename PointSequence<D, T>::int_t int_t;
    HHPointSequence(int_t _n, int_t _step):PointSequence<D, T>(_n, _step) {
/*      STATIC_ASSERT(0 <= H && H <= 1);
      STATIC_ASSERT(0 < D && D-H <= nprimes);*/
      for (int_t i = 0; i < D; ++i) {
        seeds[i] = 0;
        leaps[i] = 1;
        bases[i] = (H == 1 && i == 0)?1:primes[i - H];
        ibases[i] = 1.0 / (T) bases[i];
      }
    }
    void get_next(T *point) {
      if (PointSequence<D, T>::n <= 0) return;
      for (int_t i = 0; i < D; ++i) {
        if (H == 1 && i == 0)
          point[0] = (T) ((seeds[0] + PointSequence<D, T>::step * leaps[0]) % PointSequence<D, T>::n) / (T) PointSequence<D, T>::n;
        else {
          point[i] = 0;
          int_t v = seeds[i] + PointSequence<D, T>::step * leaps[i];
          T ib = ibases[i];
          while (v)
          {
            point[i] += (T) (v % bases[i]) * ib;
            ib *= ibases[i];
            v /= bases[i];
          }
        }
      }
      ++PointSequence<D, T>::step;
    }
  private:
    int_t seeds[D], leaps[D], bases[D];
    T ibases[D];
};

template<int D, typename T>
class HammersleyPointSequence: public HHPointSequence<D, T, 1> {
  public:
    typedef typename HHPointSequence<D, T, 1>::int_t int_t;
    HammersleyPointSequence(int_t _n, int_t _step):HHPointSequence<D, T, 1>(_n, _step) { };
};

template<int D, typename T>
class HaltonPointSequence: public HHPointSequence<D, T, 0> {
  public:
    typedef typename HHPointSequence<D, T, 0>::int_t int_t;
    HaltonPointSequence(int_t _n, int_t _step):HHPointSequence<D, T, 0>(_n, _step) { };
};


#endif /* _POINTSEQUENCE_HPP */
