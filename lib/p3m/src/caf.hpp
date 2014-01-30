/*
  Copyright (C) 2011 Olaf Lenz
  
  This file is part of ScaFaCoS.
  
  ScaFaCoS is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  ScaFaCoS is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>. 
*/
#ifndef _P3M_CAF_HPP
#define _P3M_CAF_HPP

#include "utils.hpp"

namespace P3M {

  /** The CAF. */
  class CAF {
  public:
    static const p3m_int max_cao = 7;

    /** Computes the \a i'th point of the CAF centered around \a r0 of
        the charge assignment function of order \a cao. */
    static p3m_float compute(p3m_int i, p3m_float r0, p3m_int cao);
    
    /** Computes the \a i'th point of the CAF centered around \a r0 of
        the charge assignment function of order \a cao. */
    static p3m_float computeDerivative(p3m_int i, p3m_float r0, p3m_int cao);

    /** Factory method to create the fitting CAF class. */
    static CAF* create(p3m_int cao, 
                       p3m_int n_interpol,
                       bool derivative = false);

    /** Abstract base class that caches cao values of the CAF. */
    class Cache {
    public:
      virtual void update(p3m_float r0) = 0;
      p3m_float* begin() { return pbegin; }
      p3m_float* end() { return pend; }
    protected:
      p3m_float *pbegin;
      p3m_float *pend;
    };

    class DirectCache : public Cache {
    public:
      DirectCache(CAF &_caf);
      virtual ~DirectCache();
      virtual void update(p3m_float r0);
    protected:
      CAF &caf;
    };

    class DerivativeCache : public Cache {
    public:
      DerivativeCache(CAF &_caf);
      virtual ~DerivativeCache();
      virtual void update(p3m_float r0);
    protected:
      CAF &caf;
    };

    CAF(p3m_int _cao, bool derivative = false);

    /* Factory method to create a matching cache. */
    virtual Cache *createCache();
    
  protected:
    friend class InterpolatedCAF;
    p3m_int cao;
    bool derivative;
  };

  class InterpolatedCAF: public CAF {
  public:
    class Cache: public CAF::Cache {
    public:
      Cache(InterpolatedCAF &_caf) : caf(_caf) {};
      virtual ~Cache() {};
      virtual void update(p3m_float r0);
    protected:
      InterpolatedCAF &caf;
    };

    InterpolatedCAF(p3m_int cao, p3m_int n_interpol, bool derivative = false);
    virtual ~InterpolatedCAF();

    /* Factory method to create a matching cache. */
    virtual CAF::Cache *createCache();

  protected:
    p3m_float* getBase(p3m_float r0);
    p3m_float* data;
    p3m_int n_interpol;
  };
}
#endif
