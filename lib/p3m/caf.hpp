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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

namespace ScaFaCoS {
  namespace P3M {

    /** The CAF. */
    class CAF {
    public:
      static const fcs_int max_cao = 7;

      /** Computes the \a i'th point of the CAF centered around \a r0 of
          the charge assignment function of order \a cao. */
      static fcs_float compute(fcs_int i, fcs_float r0, fcs_int cao);
    
      /** Computes the \a i'th point of the CAF centered around \a r0 of
          the charge assignment function of order \a cao. */
      static fcs_float computeDerivative(fcs_int i, fcs_float r0, fcs_int cao);

      /** Factory method to create the fitting CAF class. */
      static CAF* create(fcs_int cao, 
                         fcs_int n_interpol,
                         bool derivative = false);

      /** Abstract base class that caches cao values of the CAF. */
      class Cache {
      public:
        virtual void update(fcs_float r0) = 0;
        fcs_float* begin() { return pbegin; }
        fcs_float* end() { return pend; }
      protected:
        fcs_float *pbegin;
        fcs_float *pend;
      };

      class DirectCache : public Cache {
      public:
        DirectCache(CAF &_caf);
        virtual ~DirectCache();
        virtual void update(fcs_float r0);
      protected:
        CAF &caf;
      };

      class DerivativeCache : public Cache {
      public:
        DerivativeCache(CAF &_caf);
        virtual ~DerivativeCache();
        virtual void update(fcs_float r0);
      protected:
        CAF &caf;
      };

      CAF(fcs_int _cao, bool derivative = false);

      /* Factory method to create a matching cache. */
      virtual Cache *createCache();
    
    protected:
      friend class InterpolatedCAF;
      fcs_int cao;
      bool derivative;
    };

    class InterpolatedCAF: public CAF {
    public:
      class Cache: public CAF::Cache {
      public:
        Cache(InterpolatedCAF &_caf) : caf(_caf) {};
        virtual ~Cache() {};
        virtual void update(fcs_float r0);
      protected:
        InterpolatedCAF &caf;
      };

      InterpolatedCAF(fcs_int cao, fcs_int n_interpol, bool derivative = false);
      virtual ~InterpolatedCAF();

      /* Factory method to create a matching cache. */
      virtual CAF::Cache *createCache();

    protected:
      fcs_float* getBase(fcs_float r0);
      fcs_float* data;
      fcs_int n_interpol;
    };
  }
}
#endif
