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
 * @file   dgesv.hpp
 * @author Julian Iseringhausen <isering@ins.uni-bonn.de>
 * @date   Mon Apr 18 13:09:29 2011
 *
 * @brief  Lapack solver for general matrices.
 *
 */

#ifndef DGESV_HPP_
#define DGESV_HPP_

#ifndef HAVE_LAPACK
#error You need LAPACK in order to use MGDSYSV
#endif

#include <cassert>
#include <cstddef>

namespace VMG
{

extern "C" {

void FC_FUNC(dgesv,DGESV) (const int* N, const int* NRHS, vmg_float* A,
			   const int* LDA, int* IPIV, vmg_float* B,
			   const int* LDB, int* INFO);

void FC_FUNC(sgesv, SGESV) (const int* N, const int* NRHS, vmg_float* A,
			    const int* LDA, int* IPIV, vmg_float* B,
			    const int* LDB, int* INFO);

} /* extern "C" */

/* By default, assume fcs_float is double.  */
#if defined(FCS_FLOAT_IS_FLOAT)
#define dgesv FC_FUNC(sgesv,SGESV)
#else
#define dgesv FC_FUNC(dgesv,DGESV)
#endif

template<class T>
class DGESV : public T
{
public:
  DGESV(bool register_ = true) :
    T(register_)
  {Init();}

  DGESV(int size, bool register_ = true) :
    T(size, register_)
  {Init();}

  virtual ~DGESV();

protected:
  void Compute();

private:
  void Init();
  void Realloc();

  int la_INFO, la_LDA, la_LDB, la_N, la_NRHS;
  int *la_IPIV;
  vmg_float *la_A, *la_B;

  int la_cur_size, la_max_size;
};

template<class T>
void DGESV<T>::Compute()
{
  // Adjust size of vectors
  this->Realloc();

  la_N = la_LDA = la_LDB = la_cur_size;

  // Rewrite matrix in column-major order
  for (int i=0; i<la_N; i++) {
    la_B[i] = this->Rhs(i);
    for (int j=0; j<la_N; j++)
      la_A[i + j*la_N] = this->Mat(i,j);
  }

  // Solve system of equations
  dgesv(&la_N, &la_NRHS, la_A, &la_LDA, la_IPIV, la_B, &la_LDB, &la_INFO);

  // Assert successful exit of solver routine
  assert(la_INFO == 0);

  // Write solution back
  for (int i=0; i<la_N; i++)
    this->Sol(i) = la_B[i];
}

template<class T>
void DGESV<T>::Init()
{
  la_cur_size = 0;
  la_max_size = 0;
  la_NRHS = 1;

  la_IPIV = NULL;
  la_A = NULL;
  la_B = NULL;
}

template<class T>
void DGESV<T>::Realloc()
{
  la_cur_size = this->Size();

  if (la_cur_size > la_max_size) {

    delete [] la_IPIV;
    delete [] la_A;
    delete [] la_B;

    la_IPIV = new int[la_cur_size];
    la_A = new vmg_float[la_cur_size*la_cur_size];
    la_B = new vmg_float[la_cur_size];

    la_max_size = la_cur_size;

  }
}

template<class T>
DGESV<T>::~DGESV()
{
  delete [] la_IPIV;
  delete [] la_A;
  delete [] la_B;
}

}

#endif /* DGESV_HPP_ */
