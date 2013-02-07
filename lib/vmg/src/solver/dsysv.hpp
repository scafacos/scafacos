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
 * @file   dsysv.hpp
 * @author Julian Iseringhausen <isering@ins.uni-bonn.de>
 * @date   Mon Apr 18 13:09:54 2011
 *
 * @brief  Lapack solver for symmetric matrices
 *
 */

#ifndef DSYSV_HPP_
#define DSYSV_HPP_

#ifndef HAVE_LAPACK
#error You need LAPACK in order to use MGDSYSV
#endif

#include <cassert>
#include <cstddef>

#include "base/defs.hpp"

namespace VMG
{

extern "C"
{

void FC_FUNC(dsysv, DSYSV) (const char* UPLO, const int* N, const int* NRHS,
			    double* A, const int* LDA, int* IPIV, double* B,
			    const int* LDB, double* WORK, const int* LWORK,
			    int* INFO);

void FC_FUNC(dsyrfs, DSYRFS) (const char* UPLO, const int* N, const int* NRHS,
			      const vmg_float* A, const int* LDA, double* AF,
			      const int* LDAF, const int* IPIV, const vmg_float* B,
			      const int* LDB, double* X, const int* LDX,
			      double* FERR, double* BERR, double* WORK, int* IWORK,
			      int* INFO);

void FC_FUNC(ssysv, SSYSV) (const char* UPLO, const int* N, const int* NRHS, float* A,
			    const int* LDA, int* IPIV, float* B, const int* LDB,
			    float* WORK, const int* LWORK, int* INFO);

void FC_FUNC(ssyrfs, SSYRFS) (const char* UPLO, const int* N, const int* NRHS,
			      const vmg_float* A, const int* LDA, float* AF,
			      const int* LDAF, const int* IPIV, const vmg_float* B,
			      const int* LDB, float* X, const int* LDX, float* FERR,
			      float* BERR, double* WORK, int* IWORK, int* INFO);

} /* extern "C" */

/* By default, assume fcs_float is double.  */
#if defined(FCS_FLOAT_IS_FLOAT)
#define dsysv FC_FUNC(ssysv, SSYSV)
#define dsyrfs FC_FUNC(ssyrfs, SSYRFS)
#else
#define dsysv FC_FUNC(dsysv, DSYSV)
#define dsyrfs FC_FUNC(dsyrfs, DSYRFS)
#endif

template<class T>
class DSYSV : public T
{
public:
  DSYSV(bool register_ = true) :
    T(register_)
  {Init();}

  DSYSV(int size, bool register_ = true) :
    T(size, register_)
  {Init();}

  virtual ~DSYSV();

protected:
  void Compute();

private:
  void Init();
  void Realloc();

  char la_uplo;
  int la_n, la_nrhs, *la_ipiv, la_lwork, la_info, *la_iwork;
  vmg_float *la_A, *la_A_orig, *la_B, *la_B_orig, *la_work, *la_work2;

  int cur_lapack_size, max_lapack_size;
};

template<class T>
void DSYSV<T>::Compute()
{
  vmg_float ferr, berr;
  vmg_float opt_work;

  this->Realloc();

  for (int i=0; i<this->cur_lapack_size; i++) {
    la_B[i] = la_B_orig[i] = this->Rhs(i);
    for (int j=i; j<this->cur_lapack_size; j++)
      la_A[j + i*this->cur_lapack_size] = la_A_orig[j + i*this->cur_lapack_size] = this->Mat(i,j);
  }

  // Determine optimal size of working space
  this->la_lwork = -1;
  dsysv (&la_uplo, &la_n, &la_nrhs, la_A, &la_n, la_ipiv, la_B, &la_n, &opt_work, &la_lwork, &la_info);

  // Resize working space
  this->la_lwork = static_cast<int>(opt_work);
  delete [] this->la_work;
  this->la_work = new vmg_float[this->la_lwork];

  // Solve system
  dsysv (&la_uplo, &la_n, &la_nrhs, la_A, &la_n, la_ipiv, la_B, &la_n, la_work, &la_lwork, &la_info);

  // Improve computed solution
  dsyrfs (&la_uplo, &la_n, &la_nrhs, la_A_orig, &la_n, la_A, &la_n, la_ipiv, la_B_orig, &la_n, la_B, &la_n, &ferr, &berr, la_work2, la_iwork, &la_info);

  // Write solution back
  for (int i=0; i<this->cur_lapack_size; i++)
    this->Sol(i) = this->la_B[i];
}

template<class T>
void DSYSV<T>::Realloc()
{
  this->cur_lapack_size = this->Size();

  if (this->cur_lapack_size > this->max_lapack_size) {

    delete [] la_A;
    delete [] la_A_orig;
    delete [] la_B;
    delete [] la_B_orig;
    delete [] la_ipiv;
    delete [] la_work2;
    delete [] la_iwork;

    this->la_A = new vmg_float[this->cur_lapack_size * this->cur_lapack_size];
    this->la_A_orig = new vmg_float[this->cur_lapack_size * this->cur_lapack_size];
    this->la_B = new vmg_float[this->cur_lapack_size];
    this->la_B_orig = new vmg_float[this->cur_lapack_size];
    this->la_ipiv = new int[this->cur_lapack_size];
    this->la_work2 = new vmg_float[3 * this->cur_lapack_size];
    this->la_iwork = new int[this->cur_lapack_size];

    this->la_n = this->cur_lapack_size;

    this->max_lapack_size = this->cur_lapack_size;

  }
}

template<class T>
void DSYSV<T>::Init()
{
  this->cur_lapack_size = 0;
  this->max_lapack_size = 0;

  this->la_A = NULL;
  this->la_A_orig = NULL;
  this->la_B = NULL;
  this->la_B_orig = NULL;
  this->la_work = NULL;
  this->la_ipiv = NULL;
  this->la_work2 = NULL;
  this->la_iwork = NULL;

  this->la_nrhs = 1;
  this->la_uplo = 'L';
}

template<class T>
DSYSV<T>::~DSYSV()
{
  delete [] this->la_A;
  delete [] this->la_A_orig;
  delete [] this->la_B;
  delete [] this->la_B_orig;
  delete [] this->la_work;
  delete [] this->la_ipiv;
  delete [] this->la_work2;
  delete [] this->la_iwork;
}

}

#endif /* DSYSV_HPP_ */
