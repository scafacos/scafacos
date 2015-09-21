/*
  Copyright (C) 2011, 2012, 2013, 2014, 2015 Olaf Lenz, Michael Hofmann
  
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

#ifndef __COMMON_HPP__
#define __COMMON_HPP__


#include <stdexcept>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <math.h>


#include <mpi.h>

#include "fcs.h"

using namespace std;


#define SCAFACOS_TEST_WITH_DIPOLES  1


#define z_max(_a_, _b_)           (((_a_)>(_b_))?(_a_):(_b_))
#define z_min(_a_, _b_)           (((_a_)<(_b_))?(_a_):(_b_))

template<typename T, int N>
static inline void swap(T *x, T *y)
{
  for (int i = 0; i < N; ++i)
  {
    T t = x[i];
    x[i] = y[i];
    y[i] = t;
  }
}


template<typename T, int N>
static inline void values_swap(T *x, T *y)
{
  swap<T, N>(x, y);
}


template<typename T>
static inline void values_set(T *x, T v, fcs_int M = 1)
{
  for (fcs_int i = 0; i < M; ++i) x[i] = v;
}


template<typename T, int N>
static inline void values_set(T *x, T v, fcs_int M = 1)
{
  for (fcs_int i = 0; i < N * M; ++i) x[i] = v;
}


template<typename T, int N>
static inline void values_copy(T *x, T *y, fcs_int M)
{
  memcpy(x, y, sizeof(T) * M * N);
}


template<int N>
static inline bool values_isnan(fcs_float *x)
{
  for (int i = 0; i < N; ++i) if (isnan(x[i])) return true;
  return false;
}


template<bool> class StaticAssert;
template<> class StaticAssert<true> {};
#define STATIC_ASSERT(_a_)  {  const bool test = (_a_); StaticAssert<test> t; }

/* hack around missing MPI_OFFSET datatype prior to MPI v2.2 */
#ifndef MPI_OFFSET
#define MPI_OFFSET  ( \
  (sizeof(MPI_Offset) == sizeof(short))?MPI_SHORT:( \
  (sizeof(MPI_Offset) == sizeof(int))?MPI_INT:( \
  (sizeof(MPI_Offset) == sizeof(long))?MPI_LONG:( \
  (sizeof(MPI_Offset) == sizeof(long long))?MPI_LONG_LONG:( \
  MPI_DATATYPE_NULL)))))
#endif

/*#define FCS_ENABLE_INFO  1*/
/*#define FCS_ENABLE_DEBUG  1*/

#define MASTER(cmd)  do { if (comm_rank == MASTER_RANK) { cmd; } } while (0)

#ifdef FCS_ENABLE_INFO
#define INFO(cmd)  do { cmd; } while (0)
#else
#define INFO(cmd)  do { } while (0)
#endif

#define INFO_MASTER(cmd)  MASTER(INFO(cmd))

#ifdef FCS_ENABLE_DEBUG
#define DEBUG(cmd)  do { cmd; } while (0)
#else
#define DEBUG(cmd)  do { } while (0)
#endif

#define DEBUG_MASTER(cmd)  MASTER(DEBUG(cmd))

// Maximum string lengths
#define MAX_FILENAME_LENGTH  256
#define MAX_METHOD_LENGTH    20
#define MAX_CONF_LENGTH      1000

// How the system is decomposed
#define DECOMPOSE_ALL_ON_MASTER         0  // all particles on master node
#define DECOMPOSE_ALMOST_ALL_ON_MASTER  1  // all particles on master node BUT at least one particle on each not
#define DECOMPOSE_ATOMISTIC             2  // continuously (and equally) distributed particles
#define DECOMPOSE_RANDOM                3  // randomly distributed particles
#define DECOMPOSE_RANDOM_EQUAL          4  // randomly AND equally distributed particles
#define DECOMPOSE_DOMAIN                5  // domain decomposition

// MPI parameters
extern MPI_Comm communicator;
extern int comm_rank, comm_size;

#define MASTER_RANK  0


class ParserError : public runtime_error {
public:
  ParserError(const string &what)
  : runtime_error(what) {
  }
};


template<typename T>
static bool parse_value(string s, T &r, char *c = 0) {
  istringstream is(s);
  is.exceptions(istream::failbit | istream::badbit);
  bool r_fail = false;
  bool c_fail = false;
  T r_old = r;
  try { is >> r; } catch (istringstream::failure) { r_fail = true; r = r_old; is.clear(); }
  if (c)
  {
    char c_old = *c;
    try { is >> *c; } catch (istringstream::failure) { c_fail = true; *c = c_old; }
  }
  return (r_fail || c_fail);
}

bool parse_value(string s, fcs_int &r);
bool parse_value(string s, fcs_int &r, char &c);
bool parse_value(string s, fcs_float &r);
bool parse_value(string s, fcs_float &r, char &c);

template<typename T>
static fcs_int parse_sequence(string s, fcs_int nmax, T *rv, char *cv) {
  istringstream is(s);
  is.exceptions(istream::failbit | istream::badbit);
  fcs_int i = 0;
  string ss;
  while (i < nmax)
  {
    try { is >> ss; } catch (istringstream::failure) { break; }
    parse_value<T>(ss, rv[i], (cv)?&cv[i]:0);
    ++i;
  }
  return i;
}

fcs_int parse_sequence(string s, fcs_int nmax, fcs_int *rv, char *cv = 0);
fcs_int parse_sequence(string s, fcs_int nmax, fcs_float *rv, char *cv = 0);


template<typename T>
static void print_sequence(fcs_int n, T *v, string &s) {
  ostringstream os;
  os.precision(16);
  for (fcs_int i = 0; i < n; ++i) os << ((i > 0)?" ":"") << v[i];
  s = os.str();
}


int get_equal_distribution_count(fcs_int total_count, int size, int rank);
void get_equal_distribution(fcs_int total_count, int size, int *counts, int *displs);

void make_equal_counts_and_displs(fcs_int total_count, fcs_int ncounts, int *counts, int *displs, int *counts3, int *displs3);

typedef struct _errors_t
{
  // sums of the squared potentials/field error
  fcs_float sum_potential_error_sqr;
  fcs_float sum_energy_error_sqr;
  fcs_float sum_field_error_sqr;
  fcs_float sum_force_error_sqr;

  // maximum of the squared potentials/field error
  fcs_float max_potential_error_sqr;
  fcs_float max_energy_error_sqr;
  fcs_float max_field_error_sqr;
  fcs_float max_force_error_sqr;

  // pid with the maximal potentials/field error
  fcs_int max_potential_error_pid;
  fcs_int max_energy_error_pid;
  fcs_int max_field_error_pid;
  fcs_int max_force_error_pid;

  // absolute rms errors
  fcs_float abs_rms_potential_error;
  fcs_float abs_rms_energy_error;
  fcs_float abs_rms_field_error;
  fcs_float abs_rms_force_error;

  // relative rms errors
  fcs_float rel_rms_potential_error;
  fcs_float rel_rms_energy_error;
  fcs_float rel_rms_field_error;
  fcs_float rel_rms_force_error;

  // absolute maximum errors
  fcs_float abs_max_potential_error;
  fcs_float abs_max_energy_error;
  fcs_float abs_max_field_error;
  fcs_float abs_max_force_error;

  // relative maximum errors
  fcs_float rel_max_potential_error;
  fcs_float rel_max_energy_error;
  fcs_float rel_max_field_error;
  fcs_float rel_max_force_error;

  // energy error
  fcs_float abs_total_energy_error;
  fcs_float rel_total_energy_error;

  // total energy
  fcs_float total_energy;
  fcs_float total_energy_ref;

  bool have_potential_errors, have_field_errors;

  fcs_int total_nparticles, valid_potentials_errors, valid_field_errors;

} errors_t;

void compute_errors(errors_t *e, fcs_int nparticles, 
		    fcs_float *positions, fcs_float *charges, 
		    fcs_float *reference_potentials, fcs_float *reference_field, 
		    fcs_float *result_potentials, fcs_float *result_field, 
		    fcs_float *field_correction, fcs_float energy_correction,
		    MPI_Comm comm);
void print_errors(errors_t *e, const char *prefix = "");


#endif /* __COMMON_HPP__ */
