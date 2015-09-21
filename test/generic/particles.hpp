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

#ifndef __PARTICLES_HPP__
#define __PARTICLES_HPP__


#include "fcs.h"

#include "common.hpp"

using namespace std;


enum particle_data_type_t {
  PDT_CHARGE_POSITIONS = 0,
  PDT_CHARGE_CHARGES = 1,
  PDT_CHARGE_POTENTIALS = 2,
  PDT_CHARGE_FIELD = 3,
#if SCAFACOS_TEST_WITH_DIPOLES
  PDT_DIPOLE_POSITIONS = 4,
  PDT_DIPOLE_MOMENTS = 5,
  PDT_DIPOLE_POTENTIALS = 6,
  PDT_DIPOLE_FIELD = 7,
#endif /* SCAFACOS_TEST_WITH_DIPOLES */
  PDT_LAST = 8,
};


struct CHARGES
{
  static const fcs_int ID = 1;

  static const fcs_int POSITION_SIZE = 3;
  static const fcs_int PROP_SIZE = 1;
  static const fcs_int POTENTIAL_SIZE = 1;
  static const fcs_int FIELD_SIZE = 3;

  static const char *POSITION_TAG;
  static const char *PROP_TAG;
  static const char *POTENTIAL_TAG;
  static const char *FIELD_TAG;
};

#if MAIN_COMPILATION_UNIT
const char *CHARGES::POSITION_TAG = "position";
const char *CHARGES::PROP_TAG = "charge";
const char *CHARGES::POTENTIAL_TAG = "potential";
const char *CHARGES::FIELD_TAG = "field";
#endif /* MAIN_COMPILATION_UNIT */

#if SCAFACOS_TEST_WITH_DIPOLES

struct DIPOLES
{
  static const fcs_int ID = 2;

  static const fcs_int POSITION_SIZE = 3;
  static const fcs_int PROP_SIZE = 3;
  static const fcs_int POTENTIAL_SIZE = 3;
  static const fcs_int FIELD_SIZE = 6;

  static const char *POSITION_TAG;
  static const char *PROP_TAG;
  static const char *POTENTIAL_TAG;
  static const char *FIELD_TAG;
};

#if MAIN_COMPILATION_UNIT
const char *DIPOLES::POSITION_TAG = "position";
const char *DIPOLES::PROP_TAG = "moment";
const char *DIPOLES::POTENTIAL_TAG = "potential";
const char *DIPOLES::FIELD_TAG = "field";
#endif /* MAIN_COMPILATION_UNIT */

#endif /* SCAFACOS_TEST_WITH_DIPOLES */


template<typename P>
struct generic_particle_t: P
{
  fcs_float position[P::POSITION_SIZE];
  fcs_float prop[P::PROP_SIZE];
  fcs_float potential[P::POTENTIAL_SIZE];
  fcs_float field[P::FIELD_SIZE];
};


template<typename P>
struct generic_particle_data_t: P
{
  fcs_int n, max_n;
  fcs_float *positions;
  fcs_float *props;
  fcs_float *potentials;
  fcs_float *field;

  generic_particle_data_t()
    :n(0), max_n(0), positions(0), props(0), potentials(0), field(0) { }

  fcs_float *positions_at(fcs_int at) { return positions + at * P::POSITION_SIZE; }
  fcs_float *props_at(fcs_int at) { return props + at * P::PROP_SIZE; }
  fcs_float *potentials_at(fcs_int at) { return potentials + at * P::POTENTIAL_SIZE; }
  fcs_float *field_at(fcs_int at) { return field + at * P::FIELD_SIZE; }

  void clear()
  {
    *this = generic_particle_data_t<P>();
  }

  bool alloc(fcs_int n_)
  {
    n = 0;
    max_n = n_;

    positions = static_cast<fcs_float *>(::malloc(max_n * P::POSITION_SIZE * sizeof(fcs_float)));
    props = static_cast<fcs_float *>(::malloc(max_n * P::PROP_SIZE * sizeof(fcs_float)));
    potentials = static_cast<fcs_float *>(::malloc(max_n * P::POTENTIAL_SIZE * sizeof(fcs_float)));
    field = static_cast<fcs_float *>(::malloc(max_n * P::FIELD_SIZE * sizeof(fcs_float)));

    return (positions && props && potentials && field);
  }

  bool realloc(fcs_int n_)
  {
    const fcs_int min_realloc = 0;

    if (n_ <= max_n) return true;

    max_n = z_max(n_, max_n + min_realloc);

    positions = static_cast<fcs_float *>(::realloc(positions, max_n * P::POSITION_SIZE * sizeof(fcs_float)));
    props = static_cast<fcs_float *>(::realloc(props, max_n * P::PROP_SIZE * sizeof(fcs_float)));
    potentials = static_cast<fcs_float *>(::realloc(potentials, max_n * P::POTENTIAL_SIZE * sizeof(fcs_float)));
    field = static_cast<fcs_float *>(::realloc(field, max_n * P::FIELD_SIZE * sizeof(fcs_float)));

    return (positions && props && potentials && field);
  }

  void free()
  {
    if (positions) ::free(positions);
    if (props) ::free(props);
    if (potentials) ::free(potentials);
    if (field) ::free(field);

    clear();
  }

  fcs_int add(fcs_int add_n)
  {
    realloc(n + add_n);

    unset(n, add_n);

    return (n += add_n) - add_n;
  }

  fcs_int append(generic_particle_data_t *particles)
  {
    fcs_int add_n = particles->n;

    realloc(n + add_n);

    values_copy<fcs_float, P::POSITION_SIZE>(positions_at(n), particles->positions, add_n);
    values_copy<fcs_float, P::PROP_SIZE>(props_at(n), particles->props, add_n);
    values_copy<fcs_float, P::POTENTIAL_SIZE>(potentials_at(n), particles->potentials, add_n);
    values_copy<fcs_float, P::FIELD_SIZE>(field_at(n), particles->field, add_n);

    return (n += add_n) - add_n;
  }

  void unset(fcs_int offset, fcs_int size)
  {
    for (fcs_int i = offset; i < offset + size; ++i)
    {
      values_set<fcs_float, P::POSITION_SIZE>(positions_at(i), 0.0);
      values_set<fcs_float, P::PROP_SIZE>(props_at(i), 0.0);
      values_set<fcs_float, P::POTENTIAL_SIZE>(potentials_at(i), NAN);
      values_set<fcs_float, P::FIELD_SIZE>(field_at(i), NAN);
    }
  }

  static void copy(generic_particle_data_t *src, fcs_int src_at, generic_particle_data_t *dst, fcs_int dst_at)
  {
    values_copy<fcs_float, P::POSITION_SIZE>(dst->positions_at(dst_at), src->positions_at(src_at), 1);
    values_copy<fcs_float, P::PROP_SIZE>(dst->props_at(dst_at), src->props_at(src_at), 1);
    values_copy<fcs_float, P::POTENTIAL_SIZE>(dst->potentials_at(dst_at), src->potentials_at(src_at), 1);
    values_copy<fcs_float, P::FIELD_SIZE>(dst->field_at(dst_at), src->field_at(src_at), 1);
  }

  void copy(fcs_int from, fcs_int to)
  {
    if (from == to) return;

    copy(this, from, this, to);
  }

  static void swap(generic_particle_data_t *p0, fcs_int at0, generic_particle_data_t *p1, fcs_int at1)
  {
    values_swap<fcs_float, P::POSITION_SIZE>(p0->positions_at(at0), p1->positions_at(at1));
    values_swap<fcs_float, P::PROP_SIZE>(p0->props_at(at0), p1->props_at(at1));
    values_swap<fcs_float, P::POTENTIAL_SIZE>(p0->potentials_at(at0), p1->potentials_at(at1));
    values_swap<fcs_float, P::FIELD_SIZE>(p0->field_at(at0), p1->field_at(at1));
  }

  void swap(fcs_int at0, fcs_int at1)
  {
    if (at0 == at1) return;

    copy(this, at0, this, at1);
  }

  static void print(generic_particle_data_t *particle_data, fcs_int x)
  {
    string s;

    print_sequence(P::POSITION_SIZE, particle_data->positions_at(x), s);
    cout << "[" << s << "]";

    cout << " ";

    print_sequence(P::PROP_SIZE, particle_data->props_at(x), s);
    cout << "[" << s << "]";

    cout << " ";

    print_sequence(P::POTENTIAL_SIZE, particle_data->potentials_at(x), s);
    cout << "[" << s << "]";

    cout << " ";

    print_sequence(P::FIELD_SIZE, particle_data->field_at(x), s);
    cout << "[" << s << "]";
  }

  void print(fcs_int x)
  {
    print(this, x);
  }
};


typedef generic_particle_data_t<CHARGES> particle_data_t, charge_particle_data_t;

#if SCAFACOS_TEST_WITH_DIPOLES
typedef generic_particle_data_t<DIPOLES> dipole_particle_data_t;
#endif /* SCAFACOS_TEST_WITH_DIPOLES */


struct particles_t
{
  fcs_int total_nparticles;
  particle_data_t particles;
#if SCAFACOS_TEST_WITH_DIPOLES
  fcs_int dipole_total_nparticles;
  dipole_particle_data_t dipole_particles;
#endif /* SCAFACOS_TEST_WITH_DIPOLES */

  particles_t()
    :total_nparticles(0), particles()
#if SCAFACOS_TEST_WITH_DIPOLES
    , dipole_total_nparticles(0), dipole_particles()
#endif /* SCAFACOS_TEST_WITH_DIPOLES */
  {
  }

  void clear()
  {
    total_nparticles = 0;
    particles.clear();
#if SCAFACOS_TEST_WITH_DIPOLES
    dipole_total_nparticles = 0;
    dipole_particles.clear();
#endif /* SCAFACOS_TEST_WITH_DIPOLES */
  }

  void free()
  {
    total_nparticles = 0;
    particles.free();
#if SCAFACOS_TEST_WITH_DIPOLES
    dipole_total_nparticles = 0;
    dipole_particles.free();
#endif /* SCAFACOS_TEST_WITH_DIPOLES */
  }
};


#endif /* __PARTICLES_HPP__ */
