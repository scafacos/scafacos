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

#ifndef _GENERATOR_HPP
#define _GENERATOR_HPP

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "fcs.h"

#include "rapidxml/rapidxml.hpp"

#include "common.hpp"
#include "particles.hpp"

#define POINTSEQUENCEINT_T  fcs_int
#include "PointSequence.hpp"

using namespace std;
using namespace rapidxml;


static const char GENERATE_TAG[] = "generate";
static const char BINARY_TAG[] = "binary";
static const char PORTABLE_TAG[] = "portable";
static const char REFERENCES_TAG[] = "references";
static const char PARTICLE_TAG[] = "particle";


class ParticleSource
{
public:
  ParticleSource();

  void set_decomposition(fcs_int decomposition_) { base_params.decomposition = decomposition_; }

  virtual bool read_config(xml_node<> *config_node, const char *basename, int comm_size, int comm_rank, MPI_Comm comm) = 0;
  virtual void broadcast_config(int root, int comm_size, int comm_rank, MPI_Comm comm);
  virtual void print_config(const char *prefix, int comm_size, int comm_rank, MPI_Comm comm);
  virtual bool write_config(xml_document<> *doc, xml_node<> *config_node, int comm_size, int comm_rank, MPI_Comm comm) = 0;

  bool all_on_master() { return (base_params.decomposition == DECOMPOSE_ALL_ON_MASTER || base_params.decomposition == DECOMPOSE_ALMOST_ALL_ON_MASTER); }

  bool have(particle_data_type_t pdt) { return base_params.haves[pdt]; };

  fcs_int get_total_nparticles() { return base_params.total_nparticles; }
  virtual fcs_int get_local_nparticles(int comm_size, int comm_rank, MPI_Comm comm) = 0;
  virtual bool make_local_particles(particle_data_t *particle_data, int comm_size, int comm_rank, MPI_Comm comm) = 0;

#if SCAFACOS_TEST_WITH_DIPOLES
  fcs_int get_dipole_total_nparticles() { return base_params.dipole_total_nparticles; }
  virtual fcs_int get_dipole_local_nparticles(int comm_size, int comm_rank, MPI_Comm comm) = 0;
  virtual bool make_dipole_local_particles(dipole_particle_data_t *particle_data, int comm_size, int comm_rank, MPI_Comm comm) = 0;
#endif /* SCAFACOS_TEST_WITH_DIPOLES */

protected:
  struct {
    fcs_int decomposition;

    fcs_int total_nparticles;
#if SCAFACOS_TEST_WITH_DIPOLES
    fcs_int dipole_total_nparticles;
#endif /* SCAFACOS_TEST_WITH_DIPOLES */

    bool haves[PDT_LAST];

  } base_params;
};


class generator_type
{
  public:
    typedef enum {
      TYPE_NONE,
      TYPE_RANDOM,
      TYPE_HAMMERSLEY,
      TYPE_HALTON,
      TYPE_CONST,
      TYPE_ALTERNATE,
      TYPE_GRID
    } type_t;

    static const char *tostr(generator_type::type_t type) {
      switch (type) {
        case TYPE_NONE: return "none";
        case TYPE_RANDOM: return "random";
        case TYPE_HAMMERSLEY: return "hammersley";
        case TYPE_HALTON: return "halton";
        case TYPE_CONST: return "const";
        case TYPE_ALTERNATE: return "alternate";
        case TYPE_GRID: return "grid";
      }
      return "unknown";
    }
};

class generator_shape
{
  public:
    typedef enum {
      SHAPE_NONE,
      SHAPE_BOX,
      SHAPE_BALL,
      SHAPE_SPHERE,
      SHAPE_PLUMMER_BALL,
      SHAPE_PLUMMER
    } type_t;

    static const char *tostr(fcs_int shape) {
      switch (shape) {
        case SHAPE_NONE: return "none";
        case SHAPE_BOX: return "box";
        case SHAPE_BALL: return "ball";
        case SHAPE_SPHERE: return "sphere";
        case SHAPE_PLUMMER_BALL: return "plummer_ball";
        case SHAPE_PLUMMER: return "plummer";
      }
      return "unknown";
    }
};


struct simple_generator_params {
  generator_type::type_t type;
  fcs_int nvalues;
  fcs_float values[6];

  simple_generator_params():type(generator_type::TYPE_NONE), nvalues(0) { for (fcs_int i = 0; i < (fcs_int) (sizeof(values) / sizeof(values[0])); ++i) values[i] = 0; }
};


class Generator: public ParticleSource
{
public:
  Generator();

  virtual bool read_config(xml_node<> *config_node, const char *basename, int comm_size, int comm_rank, MPI_Comm comm);
  virtual void broadcast_config(int root, int comm_size, int comm_rank, MPI_Comm comm);
  virtual void print_config(const char *prefix, int comm_size, int comm_rank, MPI_Comm comm);
  virtual bool write_config(xml_document<> *doc, xml_node<> *config_node, int comm_size, int comm_rank, MPI_Comm comm);

private:
  bool get_numbers(fcs_int *ntotal, fcs_int *low, fcs_int *high, fcs_int *gridsize, int comm_size, int comm_rank, MPI_Comm comm);
  template<typename P>
  bool make_local_particles(generic_particle_data_t<P> *particle_data, int comm_size, int comm_rank, MPI_Comm comm);

public:
  virtual fcs_int get_local_nparticles(int comm_size, int comm_rank, MPI_Comm comm);
  virtual bool make_local_particles(particle_data_t *particle_data, int comm_size, int comm_rank, MPI_Comm comm);

#if SCAFACOS_TEST_WITH_DIPOLES
  virtual fcs_int get_dipole_local_nparticles(int comm_size, int comm_rank, MPI_Comm comm);
  virtual bool make_dipole_local_particles(dipole_particle_data_t *particle_data, int comm_size, int comm_rank, MPI_Comm comm);
#endif /* SCAFACOS_TEST_WITH_DIPOLES */

  void set_box(fcs_float *box_base, fcs_float *box_a, fcs_float *box_b, fcs_float *box_c);

private:
  struct {
    fcs_int type;
    fcs_int nlocal;
    fcs_int nntotals, ntotals[3];
    bool mult_ntotals[3];

    generator_type::type_t positions_type;
    generator_shape::type_t positions_shape;

    simple_generator_params props, potentials, field;

  } params;

  fcs_float *box_base, *box_a, *box_b, *box_c;

  void read_positions(xml_node<> *node);
  void read_simple(simple_generator_params *p, xml_node<> *node);
  void print_simple(simple_generator_params *p, const char *name, const char *prefix = "");

  template<template<int, typename> class S>
  void make_box(fcs_int n, fcs_int offset, fcs_float *v);

  template<template<int, typename> class S>
  void make_ball(fcs_int n, fcs_int offset, fcs_float *v);

  template<template<int, typename> class S>
  void make_sphere(fcs_int n, fcs_int offset, fcs_float *v);

  template<template<int, typename> class S>
  void make_plummer_ball(fcs_int n, fcs_int offset, fcs_float *v);

  template<template<int, typename> class S>
  void make_plummer(fcs_int n, fcs_int offset, fcs_float *v);

  void make_grid(fcs_int n, fcs_int offset, fcs_float *v, fcs_int *gridsize);

  template<int D, typename T>
  void make_const(fcs_int n, fcs_int offset, T *v, T *c);

  template<int D, typename T>
  void make_alternate(fcs_int n, fcs_int offset, T *v, fcs_int nc, T *c, fcs_int *gridsize = 0);

  template<int D, typename T>
  bool make_simple(simple_generator_params *p, fcs_int n, fcs_int offset, T *v, fcs_int *gridsize = 0);
};


class PlainParticles: public ParticleSource
{
public:
  PlainParticles();
  ~PlainParticles();

  virtual bool read_config(xml_node<> *config_node, const char *basename, int comm_size, int comm_rank, MPI_Comm comm);
  virtual void broadcast_config(int root, int comm_size, int comm_rank, MPI_Comm comm);
  virtual void print_config(const char *prefix, int comm_size, int comm_rank, MPI_Comm comm);
  virtual bool write_config(xml_document<> *doc, xml_node<> *config_node, int comm_size, int comm_rank, MPI_Comm comm);

  static void write_config(xml_document<> *doc, xml_node<> *config_node, const char *filename, particles_t *particles, int comm_size, int comm_rank, MPI_Comm comm);

  virtual fcs_int get_local_nparticles(int comm_size, int comm_rank, MPI_Comm comm);
  virtual bool make_local_particles(particle_data_t *particle_data, int comm_size, int comm_rank, MPI_Comm comm);

  particle_data_t &get_local_particles() { return local_particles; };

#if SCAFACOS_TEST_WITH_DIPOLES
  virtual fcs_int get_dipole_local_nparticles(int comm_size, int comm_rank, MPI_Comm comm);
  virtual bool make_dipole_local_particles(dipole_particle_data_t *particle_data, int comm_size, int comm_rank, MPI_Comm comm);

  dipole_particle_data_t &get_dipole_local_particles() { return dipole_local_particles; };
#endif /* SCAFACOS_TEST_WITH_DIPOLES */

private:
  struct {
    fcs_int total_nparticles;
#if SCAFACOS_TEST_WITH_DIPOLES
    fcs_int dipole_total_nparticles;
#endif /* SCAFACOS_TEST_WITH_DIPOLES */

  } params;

  particle_data_t local_particles;
#if SCAFACOS_TEST_WITH_DIPOLES
  dipole_particle_data_t dipole_local_particles;
#endif /* SCAFACOS_TEST_WITH_DIPOLES */
};

typedef long long sparse_int_t;

class FormatBinary
{
  public:
    static long long format_id, int_size, float_size;

    static char *read_full(char *buf, fcs_float *data, int s);
    static char *write_full(char *buf, fcs_float *data, int s);

    static char *read_sparse(char *buf, sparse_int_t *pid, fcs_float *data, fcs_int s);
    static char *write_sparse(char *buf, sparse_int_t pid, fcs_float *data, fcs_int s);
};

class FormatPortable
{
  public:
    static long long format_id, int_size, float_size;

    static char *read_full(char *buf, fcs_float *data, int s);
    static char *write_full(char *buf, fcs_float *data, int s);

    static char *read_sparse(char *buf, sparse_int_t *pid, fcs_float *data, fcs_int s);
    static char *write_sparse(char *buf, sparse_int_t pid, fcs_float *data, fcs_int s);
};

class FileParticles: public ParticleSource
{
public:
  FileParticles();

  virtual bool read_config(xml_node<> *config_node, const char *basename, int comm_size, int comm_rank, MPI_Comm comm);
  virtual void broadcast_config(int root, int comm_size, int comm_rank, MPI_Comm comm);
  virtual void print_config(const char *prefix, int comm_size, int comm_rank, MPI_Comm comm);
  virtual bool write_config(xml_document<> *doc, xml_node<> *config_node, int comm_size, int comm_rank, MPI_Comm comm);

private:
  template<class F>
  static void write_config(xml_document<> *doc, xml_node<> *parent_node, const char *node_name, const char *filename, particles_t *particles, int comm_size, int comm_rank, MPI_Comm comm);

public:
  static void write_config_binary(xml_document<> *doc, xml_node<> *parent_node, const char *node_name, const char *filename, particles_t *particles, int comm_size, int comm_rank, MPI_Comm comm);
  static void write_config_portable(xml_document<> *doc, xml_node<> *parent_node, const char *node_name, const char *filename, particles_t *particles, int comm_size, int comm_rank, MPI_Comm comm);

private:
  template<typename P>
  bool make_local_particles(generic_particle_data_t<P> *particle_data, fcs_int nlocal, int comm_size, int comm_rank, MPI_Comm comm);

public:
  virtual fcs_int get_local_nparticles(int comm_size, int comm_rank, MPI_Comm comm);
  virtual bool make_local_particles(particle_data_t *particle_data, int comm_size, int comm_rank, MPI_Comm comm);
  bool make_local_particles(particle_data_t *particle_data, fcs_int nlocal, int comm_size, int comm_rank, MPI_Comm comm);

#if SCAFACOS_TEST_WITH_DIPOLES
  virtual fcs_int get_dipole_local_nparticles(int comm_size, int comm_rank, MPI_Comm comm);
  virtual bool make_dipole_local_particles(dipole_particle_data_t *particle_data, int comm_size, int comm_rank, MPI_Comm comm);
  bool make_dipole_local_particles(dipole_particle_data_t *particle_data, fcs_int nlocal, int comm_size, int comm_rank, MPI_Comm comm);
#endif /* SCAFACOS_TEST_WITH_DIPOLES */

private:
  struct {
    fcs_int type, format;
    long long offset;
    char filename[MAX_FILENAME_LENGTH];

    fcs_int nsparse_potentials, nsparse_field;

  } params;
};


class Duplicate
{
public:
  Duplicate();

  void read_config(xml_node<> *config_node, const char *basename);
  void broadcast_config(int root, MPI_Comm comm);
  void print_config(const char *prefix = "");
  xml_node<> *write_config(xml_document<> *doc, xml_node<> *config_node);

  struct {
    fcs_int times[3];
    fcs_int rescale;

  } params;
};


#endif
