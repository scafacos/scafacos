/*
  Copyright (C) 2011,2012 Olaf Lenz, Michael Hofmann
  
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

#define POINTSEQUENCEINT_T  fcs_int
#include "PointSequence.hpp"

using namespace std;
using namespace rapidxml;


static const char GENERATE_TAG[] = "generate";
static const char BINARY_TAG[] = "binary";
static const char PORTABLE_TAG[] = "portable";
static const char REFERENCES_TAG[] = "references";


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
      SHAPE_PLUMMER,
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


typedef struct _simple_generator_params {
  generator_type::type_t type;
  fcs_int nvalues;
  fcs_float values[6];

  _simple_generator_params():type(generator_type::TYPE_NONE), nvalues(0) { for (fcs_int i = 0; i < (fcs_int) (sizeof(values) / sizeof(values[0])); ++i) values[i] = 0; }

} simple_generator_params;


class Generator {
public:
  Generator();

  void read_config(xml_node<> *config_node);
  void broadcast_config(int root, MPI_Comm comm);
  void print_config(const char *prefix = "");
//  static void write_config(xml_document<> *doc, xml_node<> *config_node);

  fcs_int get_local_nparticles(bool all_on_master, int comm_size, int comm_rank, MPI_Comm comm);
  fcs_int get_local_particles(fcs_float *positions, fcs_float *charges, fcs_float *potentials, fcs_float *field, bool all_on_master, int comm_size, int comm_rank, MPI_Comm comm);

  bool have_positions() { return params.have_positions; }
  bool have_charges() { return params.have_charges; }
  bool have_potentials() { return params.have_potentials; }
  bool have_field() { return params.have_field; }

  void set_box(fcs_float *box_base, fcs_float *box_a, fcs_float *box_b, fcs_float *box_c);

private:
  struct {
    fcs_int nlocal;
    fcs_int nntotals, ntotals[3];
    bool mult_ntotals[3];

    generator_type::type_t positions_type;
    generator_shape::type_t positions_shape;

    simple_generator_params charges, potentials, field;

    bool have_positions, have_charges, have_potentials, have_field;

  } params;

  fcs_float *box_base, *box_a, *box_b, *box_c;

  void read_positions(xml_node<> *config_node);
  void read_simple(simple_generator_params *p, xml_node<> *config_node);
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


class PlainParticles
{
public:
  fcs_int nparticles, allocated_nparticles;
  fcs_float *positions, *charges, *potentials, *field;

  PlainParticles();

  void read_config(xml_node<> *node, const char *basename);
  void broadcast_config(int root, MPI_Comm comm, bool particles);
  void print_config(const char *prefix = "");

  static void write_config(xml_document<> *doc, xml_node<> *config_node, const char *filename, fcs_int ntotal, fcs_int nparticles, fcs_float *positions, fcs_float *charges, fcs_float *potentials, fcs_float *field, int comm_size, int comm_rank, MPI_Comm comm);

  fcs_int get_local_nparticles(bool all_on_master, int comm_size, int comm_rank, MPI_Comm comm);

  bool have_positions() { return params.have_positions; }
  bool have_charges() { return params.have_charges; }
  bool have_potentials() { return params.have_potentials; }
  bool have_field() { return params.have_field; }
  fcs_int get_total_nparticles() { return params.total_nparticles; }

private:
  struct {
    fcs_int total_nparticles;
    bool have_positions, have_charges, have_potentials, have_field;

  } params;

  fcs_int add_particles(fcs_int add_nparticles);
  void free_input_particles();
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

class FileParticles
{
public:
  FileParticles();

  void read_config(xml_node<> *config_node, const char *basename);
  void broadcast_config(int root, MPI_Comm comm);
  void print_config(const char *prefix = "");

  template<class F>
  static void write_config(xml_document<> *doc, xml_node<> *parent_node, const char *node_name, const char *filename, fcs_int ntotal, fcs_int nparticles, fcs_float *positions, fcs_float *charges, fcs_float *potentials, fcs_float *field, int comm_size, int comm_rank, MPI_Comm comm);

  fcs_int get_local_nparticles(bool all_on_master, int comm_size, int comm_rank, MPI_Comm comm);
  fcs_int get_local_particles(fcs_float *positions, fcs_float *charges, fcs_float *potentials, fcs_float *field, bool all_on_master, int comm_size, int comm_rank, MPI_Comm comm);
  fcs_int get_local_particles(fcs_float *positions, fcs_float *charges, fcs_float *potentials, fcs_float *field, fcs_int nlocal, int comm_size, int comm_rank, MPI_Comm comm);

  bool have_positions() { return params.have_positions; }
  bool have_charges() { return params.have_charges; }
  bool have_potentials() { return params.have_potentials; }
  bool have_field() { return params.have_field; }
  fcs_int get_total_nparticles() { return params.total_nparticles; }

private:
  struct {
    fcs_int format;
    long long offset;
    char filename[MAX_FILENAME_LENGTH];

    fcs_int total_nparticles;

    bool have_positions, have_charges, have_potentials, have_field;

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
  xml_node<> *write_config(xml_document<> *doc, xml_node<> *node);

  struct {
    fcs_int times[3];
    fcs_int rescale;

  } params;
};


#endif
