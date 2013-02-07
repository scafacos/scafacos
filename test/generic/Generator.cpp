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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <sstream>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <cfloat>
#include <mpi.h>

#include "Generator.hpp"


using namespace std;


Generator::Generator()
{
  params.nlocal = -1;
  params.nntotals = 0;
  params.ntotals[0] = params.ntotals[1] = params.ntotals[2] = -1;
  params.mult_ntotals[0] = params.mult_ntotals[1] = params.mult_ntotals[2] = false;
  
  
  params.positions_type = generator_type::TYPE_NONE;
  params.positions_shape = generator_shape::SHAPE_NONE;

  params.have_positions = params.have_charges = params.have_potentials = params.have_field = false;
}


void Generator::read_config(xml_node<> *config_node)
{
  for (xml_attribute<> *attr = config_node->first_attribute(); attr; attr = attr->next_attribute()) {
    string aname = attr->name();
    if (aname == "nlocal") {
      parse_value(attr->value(), params.nlocal);
      params.nntotals = 0;
    } else if (aname == "ntotal") {
      char cv[] = { '_', '_', '_' };
      params.nntotals = parse_sequence(attr->value(), 3, params.ntotals, cv);
      for (fcs_int i = 0; i < params.nntotals; ++i) params.mult_ntotals[i] = (cv[i] == 'p');
      params.nlocal = -1;
    }
  }

  for (xml_node<> *node = config_node->first_node(); node; node = node->next_sibling()) {
    string nname = node->name();
    if (nname == "positions") read_positions(node);
    else if (nname == "charges") read_simple(&params.charges, node);
    else if (nname == "potentials") read_simple(&params.potentials, node);
    else if (nname == "field") read_simple(&params.field, node);
  }
}


void Generator::broadcast_config(int root, MPI_Comm comm)
{
  MPI_Bcast(&params, sizeof(params), MPI_BYTE, root, comm);
}


void Generator::print_config(const char *prefix)
{
  cout << prefix << "nlocal = " << params.nlocal << ", ntotal =";
  for (fcs_int i = 0; i < params.nntotals; ++i) cout << " " << params.ntotals[i] << (params.mult_ntotals[i]?"p":"");
  if (params.nntotals < 1) cout << " (none)";
  cout << endl;

  cout << prefix << "positions: type: " << generator_type::tostr(params.positions_type) << ", shape: " << generator_shape::tostr(params.positions_shape) << endl;

  print_simple(&params.charges, "charges", prefix);
  print_simple(&params.potentials, "potentials", prefix);
  print_simple(&params.field, "field", prefix);
}


void Generator::read_positions(xml_node<> *config_node)
{
  for (xml_attribute<> *attr = config_node->first_attribute(); attr; attr = attr->next_attribute()) {
    string aname = attr->name();
    string aval = attr->value();
    if (aname == "type") {
      if (aval == "random") params.positions_type = generator_type::TYPE_RANDOM;
      else if (aval == "hammersley") params.positions_type = generator_type::TYPE_HAMMERSLEY;
      else if (aval == "halton") params.positions_type = generator_type::TYPE_HALTON;
      else if (aval == "grid") params.positions_type = generator_type::TYPE_GRID;
      else params.positions_type = generator_type::TYPE_NONE;
    } else if (aname == "shape") {
      if (aval == "box") params.positions_shape = generator_shape::SHAPE_BOX;
      else if (aval == "ball") params.positions_shape = generator_shape::SHAPE_BALL;
      else if (aval == "sphere") params.positions_shape = generator_shape::SHAPE_SPHERE;
      else if (aval == "plummer_ball") params.positions_shape = generator_shape::SHAPE_PLUMMER_BALL;
      else if (aval == "plummer") params.positions_shape = generator_shape::SHAPE_PLUMMER;
      else params.positions_shape = generator_shape::SHAPE_NONE;
    }
  }
}


fcs_int Generator::get_local_nparticles(bool all_on_master, int comm_size, int comm_rank, MPI_Comm comm)
{
  return get_local_particles(NULL, NULL, NULL, NULL, all_on_master, comm_size, comm_rank, comm);
}

fcs_int Generator::get_local_particles(fcs_float *positions, fcs_float *charges, fcs_float *potentials, fcs_float *field, bool all_on_master, int comm_size, int comm_rank, MPI_Comm comm)
{
  fcs_int low = 0, high = 0;
  fcs_int gridsize[3] = { 1, 1, 1 };

//  print_config();

  if (params.nlocal >= 0)
  {
    if (all_on_master)
    {
      if (comm_rank == MASTER_RANK) { low = 0; high = params.nlocal * comm_size; }

    } else {
      low = comm_rank * params.nlocal;
      high = (comm_rank + 1) * params.nlocal;
    }

    /* FIXME: make a better grid when nlocal is given (currently: 2d grid: nlocal x nproc) */
    gridsize[0] = params.nlocal;
    gridsize[1] = comm_size;

  } else if (params.nntotals >= 0)
  {
    fcs_int ntotals = 1;

    for (fcs_int i = 0; i < params.nntotals; ++i) ntotals *= (gridsize[i] = params.ntotals[i] * (params.mult_ntotals[i]?comm_size:1));
    
    if (all_on_master)
    {
      if (comm_rank == MASTER_RANK) { low = 0; high = ntotals; }

    } else {
      low = (fcs_int) (((fcs_float) ntotals * (fcs_float) comm_rank) / (fcs_float) comm_size);
      high = (fcs_int) (((fcs_float) ntotals * (fcs_float) (comm_rank + 1)) / (fcs_float) comm_size);
    }

  } else cout << "ERROR: neither valid local nor total numbers of particles available" << endl;

/*  cout << "generating " << high - low << " particles (" << low << "-" << high << ") on node " << comm_rank << endl;*/

  if (positions == NULL) return high - low;

  if (params.positions_shape == generator_shape::SHAPE_BOX)
  switch (params.positions_type)
  {
    case generator_type::TYPE_RANDOM:
      make_box<RandomPointSequence>(high - low, low, positions);
      break;
    case generator_type::TYPE_HAMMERSLEY:
      make_box<HammersleyPointSequence>(high - low, low, positions);
      break;
    case generator_type::TYPE_HALTON:
      make_box<HaltonPointSequence>(high - low, low, positions);
      break;
    default:
      break;
  }

  if (params.positions_shape == generator_shape::SHAPE_BALL)
  switch (params.positions_type)
  {
    case generator_type::TYPE_RANDOM:
      make_ball<RandomPointSequence>(high - low, low, positions);
      break;
    case generator_type::TYPE_HAMMERSLEY:
      make_ball<HammersleyPointSequence>(high - low, low, positions);
      break;
    case generator_type::TYPE_HALTON:
      make_ball<HaltonPointSequence>(high - low, low, positions);
      break;
    default:
      break;
  }

  if (params.positions_shape == generator_shape::SHAPE_SPHERE)
  switch (params.positions_type)
  {
    case generator_type::TYPE_RANDOM:
      make_sphere<RandomPointSequence>(high - low, low, positions);
      break;
    case generator_type::TYPE_HAMMERSLEY:
      make_sphere<HammersleyPointSequence>(high - low, low, positions);
      break;
    case generator_type::TYPE_HALTON:
      make_sphere<HaltonPointSequence>(high - low, low, positions);
      break;
    default:
      break;
  }

  if (params.positions_shape == generator_shape::SHAPE_PLUMMER_BALL)
  switch (params.positions_type)
  {
    case generator_type::TYPE_RANDOM:
      make_plummer_ball<RandomPointSequence>(high - low, low, positions);
      break;
    case generator_type::TYPE_HAMMERSLEY:
      make_plummer_ball<HammersleyPointSequence>(high - low, low, positions);
      break;
    case generator_type::TYPE_HALTON:
      make_plummer_ball<HaltonPointSequence>(high - low, low, positions);
      break;
    default:
      break;
  }

  if (params.positions_shape == generator_shape::SHAPE_PLUMMER)
  switch (params.positions_type)
  {
    case generator_type::TYPE_RANDOM:
      make_plummer<RandomPointSequence>(high - low, low, positions);
      break;
    case generator_type::TYPE_HAMMERSLEY:
      make_plummer<HammersleyPointSequence>(high - low, low, positions);
      break;
    case generator_type::TYPE_HALTON:
      make_plummer<HaltonPointSequence>(high - low, low, positions);
      break;
    default:
      break;
  }

  if (params.positions_type == generator_type::TYPE_GRID) make_grid(high - low, low, positions, gridsize);

  make_simple<1, fcs_float>(&params.charges, high - low, low, charges, (params.positions_type == generator_type::TYPE_GRID)?gridsize:0);
  params.have_potentials = make_simple<1, fcs_float>(&params.potentials, high - low, low, potentials, (params.positions_type == generator_type::TYPE_GRID)?gridsize:0);
  params.have_field = make_simple<3, fcs_float>(&params.field, high - low, low, field, (params.positions_type == generator_type::TYPE_GRID)?gridsize:0);

  return high - low;
}


void Generator::set_box(fcs_float *box_base, fcs_float *box_a, fcs_float *box_b, fcs_float *box_c)
{
  Generator::box_base = box_base;
  Generator::box_a = box_a;
  Generator::box_b = box_b;
  Generator::box_c = box_c;
}


void Generator::read_simple(simple_generator_params *p, xml_node<> *config_node)
{
  for (xml_attribute<> *attr = config_node->first_attribute(); attr; attr = attr->next_attribute()) {
    string aname = attr->name();
    string aval = attr->value();
    if (aname == "type") {
      if (aval == "const") p->type = generator_type::TYPE_CONST;
      else if (aval == "alternate") p->type = generator_type::TYPE_ALTERNATE;
      else p->type = generator_type::TYPE_NONE;
    } else if (aname == "value") {
      p->nvalues = parse_sequence(aval, sizeof(p->values) / sizeof(p->values[0]), p->values);
    }
  }
}


void Generator::print_simple(simple_generator_params *p, const char *name, const char *prefix)
{
  cout << prefix << name << ": type: " << generator_type::tostr(p->type) << ", values =";
  for (fcs_int i = 0; i < p->nvalues; ++i) cout << " " << p->values[i];
  if (p->nvalues < 1) cout << " (none)";
  cout << endl;
}


template<template<int, typename> class S>
void Generator::make_box(fcs_int n, fcs_int offset, fcs_float *v)
{
  S<3, fcs_float> seq(n, offset);
  fcs_float rnd[3];


  for (fcs_int i = 0; i < n; ++i)
  {
    seq.get_next(rnd);

    v[3 * i + 0] = box_base[0] + box_a[0] * rnd[0] + box_b[0] * rnd[1] + box_c[0] * rnd[2];
    v[3 * i + 1] = box_base[1] + box_a[1] * rnd[0] + box_b[1] * rnd[1] + box_c[1] * rnd[2];
    v[3 * i + 2] = box_base[2] + box_a[2] * rnd[0] + box_b[2] * rnd[1] + box_c[2] * rnd[2];
  }
}


template<template <int, typename> class S>
void Generator::make_ball(fcs_int n, fcs_int offset, fcs_float *v)
{
  S<3, fcs_float> seq(n, offset);
  fcs_float rnd[3], base[3], m;
  fcs_float phi, theta, r;


  m = z_min(box_a[0], z_min(box_b[1], box_c[2])) / 2.0;
  base[0] = (box_base[0] + box_a[0] / 2.0);
  base[1] = (box_base[1] + box_b[1] / 2.0);
  base[2] = (box_base[2] + box_c[2] / 2.0);

  for (fcs_int i = 0; i < n; ++i)
  {
    seq.get_next(rnd);

    phi = acos(2.0 * rnd[0] - 1.0);
    theta = 2.0 * FCS_PI * rnd[1];
    r = pow(rnd[2], 1.0 / 3.0) * m;

//    cout << i << ": " << r[0] << ", " << r[1] << ", " << r[2] << endl;

    v[3 * i + 0] = base[0] + (r * cos (theta) * sin(phi));
    v[3 * i + 1] = base[1] + (r * sin (theta) * sin(phi));
    v[3 * i + 2] = base[2] + (r * cos (phi));
  }
}


template<template<int, typename> class S>
void Generator::make_sphere(fcs_int n, fcs_int offset, fcs_float *v)
{
  S<2, fcs_float> seq(n, offset);
  fcs_float rnd[2], base[3], m;
  fcs_float phi, theta, r;


  m = z_min(box_a[0], z_min(box_b[1], box_c[2])) / 2.0;
  base[0] = (box_base[0] + box_a[0] / 2.0);
  base[1] = (box_base[1] + box_b[1] / 2.0);
  base[2] = (box_base[2] + box_c[2] / 2.0);

  for (fcs_int i = 0; i < n; ++i)
  {
    seq.get_next(rnd);

    phi = acos(2.0 * rnd[0] - 1.0);
    theta = 2.0 * FCS_PI * rnd[1];
    r = m;

//    cout << i << ": " << r[0] << ", " << r[1] << ", " << r[2] << endl;

    v[3 * i + 0] = base[0] + (r * cos (theta) * sin(phi));
    v[3 * i + 1] = base[1] + (r * sin (theta) * sin(phi));
    v[3 * i + 2] = base[2] + (r * cos (phi));
  }
}


#define PLUM_NORM 0.033510321638291128

static fcs_float int_plummer(fcs_float r) {
  return 4 * FCS_PI * pow(r,3) / (3 * pow(1.0 + 25 * r * r, 1.5)) / PLUM_NORM;
}

template<template<int, typename> class S>
void Generator::make_plummer_ball(fcs_int n, fcs_int offset, fcs_float *v)
{
  fcs_float step = 1.0 / n;
  fcs_float b_left, b_right, integral_left, integral_right, integral_val, b_m;
  fcs_float radius, phi, theta;
  fcs_float r_max = 1.0e8;
  fcs_float b_eps = 1.0 / n * 1.0e-3;

  S<2, fcs_float> seq(n, offset);
  fcs_float rnd[2], base[3], m;


  m = z_min(box_a[0], z_min(box_b[1], box_c[2])) / 2.0;
  base[0] = (box_base[0] + box_a[0] / 2.0);
  base[1] = (box_base[1] + box_b[1] / 2.0);
  base[2] = (box_base[2] + box_c[2] / 2.0);

  if (b_eps < DBL_EPSILON) b_eps = DBL_EPSILON;

  b_left = 0.0;
  b_right = r_max;
  integral_right = 1.0 - step;
  while (b_right - b_left > DBL_EPSILON)
  {
    integral_val = int_plummer(b_right);
    if (integral_val < integral_right) break;
    b_right = 0.5 * (b_right - b_left);
  }
  /* Get upper boundary of search radius */
  r_max = b_right;

  integral_val = 0.0;

  for (fcs_int i = 0; i < n; ++i)
  {
    integral_left = step * i;
    integral_right = step * (i + 1);

    b_left = 0.0;
    b_right = r_max;

    /* Get radius by binary search */    
    while (b_right - b_left > b_eps)
    {
      b_m = 0.5 * (b_left + b_right);
      integral_val = int_plummer(b_m);
      if (integral_val < integral_left) b_left = b_m;
      else if (integral_val > integral_right) b_right = b_m;
      else
      {
        b_left = b_m - 0.5 * (b_m - b_left);
        b_right = b_m + 0.5 * (b_m - b_left);
      }
    }
    radius = b_m;

    seq.get_next(rnd);

    /* Random phi in [0, 2*PI] */
    phi = 2 * FCS_PI * rnd[0];

    /* Random theta in [0, PI] */
    theta = acos(1.0 - 2 * rnd[1]);

    v[3 * i + 0] = base[0] + m * radius * sin(theta) * cos(phi);
    v[3 * i + 1] = base[1] + m * radius * sin(theta) * sin(phi);
    v[3 * i + 2] = base[2] + m * radius * cos(theta);
  }
}


template<template<int, typename> class S>
void Generator::make_plummer(fcs_int n, fcs_int offset, fcs_float *v)
{
  S<3, fcs_float> seq(n, offset);
  fcs_float rnd[3], base[3], m;
  fcs_float r;


  m = z_min(box_a[0], z_min(box_b[1], box_c[2])) / 2.0;
  base[0] = (box_base[0] + box_a[0] / 2.0);
  base[1] = (box_base[1] + box_b[1] / 2.0);
  base[2] = (box_base[2] + box_c[2] / 2.0);

  for (fcs_int i = 0; i < n; ++i)
  {
    do
    {
      seq.get_next(rnd);

      r = 1.0 / sqrt(pow(rnd[0], -2.0 / 3.0) - 1.0);

    } while (r <= 0.1 || r >= 3.0);

    v[3 * i + 2] = (2.0 * rnd[1] - 1.0) * r;
    v[3 * i + 0] = sqrt(r * r - v[3 * i + 2] * v[3 * i + 2]) * cos(2.0 * FCS_PI * rnd[2]);
    v[3 * i + 1] = sqrt(r * r - v[3 * i + 2] * v[3 * i + 2]) * sin(2.0 * FCS_PI * rnd[2]);

    v[3 * i + 0] = (v[3 * i + 0] / 3.0 * m) + base[0];
    v[3 * i + 1] = (v[3 * i + 1] / 3.0 * m) + base[1];
    v[3 * i + 2] = (v[3 * i + 2] / 3.0 * m) + base[2];
  }
}


void Generator::make_grid(fcs_int n, fcs_int offset, fcs_float *v, fcs_int *gridsize)
{
  fcs_float base[3];

  base[0] = box_base[0] + 0.5 * (box_a[0] / gridsize[0] + box_b[0] / gridsize[1] + box_c[0] / gridsize[2]);
  base[1] = box_base[1] + 0.5 * (box_a[1] / gridsize[0] + box_b[1] / gridsize[1] + box_c[1] / gridsize[2]);
  base[2] = box_base[2] + 0.5 * (box_a[2] / gridsize[0] + box_b[2] / gridsize[1] + box_c[2] / gridsize[2]);

  for (fcs_int i = 0; i < n; ++i)
  {
    fcs_int k, g[3];

    k = i + offset;
    g[0] = k % gridsize[0]; k /= gridsize[0];
    g[1] = k % gridsize[1]; k /= gridsize[1];
    g[2] = k;

    v[3 * i + 0] = base[0] + (g[0] * box_a[0] / gridsize[0]) + (g[1] * box_b[0] / gridsize[1]) + (g[2] * box_c[0] / gridsize[2]);
    v[3 * i + 1] = base[1] + (g[0] * box_a[1] / gridsize[0]) + (g[1] * box_b[1] / gridsize[1]) + (g[2] * box_c[1] / gridsize[2]);
    v[3 * i + 2] = base[2] + (g[0] * box_a[2] / gridsize[0]) + (g[1] * box_b[2] / gridsize[1]) + (g[2] * box_c[2] / gridsize[2]);
  }
}


template<int D, typename T>
void Generator::make_const(fcs_int n, fcs_int offset, T *v, T *c)
{
  for (fcs_int i = 0; i < n; ++i)
  for (fcs_int j = 0; j < D; ++j) v[D * i + j] = c[j];
}


template<int D, typename T>
void Generator::make_alternate(fcs_int n, fcs_int offset, T *v, fcs_int nc, T *c, fcs_int *gridsize)
{
  if (gridsize)
  {
    for (fcs_int i = 0; i < n; ++i)
    for (fcs_int j = 0; j < D; ++j)
    {
      fcs_int k, g[3];

      k = i + offset;
      g[0] = k % gridsize[0]; k /= gridsize[0];
      g[1] = k % gridsize[1]; k /= gridsize[1];
      g[2] = k;

      v[D * i + j] = c[D * ((g[0] + g[1] + g[2]) % nc) + j];
    }

  } else {
    for (fcs_int i = 0; i < n; ++i)
    for (fcs_int j = 0; j < D; ++j) v[D * i + j] = c[D * ((offset + i) % nc) + j];
  }
}


template<int D, typename T>
bool Generator::make_simple(simple_generator_params *p, fcs_int n, fcs_int offset, T *v, fcs_int *gridsize)
{
  bool r = false;

  switch (p->type)
  {
    case generator_type::TYPE_CONST:
      make_const<D, T>(n, offset, v, p->values);
      r = true;
      break;
    case generator_type::TYPE_ALTERNATE:
      make_alternate<D, T>(n, offset, v, p->nvalues, p->values, gridsize);
      r = true;
      break;
    default:
      break;
  }

  return r;
}


PlainParticles::PlainParticles()
{
  nparticles = allocated_nparticles = 0;
  positions = 0;
  charges = 0;
  field = 0;
  potentials = 0;

  params.total_nparticles = 0;
  params.have_potentials = params.have_field = false;
}

static const char PARTICLE_TAG[] = "particle";

void PlainParticles::read_config(xml_node<> *node, const char *basename)
{
  xml_attribute<> *attr;
  string aname;

  for (xml_node<> *particle_node = node->first_node(PARTICLE_TAG); particle_node; particle_node = particle_node->next_sibling(PARTICLE_TAG)) {
    fcs_int pid = add_particles(1);

    for (attr = particle_node->first_attribute(); attr; attr = attr->next_attribute()) {
      aname = attr->name();
      if (aname == "position") parse_sequence(attr->value(), 3, &positions[3 * pid]);
      else if (aname == "q") parse_value(attr->value(), charges[pid]);
      else if (aname == "potential") {
        parse_value(attr->value(), potentials[pid]);
        params.have_potentials = true;
      } else if (aname == "field") {
        parse_sequence(attr->value(), 3, &field[3 * pid]);
        params.have_field = true;
      }
    }
  }
  
  params.total_nparticles = nparticles;
}

void PlainParticles::broadcast_config(int root, MPI_Comm comm, bool particles)
{
  MPI_Bcast(&params, sizeof(params), MPI_BYTE, root, comm);

  if (particles)
  {
    if (comm_rank != root)
    {
      nparticles = 0;
      add_particles(params.total_nparticles);
    }
    MPI_Bcast(positions, 3 * nparticles, FCS_MPI_FLOAT, root, comm);
    MPI_Bcast(charges, nparticles, FCS_MPI_FLOAT, root, comm);
    if (params.have_potentials) MPI_Bcast(potentials, nparticles, FCS_MPI_FLOAT, root, comm);
    if (params.have_field) MPI_Bcast(field, 3 * nparticles, FCS_MPI_FLOAT, root, comm);
  }
}

void PlainParticles::print_config(const char *prefix)
{
}

static void write_plain_full(xml_document<> *doc, xml_node<> *config_node, fcs_int nparticles, fcs_float *positions, fcs_float *charges, fcs_float *potentials, fcs_float *field, int comm_size, int comm_rank, MPI_Comm comm)
{
  fcs_int max_nparticles, my_nparticles;
  
  if (doc && config_node) my_nparticles = 0;
  else my_nparticles = nparticles;

  MPI_Reduce(&my_nparticles, &max_nparticles, 1, FCS_MPI_INT, MPI_MAX, MASTER_RANK, comm);

  fcs_int all_nparticles[comm_size];
  MPI_Gather(&nparticles, 1, FCS_MPI_INT, all_nparticles, 1, FCS_MPI_INT, MASTER_RANK, comm);

  if (doc && config_node)
  {
    xml_attribute<> *attr;
    ostringstream os;
    char* s;

    os.precision(16);

    fcs_float *tmp_positions = new fcs_float[3*max_nparticles];
    fcs_float *tmp_charges = new fcs_float[max_nparticles];
    fcs_float *tmp_potentials = new fcs_float[max_nparticles];
    fcs_float *tmp_field = new fcs_float[3*max_nparticles];

    fcs_float *pos, *c, *pot, *f;

    for (int i = 0; i < comm_size; ++i)
    if (all_nparticles[i] > 0)
    {
      if (i == MASTER_RANK)
      {
        pos = positions;
        c = charges;
        pot = potentials;
        f = field;

      } else {
        MPI_Status status;
        MPI_Recv(tmp_positions, 3 * all_nparticles[i], FCS_MPI_FLOAT, i, 0, comm, &status);
        MPI_Recv(tmp_charges, all_nparticles[i], FCS_MPI_FLOAT, i, 0, comm, &status);
        if (potentials) MPI_Recv(tmp_potentials, all_nparticles[i], FCS_MPI_FLOAT, i, 0, comm, &status);
        if (field) MPI_Recv(tmp_field, 3 * all_nparticles[i], FCS_MPI_FLOAT, i, 0, comm, &status);

        pos = tmp_positions;
        c = tmp_charges;
        pot = tmp_potentials;
        f = tmp_field;
      }

      for (fcs_int pid=0; pid < all_nparticles[i]; pid++)
      {
        // Create particle node
        xml_node<> *particle_node = doc->allocate_node(node_element, PARTICLE_TAG);
        config_node->append_node(particle_node);

        os.str("");
        os << pos[3*pid] << " " << pos[3*pid+1] << " " << pos[3*pid+2];
        s = doc->allocate_string(os.str().c_str());
        attr = doc->allocate_attribute("position", s);
        particle_node->append_attribute(attr);

        os.str("");
        os << c[pid];
        s = doc->allocate_string(os.str().c_str());
        attr = doc->allocate_attribute("q", s);
        particle_node->append_attribute(attr);

        if (potentials && pot[pid] == pot[pid])
        {
          os.str("");
          os << pot[pid];
          s = doc->allocate_string(os.str().c_str());
          attr = doc->allocate_attribute("potential", s);
          particle_node->append_attribute(attr);
        }

        if (field && f[3*pid] == f[3*pid] && f[3*pid+1] == f[3*pid+1] && f[3*pid+2] == f[3*pid+2])
        {
          os.str("");
          os << f[3*pid] << " " << f[3*pid+1] << " " << f[3*pid+2];
          s = doc->allocate_string(os.str().c_str());
          attr = doc->allocate_attribute("field", s);
          particle_node->append_attribute(attr);
        }
      }
    }
  
    delete[] tmp_positions;
    delete[] tmp_charges;
    delete[] tmp_potentials;
    delete[] tmp_field;

  } else
  {
    if (nparticles > 0)
    {
      MPI_Send(positions, 3 * nparticles, FCS_MPI_FLOAT, MASTER_RANK, 0, comm);
      MPI_Send(charges, nparticles, FCS_MPI_FLOAT, MASTER_RANK, 0, comm);
      if (potentials) MPI_Send(potentials, nparticles, FCS_MPI_FLOAT, MASTER_RANK, 0, comm);
      if (field) MPI_Send(field, 3 * nparticles, FCS_MPI_FLOAT, MASTER_RANK, 0, comm);
    }
  }
}

void PlainParticles::write_config(xml_document<> *doc, xml_node<> *config_node, const char *filename, fcs_int ntotal, fcs_int nparticles, fcs_float *positions, fcs_float *charges, fcs_float *potentials, fcs_float *field, int comm_size, int comm_rank, MPI_Comm comm)
{
  write_plain_full(doc, config_node, nparticles, positions, charges, potentials, field, comm_size, comm_rank, comm);
}

fcs_int PlainParticles::get_local_nparticles(bool all_on_master, int comm_size, int comm_rank, MPI_Comm comm)
{
  return nparticles;
}

fcs_int PlainParticles::add_particles(fcs_int add_nparticles)
{
  const fcs_int min_alloc_step = 1;
  fcs_int i;


  if (nparticles + add_nparticles > allocated_nparticles)
  {
    allocated_nparticles = z_max(nparticles + add_nparticles, allocated_nparticles + min_alloc_step);

    positions = (fcs_float *) realloc(positions, allocated_nparticles * 3 * sizeof(fcs_float));
    charges = (fcs_float *) realloc(charges, allocated_nparticles * sizeof(fcs_float));
    potentials = (fcs_float *) realloc(potentials, allocated_nparticles * sizeof(fcs_float));
    field = (fcs_float *) realloc(field, allocated_nparticles * 3 * sizeof(fcs_float));
  }

  for (i = nparticles; i < nparticles + add_nparticles; ++i)
  {
    positions[3 * i + 0] = 0.0;
    positions[3 * i + 1] = 0.0;
    positions[3 * i + 2] = 0.0;

    charges[i] = 0.0;

    potentials[i] = NAN;

    field[3 * i + 0] = NAN;
    field[3 * i + 1] = NAN;
    field[3 * i + 2] = NAN;
  }

  return (nparticles += add_nparticles) - add_nparticles;
}


void PlainParticles::free_input_particles()
{
  if (positions) free(positions);
  if (charges) free(charges);
  if (potentials) free(potentials);
  if (field) free(field);

  nparticles = allocated_nparticles = 0;
  positions = 0;
  charges = 0;
  field = 0;
  potentials = 0;
}


long long FormatBinary::format_id = 0;

long long FormatBinary::int_size = sizeof(sparse_int_t);
long long FormatBinary::float_size = sizeof(double);

char *FormatBinary::read_full(char *buf, fcs_float *data, int s)
{
  for (int i = 0; i < s; ++i)
  {
    if (data) data[i] = *((double *) buf);
    buf += sizeof(double);
  }

  return buf;
}

char *FormatBinary::write_full(char *buf, fcs_float *data, int s)
{
  for (int i = 0; i < s; ++i)
  {
    *((double *) buf) = data[i];
    buf += sizeof(double);
  }

  return buf;
}

char *FormatBinary::read_sparse(char *buf, sparse_int_t *pid, fcs_float *data, fcs_int s)
{
  if (pid) *pid = *((sparse_int_t *) buf);
  buf += sizeof(sparse_int_t);

  buf = read_full(buf, data, s);

  return buf;
}

char *FormatBinary::write_sparse(char *buf, sparse_int_t pid, fcs_float *data, fcs_int s)
{
  *((sparse_int_t *) buf) = pid;
  buf += sizeof(sparse_int_t);

  buf = write_full(buf, data, s);
  
  return buf;
}


#define STRINGIFY(_s_)   XSTRINGIFY(_s_)
#define XSTRINGIFY(_s_)  #_s_

#define PORTABLE_DEFAULT_FORMAT  1

long long FormatPortable::format_id = PORTABLE_DEFAULT_FORMAT;

#define PORTABLE_INT_LENGTH  21
long long FormatPortable::int_size = PORTABLE_INT_LENGTH + 1;
#define PORTABLE_FLOAT_PREC    16
#define PORTABLE_FLOAT_LENGTH  24
long long FormatPortable::float_size = PORTABLE_FLOAT_LENGTH + 1;

#define PORTABLE_SEPARATOR     '\n'

typedef long long portable_int_t;
typedef long double portable_float_t;

#define PORTABLE_INT_PRNFMT     "%" STRINGIFY(PORTABLE_INT_LENGTH) "lld"
#define PORTABLE_FLOAT_PRNFMT   "%" STRINGIFY(PORTABLE_FLOAT_LENGTH) "." STRINGIFY(PORTABLE_FLOAT_PREC) "Le"

#define PORTABLE_INT_SCNFMT     "%lld"
#define PORTABLE_FLOAT_SCNFMT   "%Le"

char *FormatPortable::read_full(char *buf, fcs_float *data, int s)
{
  char *p;

  for (int i = 0; i < s; ++i)
  {
    portable_float_t x;

    p = strchr(buf, PORTABLE_SEPARATOR);
    if (!p) return NULL;

    if (data)
    {
      *p = 0; sscanf(buf, PORTABLE_FLOAT_SCNFMT, &x); *p = PORTABLE_SEPARATOR;
      data[i] = (fcs_float) x;
    }

    buf = p + 1;
  }

  return buf;
}

char *FormatPortable::write_full(char *buf, fcs_float *data, int s)
{
  int n;

  for (int i = 0; i < s; ++i)
  {
    portable_float_t x = (portable_float_t) data[i];

    n = sprintf(buf, PORTABLE_FLOAT_PRNFMT, x);
    buf += n;
    *buf = PORTABLE_SEPARATOR;
    ++buf;
  }

  return buf;
}

char *FormatPortable::read_sparse(char *buf, sparse_int_t *pid, fcs_float *data, fcs_int s)
{
  char *p;

  p = strchr(buf, PORTABLE_SEPARATOR);
  if (!p) return NULL;

  if (pid)
  {
    portable_int_t x;

    *p = 0; sscanf(buf, PORTABLE_INT_SCNFMT, &x); *p = PORTABLE_SEPARATOR;

    *pid = (sparse_int_t) x;
  }

  buf = p + 1;

  buf = read_full(buf, data, s);

  return buf;
}

char *FormatPortable::write_sparse(char *buf, sparse_int_t pid, fcs_float *data, fcs_int s)
{
  int n;
  portable_int_t x = (portable_int_t) pid;

  n = sprintf(buf, PORTABLE_INT_PRNFMT, x);
  buf += n;
  *buf = PORTABLE_SEPARATOR;
  ++buf;

  buf = write_full(buf, data, s);

  return buf;
}


FileParticles::FileParticles()
{
  params.format = 0;

  params.offset = -1;
  strcpy(params.filename, "");

  params.total_nparticles = 0;

  params.have_positions = params.have_charges = params.have_potentials = params.have_field = false;
  
  params.nsparse_potentials = params.nsparse_field = -1;
}

void FileParticles::read_config(xml_node<> *config_node, const char *basename)
{
  if (config_node->name() == string("portable")) params.format = 1;

  for (xml_attribute<> *attr = config_node->first_attribute(); attr; attr = attr->next_attribute()) {
    string aname = attr->name();
    string aval = attr->value();
    if (aname == "file") {
      if (aval.c_str()[0] != '/')
        snprintf(params.filename, MAX_FILENAME_LENGTH, "%s%s", basename, aval.c_str());
      else
        strncpy(params.filename, aval.c_str(), MAX_FILENAME_LENGTH);
    } else if (aname == "format") {
      params.format = atoll(aval.c_str());
    } else if (aname == "offset") {
      params.offset = atoll(aval.c_str());
    } else if (aname == "ntotal") {
      parse_value(aval, params.total_nparticles);
    }
  }

  for (xml_node<> *node = config_node->first_node(); node; node = node->next_sibling()) {
    string nname = node->name();
    if (nname == "positions") params.have_positions = true;
    else if (nname == "charges") params.have_charges = true;
    else if (nname == "potentials")
    {
      params.have_potentials = true;
      
      xml_attribute<> *attr = node->first_attribute("nsparse");
      if (attr) params.nsparse_potentials = atol(attr->value());

    } else if (nname == "field")
    {
      params.have_field = true;

      xml_attribute<> *attr = node->first_attribute("nsparse");
      if (attr) params.nsparse_field = atol(attr->value());
    }
  }
}

void FileParticles::broadcast_config(int root, MPI_Comm comm)
{
  MPI_Bcast(&params, sizeof(params), MPI_BYTE, root, comm);
}

void FileParticles::print_config(const char *prefix)
{
  cout << prefix << "file: " << params.filename << endl;
  cout << prefix << "format: " << params.format << endl;
  cout << prefix << "offset: " << params.offset << endl;
  cout << prefix << "ntotal: " << params.total_nparticles << endl;
}

template<class F>
static long long write_data_full(MPI_File file, fcs_float *data, fcs_int s, fcs_int nlocal, int comm_size, int comm_rank, MPI_Comm comm)
{
  const long long single_size = s * F::float_size;
  const long long write_size = single_size * nlocal;

  if (file == MPI_FILE_NULL) return write_size;

  fcs_int pid_offset = 0;
  MPI_Exscan(&nlocal, &pid_offset, 1, FCS_MPI_INT, MPI_SUM, comm);

  const long long buffer_count = 10000;
  const long long buffer_size = single_size * buffer_count;

  fcs_int local_nrounds = (nlocal / buffer_count) + ((nlocal % buffer_count)?1:0);
  fcs_int global_nrounds;
  MPI_Allreduce(&local_nrounds, &global_nrounds, 1, FCS_MPI_INT, MPI_MAX, comm);

  char *buf, *p;
  
  buf = new char[buffer_size];

  fcs_int r, i, n;

  MPI_Offset displ = pid_offset * single_size;

  MPI_File_seek(file, displ, MPI_SEEK_CUR);

  i = 0;
  for (r = 0; r < global_nrounds; ++r)
  {
    n = z_min(nlocal - i, buffer_count);

    p = buf;

    for (fcs_int j = i; j < i + n; ++j) p = F::write_full(p, &data[s * j], s);
    
    MPI_Status status;
    MPI_File_write_all(file, buf, n * single_size, MPI_BYTE, &status);

    i += n;
  }

  delete[] buf;

  MPI_File_get_position(file, &displ);
  MPI_Bcast(&displ, 1, MPI_OFFSET, comm_size - 1, comm);
  MPI_File_seek(file, displ, MPI_SEEK_SET);

  return write_size;
}

template<class F>
static void read_data_full(MPI_File file, fcs_float *data, fcs_int s, fcs_int nlocal, int comm_size, int comm_rank, MPI_Comm comm)
{
  const long long single_size = s * F::float_size;

  fcs_int pid_offset = 0;
  MPI_Exscan(&nlocal, &pid_offset, 1, FCS_MPI_INT, MPI_SUM, comm);

  const long long buffer_count = 10000;
  const long long buffer_size = single_size * buffer_count;

  fcs_int local_nrounds = (nlocal / buffer_count) + ((nlocal % buffer_count)?1:0);
  fcs_int global_nrounds;
  MPI_Allreduce(&local_nrounds, &global_nrounds, 1, FCS_MPI_INT, MPI_MAX, comm);

  char *buf, *p;
  
  buf = new char[buffer_size];

  fcs_int r, i, n;
  
  MPI_Offset displ = pid_offset * single_size;

  MPI_File_seek(file, displ, MPI_SEEK_CUR);

  i = 0;
  for (r = 0; r < global_nrounds; ++r)
  {
    n = z_min(nlocal - i, buffer_count);

    MPI_Status status;
    MPI_File_read_all(file, buf, n * single_size, MPI_BYTE, &status);

    p = buf;

    for (fcs_int j = i; j < i + n; ++j) p = F::read_full(p, &data[s * j], s);
    
    i += n;
  }

  delete[] buf;

  MPI_File_get_position(file, &displ);
  MPI_Bcast(&displ, 1, MPI_OFFSET, comm_size - 1, comm);
  MPI_File_seek(file, displ, MPI_SEEK_SET);
}

template<class F>
static long long write_data_sparse(MPI_File file, fcs_float *data, fcs_int s, fcs_int nlocal, fcs_int &nsparse, int comm_size, int comm_rank, MPI_Comm comm)
{
  const long long single_size = F::int_size + s * F::float_size;

  fcs_int i, j;

  if (nsparse < 0)
  {
    nsparse = 0;
    for (i = 0; i < nlocal; ++i)
    {
      for (j = 0; j < s; ++j) if (data[s * i + j] != data[s * i + j]) break;
      if (j >= s) ++nsparse;
    }
  }

  long long write_size = single_size * nsparse;

  if (file == MPI_FILE_NULL) return write_size;

  fcs_int pid_offset = 0;
  MPI_Exscan(&nlocal, &pid_offset, 1, FCS_MPI_INT, MPI_SUM, comm);

  fcs_int sparse_offset = 0;
  MPI_Exscan(&nsparse, &sparse_offset, 1, FCS_MPI_INT, MPI_SUM, comm);

  char *buf, *p;

  buf = p = new char[write_size];

  for (i = 0; i < nlocal; ++i)
  {
    for (j = 0; j < s; ++j) if (data[s * i + j] != data[s * i + j]) break;

    if (j >= s) p = F::write_sparse(p, pid_offset + i, &data[s * i], s);
  }

  MPI_Offset offset = sparse_offset * single_size;

  MPI_File_seek(file, offset, MPI_SEEK_CUR);

  MPI_Status status;
  MPI_File_write_all(file, buf, nsparse * single_size, MPI_BYTE, &status);
  
  delete[] buf;

  MPI_File_get_position(file, &offset);
  MPI_Bcast(&offset, 1, MPI_OFFSET, comm_size - 1, comm);
  MPI_File_seek(file, offset, MPI_SEEK_SET);
  
  return write_size;
}


static int pid2rank(fcs_int pid, fcs_int noffsets, fcs_int *offsets)
{
  fcs_int l, h, m;
  
  l = 0;
  h = noffsets - 1;

  while (l <= h)
  {
    m = (l + h) / 2;
    if (pid < offsets[m]) h = m - 1;
    else l = m + 1;
  }

  return l - 1;
}


template<class F>
static void read_data_sparse(MPI_File file, fcs_float *data, fcs_int s, fcs_int nlocal, fcs_int total_nsparse, int comm_size, int comm_rank, MPI_Comm comm)
{
  const long long single_size = F::int_size + s * F::float_size;

  fcs_int offset = 0;
  MPI_Exscan(&nlocal, &offset, 1, FCS_MPI_INT, MPI_SUM, comm);

  fcs_int offsets[comm_size];
  MPI_Allgather(&offset, 1, FCS_MPI_INT, offsets, 1, FCS_MPI_INT, comm);

  fcs_int low = (fcs_int) (((fcs_float) total_nsparse * (fcs_float) comm_rank) / (fcs_float) comm_size);
  fcs_int high = (fcs_int) (((fcs_float) total_nsparse * (fcs_float) (comm_rank + 1)) / (fcs_float) comm_size);

  fcs_int local_nsparse = high - low;

  long long read_size = single_size * local_nsparse;

  char *read_buf = new char[read_size];

  MPI_Offset displ = low * single_size;

  MPI_File_seek(file, displ, MPI_SEEK_CUR);

  MPI_Status status;
  MPI_File_read_all(file, read_buf, local_nsparse * single_size, MPI_BYTE, &status);

  MPI_File_get_position(file, &displ);
  MPI_Bcast(&displ, 1, MPI_OFFSET, comm_size - 1, comm);
  MPI_File_seek(file, displ, MPI_SEEK_SET);

  int r, scounts[comm_size], sdispls[comm_size], rcounts[comm_size], rdispls[comm_size];
  fcs_int i, j;
  char *p, *p_new;
  sparse_int_t pid;

  for (i = 0; i < comm_size; ++i) scounts[i] = 0;

  p = read_buf;
  for (i = 0; i < local_nsparse; ++i)
  {
    p = F::read_sparse(p, &pid, NULL, s);

    r = pid2rank(pid, comm_size, offsets);
    ++scounts[r];
  }

  char *send_buf = new char[read_size];

  sdispls[0] = 0;
  for (i = 1; i < comm_size; ++i) sdispls[i] = sdispls[i - 1] + scounts[i - 1];
  
  p = read_buf;
  for (i = 0; i < local_nsparse; ++i)
  {
    p_new = F::read_sparse(p, &pid, NULL, s);

    r = pid2rank(pid, comm_size, offsets);
    memcpy(send_buf + sdispls[r] * single_size, p, single_size);
    p = p_new;

    ++sdispls[r];
  }

  delete[] read_buf;

  MPI_Alltoall(scounts, 1, MPI_INT, rcounts, 1, MPI_INT, comm);

  for (i = 0; i < comm_size; ++i)
  {
    scounts[i] *= single_size;
    rcounts[i] *= single_size;
  }
  
  sdispls[0] = rdispls[0] = 0;
  for (i = 1; i < comm_size; ++i)
  {
    sdispls[i] = sdispls[i - 1] + scounts[i - 1];
    rdispls[i] = rdispls[i - 1] + rcounts[i - 1];
  }

  
  fcs_int recv_size = rdispls[comm_size - 1] + rcounts[comm_size - 1];

  char *recv_buf = new char[recv_size];

  MPI_Alltoallv(send_buf, scounts, sdispls, MPI_BYTE, recv_buf, rcounts, rdispls, MPI_BYTE, comm);

  local_nsparse = recv_size / single_size;

  fcs_float d[s];

  p = recv_buf;
  for (i = 0; i < local_nsparse; ++i)
  {
    p = F::read_sparse(p, &pid, d, s);
  
    pid -= offsets[comm_rank];
    for (j = 0; j < s; ++j) data[s * pid + j] = d[j];
  }

  delete[] send_buf;
  delete[] recv_buf;
}


template<class F>
static long long write_data(MPI_File file, fcs_float *data, fcs_int s, fcs_int nlocal, fcs_int nsparse, int comm_size, int comm_rank, MPI_Comm comm)
{
  if (nsparse >= 0) return write_data_sparse<F>(file, data, s, nlocal, nsparse, comm_size, comm_rank, comm);
  else return write_data_full<F>(file, data, s, nlocal, comm_size, comm_rank, comm);
}


static void read_data(int format, MPI_File file, fcs_float *data, fcs_int s, fcs_int nlocal, fcs_int total_nsparse, int comm_size, int comm_rank, MPI_Comm comm)
{
  if (total_nsparse >= 0)
  {
    if (format == 0) read_data_sparse<FormatBinary>(file, data, s, nlocal, total_nsparse, comm_size, comm_rank, comm);
    else read_data_sparse<FormatPortable>(file, data, s, nlocal, total_nsparse, comm_size, comm_rank, comm);

  } else
  {
    if (format == 0) read_data_full<FormatBinary>(file, data, s, nlocal, comm_size, comm_rank, comm);
    else read_data_full<FormatPortable>(file, data, s, nlocal, comm_size, comm_rank, comm);
  }
}

template<class F>
void FileParticles::write_config(xml_document<> *doc, xml_node<> *parent_node, const char *node_name, const char *filename, fcs_int ntotal, fcs_int nparticles, fcs_float *positions, fcs_float *charges, fcs_float *potentials, fcs_float *field, int comm_size, int comm_rank, MPI_Comm comm)
{
  MPI_File file;
  MPI_Offset offset;

  fcs_int local_nsparse_potentials = -1, total_nsparse_potentials = -1;

#define SPARSE_THRESHOLD  0.75

  if (potentials)
  {
    long long local_sums[3], total_sums[3];

    local_sums[0] = write_data_full<F>(MPI_FILE_NULL, potentials, 1, nparticles, comm_size, comm_rank, comm);
    local_sums[1] = write_data_sparse<F>(MPI_FILE_NULL, potentials, 1, nparticles, local_nsparse_potentials, comm_size, comm_rank, comm);
    local_sums[2] = local_nsparse_potentials;

    MPI_Allreduce(&local_sums, &total_sums, 3, MPI_LONG_LONG, MPI_SUM, comm);

    if (total_sums[0] * SPARSE_THRESHOLD >= total_sums[1]) total_nsparse_potentials = total_sums[2]; /* sparse mode */
    else local_nsparse_potentials = total_nsparse_potentials = -1; /* full mode */
  }

  fcs_int local_nsparse_field = -1, total_nsparse_field = -1;

  if (field)
  {
    long long local_sums[3], total_sums[3];

    local_sums[0] = write_data_full<F>(MPI_FILE_NULL, field, 3, nparticles, comm_size, comm_rank, comm);
    local_sums[1] = write_data_sparse<F>(MPI_FILE_NULL, field, 3, nparticles, local_nsparse_field, comm_size, comm_rank, comm);
    local_sums[2] = local_nsparse_field;

    MPI_Allreduce(&local_sums, &total_sums, 3, MPI_LONG_LONG, MPI_SUM, comm);

    if (total_sums[0] * SPARSE_THRESHOLD >= total_sums[1]) total_nsparse_field = total_sums[2]; /* sparse mode */
    else local_nsparse_field = total_nsparse_field = -1; /* full mode */
  }

  MPI_File_open(comm, (char *) filename, MPI_MODE_CREATE|MPI_MODE_APPEND|MPI_MODE_WRONLY, MPI_INFO_NULL, &file);

  MPI_File_get_position(file, &offset);

  if (doc != NULL)
  {
    xml_attribute<> *attr;
    ostringstream os;
    char* s;

    xml_node<> *file_node = doc->allocate_node(node_element, node_name);
    parent_node->append_node(file_node);

    const char *filenamestr = strrchr(filename, '/');
    if (filenamestr == NULL) filenamestr = filename;
    else ++filenamestr;

    s = doc->allocate_string(filenamestr);
    attr = doc->allocate_attribute("file", s);
    file_node->append_attribute(attr);
    
    os.str("");
    os << F::format_id;
    s = doc->allocate_string(os.str().c_str());
    attr = doc->allocate_attribute("format", s);
    file_node->append_attribute(attr);
    
    os.str("");
    os << offset;
    s = doc->allocate_string(os.str().c_str());
    attr = doc->allocate_attribute("offset", s);
    file_node->append_attribute(attr);
    
    os.str("");
    os << ntotal;
    s = doc->allocate_string(os.str().c_str());
    attr = doc->allocate_attribute("ntotal", s);
    file_node->append_attribute(attr);

    if (positions) file_node->append_node(doc->allocate_node(node_element, "positions"));
    
    if (charges) file_node->append_node(doc->allocate_node(node_element, "charges"));

    if (potentials)
    {
      xml_node<> *node = doc->allocate_node(node_element, "potentials");
      file_node->append_node(node);

      if (total_nsparse_potentials >= 0)
      {
        os.str("");
        os << total_nsparse_potentials;
        s = doc->allocate_string(os.str().c_str());
        attr = doc->allocate_attribute("nsparse", s);
        node->append_attribute(attr);
      }
    }

    if (field)
    {
      xml_node<> *node = doc->allocate_node(node_element, "field");
      file_node->append_node(node);

      if (total_nsparse_field >= 0)
      {
        os.str("");
        os << total_nsparse_field;
        s = doc->allocate_string(os.str().c_str());
        attr = doc->allocate_attribute("nsparse", s);
        node->append_attribute(attr);
      }
    }
  }

  // write positions
  if (positions) write_data<F>(file, positions, 3, nparticles, -1, comm_size, comm_rank, comm);

  // write charges
  if (charges) write_data<F>(file, charges, 1, nparticles, -1, comm_size, comm_rank, comm);

  // write potentials
  if (potentials) write_data<F>(file, potentials, 1, nparticles, local_nsparse_potentials, comm_size, comm_rank, comm);

  // write field
  if (field) write_data<F>(file, field, 3, nparticles, local_nsparse_field, comm_size, comm_rank, comm);
  
  MPI_File_close(&file);
}

template void FileParticles::write_config<FormatBinary>(xml_document<> *doc, xml_node<> *parent_node, const char *node_name, const char *filename, fcs_int ntotal, fcs_int nparticles, fcs_float *positions, fcs_float *charges, fcs_float *potentials, fcs_float *field, int comm_size, int comm_rank, MPI_Comm comm);

template void FileParticles::write_config<FormatPortable>(xml_document<> *doc, xml_node<> *parent_node, const char *node_name, const char *filename, fcs_int ntotal, fcs_int nparticles, fcs_float *positions, fcs_float *charges, fcs_float *potentials, fcs_float *field, int comm_size, int comm_rank, MPI_Comm comm);

fcs_int FileParticles::get_local_nparticles(bool all_on_master, int comm_size, int comm_rank, MPI_Comm comm)
{
  return get_local_particles(NULL, NULL, NULL, NULL, all_on_master, comm_size, comm_rank, comm);
}

fcs_int FileParticles::get_local_particles(fcs_float *positions, fcs_float *charges, fcs_float *potentials, fcs_float *field, bool all_on_master, int comm_size, int comm_rank, MPI_Comm comm)
{
  fcs_int low = 0, high = 0;

  if (all_on_master)
  {
    if (comm_rank == MASTER_RANK) { low = 0; high = params.total_nparticles; }

  } else {

    low = (fcs_int) (((fcs_float) params.total_nparticles * (fcs_float) comm_rank) / (fcs_float) comm_size);
    high = (fcs_int) (((fcs_float) params.total_nparticles * (fcs_float) (comm_rank + 1)) / (fcs_float) comm_size);
  }

  return get_local_particles(positions, charges, potentials, field, high - low, comm_size, comm_rank, comm);
}

fcs_int FileParticles::get_local_particles(fcs_float *positions, fcs_float *charges, fcs_float *potentials, fcs_float *field, fcs_int nlocal, int comm_size, int comm_rank, MPI_Comm comm)
{
  if (positions != NULL)
  {
    MPI_File file;
  
    MPI_File_open(comm, params.filename, MPI_MODE_RDONLY, MPI_INFO_NULL, &file);

    if (params.offset >= 0) MPI_File_seek(file, params.offset, MPI_SEEK_SET);

    // read positions
    if (params.have_positions) read_data(params.format, file, positions, 3, nlocal, -1, comm_size, comm_rank, comm);

    // read charges
    if (params.have_charges) read_data(params.format, file, charges, 1, nlocal, -1, comm_size, comm_rank, comm);

    // read potentials
    if (params.have_potentials) read_data(params.format, file, potentials, 1, nlocal, params.nsparse_potentials, comm_size, comm_rank, comm);

    // read field
    if (params.have_field) read_data(params.format, file, field, 3, nlocal, params.nsparse_field, comm_size, comm_rank, comm);

    MPI_File_close(&file);
  }

  return nlocal;
}


Duplicate::Duplicate()
{
  params.times[0] = params.times[1] = params.times[2] = 1;
  params.rescale = 0;
}

void Duplicate::read_config(xml_node<> *node, const char *basename)
{
  xml_attribute<> *attr;
  string aname;

  for (attr = node->first_attribute(); attr; attr = attr->next_attribute()) {
    aname = attr->name();
    if (aname == "times") parse_sequence(attr->value(), 3, params.times);
    else if (aname == "rescale") parse_value(attr->value(), params.rescale);
  }
}

void Duplicate::broadcast_config(int root, MPI_Comm comm)
{
  MPI_Bcast(&params, sizeof(params), MPI_BYTE, root, comm);
}

void Duplicate::print_config(const char *prefix)
{
  cout << prefix << params.times[0] << "x" << params.times[1] << "x" << params.times[2]
       << " " << ((params.rescale)?"with":"without") << " rescaling" << endl;
}

xml_node<> *Duplicate::write_config(xml_document<> *doc, xml_node<> *node)
{
  if (params.times[0] == 1 && params.times[1] == 1 && params.times[2] == 1) return node;

  xml_node<> *new_node = NULL;

  if (doc != NULL)
  {
    xml_attribute<> *attr;
    ostringstream os;
    char* s;

    new_node = doc->allocate_node(node_element, "duplicate");
    node->append_node(new_node);

    os.str("");
    os << params.times[0] << " " << params.times[1] << " " << params.times[2];
    s = doc->allocate_string(os.str().c_str());
    attr = doc->allocate_attribute("times", s);
    new_node->append_attribute(attr);
    
    os.str("");
    os << params.rescale;
    s = doc->allocate_string(os.str().c_str());
    attr = doc->allocate_attribute("rescale", s);
    new_node->append_attribute(attr);
  }

  return new_node;
}
