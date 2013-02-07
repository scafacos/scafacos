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
#include <cerrno>
#ifdef HAVE_ZLIB_H
# include <zlib.h>
# define ZLIB_IFELSE(_if_, _else_) _if_
#else
# define ZLIB_IFELSE(_if_, _else_) _else_
#endif
#include <mpi.h>

#include "rapidxml/rapidxml.hpp"
#include "rapidxml/rapidxml_utils.hpp"
#include "rapidxml/rapidxml_print.hpp"

#include "common.hpp"
#include "Generator.hpp"
#include "Testcase.hpp"

using namespace std;
using namespace rapidxml;


static const char TESTCASE_TAG[] = "scafacos_test";
static const char METHOD_TAG[] = "method";
static const char PARAM_TAG[] = "param";
static const char CONFIGURATION_TAG[] = "configuration";
/*static const char PARTICLE_TAG[] = "particle";*/
static const char DUPLICATE_TAG[] = "duplicate";


Configuration::Configuration() {

  params.offset[0] = params.offset[1] = params.offset[2] = 0.0;
  params.box_a[0] = params.box_a[1] = params.box_a[2] = 0.0;
  params.box_b[0] = params.box_b[1] = params.box_b[2] = 0.0;
  params.box_c[0] = params.box_c[1] = params.box_c[2] = 0.0;
  params.periodicity[0] = params.periodicity[1] = params.periodicity[2] = 0;
  params.epsilon = 0.0;
  params.decomposition = DECOMPOSE_ATOMISTIC;

  params.periodic_duplications[0] = params.periodic_duplications[1] = params.periodic_duplications[2] = 1;

  total_duplications[0] = total_duplications[1] = total_duplications[2] = total_duplication = 0;
  
  unscale_box[0] = unscale_box[1] = unscale_box[2] = 1.0;

  input_plain_nparticles = 0;

  input_file_nparticles = 0;

  input_generator_nparticles = 0;
  
  input_ref_nparticles = 0;

  dup_input_total_nparticles = 0;
  dup_input_nparticles = dup_input_nparticles_allocated = 0;
  dup_input_positions = 0;
  dup_input_charges = 0;
  dup_input_field = 0;
  dup_input_potentials = 0;

  have_reference_values[0] = have_reference_values[1] = 0;
  have_result_values[0] = have_result_values[1] = 0;

  decomp_total_nparticles = 0;
  decomp_nparticles = 0;
  decomp_positions = 0;
  decomp_charges = 0;
  decomp_potentials = 0;
  decomp_field = 0;
  decomp_comm = MPI_COMM_NULL;

  result_potentials = 0;
  result_field = 0;

  cart_comm = MPI_COMM_NULL;
}

Configuration::~Configuration() {

  free_dup_input_particles();
  
  free_decomp_particles(true);
}

fcs_int Configuration::add_dup_input_particles(fcs_int add_nparticles)
{
  const fcs_int min_alloc_step = 1;
  fcs_int i;


  if (dup_input_nparticles + add_nparticles > dup_input_nparticles_allocated) {

    dup_input_nparticles_allocated = z_max(dup_input_nparticles + add_nparticles, dup_input_nparticles_allocated + min_alloc_step);

    dup_input_positions = (fcs_float*) realloc(dup_input_positions, dup_input_nparticles_allocated*3*sizeof(fcs_float));
    dup_input_charges = (fcs_float*) realloc(dup_input_charges, dup_input_nparticles_allocated*sizeof(fcs_float));
    dup_input_potentials = (fcs_float*) realloc(dup_input_potentials, dup_input_nparticles_allocated*sizeof(fcs_float));
    dup_input_field = (fcs_float*) realloc(dup_input_field, dup_input_nparticles_allocated*3*sizeof(fcs_float));
  }

  for (i = dup_input_nparticles; i < dup_input_nparticles + add_nparticles; ++i)
  {
    dup_input_positions[3 * i + 0] = 0.0;
    dup_input_positions[3 * i + 1] = 0.0;
    dup_input_positions[3 * i + 2] = 0.0;

    dup_input_charges[i] = 0.0;

    dup_input_potentials[i] = NAN;

    dup_input_field[3 * i + 0] = NAN;
    dup_input_field[3 * i + 1] = NAN;
    dup_input_field[3 * i + 2] = NAN;
  }

  return (dup_input_nparticles += add_nparticles) - add_nparticles;
}

void Configuration::free_dup_input_particles()
{
  if (dup_input_positions != input_plain.positions) free(dup_input_positions);
  if (dup_input_charges != input_plain.charges) free(dup_input_charges);
  if (dup_input_potentials != input_plain.potentials) free(dup_input_potentials);
  if (dup_input_field != input_plain.field) free(dup_input_field);

  dup_input_nparticles = 0;
  dup_input_positions = 0;
  dup_input_charges = 0;
  dup_input_field = 0;
  dup_input_potentials = 0;
}

void Configuration::read_config(xml_node<> *config_node, const char *basename)
{
  xml_attribute<> *attr;
  string aname;

  // Read configuration data
  for (attr = config_node->first_attribute(); attr; attr = attr->next_attribute()) {
    aname = attr->name();
    if (aname == "offset") parse_sequence(attr->value(), 3, params.offset);
    else if (aname == "box_a") parse_sequence(attr->value(), 3, params.box_a);
    else if (aname == "box_b") parse_sequence(attr->value(), 3, params.box_b);
    else if (aname == "box_c") parse_sequence(attr->value(), 3, params.box_c);
    else if (aname == "periodicity") parse_sequence(attr->value(), 3, params.periodicity);
    else if (aname == "epsilon") {
      if (strcmp(attr->value(), "metallic") == 0)
        params.epsilon = EPSILON_METALLIC;
      else if (strcmp(attr->value(), "vacuum") == 0)
        params.epsilon = 0.0;
      else parse_value(attr->value(), params.epsilon);
    } else if (aname == "decomposition") {
      if (strcmp(attr->value(), "all_on_master") == 0)
        params.decomposition = DECOMPOSE_ALL_ON_MASTER;
      else if (strcmp(attr->value(), "atomistic") == 0)
        params.decomposition = DECOMPOSE_ATOMISTIC;
      else if (strcmp(attr->value(), "random") == 0)
        params.decomposition = DECOMPOSE_RANDOM;
      else if (strcmp(attr->value(), "domain") == 0)
        params.decomposition = DECOMPOSE_DOMAIN;
    }
  }

  xml_node<> *particles_node = config_node;

  // Read duplication information
  xml_node<> *dup_node = config_node->first_node(DUPLICATE_TAG);
  if (dup_node)
  {
    input_duplication.read_config(dup_node, basename);
    particles_node = dup_node;
  }

  determine_total_duplication();

  // Read plain particles
  input_plain.read_config(particles_node, basename);
  input_plain_nparticles += input_plain.get_local_nparticles(true, comm_size, comm_rank, communicator);

  // Loop over binary particle input files
  for (xml_node<> *binary_node = particles_node->first_node(BINARY_TAG); binary_node; binary_node = binary_node->next_sibling(BINARY_TAG)) {
    input_files.push_back(FileParticles());
    input_files.back().read_config(binary_node, basename);
    input_file_nparticles += input_files.back().get_local_nparticles(true, comm_size, comm_rank, communicator);
  }

  // Loop over portable particle input files
  for (xml_node<> *portable_node = particles_node->first_node(PORTABLE_TAG); portable_node; portable_node = portable_node->next_sibling(PORTABLE_TAG)) {
    input_files.push_back(FileParticles());
    input_files.back().read_config(portable_node, basename);
    input_file_nparticles += input_files.back().get_local_nparticles(true, comm_size, comm_rank, communicator);
  }

  // Loop over particle generators
  for (xml_node<> *generator_node = particles_node->first_node(GENERATE_TAG); generator_node; generator_node = generator_node->next_sibling(GENERATE_TAG)) {
    input_generators.push_back(Generator());
    input_generators.back().read_config(generator_node);
    input_generator_nparticles += input_generators.back().get_local_nparticles(true, comm_size, comm_rank, communicator);
  }

  // Read extra references
  xml_node<> *ref_node = config_node->first_node(REFERENCES_TAG);
  if (ref_node)
  {
    input_ref.read_config(ref_node, basename);
    input_ref_nparticles += input_ref.get_local_nparticles(true, comm_size, comm_rank, communicator);
  }
}

void Configuration::write_config(xml_document<> *doc, xml_node<> *config_node, const char *binfilename, const char* portable_filename, bool keep_dupgen)
{
  if (doc && config_node)
  {
    xml_attribute<> *attr;
    ostringstream os;
    char* s;

    fcs_float out_box_a[3] = { params.box_a[0], params.box_a[1], params.box_a[2] };
    fcs_float out_box_b[3] = { params.box_b[0], params.box_b[1], params.box_b[2] };
    fcs_float out_box_c[3] = { params.box_c[0], params.box_c[1], params.box_c[2] };

    if (keep_dupgen)
    {
      out_box_a[0] *= unscale_box[0]; out_box_a[1] *= unscale_box[0]; out_box_a[2] *= unscale_box[0];
      out_box_b[0] *= unscale_box[1]; out_box_b[1] *= unscale_box[1]; out_box_b[2] *= unscale_box[1];
      out_box_c[0] *= unscale_box[2]; out_box_c[1] *= unscale_box[2]; out_box_c[2] *= unscale_box[2];
    }

    os.precision(16);

    os.str("");
    os << params.offset[0] << " " << params.offset[1] << " " << params.offset[2];
    s = doc->allocate_string(os.str().c_str());
    attr = doc->allocate_attribute("offset", s);
    config_node->append_attribute(attr);

    os.str("");
    os << out_box_a[0] << " " << out_box_a[1] << " " << out_box_a[2];
    s = doc->allocate_string(os.str().c_str());
    attr = doc->allocate_attribute("box_a", s);
    config_node->append_attribute(attr);

    os.str("");
    os << out_box_b[0] << " " << out_box_b[1] << " " << out_box_b[2];
    s = doc->allocate_string(os.str().c_str());
    attr = doc->allocate_attribute("box_b", s);
    config_node->append_attribute(attr);

    os.str("");
    os << out_box_c[0] << " " << out_box_c[1] << " " << out_box_c[2];
    s = doc->allocate_string(os.str().c_str());
    attr = doc->allocate_attribute("box_c", s);
    config_node->append_attribute(attr);

    os.str("");
    os << params.periodicity[0] << " " << params.periodicity[1] << " " << params.periodicity[2];
    s = doc->allocate_string(os.str().c_str());
    attr = doc->allocate_attribute("periodicity", s);
    config_node->append_attribute(attr);

    if (params.epsilon == EPSILON_METALLIC)
      attr = doc->allocate_attribute("epsilon", "metallic");
    else if (params.epsilon == 0.0)
      attr = doc->allocate_attribute("epsilon", "vacuum");
    else {
      os.str("");
      os << params.epsilon;
      s = doc->allocate_string(os.str().c_str());
      attr = doc->allocate_attribute("epsilon", s);
    }
    config_node->append_attribute(attr);
  }

  if (keep_dupgen)
  {
    xml_node<> *particles_node = input_duplication.write_config(doc, config_node);

    if (comm_rank == MASTER_RANK)
      PlainParticles::write_config(doc, particles_node, NULL, input_plain.get_total_nparticles(), input_plain.nparticles,
        input_plain.positions, input_plain.charges, NULL, NULL, comm_size, comm_rank, communicator);
    else
      PlainParticles::write_config(NULL, NULL, NULL, input_plain.get_total_nparticles(), 0,
        (fcs_float *) 1, (fcs_float *) 1, NULL, NULL, comm_size, comm_rank, communicator);

    if (binfilename)
    {
      FileParticles::write_config<FormatBinary>(doc, config_node, REFERENCES_TAG, binfilename, dup_input_total_nparticles, dup_input_nparticles,
        NULL, NULL,
        (have_reference_values[0] || have_result_values[0])?(dup_input_potentials?dup_input_potentials:(fcs_float *) 1):NULL,
        (have_reference_values[1] || have_result_values[1])?(dup_input_field?dup_input_field:(fcs_float *) 1):NULL,
        comm_size, comm_rank, communicator);

    } else if (portable_filename)
    {
      FileParticles::write_config<FormatPortable>(doc, config_node, REFERENCES_TAG, portable_filename, dup_input_total_nparticles, dup_input_nparticles,
        NULL, NULL,
        (have_reference_values[0] || have_result_values[0])?(dup_input_potentials?dup_input_potentials:(fcs_float *) 1):NULL,
        (have_reference_values[1] || have_result_values[1])?(dup_input_field?dup_input_field:(fcs_float *) 1):NULL,
        comm_size, comm_rank, communicator);

    } else {
      /* TODO: write references plain? */
      MASTER(cout << "WARNING: writing references to XML file when keeping duplication/generation information not implemented, yet!" << endl);
    }

  } else {

    if (binfilename)
      FileParticles::write_config<FormatBinary>(doc, config_node, BINARY_TAG, binfilename, dup_input_total_nparticles, dup_input_nparticles,
        dup_input_positions?dup_input_positions:(fcs_float *) 1,
        dup_input_charges?dup_input_charges:(fcs_float *) 1,
        (have_reference_values[0] || have_result_values[0])?(dup_input_potentials?dup_input_potentials:(fcs_float *) 1):NULL,
        (have_reference_values[1] || have_result_values[1])?(dup_input_field?dup_input_field:(fcs_float *) 1):NULL,
        comm_size, comm_rank, communicator);
    else if (portable_filename)
      FileParticles::write_config<FormatPortable>(doc, config_node, PORTABLE_TAG, portable_filename, dup_input_total_nparticles, dup_input_nparticles,
        dup_input_positions?dup_input_positions:(fcs_float *) 1,
        dup_input_charges?dup_input_charges:(fcs_float *) 1,
        (have_reference_values[0] || have_result_values[0])?(dup_input_potentials?dup_input_potentials:(fcs_float *) 1):NULL,
        (have_reference_values[1] || have_result_values[1])?(dup_input_field?dup_input_field:(fcs_float *) 1):NULL,
        comm_size, comm_rank, communicator);
    else
      PlainParticles::write_config(doc, config_node, NULL, dup_input_total_nparticles, dup_input_nparticles,
        dup_input_positions, dup_input_charges,
        (have_reference_values[0] || have_result_values[0])?(dup_input_potentials?dup_input_potentials:(fcs_float *) 1):NULL,
        (have_reference_values[1] || have_result_values[1])?(dup_input_field?dup_input_field:(fcs_float *) 1):NULL,
        comm_size, comm_rank, communicator);
  }
}

void Configuration::print_config(fcs_int n)
{
  cout << "  config " << n << ":" << endl;
  cout << "    system:" << endl;
  cout << "      box: offset: [" << params.offset[0] << "," << params.offset[1] << "," << params.offset[2] << "], size: "
       << "[" << params.box_a[0] << "," << params.box_a[1] << "," << params.box_a[2] << "]x"
       << "[" << params.box_b[0] << "," << params.box_b[1] << "," << params.box_b[2] << "]x"
       << "[" << params.box_c[0] << "," << params.box_c[1] << "," << params.box_c[2] << "]" << endl;
  cout << "      periodicity = (" << params.periodicity[0] << ", " << params.periodicity[1] << ", " << params.periodicity[2] << ")" << endl;

  cout << "    input particles: " << input_plain_nparticles << endl;
  cout << "    input file: " << input_file_nparticles << endl;
  fcs_int i = 0;
  for (vector<Generator>::iterator g = input_generators.begin(); g != input_generators.end(); ++g)
  {
    cout << "    particle generator " << i << ":" << endl;
    g->print_config("      ");
    ++i;
  }
  cout << "    total generated particles: " << input_generator_nparticles << endl;

  input_duplication.print_config("    general duplication: ");

  cout << "    input references: " << input_ref_nparticles << endl;
  input_ref.print_config("      ");

  cout << "    periodic duplication: " << params.periodic_duplications[0] << "x" << params.periodic_duplications[1] << "x" << params.periodic_duplications[2] << " (without rescaling)" << endl;
  cout << "    decomposition: ";
  switch (params.decomposition)
  {
    case DECOMPOSE_ALL_ON_MASTER: cout << "all-on-master" << endl; break;
    case DECOMPOSE_ATOMISTIC: cout << "atomistic" << endl; break;
    case DECOMPOSE_RANDOM: cout << "random" << endl; break;
    case DECOMPOSE_DOMAIN: cout << "domain" << endl; break;
  }

  cout << "    total particles (incl. duplications): (" << input_plain_nparticles << " + " << input_file_nparticles << " + " << input_generator_nparticles << ") * "<< total_duplication
       << " = " << (input_plain_nparticles + input_file_nparticles + input_generator_nparticles) * total_duplication << endl;
}

void Configuration::broadcast_config()
{
  INFO_MASTER(cout << "Broadcast config ..." << endl);

  MPI_Bcast(&params, sizeof(params), MPI_BYTE, MASTER_RANK, communicator);
}

void Configuration::broadcast_input()
{
  INFO_MASTER(cout << "Broadcast input ..." << endl);

  input_duplication.broadcast_config(MASTER_RANK, communicator);

  determine_total_duplication();

  // Broadcast plain particle information
  input_plain.broadcast_config(MASTER_RANK, communicator, (total_duplication > 1 && params.decomposition != DECOMPOSE_ALL_ON_MASTER));

  fcs_int n;

  // Broadcast input file information
  n = input_files.size();
  MPI_Bcast(&n, 1, FCS_MPI_INT, MASTER_RANK, communicator);
  while ((fcs_int) input_files.size() < n) input_files.push_back(FileParticles());
  for (vector<FileParticles>::iterator b = input_files.begin(); b != input_files.end(); ++b)
    b->broadcast_config(MASTER_RANK, communicator);

  // Broadcast generator information
  n = input_generators.size();
  MPI_Bcast(&n, 1, FCS_MPI_INT, MASTER_RANK, communicator);
  while ((fcs_int) input_generators.size() < n) input_generators.push_back(Generator());
  for (vector<Generator>::iterator g = input_generators.begin(); g != input_generators.end(); ++g)
    g->broadcast_config(MASTER_RANK, communicator);

  // Broadcast extra reference information
  input_ref.broadcast_config(MASTER_RANK, communicator);
}

void Configuration::generate_input_particles()
{
  broadcast_input();

  for (fcs_int i = -2; i < (fcs_int) input_generators.size(); ++i)
  {
    fcs_int my_first, my_count, my_modulo;

    if (i == -2)
    {
      if (input_plain.get_total_nparticles() <= 0) continue;

      if (total_duplication == 1 && params.decomposition != DECOMPOSE_ALL_ON_MASTER)
      {
        int scounts[comm_size], sdispls[comm_size], scounts3[comm_size], sdispls3[comm_size];

        make_equal_counts_and_displs(input_plain.get_total_nparticles(), comm_size, scounts, sdispls, scounts3, sdispls3);

        DEBUG(cout << comm_rank << ": scattering input: first: " << sdispls[comm_rank] << ", count: " << scounts[comm_rank] << endl);

        fcs_int pid = add_dup_input_particles(scounts[comm_rank]);

        MPI_Scatterv(input_plain.positions, scounts3, sdispls3, FCS_MPI_FLOAT, dup_input_positions + 3 * pid, scounts3[comm_rank], FCS_MPI_FLOAT, MASTER_RANK, communicator);
        MPI_Scatterv(input_plain.charges, scounts, sdispls, FCS_MPI_FLOAT, dup_input_charges + pid, scounts[comm_rank], FCS_MPI_FLOAT, MASTER_RANK, communicator);
        MPI_Scatterv(input_plain.potentials, scounts, sdispls, FCS_MPI_FLOAT, dup_input_potentials + pid, scounts[comm_rank], FCS_MPI_FLOAT, MASTER_RANK, communicator);
        MPI_Scatterv(input_plain.field, scounts3, sdispls3, FCS_MPI_FLOAT, dup_input_field + 3 * pid, scounts3[comm_rank], FCS_MPI_FLOAT, MASTER_RANK, communicator);

        if (input_plain.have_potentials()) have_reference_values[0] = 1;
        if (input_plain.have_field()) have_reference_values[1] = 1;

        continue;
      }

      my_modulo = input_plain.get_total_nparticles();

      fcs_int total_count = total_duplication * input_plain.get_total_nparticles();

      if (params.decomposition == DECOMPOSE_ALL_ON_MASTER)
      {
        my_first = 0;
        my_count = (comm_rank == MASTER_RANK)?total_count:0;

      } else
      {
        my_first = (fcs_int) (((fcs_float) (comm_rank + 0) * (fcs_float) total_count / (fcs_float) comm_size) + 0.5);
        my_count = (fcs_int) (((fcs_float) (comm_rank + 1) * (fcs_float) total_count / (fcs_float) comm_size) + 0.5) - my_first;
      }

      DEBUG(cout << comm_rank << ": get from input: first: " << my_first << ", count: " << my_count << ", total: " << total_count << endl);

    } else if (i == -1)
    {
      if (input_files.size() < 1) continue;

      my_modulo = input_files[0].get_local_nparticles((params.decomposition == DECOMPOSE_ALL_ON_MASTER), comm_size, comm_rank, communicator);

      my_first = 0;
      my_count = total_duplication * my_modulo;

      DEBUG(cout << comm_rank << ": get from file #" << i << ": first: " << my_first << ", count: " << my_count << endl);

    } else
    {
      input_generators[i].set_box(params.offset, params.box_a, params.box_b, params.box_c);

      my_modulo = input_generators[i].get_local_nparticles((params.decomposition == DECOMPOSE_ALL_ON_MASTER), comm_size, comm_rank, communicator);

      my_first = 0;
      my_count = total_duplication * my_modulo;

      DEBUG(cout << comm_rank << ": get from generator #" << i << ": first: " << my_first << ", count: " << my_count << endl);
    }

    fcs_int plow, dlow[4], phigh, dhigh[4], x;

    x = my_first;
    plow = x % my_modulo; x /= my_modulo;
    dlow[0] = x % total_duplications[0]; x /= total_duplications[0];
    dlow[1] = x % total_duplications[1]; x /= total_duplications[1];
    dlow[2] = x % total_duplications[2]; x /= total_duplications[2];
    dlow[3] = x;

    x = my_first + my_count;
    phigh = x % my_modulo; x /= my_modulo;
    dhigh[0] = x % total_duplications[0]; x /= total_duplications[0];
    dhigh[1] = x % total_duplications[1]; x /= total_duplications[1];
    dhigh[2] = x % total_duplications[2]; x /= total_duplications[2];
    dhigh[3] = x;

    DEBUG(cout << comm_rank << ": low: " << plow << " @ " << dlow[0] << "," << dlow[1] << "," << dlow[2] << "," << dlow[3] << endl);
    DEBUG(cout << comm_rank << ": high: " << phigh << " @ " << dhigh[0] << "," << dhigh[1] << "," << dhigh[2] << "," << dhigh[3] << endl);

    fcs_float scale[3] = { 1.0, 1.0, 1.0 };
    if (input_duplication.params.rescale)
    {
      scale[0] = 1.0 / input_duplication.params.times[0];
      scale[1] = 1.0 / input_duplication.params.times[1];
      scale[2] = 1.0 / input_duplication.params.times[2];
    }

    fcs_int pid = add_dup_input_particles(my_count);

    bool first = true;
    fcs_int first_pid;
    
    for (fcs_int d3 = dlow[3]; d3 <= dhigh[3]; ++d3)
    {
//      cout << "d3 = " << d3 << endl;
      for (fcs_int d2 = dlow[2]; d2 <= ((d3 == dhigh[3])?dhigh[2]:(total_duplications[2] - 1)); ++d2)
      {
//        cout << "d2 = " << d2 << endl;
        for (fcs_int d1 = dlow[1]; d1 <= ((d3 == dhigh[3] && d2 == dhigh[2])?dhigh[1]:(total_duplications[1] - 1)); ++d1)
        {
//          cout << "d1 = " << d1 << endl;
          for (fcs_int d0 = dlow[0]; d0 <= ((d3 == dhigh[3] && d2 == dhigh[2] && d1 == dhigh[1])?dhigh[0]:(total_duplications[0] - 1)); ++d0)
          {
//            cout << "d0 = " << d0 << endl;
            fcs_int x_low = plow;
            fcs_int x_high = (d3 == dhigh[3] && d2 == dhigh[2] && d1 == dhigh[1] && d0 == dhigh[0])?phigh:my_modulo;

//            cout << "x_low = " << x_low << ", x_high = " << x_high << endl;

            if (i == -2)
            {
              for (x = x_low; x < x_high; ++x)
              {
//               cout << "duplicate particle #" << x << endl;

                // Make duplicated particle positions and scale if necessary
                dup_input_positions[3 * pid + 0] = params.offset[0] * (1.0 - scale[0]) +
                  scale[0] * (input_plain.positions[3 * x + 0] + d0 * params.box_a[0] +  d1 * params.box_b[0] +  d2 * params.box_c[0]);
                dup_input_positions[3 * pid + 1] = params.offset[1] * (1.0 - scale[1]) +
                  scale[1] * (input_plain.positions[3 * x + 1] + d0 * params.box_a[1] +  d1 * params.box_b[1] +  d2 * params.box_c[1]);
                dup_input_positions[3 * pid + 2] = params.offset[2] * (1.0 - scale[2]) +
                  scale[2] * (input_plain.positions[3 * x + 2] + d0 * params.box_a[2] +  d1 * params.box_b[2] +  d2 * params.box_c[2]);

                // Duplicated particles have the same charge as their orignals
                dup_input_charges[pid] = input_plain.charges[x];

                // Use same reference potential and field values
                dup_input_potentials[pid] = input_plain.potentials[x];
                dup_input_field[3 * pid + 0] = input_plain.field[3 * x + 0];
                dup_input_field[3 * pid + 1] = input_plain.field[3 * x + 1];
                dup_input_field[3 * pid + 2] = input_plain.field[3 * x + 2];

/*                cout << "duplicated particle #" << pid << ": " << dup_input_positions[3 * pid + 0] << ","
                                                << dup_input_positions[3 * pid + 1] << ","
                                                << dup_input_positions[3 * pid + 2] << " " << dup_input_charges[pid] << endl;*/

                ++pid;
              }

              if (input_plain.have_potentials()) have_reference_values[0] = 1;
              if (input_plain.have_field()) have_reference_values[1] = 1;

            } else
            {
              if (first)
              {
                first = false;
                first_pid = pid;

                if (i < 0)
                {
                  input_files[0].get_local_particles(&dup_input_positions[3 * first_pid], &dup_input_charges[first_pid], &dup_input_potentials[first_pid], &dup_input_field[3 * first_pid],
                    (params.decomposition == DECOMPOSE_ALL_ON_MASTER), comm_size, comm_rank, communicator);

                  if (input_files[0].have_potentials()) have_reference_values[0] = 1;
                  if (input_files[0].have_field()) have_reference_values[1] = 1;

                } else
                {
                  input_generators[i].get_local_particles(&dup_input_positions[3 * first_pid], &dup_input_charges[first_pid], &dup_input_potentials[first_pid], &dup_input_field[3 * first_pid],
                    (params.decomposition == DECOMPOSE_ALL_ON_MASTER), comm_size, comm_rank, communicator);

                  if (input_generators[i].have_potentials()) have_reference_values[0] = 1;
                  if (input_generators[i].have_field()) have_reference_values[1] = 1;
                }
              }

              for (x = x_low; x < x_high; ++x)
              {
//                cout << "duplicate particle #" << x << endl;

                dup_input_positions[3 * pid + 0] = params.offset[0] * (1.0 - scale[0]) +
                  scale[0] * (dup_input_positions[3 * (first_pid + x) + 0] + d0 * params.box_a[0] +  d1 * params.box_b[0] +  d2 * params.box_c[0]);
                dup_input_positions[3 * pid + 1] = params.offset[1] * (1.0 - scale[1]) +
                  scale[1] * (dup_input_positions[3 * (first_pid + x) + 1] + d0 * params.box_a[1] +  d1 * params.box_b[1] +  d2 * params.box_c[1]);
                dup_input_positions[3 * pid + 2] = params.offset[2] * (1.0 - scale[2]) +
                  scale[2] * (dup_input_positions[3 * (first_pid + x) + 2] + d0 * params.box_a[2] +  d1 * params.box_b[2] +  d2 * params.box_c[2]);

                dup_input_charges[pid] = dup_input_charges[first_pid + x];

                dup_input_potentials[pid] = dup_input_potentials[first_pid + x];

                dup_input_field[3 * pid + 0] = dup_input_field[3 * (first_pid + x) + 0];
                dup_input_field[3 * pid + 1] = dup_input_field[3 * (first_pid + x) + 1];
                dup_input_field[3 * pid + 2] = dup_input_field[3 * (first_pid + x) + 2];

/*                cout << "duplicated particle #" << pid << ": "
                                                << dup_input_positions[3 * pid + 0] << "," << dup_input_positions[3 * pid + 1] << "," << dup_input_positions[3 * pid + 2] << " "
                                                << dup_input_charges[pid] << " "
                                                << dup_input_potentials[pid] << " "
                                                << dup_input_field[3 * pid + 0] << "," << dup_input_field[3 * pid + 1] << "," << dup_input_field[3 * pid + 2] << endl;*/

                ++pid;
              }
            }

            plow = 0;
          }
          dlow[0] = 0;
        }
        dlow[1] = 0;
      }
      dlow[2] = 0;
    }
    dlow[3] = 0;
  }

  params.box_a[0] *= total_duplications[0]; params.box_a[1] *= total_duplications[0]; params.box_a[2] *= total_duplications[0];
  params.box_b[0] *= total_duplications[1]; params.box_b[1] *= total_duplications[1]; params.box_b[2] *= total_duplications[1];
  params.box_c[0] *= total_duplications[2]; params.box_c[1] *= total_duplications[2]; params.box_c[2] *= total_duplications[2];

  unscale_box[0] /= total_duplications[0];
  unscale_box[1] /= total_duplications[1];
  unscale_box[2] /= total_duplications[2];

  if (input_duplication.params.rescale)
  {
    params.box_a[0] /= input_duplication.params.times[0]; params.box_a[1] /= input_duplication.params.times[0]; params.box_a[2] /= input_duplication.params.times[0];
    params.box_b[0] /= input_duplication.params.times[1]; params.box_b[1] /= input_duplication.params.times[1]; params.box_b[2] /= input_duplication.params.times[1];
    params.box_c[0] /= input_duplication.params.times[2]; params.box_c[1] /= input_duplication.params.times[2]; params.box_c[2] /= input_duplication.params.times[2];

    unscale_box[0] *= input_duplication.params.times[0];
    unscale_box[1] *= input_duplication.params.times[1];
    unscale_box[2] *= input_duplication.params.times[2];
  }

  DEBUG(cout << comm_rank << ": local input particles: " << dup_input_nparticles << endl);

  MPI_Allreduce(&dup_input_nparticles, &dup_input_total_nparticles, 1, FCS_MPI_INT, MPI_SUM, communicator);

  if (input_ref.get_total_nparticles() > 0)
  {
    input_ref.get_local_particles(dup_input_positions, dup_input_charges, dup_input_potentials, dup_input_field,
      dup_input_nparticles, comm_size, comm_rank, communicator);

    if (input_ref.have_potentials()) have_reference_values[0] = 1;
    if (input_ref.have_field()) have_reference_values[1] = 1;
  }
}

void Configuration::create_cart_comm()
{
  int cart_dims[3], cart_periods[3];

  cart_dims[0] = cart_dims[1] = cart_dims[2] = 0;
  MPI_Dims_create(comm_size, 3, cart_dims);

  cart_periods[0] = params.periodicity[0];
  cart_periods[1] = params.periodicity[1];
  cart_periods[2] = params.periodicity[2];

  INFO_MASTER(
    cout << "Creating " << cart_dims[0] << "x" << cart_dims[1] << "x" << cart_dims[2]
         << " Cartesian communicator (periodicity: " << cart_periods[0] << ","  << cart_periods[1] << "," << cart_periods[2] << ")" << endl;
  );

  MPI_Cart_create(communicator, 3, cart_dims, cart_periods, 0, &cart_comm);
}

void Configuration::destroy_cart_comm()
{
  if (cart_comm != MPI_COMM_NULL) MPI_Comm_free(&cart_comm);
}

void Configuration::decompose_particles()
{
  generate_input_particles();
  
  MPI_Bcast(have_reference_values, 2, FCS_MPI_INT, MASTER_RANK, communicator);

  decomp_comm = communicator;
  
  decomp_total_nparticles = dup_input_total_nparticles;

  switch (params.decomposition)
  {
    case DECOMPOSE_ALL_ON_MASTER:
      INFO_MASTER(cout << "Decomposing system (all-on-master)..." << endl);
      if (comm_rank == 0) {
        decomp_nparticles = dup_input_nparticles;
        decomp_positions = dup_input_positions;
        decomp_charges = dup_input_charges;
      } else {
        decomp_nparticles = 0;
        decomp_positions = 0;
        decomp_charges = 0;
      }
      break;
    case DECOMPOSE_ATOMISTIC:
      INFO_MASTER(cout << "Decomposing system (atomistic)..." << endl);
      decomp_nparticles = dup_input_nparticles;
      decomp_positions = dup_input_positions;
      decomp_charges = dup_input_charges;
      break;
    case DECOMPOSE_RANDOM:
      INFO_MASTER(cout << "Decomposing system (random)..." << endl);
      fcs_gridsort_create(&gridsort);
      fcs_gridsort_set_particles(&gridsort, dup_input_nparticles, dup_input_positions, dup_input_charges);
      fcs_gridsort_sort_random(&gridsort, communicator);
      fcs_gridsort_get_sorted_particles(&gridsort, &decomp_nparticles, &decomp_positions, &decomp_charges, NULL);
      break;
    case DECOMPOSE_DOMAIN:
      INFO_MASTER(cout << "Decomposing system (domain)..." << endl);
      create_cart_comm();
      fcs_gridsort_create(&gridsort);
      fcs_gridsort_set_system(&gridsort, params.offset, params.box_a, params.box_b, params.box_c, NULL);
      fcs_gridsort_set_particles(&gridsort, dup_input_nparticles, dup_input_positions, dup_input_charges);
      fcs_gridsort_sort_forward(&gridsort, 0.0, cart_comm);
      fcs_gridsort_get_sorted_particles(&gridsort, &decomp_nparticles, &decomp_positions, &decomp_charges, NULL);
      decomp_comm = cart_comm;
      break;
    default:
      break;
  }

  if (decomp_nparticles > 0)
  {
    decomp_potentials = new fcs_float[decomp_nparticles];
    decomp_field = new fcs_float[3*decomp_nparticles];
    for (fcs_int i = 0; i < decomp_nparticles; ++i)
    {
      decomp_potentials[i] = NAN;
      decomp_field[3 * i + 0] = decomp_field[3 * i + 1] = decomp_field[3 * i + 2] = NAN;
    }
  } else
  {
    decomp_potentials = (fcs_float *) 1;
    decomp_field = (fcs_float *) 1;
  }
}

void Configuration::gather_particles()
{
  result_potentials = 0;
  result_field = 0;

  switch (params.decomposition)
  {
    case DECOMPOSE_ALL_ON_MASTER:
      INFO_MASTER(cout << "Gathering results (all-on-master)..." << endl);
      if (comm_rank == 0) {
        result_potentials = decomp_potentials;
        result_field = decomp_field;
      }
      break;
    case DECOMPOSE_ATOMISTIC:
      INFO_MASTER(cout << "Gathering results (atomistic)..." << endl);
      result_potentials = decomp_potentials;
      result_field = decomp_field;
      break;
    case DECOMPOSE_RANDOM:
      INFO_MASTER(cout << "Gathering results (random)..." << endl);
      if (dup_input_nparticles > 0) {
        result_potentials = new fcs_float[dup_input_nparticles];
        result_field = new fcs_float[3*dup_input_nparticles];
      }
      fcs_gridsort_sort_backward(&gridsort, decomp_field, decomp_potentials,
        (result_field)?result_field:(fcs_float *)1, (result_potentials)?result_potentials:(fcs_float *)1, 1, decomp_comm);
      break;
    case DECOMPOSE_DOMAIN:
      INFO_MASTER(cout << "Gathering results (domain)..." << endl);
      if (dup_input_nparticles > 0) {
        result_potentials = new fcs_float[dup_input_nparticles];
        result_field = new fcs_float[3*dup_input_nparticles];
      }
      fcs_gridsort_sort_backward(&gridsort, decomp_field, decomp_potentials,
        (result_field)?result_field:(fcs_float *)1, (result_potentials)?result_potentials:(fcs_float *)1, 1, decomp_comm);
      break;
    default:
      break;
  }
}

bool Configuration::compute_errors(errors_t *e) {

  ::compute_errors(e, dup_input_nparticles, dup_input_positions, dup_input_charges,
    have_reference_values[0]?dup_input_potentials:NULL, have_reference_values[1]?dup_input_field:NULL,
    have_result_values[0]?result_potentials:NULL, have_result_values[1]?result_field:NULL, decomp_comm);

  return true;
}

void Configuration::free_decomp_particles(bool quiet)
{
  switch (params.decomposition)
  {
    case DECOMPOSE_ALL_ON_MASTER:
      if (!quiet) INFO_MASTER(cout << "Freeing data (all-on-master)..." << endl);
      break;
    case DECOMPOSE_ATOMISTIC:
      if (!quiet) INFO_MASTER(cout << "Freeing data (atomistic)..." << endl);
      break;
    case DECOMPOSE_RANDOM:
      if (!quiet) INFO_MASTER(cout << "Freeing data (random)..." << endl);
      fcs_gridsort_free(&gridsort);
      fcs_gridsort_destroy(&gridsort);
      delete[] result_potentials;
      delete[] result_field;
      break;
    case DECOMPOSE_DOMAIN:
      if (!quiet) INFO_MASTER(cout << "Freeing data (domain)..." << endl);
      fcs_gridsort_free(&gridsort);
      fcs_gridsort_destroy(&gridsort);
      destroy_cart_comm();
      delete[] result_potentials;
      delete[] result_field;
      break;
    default:
      break;
  }

  if (decomp_nparticles > 0)
  {
    delete[] decomp_potentials;
    delete[] decomp_field;
  }

  decomp_nparticles = 0;
  decomp_positions = 0;
  decomp_charges = 0;
  decomp_potentials = 0;
  decomp_field = 0;

  decomp_comm = MPI_COMM_NULL;

  result_potentials = 0;
  result_field = 0;
}

void Configuration::determine_total_duplication()
{
  total_duplications[0] = ((params.periodicity[0])?params.periodic_duplications[0]:1);
  total_duplications[1] = ((params.periodicity[1])?params.periodic_duplications[0]:1);
  total_duplications[2] = ((params.periodicity[2])?params.periodic_duplications[0]:1);

  total_duplications[0] *= input_duplication.params.times[0];
  total_duplications[1] *= input_duplication.params.times[1];
  total_duplications[2] *= input_duplication.params.times[2];
  
  total_duplication = total_duplications[0] * total_duplications[1] * total_duplications[2];
}

Testcase::Testcase()
  : name(""), description(""), reference_method(""),
    error_potential(1.0), error_field(1.0),
    configurations()
{}

void
Testcase::read_file(const char* filename, fcs_int *periodic_duplications, fcs_int decomposition) {
  static const fcs_int READSIZE = 8192;
  ZLIB_IFELSE(gzFile, FILE*) inputfile;
  fcs_int inputsize;
  char* inputdata;
  char* datap;
  int dataread;
  xml_document<> doc;
  xml_attribute<> *attr;
  string aname;
  fcs_int config_count;

  cout << "Reading testcase file " << filename << "..." << endl;
  inputfile = ZLIB_IFELSE(gzopen,fopen)(filename, ZLIB_IFELSE("rb","r"));
  if (!inputfile) {
    string msg = filename;
    msg += ": ";
    msg += strerror(errno);
    throw ParserError(msg);
  }

  // read in blocks
  inputsize = READSIZE + 1;
  inputdata = (char*)malloc(inputsize);
  datap = inputdata;
  while ((dataread = ZLIB_IFELSE(gzread(inputfile, datap, READSIZE),fread(datap, 1, READSIZE, inputfile))) == READSIZE) {
    inputsize += READSIZE;
    inputdata = (char*)realloc(inputdata, inputsize);
    datap = inputdata + inputsize - READSIZE - 1;
  }
  datap += dataread;
  *datap = 0;

  ZLIB_IFELSE(gzclose,fclose)(inputfile);
  
  cout << "Parsing file..." << endl;
  try {
    doc.parse<0>(inputdata);
  } catch (parse_error &er) {
    string msg = filename;
    msg += ": ";
    msg += er.what();
    msg += " instead of ";
    msg += er.where<char>()[0];
    throw ParserError(msg);
  }

  // Parse system data
  xml_node<> *system_node = doc.first_node(TESTCASE_TAG);
  if (!system_node) {
    string msg = filename;
    msg += ": did not find top level element ";
    msg += TESTCASE_TAG;
    throw ParserError(msg);
  }

  for (attr = system_node->first_attribute();
       attr; attr = attr->next_attribute()) {
    aname = attr->name();

    if (aname == "name")
      this->name = attr->value();
    else if (aname == "description")
      this->description = attr->value();
    else if (aname == "reference_method")
      this->reference_method = attr->value();
    else if (aname == "error_potential") {
      istringstream is(attr->value());
      is.exceptions(istream::failbit | istream::badbit);
      is >> this->error_potential;
    } else if (aname == "error_field") {
      istringstream is(attr->value());
      is.exceptions(istream::failbit | istream::badbit);
      is >> this->error_field;
    }
  }

  // Output system data
  cout << "Read testcase \"" << this->name << "\"" << endl;
  cout << "  \"" << this->description << '"' << endl;
  cout << "  Reference: " 
       << this->reference_method
       << " (error_potential=" << this->error_potential
       << " error_field=" << this->error_field
       << ")" << endl;

  for (xml_node<> *config_node = system_node->first_node(METHOD_TAG); config_node; config_node = config_node->next_sibling(METHOD_TAG)) {
    read_method_config(config_node);
  }

  config_count = 0;
  
  string basename = filename;
  const char *t = strrchr(filename, '/');
  if (t == NULL) basename.resize(0);
  else basename.resize(t - filename + 1);

  // Loop over all configuration nodes
  for (xml_node<> *config_node = system_node->first_node(CONFIGURATION_TAG);
    config_node; config_node = config_node->next_sibling(CONFIGURATION_TAG)) {
    this->configurations.push_back(new Configuration());
    Configuration *config = this->configurations.back();

    config->params.periodic_duplications[0] = periodic_duplications[0];
    config->params.periodic_duplications[1] = periodic_duplications[1];
    config->params.periodic_duplications[2] = periodic_duplications[2];

    config->read_config(config_node, basename.c_str());

    // override configured decomposition mode
    if (decomposition >= 0) config->params.decomposition = decomposition;

    config->print_config(config_count);

    config_count++;
  }
  
  cout << "  Got " << this->configurations.size() 
       << " configurations." << endl;

  free(inputdata);
}

void Testcase::write_file(const char* outfilename, const char* binfilename, const char* portable_filename, bool keep_dupgen) {
  xml_document<> doc;
  xml_attribute<> *attr;
  xml_node<> *system_node;
  ostringstream os;
  char* s;

  os.precision(16);

  if (comm_rank == MASTER_RANK)
  {
    xml_node<> *declaration_node =
      doc.allocate_node(node_declaration);
    doc.append_node(declaration_node);

    attr = doc.allocate_attribute("version", "1.0");
    declaration_node->append_attribute(attr);

    attr = doc.allocate_attribute("encoding", "UTF-8");
    declaration_node->append_attribute(attr);

    xml_node<> *doctype_node =
      doc.allocate_node(node_doctype);
    const char* doctype = "scafacos_test SYSTEM 'scafacos_test.dtd'";
    doctype_node->value(doctype, strlen(doctype));
    doc.append_node(doctype_node);

    // Create system node
    system_node = 
      doc.allocate_node(node_element, TESTCASE_TAG);
    doc.append_node(system_node);

    attr = doc.allocate_attribute("name", 
      this->name.c_str());
    system_node->append_attribute(attr);

    attr = doc.allocate_attribute("description", 
      this->description.c_str());
    system_node->append_attribute(attr);

    attr = doc.allocate_attribute("reference_method", 
      this->reference_method.c_str());
    system_node->append_attribute(attr);

    os.str("");
    os << this->error_potential;
    s = doc.allocate_string(os.str().c_str());
    attr = doc.allocate_attribute("error_potential", s);
    system_node->append_attribute(attr);

    os.str("");
    os << this->error_field;
    s = doc.allocate_string(os.str().c_str());
    attr = doc.allocate_attribute("error_field", s);
    system_node->append_attribute(attr);
  }

  if (binfilename != NULL)
  {
    if (comm_rank == MASTER_RANK)
      cout << "Writing particle data to file " << binfilename << "..." << endl;
  }

  for (vector<Configuration*>::iterator config = this->configurations.begin();
       config !=  this->configurations.end(); config++) {

    if (comm_rank == MASTER_RANK) {
      // Create config node
      xml_node<> *config_node = doc.allocate_node(node_element, CONFIGURATION_TAG);
      system_node->append_node(config_node);

      (*config)->write_config(&doc, config_node, binfilename, portable_filename, keep_dupgen);

    } else (*config)->write_config(NULL, NULL, binfilename, portable_filename, keep_dupgen);
  }

  if (comm_rank == MASTER_RANK) {
    cout << "Writing testcase data to file " << outfilename << "..." << endl;
    ZLIB_IFELSE(gzFile,FILE*) outputfile = ZLIB_IFELSE(gzopen,fopen)(outfilename, ZLIB_IFELSE("wb","w"));

    // Print to string using output iterator
    std::string xmlstring;
    print(std::back_inserter(xmlstring), doc, 0);
    ZLIB_IFELSE(gzwrite(outputfile, xmlstring.c_str(), strlen(xmlstring.c_str())),fwrite(xmlstring.c_str(), 1, strlen(xmlstring.c_str()), outputfile));
    ZLIB_IFELSE(gzclose,fclose)(outputfile);
  }
}

void Testcase::broadcast_config(int root, MPI_Comm comm)
{
  fcs_int n;
  if (root == comm_rank) n = method_config.length();
  MPI_Bcast(&n, 1, FCS_MPI_INT, root, comm);

  if (root != comm_rank) method_config.resize(n);

  MPI_Bcast((void *) method_config.c_str(), n, MPI_CHAR, root, comm);
}

const char *Testcase::get_method_config()
{
  return method_config.c_str();
}

void Testcase::read_method_config(xml_node<> *config_node)
{
  for (xml_node<> *param_node = config_node->first_node(PARAM_TAG); param_node; param_node = param_node->next_sibling(PARAM_TAG))
  {
    xml_attribute<> *name_attr = param_node->first_attribute("name");

    if (name_attr == 0) continue;

    vector<xml_attribute<> *> value_attrs;
    xml_attribute<> *attr;

    attr = param_node->first_attribute("value");
    if (attr) value_attrs.push_back(attr);

    for (int i = 0; i < 3; ++i)
    {
      char valstr[16];
      sprintf(valstr, "value%d", i);
      attr = param_node->first_attribute(valstr);
      if (attr) value_attrs.push_back(attr);
    }

    INFO_MASTER(cout << "    set parameter '" << name_attr->value() << "' to value(s)";);
    if (method_config.size() > 0) method_config += ",";
    method_config += name_attr->value();
    for (fcs_int i = 0; i < (fcs_int) value_attrs.size(); ++i)
    {
      INFO_MASTER(cout << " '" << value_attrs[i]->value() << "'";);
      method_config += ",";
      method_config += value_attrs[i]->value();
    }
    INFO_MASTER(cout << endl;);
  }
}
