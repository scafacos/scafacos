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

#ifndef _TESTCASE_HPP
#define _TESTCASE_HPP

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <vector>
#include <string>

#include "fcs.h"

#include "common/gridsort/gridsort.h"

#include "common.hpp"
#include "Generator.hpp"

#include "rapidxml/rapidxml.hpp"

using namespace std;
using namespace rapidxml;

static const fcs_float EPSILON_METALLIC = -1.0;


class Configuration {
public:
  struct {
    fcs_float offset[3];
    fcs_float box_a[3];
    fcs_float box_b[3];
    fcs_float box_c[3];
    fcs_int periodicity[3];
    fcs_float epsilon;
    fcs_int decomposition;
    fcs_int periodic_duplications[3];

    fcs_float result_density;

  } params;

  fcs_int total_duplications[3], total_duplication;

  fcs_float unscale_box[3];

  PlainParticles input_plain;
  fcs_int input_plain_nparticles;

  vector<FileParticles> input_files;
  fcs_int input_file_nparticles;

  vector<Generator> input_generators;
  fcs_int input_generator_nparticles;
  
  Duplicate input_duplication;

  FileParticles input_ref;
  fcs_int input_ref_nparticles;

  fcs_int dup_input_total_nparticles;
  fcs_int dup_input_nparticles, dup_input_nparticles_allocated;
  fcs_float dup_input_overalloc;
  fcs_float *dup_input_positions;
  fcs_float *dup_input_charges;
  fcs_float *dup_input_potentials;
  fcs_float *dup_input_field;

  fcs_int have_reference_values[2], have_result_values[2];

  fcs_int decomp_total_nparticles;
  fcs_int decomp_nparticles, decomp_max_nparticles;
  fcs_float *decomp_positions;
  fcs_float *decomp_charges;
  fcs_float *decomp_potentials;
  fcs_float *decomp_field;
  MPI_Comm decomp_comm;

  fcs_float *reference_potentials, *reference_field;
  
  fcs_float field_correction[3], energy_correction;

  MPI_Comm cart_comm;
  fcs_gridsort_t gridsort;

  Configuration();
  ~Configuration();

  fcs_int add_dup_input_particles(fcs_int add_nparticles);
  void free_dup_input_particles();

  void read_config(xml_node<> *config_node, const char *basename);
  void write_config(xml_document<> *doc, xml_node<> *config_node, const char *binfilename = NULL, const char* portable_filename = NULL, bool keep_dupgen = false);

  void write_particles(xml_document<> *doc, xml_node<> *config_node);

  void print_config(fcs_int n);

  void broadcast_config();

  void broadcast_input();
  void generate_input_particles();

  void create_cart_comm();
  void destroy_cart_comm();

  void decompose_particles(fcs_float overalloc = 0);
  void gather_particles();
  bool compute_errors(errors_t *e);
  void free_decomp_particles(bool quiet = false);

private:
  void determine_total_duplication();
};


class Testcase {
public:
  string name;
  string description;
  string reference_method;
  fcs_float error_potential;
  fcs_float error_field;
  vector<Configuration*> configurations;

  Testcase();

  void read_file(const char* filename, fcs_int *periodic_duplications, fcs_int decomposition = -1);
  void write_file(const char* outfilename, const char* binfilename = NULL, const char* portable_filename = NULL, bool keep_dupgen = false);

  void broadcast_config(int root, MPI_Comm comm);

  const char *get_method_config();

private:
  string method_config;

  void read_method_config(xml_node<> *config_node);
};

#endif
