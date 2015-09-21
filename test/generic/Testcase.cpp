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

#include "common/fcs-common/FCSCommon.h"

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


Configuration::Configuration()
{
  params.box_origin[0] = params.box_origin[1] = params.box_origin[2] = 0.0;
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

  input_particles_allocated = false;

  have_reference_values[0] = have_reference_values[1] = 0;
  have_result_values[0] = have_result_values[1] = 0;

#if SCAFACOS_TEST_WITH_DIPOLES
  dipole_have_reference_values[0] = dipole_have_reference_values[1] = 0;
  dipole_have_result_values[0] = dipole_have_result_values[1] = 0;
#endif /* SCAFACOS_TEST_WITH_DIPOLES */

  decomp_comm = MPI_COMM_NULL;

  reference_potentials = 0;
  reference_field = 0;

#if SCAFACOS_TEST_WITH_DIPOLES
  dipole_reference_potentials = 0;
  dipole_reference_field = 0;
#endif /* SCAFACOS_TEST_WITH_DIPOLES */

  field_correction[0] = field_correction[1] = field_correction[2] = 0;
  energy_correction = 0;

  cart_comm = MPI_COMM_NULL;
}

Configuration::~Configuration()
{
  input_particles.particles.free();
  input_particles.particles.free();
  
  free_decomp_particles(true);
}


void Configuration::read_config(xml_node<> *config_node, const char *basename)
{
  xml_attribute<> *attr;
  string aname;

  // Read configuration data
  for (attr = config_node->first_attribute(); attr; attr = attr->next_attribute()) {
    aname = attr->name();
    if (aname == "offset") parse_sequence(attr->value(), 3, params.box_origin);
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

    os.str(string(""));
    os << params.box_origin[0] << " " << params.box_origin[1] << " " << params.box_origin[2];
    s = doc->allocate_string(os.str().c_str());
    attr = doc->allocate_attribute("offset", s);
    config_node->append_attribute(attr);

    os.str(string(""));
    os << out_box_a[0] << " " << out_box_a[1] << " " << out_box_a[2];
    s = doc->allocate_string(os.str().c_str());
    attr = doc->allocate_attribute("box_a", s);
    config_node->append_attribute(attr);

    os.str(string(""));
    os << out_box_b[0] << " " << out_box_b[1] << " " << out_box_b[2];
    s = doc->allocate_string(os.str().c_str());
    attr = doc->allocate_attribute("box_b", s);
    config_node->append_attribute(attr);

    os.str(string(""));
    os << out_box_c[0] << " " << out_box_c[1] << " " << out_box_c[2];
    s = doc->allocate_string(os.str().c_str());
    attr = doc->allocate_attribute("box_c", s);
    config_node->append_attribute(attr);

    os.str(string(""));
    os << params.periodicity[0] << " " << params.periodicity[1] << " " << params.periodicity[2];
    s = doc->allocate_string(os.str().c_str());
    attr = doc->allocate_attribute("periodicity", s);
    config_node->append_attribute(attr);

    if (params.epsilon == EPSILON_METALLIC)
      attr = doc->allocate_attribute("epsilon", "metallic");
    else if (params.epsilon == 0.0)
      attr = doc->allocate_attribute("epsilon", "vacuum");
    else {
      os.str(string(""));
      os << params.epsilon;
      s = doc->allocate_string(os.str().c_str());
      attr = doc->allocate_attribute("epsilon", s);
    }
    config_node->append_attribute(attr);
  }

  if (keep_dupgen)
  {
    xml_node<> *particles_node = input_duplication.write_config(doc, config_node);

    particles_t particles;

    particles.total_nparticles = input_plain.get_total_nparticles();
    particles.particles = input_plain.local_particles;
#if SCAFACOS_TEST_WITH_DIPOLES
    particles.dipole_total_nparticles = input_plain.get_dipole_total_nparticles();
    particles.dipole_particles = input_plain.dipole_local_particles;
#endif /* SCAFACOS_TEST_WITH_DIPOLES */
    PlainParticles::write_config(doc, particles_node, NULL, &particles, comm_size, comm_rank, communicator);

    if (binfilename)
    {
      FileParticles::write_config<FormatBinary>(doc, config_node, REFERENCES_TAG, binfilename, input_particles.total_nparticles, input_particles.particles.n,
        NULL, NULL,
        (have_reference_values[0] || have_result_values[0])?(input_particles.particles.potentials?input_particles.particles.potentials:(fcs_float *) 1):NULL,
        (have_reference_values[1] || have_result_values[1])?(input_particles.particles.field?input_particles.particles.field:(fcs_float *) 1):NULL,
        comm_size, comm_rank, communicator);

    } else if (portable_filename)
    {
      FileParticles::write_config<FormatPortable>(doc, config_node, REFERENCES_TAG, portable_filename, input_particles.total_nparticles, input_particles.particles.n,
        NULL, NULL,
        (have_reference_values[0] || have_result_values[0])?(input_particles.particles.potentials?input_particles.particles.potentials:(fcs_float *) 1):NULL,
        (have_reference_values[1] || have_result_values[1])?(input_particles.particles.field?input_particles.particles.field:(fcs_float *) 1):NULL,
        comm_size, comm_rank, communicator);

    } else {
      /* TODO: write references plain? */
      MASTER(cout << "WARNING: writing references to XML file when keeping duplication/generation information not implemented, yet!" << endl);
    }

  } else {

    if (binfilename)
      FileParticles::write_config<FormatBinary>(doc, config_node, BINARY_TAG, binfilename, input_particles.total_nparticles, input_particles.particles.n,
        input_particles.particles.positions?input_particles.particles.positions:(fcs_float *) 1,
        input_particles.particles.props?input_particles.particles.props:(fcs_float *) 1,
        (have_reference_values[0] || have_result_values[0])?(input_particles.particles.potentials?input_particles.particles.potentials:(fcs_float *) 1):NULL,
        (have_reference_values[1] || have_result_values[1])?(input_particles.particles.field?input_particles.particles.field:(fcs_float *) 1):NULL,
        comm_size, comm_rank, communicator);
    else if (portable_filename)
      FileParticles::write_config<FormatPortable>(doc, config_node, PORTABLE_TAG, portable_filename, input_particles.total_nparticles, input_particles.particles.n,
        input_particles.particles.positions?input_particles.particles.positions:(fcs_float *) 1,
        input_particles.particles.props?input_particles.particles.props:(fcs_float *) 1,
        (have_reference_values[0] || have_result_values[0])?(input_particles.particles.potentials?input_particles.particles.potentials:(fcs_float *) 1):NULL,
        (have_reference_values[1] || have_result_values[1])?(input_particles.particles.field?input_particles.particles.field:(fcs_float *) 1):NULL,
        comm_size, comm_rank, communicator);
    else
    {
      particles_t particles = input_particles;

      if (!(have_reference_values[0] || have_result_values[0])) particles.particles.potentials = NULL;
      if (!(have_reference_values[1] || have_result_values[1])) particles.particles.field = NULL;
      
      PlainParticles::write_config(doc, config_node, NULL, &particles, comm_size, comm_rank, communicator);
    }
  }
}

void Configuration::print_config(fcs_int n)
{
  cout << "  config " << n << ":" << endl;
  cout << "    system:" << endl;
  cout << "      box: box_origin: [" << params.box_origin[0] << "," << params.box_origin[1] << "," << params.box_origin[2] << "], size: "
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
    case DECOMPOSE_ALMOST_ALL_ON_MASTER: cout << "almost-all-on-master" << endl; break;
    case DECOMPOSE_ATOMISTIC: cout << "atomistic" << endl; break;
    case DECOMPOSE_RANDOM: cout << "random" << endl; break;
    case DECOMPOSE_RANDOM_EQUAL: cout << "random-equal" << endl; break;
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
  input_plain.broadcast_config(MASTER_RANK, communicator);

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

void Configuration::generate_input_particles(fcs_float minalloc, fcs_float overalloc)
{
  broadcast_input();
  
  bool all_on_master = (params.decomposition == DECOMPOSE_ALL_ON_MASTER || params.decomposition == DECOMPOSE_ALMOST_ALL_ON_MASTER);

  for (fcs_int i = -2; i < (fcs_int) input_generators.size(); ++i)
  {
    fcs_float scale[3] = { 1.0, 1.0, 1.0 };
    if (input_duplication.params.rescale)
    {
      scale[0] = 1.0 / input_duplication.params.times[0];
      scale[1] = 1.0 / input_duplication.params.times[1];
      scale[2] = 1.0 / input_duplication.params.times[2];
    }

    fcs_int local_nparticles;
#if SCAFACOS_TEST_WITH_DIPOLES
    fcs_int dipole_local_nparticles;
#endif /* SCAFACOS_TEST_WITH_DIPOLES */

    if (i == -2)
    {
      local_nparticles = input_plain.get_local_nparticles(all_on_master, comm_size, comm_rank, communicator);
      DEBUG(cout << comm_rank << ": get from plain: " << local_nparticles << endl);
#if SCAFACOS_TEST_WITH_DIPOLES
      dipole_local_nparticles = input_plain.get_dipole_local_nparticles(all_on_master, comm_size, comm_rank, communicator);
      DEBUG(cout << comm_rank << ": get dipoles from plain: " << dipole_local_nparticles << endl);
#endif /* SCAFACOS_TEST_WITH_DIPOLES */

    } else if (i == -1)
    {
      if (input_files.size() < 1) continue;

      local_nparticles = input_files[0].get_local_nparticles(all_on_master, comm_size, comm_rank, communicator);
      DEBUG(cout << comm_rank << ": get from file #" << i << ": local_nparticles: " << local_nparticles << endl);

    } else
    {
      input_generators[i].set_box(params.box_origin, params.box_a, params.box_b, params.box_c);
      local_nparticles = input_generators[i].get_local_nparticles(all_on_master, comm_size, comm_rank, communicator);
      DEBUG(cout << comm_rank << ": get from generator #" << i << ": local_nparticles: " << local_nparticles << endl);
    }

    input_particles.particles.realloc(input_particles.particles.n + (local_nparticles * total_duplication));
    fcs_int pid = input_particles.particles.n;
    fcs_int first_pid = pid;
#if SCAFACOS_TEST_WITH_DIPOLES
    input_particles.dipole_particles.realloc(input_particles.dipole_particles.n + (dipole_local_nparticles * total_duplication));
    fcs_int dipole_pid = input_particles.dipole_particles.n;
    fcs_int dipole_first_pid = dipole_pid;
#endif /* SCAFACOS_TEST_WITH_DIPOLES */

    if (i == -2)
    {
      input_plain.get_local_particles(&input_particles.particles, all_on_master, comm_size, comm_rank, communicator);

      if (input_plain.have(PDT_CHARGE_POTENTIALS)) have_reference_values[0] = 1;
      if (input_plain.have(PDT_CHARGE_FIELD)) have_reference_values[1] = 1;

#if SCAFACOS_TEST_WITH_DIPOLES
      input_plain.get_dipole_local_particles(&input_particles.dipole_particles, all_on_master, comm_size, comm_rank, communicator);

      if (input_plain.have(PDT_DIPOLE_POTENTIALS)) dipole_have_reference_values[0] = 1;
      if (input_plain.have(PDT_DIPOLE_FIELD)) dipole_have_reference_values[1] = 1;
#endif /* SCAFACOS_TEST_WITH_DIPOLES */

    } else if (i == -1)
    {
      input_files[0].get_local_particles(&input_particles.particles, all_on_master, comm_size, comm_rank, communicator);

      if (input_files[0].have_potentials()) have_reference_values[0] = 1;
      if (input_files[0].have_field()) have_reference_values[1] = 1;

    } else
    {
      input_generators[i].get_local_particles(&input_particles.particles, all_on_master, comm_size, comm_rank, communicator);

      if (input_generators[i].have_potentials()) have_reference_values[0] = 1;
      if (input_generators[i].have_field()) have_reference_values[1] = 1;
    }

    for (fcs_int d2 = 0; d2 < total_duplications[2]; ++d2)
    {
      for (fcs_int d1 = 0; d1 < total_duplications[1]; ++d1)
      {
        for (fcs_int d0 = 0; d0 < total_duplications[0]; ++d0)
        {
/*          cout << "d0 = " << d0 << ", d1 = " << d1 << ", d2 = " << d2 << endl;*/

          fcs_float offset[3];

          offset[0] = params.box_origin[0] * (1.0 - scale[0]) + scale[0] * (d0 * params.box_a[0] +  d1 * params.box_b[0] +  d2 * params.box_c[0]);
          offset[1] = params.box_origin[1] * (1.0 - scale[1]) + scale[1] * (d0 * params.box_a[1] +  d1 * params.box_b[1] +  d2 * params.box_c[1]);
          offset[2] = params.box_origin[2] * (1.0 - scale[2]) + scale[2] * (d0 * params.box_a[2] +  d1 * params.box_b[2] +  d2 * params.box_c[2]);

          for (fcs_int x = 0; x < local_nparticles; ++x)
          {
/*            cout << "rescale/duplicate particle #" << first_pid + x << " to new #" << pid << endl;*/

            input_particles.particles.copy(first_pid + x, pid);

            fcs_float *pos = input_particles.particles.positions_at(pid);

            pos[0] = offset[0] + (pos[0] * scale[0]);
            pos[1] = offset[1] + (pos[1] * scale[1]);
            pos[2] = offset[2] + (pos[2] * scale[2]);

/*            cout << "rescaled/duplicated particle #" << pid << ": ";
            input_particles.particles.print(pid);
            cout << endl;*/

            ++pid;
          }

#if SCAFACOS_TEST_WITH_DIPOLES
          for (fcs_int x = 0; x < dipole_local_nparticles; ++x)
          {
/*            cout << "rescale/duplicate dipole particle #" << dipole_first_pid + x << " to new #" << dipole_pid << endl;*/

            input_particles.dipole_particles.copy(dipole_first_pid + x, dipole_pid);

            fcs_float *pos = input_particles.dipole_particles.positions_at(dipole_pid);

            pos[0] = offset[0] + (pos[0] * scale[0]);
            pos[1] = offset[1] + (pos[1] * scale[1]);
            pos[2] = offset[2] + (pos[2] * scale[2]);

/*            cout << "rescaled/duplicated dipole particle #" << dipole_pid << ": ";
            input_particles.dipole_particles.print(dipole_pid);
            cout << endl;*/

            ++dipole_pid;
          }
#endif /* SCAFACOS_TEST_WITH_DIPOLES */

          if (d0 + d1 + d2 > 0)
          {
            input_particles.particles.n += local_nparticles;
#if SCAFACOS_TEST_WITH_DIPOLES
            input_particles.dipole_particles.n += dipole_local_nparticles;
#endif /* SCAFACOS_TEST_WITH_DIPOLES */
          }
        }
      }
    }
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

  DEBUG(cout << comm_rank << ": local input particles: " << input_particles.particles.n << endl);

  fcs_int nlocals[2], ntotals[2];

  nlocals[0] = input_particles.particles.n;
#if SCAFACOS_TEST_WITH_DIPOLES
  nlocals[1] = input_particles.dipole_particles.n;
#endif /* SCAFACOS_TEST_WITH_DIPOLES */

  MPI_Allreduce(nlocals, ntotals, 2, FCS_MPI_INT, MPI_SUM, communicator);

  input_particles.total_nparticles = ntotals[0];
#if SCAFACOS_TEST_WITH_DIPOLES
  input_particles.dipole_total_nparticles = ntotals[1];
#endif /* SCAFACOS_TEST_WITH_DIPOLES */

  if (input_ref.get_total_nparticles() > 0)
  {
    fcs_int n = input_particles.particles.n;

    input_particles.particles.n = 0;

    input_ref.get_local_particles(&input_particles.particles, n, comm_size, comm_rank, communicator);

    if (input_ref.have_potentials()) have_reference_values[0] = 1;
    if (input_ref.have_field()) have_reference_values[1] = 1;
  }

#if SCAFACOS_TEST_WITH_DIPOLES && 0
  /* FIXME: */
  if (input_ref.get_total_nparticles() > 0)
  {
    fcs_int n = input_particles.particles.n;

    input_particles.particles.n = 0;

    input_ref.get_local_particles(&input_particles.particles, n, comm_size, comm_rank, communicator);

    if (input_ref.have_potentials()) have_reference_values[0] = 1;
    if (input_ref.have_field()) have_reference_values[1] = 1;
  }
#endif /* SCAFACOS_TEST_WITH_DIPOLES */

  /* FIXME: minalloc and overalloc */
/*  input_particles.particles.realloc();*/
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

static fcs_int absolute_minalloc(fcs_float minalloc, fcs_int decomposition, fcs_int ntotal, int comm_size)
{
  fcs_int n = 0;

  if (minalloc < 0) n = (fcs_int) fcs_ceil(ntotal / comm_size * -minalloc);
  else n = (fcs_int) fcs_ceil(minalloc);

  if (decomposition == DECOMPOSE_RANDOM_EQUAL) n = z_max(n, (fcs_int) fcs_ceil(ntotal / comm_size));

  return n;
}

void Configuration::decompose_particles(bool alloc_potentials, bool alloc_field, fcs_float minalloc, fcs_float overalloc)
{
  bool decompose_separate = (params.decomposition == DECOMPOSE_RANDOM || params.decomposition == DECOMPOSE_RANDOM_EQUAL || params.decomposition == DECOMPOSE_DOMAIN);

  generate_input_particles(decompose_separate?0:minalloc, decompose_separate?0:overalloc);
  MPI_Bcast(have_reference_values, 2, FCS_MPI_INT, MASTER_RANK, communicator);

  decomp_comm = communicator;

  decomp_particles.clear();

  fcs_int input_particles_minalloc = absolute_minalloc(minalloc, params.decomposition, input_particles.total_nparticles, comm_size);

  /* wrap particle positions of periodic dimensions */
  fcs_wrap_positions(input_particles.particles.n, input_particles.particles.positions, params.box_a, params.box_b, params.box_c, params.box_origin, params.periodicity);
  
  /* increase particle system in open dimensions to enclose all particles */
  fcs_expand_system_box(input_particles.particles.n, input_particles.particles.positions, params.box_a, params.box_b, params.box_c, params.box_origin, params.periodicity);

  decomp_particles.total_nparticles = input_particles.total_nparticles;

  reference_potentials = 0;
  reference_field = 0;

#if SCAFACOS_TEST_WITH_DIPOLES
/*  fcs_int dipole_input_particles_minalloc = absolute_minalloc(minalloc, params.decomposition, input_particles.dipole_total_nparticles, comm_size);*/

  /* wrap particle positions of periodic dimensions */
  fcs_wrap_positions(input_particles.dipole_particles.n, input_particles.dipole_particles.positions, params.box_a, params.box_b, params.box_c, params.box_origin, params.periodicity);
  
  /* increase particle system in open dimensions to enclose all particles */
  fcs_expand_system_box(input_particles.dipole_particles.n, input_particles.dipole_particles.positions, params.box_a, params.box_b, params.box_c, params.box_origin, params.periodicity);

  decomp_particles.dipole_total_nparticles = input_particles.dipole_total_nparticles;

  dipole_reference_potentials = 0;
  dipole_reference_field = 0;
#endif /* SCAFACOS_TEST_WITH_DIPOLES */

  fcs_gridsort_resort_t gridsort_resort;

  switch (params.decomposition)
  {
    case DECOMPOSE_ALL_ON_MASTER:
    case DECOMPOSE_ALMOST_ALL_ON_MASTER:
      INFO_MASTER(cout << "Decomposing system (" << ((params.decomposition == DECOMPOSE_ALMOST_ALL_ON_MASTER)?"almost-":"") << "all-on-master)..." << endl);
      decomp_particles.particles = input_particles.particles;
      if (have_reference_values[0]) reference_potentials = input_particles.particles.potentials;
      if (have_reference_values[1]) reference_field = input_particles.particles.field;
#if SCAFACOS_TEST_WITH_DIPOLES
      decomp_particles.dipole_particles = input_particles.dipole_particles;
      if (dipole_have_reference_values[0]) dipole_reference_potentials = input_particles.dipole_particles.potentials;
      if (dipole_have_reference_values[1]) dipole_reference_field = input_particles.dipole_particles.field;
#endif /* SCAFACOS_TEST_WITH_DIPOLES */
      if (params.decomposition == DECOMPOSE_ALMOST_ALL_ON_MASTER) almost_master_particles();
      break;
    case DECOMPOSE_ATOMISTIC:
      INFO_MASTER(cout << "Decomposing system (atomistic)..." << endl);
      decomp_particles.particles = input_particles.particles;
      if (have_reference_values[0]) reference_potentials = input_particles.particles.potentials;
      if (have_reference_values[1]) reference_field = input_particles.particles.field;
#if SCAFACOS_TEST_WITH_DIPOLES
      decomp_particles.dipole_particles = input_particles.dipole_particles;
      if (dipole_have_reference_values[0]) dipole_reference_potentials = input_particles.dipole_particles.potentials;
      if (dipole_have_reference_values[1]) dipole_reference_field = input_particles.dipole_particles.field;
#endif /* SCAFACOS_TEST_WITH_DIPOLES */
      break;
    case DECOMPOSE_RANDOM:
    case DECOMPOSE_RANDOM_EQUAL:
      INFO_MASTER(cout << "Decomposing system (random" << ((params.decomposition == DECOMPOSE_RANDOM_EQUAL)?"-equal":"") << ")..." << endl);
      fcs_gridsort_create(&gridsort);
      fcs_gridsort_set_particles(&gridsort, input_particles.particles.n, input_particles.particles.n, input_particles.particles.positions, input_particles.particles.props);
#if SCAFACOS_TEST_WITH_DIPOLES
      fcs_gridsort_set_dipole_particles(&gridsort, input_particles.dipole_particles.n, input_particles.dipole_particles.n, input_particles.dipole_particles.positions, input_particles.dipole_particles.props);
#endif /* SCAFACOS_TEST_WITH_DIPOLES */
      fcs_gridsort_set_minalloc(&gridsort, input_particles_minalloc);
      fcs_gridsort_set_overalloc(&gridsort, overalloc);
      fcs_gridsort_sort_random(&gridsort, communicator);
      fcs_gridsort_get_sorted_particles(&gridsort, &decomp_particles.particles.n, &decomp_particles.particles.max_n, &decomp_particles.particles.positions, &decomp_particles.particles.props, NULL);
#if SCAFACOS_TEST_WITH_DIPOLES
      fcs_gridsort_get_sorted_dipole_particles(&gridsort, &decomp_particles.dipole_particles.n, &decomp_particles.dipole_particles.max_n, &decomp_particles.dipole_particles.positions, &decomp_particles.dipole_particles.props, NULL);
#endif /* SCAFACOS_TEST_WITH_DIPOLES */
      fcs_gridsort_prepare_resort(&gridsort, communicator);

      fcs_gridsort_resort_create(&gridsort_resort, &gridsort, communicator);
      if (have_reference_values[0])
      {
        reference_potentials = new fcs_float[decomp_particles.particles.max_n * decomp_particles.particles.POTENTIAL_SIZE];
        fcs_gridsort_resort_floats(gridsort_resort, input_particles.particles.potentials, reference_potentials, decomp_particles.particles.POTENTIAL_SIZE, communicator);
      }
      if (have_reference_values[1])
      {
        reference_field = new fcs_float[decomp_particles.particles.max_n * decomp_particles.particles.FIELD_SIZE];
        fcs_gridsort_resort_floats(gridsort_resort, input_particles.particles.field, reference_field, decomp_particles.particles.FIELD_SIZE, communicator);
      }
#if SCAFACOS_TEST_WITH_DIPOLES
      if (dipole_have_reference_values[0])
      {
        dipole_reference_potentials = new fcs_float[decomp_particles.dipole_particles.max_n * decomp_particles.dipole_particles.POTENTIAL_SIZE];
        fcs_gridsort_resort_floats(gridsort_resort, input_particles.dipole_particles.potentials, dipole_reference_potentials, decomp_particles.dipole_particles.POTENTIAL_SIZE, communicator);
      }
      if (dipole_have_reference_values[1])
      {
        dipole_reference_field = new fcs_float[decomp_particles.dipole_particles.max_n * decomp_particles.dipole_particles.FIELD_SIZE];
        fcs_gridsort_resort_floats(gridsort_resort, input_particles.dipole_particles.field, dipole_reference_field, decomp_particles.dipole_particles.FIELD_SIZE, communicator);
      }
#endif /* SCAFACOS_TEST_WITH_DIPOLES */
      fcs_gridsort_resort_destroy(&gridsort_resort);
      if (params.decomposition == DECOMPOSE_RANDOM_EQUAL) equalize_particles();
      break;
    case DECOMPOSE_DOMAIN:
      INFO_MASTER(cout << "Decomposing system (domain)..." << endl);
      create_cart_comm();
      decomp_comm = cart_comm;

      fcs_gridsort_create(&gridsort);
      fcs_gridsort_set_system(&gridsort, params.box_origin, params.box_a, params.box_b, params.box_c, params.periodicity);
      fcs_gridsort_set_particles(&gridsort, input_particles.particles.n, input_particles.particles.n, input_particles.particles.positions, input_particles.particles.props);
#if SCAFACOS_TEST_WITH_DIPOLES
      fcs_gridsort_set_dipole_particles(&gridsort, input_particles.dipole_particles.n, input_particles.dipole_particles.n, input_particles.dipole_particles.positions, input_particles.dipole_particles.props);
#endif /* SCAFACOS_TEST_WITH_DIPOLES */
      fcs_gridsort_set_minalloc(&gridsort, input_particles_minalloc);
      fcs_gridsort_set_overalloc(&gridsort, overalloc);
      fcs_gridsort_sort_forward(&gridsort, 0.0, cart_comm);
      fcs_gridsort_get_sorted_particles(&gridsort, &decomp_particles.particles.n, &decomp_particles.particles.max_n, &decomp_particles.particles.positions, &decomp_particles.particles.props, NULL);
#if SCAFACOS_TEST_WITH_DIPOLES
      fcs_gridsort_get_sorted_dipole_particles(&gridsort, &decomp_particles.dipole_particles.n, &decomp_particles.dipole_particles.max_n, &decomp_particles.dipole_particles.positions, &decomp_particles.dipole_particles.props, NULL);
#endif /* SCAFACOS_TEST_WITH_DIPOLES */
      fcs_gridsort_prepare_resort(&gridsort, communicator);

      fcs_gridsort_resort_create(&gridsort_resort, &gridsort, communicator);
      if (have_reference_values[0])
      {
        reference_potentials = new fcs_float[decomp_particles.particles.max_n * decomp_particles.particles.POTENTIAL_SIZE];
        fcs_gridsort_resort_floats(gridsort_resort, input_particles.particles.potentials, reference_potentials, decomp_particles.particles.POTENTIAL_SIZE, communicator);
      }
      if (have_reference_values[1])
      {
        reference_field = new fcs_float[decomp_particles.particles.max_n * decomp_particles.particles.FIELD_SIZE];
        fcs_gridsort_resort_floats(gridsort_resort, input_particles.particles.field, reference_field, decomp_particles.particles.FIELD_SIZE, communicator);
      }
#if SCAFACOS_TEST_WITH_DIPOLES
      if (dipole_have_reference_values[0])
      {
        dipole_reference_potentials = new fcs_float[decomp_particles.dipole_particles.max_n * decomp_particles.dipole_particles.POTENTIAL_SIZE];
        fcs_gridsort_resort_floats(gridsort_resort, input_particles.dipole_particles.potentials, dipole_reference_potentials, decomp_particles.dipole_particles.POTENTIAL_SIZE, communicator);
      }
      if (dipole_have_reference_values[1])
      {
        dipole_reference_field = new fcs_float[decomp_particles.dipole_particles.max_n * decomp_particles.dipole_particles.FIELD_SIZE];
        fcs_gridsort_resort_floats(gridsort_resort, input_particles.dipole_particles.field, dipole_reference_field, decomp_particles.dipole_particles.FIELD_SIZE, communicator);
      }
#endif /* SCAFACOS_TEST_WITH_DIPOLES */
      fcs_gridsort_resort_destroy(&gridsort_resort);
      break;
    default:
      break;
  }

  if (decomp_particles.particles.max_n > 0)
  {
    if (alloc_potentials)
    {
      decomp_particles.particles.potentials = new fcs_float[decomp_particles.particles.max_n * particle_data_t::POTENTIAL_SIZE];
      values_set<fcs_float, particle_data_t::POTENTIAL_SIZE>(decomp_particles.particles.potentials, NAN, decomp_particles.particles.n);

    } else decomp_particles.particles.potentials = NULL;

    if (alloc_field)
    {
      decomp_particles.particles.field = new fcs_float[decomp_particles.particles.max_n * particle_data_t::FIELD_SIZE];
      values_set<fcs_float, particle_data_t::FIELD_SIZE>(decomp_particles.particles.field, NAN, decomp_particles.particles.n);

    } else decomp_particles.particles.field = NULL;

  } else
  {
    decomp_particles.particles.potentials = (fcs_float *) 1;
    decomp_particles.particles.field = (fcs_float *) 1;
  }
}


void Configuration::almost_master_particles()
{
  /* FIXME: only on particles (not dipoles?) */

  if (comm_rank == 0 && decomp_particles.particles.n < comm_size)
  {
    cout << "ERROR: not enough particles for almost-all-on-master setup " << endl;
    MPI_Abort(communicator, 1);
  }

  fcs_float *new_positions = new fcs_float[particle_data_t::POSITION_SIZE];
  fcs_float *new_charges = new fcs_float[particle_data_t::PROP_SIZE];
  fcs_float *new_potentials = new fcs_float[particle_data_t::POTENTIAL_SIZE];
  fcs_float *new_field = new fcs_float[particle_data_t::FIELD_SIZE];

  fcs_int offset = decomp_particles.particles.n - comm_size;

  MPI_Scatter(decomp_particles.particles.positions_at(offset), particle_data_t::POSITION_SIZE, FCS_MPI_FLOAT, new_positions, particle_data_t::POSITION_SIZE, FCS_MPI_FLOAT, MASTER_RANK, communicator);
  MPI_Scatter(decomp_particles.particles.props_at(offset), particle_data_t::PROP_SIZE, FCS_MPI_FLOAT, new_charges, particle_data_t::PROP_SIZE, FCS_MPI_FLOAT, MASTER_RANK, communicator);
  if (have_reference_values[0])
    MPI_Scatter(reference_potentials + (offset * particle_data_t::POTENTIAL_SIZE), particle_data_t::POTENTIAL_SIZE, FCS_MPI_FLOAT, new_potentials, particle_data_t::POTENTIAL_SIZE, FCS_MPI_FLOAT, MASTER_RANK, communicator);
  if (have_reference_values[1])
    MPI_Scatter(reference_field + (offset * particle_data_t::FIELD_SIZE), particle_data_t::FIELD_SIZE, FCS_MPI_FLOAT, new_field, particle_data_t::FIELD_SIZE, FCS_MPI_FLOAT, MASTER_RANK, communicator);

  if (comm_rank == 0)
  {
    decomp_particles.particles.n -= comm_size - 1;
    delete[] new_positions;
    delete[] new_charges;
    delete[] new_potentials;
    delete[] new_field;

  } else
  {
    decomp_particles.particles.n = 1;
    decomp_particles.particles.max_n = 1;
    decomp_particles.particles.positions = new_positions;
    decomp_particles.particles.props = new_charges;
    if (have_reference_values[0]) reference_potentials = new_potentials; else delete[] new_potentials;
    if (have_reference_values[1]) reference_field = new_field; else delete[] new_field;
  }
}


void Configuration::equalize_particles()
{
  fcs_int nparticles, my_diff, all_diffs[comm_size], n;
  int sr, rr;
  

  nparticles = (decomp_particles.total_nparticles / comm_size) + ((decomp_particles.total_nparticles % comm_size < comm_rank)?1:0);
  nparticles = z_min(nparticles, decomp_particles.particles.max_n);
  my_diff = decomp_particles.particles.n - nparticles;
  
  MPI_Allgather(&my_diff, 1, FCS_MPI_INT, all_diffs, 1, FCS_MPI_INT, communicator);
  
  for (sr = 0, rr = 0; rr < comm_size; ++rr)
  {
    while (all_diffs[rr] < 0)
    {
      while (sr < comm_size && all_diffs[sr] <= 0) ++sr;

      if (sr >= comm_size) break;
    
      n = z_min(-all_diffs[rr], all_diffs[sr]);

/*      if (comm_rank == 0) printf("send %" FCS_LMOD_INT "d from %d to %d\n", n, sr, rr);*/

      if (comm_rank == sr)
      {
        fcs_int offset = decomp_particles.particles.n - n;

        MPI_Send(decomp_particles.particles.positions_at(offset), particle_data_t::POSITION_SIZE * n, FCS_MPI_FLOAT, rr, 0, communicator);
        MPI_Send(decomp_particles.particles.props_at(offset), particle_data_t::PROP_SIZE * n, FCS_MPI_FLOAT, rr, 0, communicator);
        MPI_Send(reference_potentials + (offset * particle_data_t::POTENTIAL_SIZE), particle_data_t::POTENTIAL_SIZE * n, FCS_MPI_FLOAT, rr, 0, communicator);
        MPI_Send(reference_field + (offset * particle_data_t::FIELD_SIZE), particle_data_t::FIELD_SIZE * n, FCS_MPI_FLOAT, rr, 0, communicator);
        decomp_particles.particles.n -= n;

      } else if (comm_rank == rr)
      {
        fcs_int offset = decomp_particles.particles.n;

        MPI_Recv(decomp_particles.particles.positions_at(offset), particle_data_t::POSITION_SIZE * n, FCS_MPI_FLOAT, sr, 0, communicator, MPI_STATUS_IGNORE);
        MPI_Recv(decomp_particles.particles.props_at(offset), particle_data_t::PROP_SIZE * n, FCS_MPI_FLOAT, sr, 0, communicator, MPI_STATUS_IGNORE);
        MPI_Recv(reference_potentials + (offset * particle_data_t::POTENTIAL_SIZE), particle_data_t::POTENTIAL_SIZE * n, FCS_MPI_FLOAT, sr, 0, communicator, MPI_STATUS_IGNORE);
        MPI_Recv(reference_field + (offset * particle_data_t::FIELD_SIZE), particle_data_t::FIELD_SIZE * n, FCS_MPI_FLOAT, sr, 0, communicator, MPI_STATUS_IGNORE);
        decomp_particles.particles.n += n;
      }

      all_diffs[rr] += n;
      all_diffs[sr] -= n;
    }
  }
}

bool Configuration::compute_errors(errors_t *e)
{
  ::compute_errors(e, decomp_particles.particles.n, decomp_particles.particles.positions, decomp_particles.particles.props,
    have_reference_values[0]?reference_potentials:NULL, have_reference_values[1]?reference_field:NULL,
    have_result_values[0]?decomp_particles.particles.potentials:NULL, have_result_values[1]?decomp_particles.particles.field:NULL,
    field_correction, energy_correction,
    decomp_comm);

  return true;
}

void Configuration::free_decomp_particles(bool quiet)
{
  switch (params.decomposition)
  {
    case DECOMPOSE_ALL_ON_MASTER:
    case DECOMPOSE_ALMOST_ALL_ON_MASTER:
      if (!quiet) INFO_MASTER(cout << "Freeing data (" << ((params.decomposition == DECOMPOSE_ALMOST_ALL_ON_MASTER)?"almost-":"") << "all-on-master)..." << endl);
      break;
    case DECOMPOSE_ATOMISTIC:
      if (!quiet) INFO_MASTER(cout << "Freeing data (atomistic)..." << endl);
      break;
    case DECOMPOSE_RANDOM:
    case DECOMPOSE_RANDOM_EQUAL:
      if (!quiet) INFO_MASTER(cout << "Freeing data (rand" << ((params.decomposition == DECOMPOSE_RANDOM)?"om":"eq") << ")..." << endl);
      fcs_gridsort_free(&gridsort);
      fcs_gridsort_destroy(&gridsort);
      delete[] reference_potentials;
      delete[] reference_field;
      break;
    case DECOMPOSE_DOMAIN:
      if (!quiet) INFO_MASTER(cout << "Freeing data (domain)..." << endl);
      fcs_gridsort_free(&gridsort);
      fcs_gridsort_destroy(&gridsort);
      destroy_cart_comm();
      delete[] reference_potentials;
      delete[] reference_field;
      break;
    default:
      break;
  }

  if (decomp_particles.particles.max_n > 0)
  {
    delete[] decomp_particles.particles.potentials;
    delete[] decomp_particles.particles.field;
  }

  decomp_particles.clear();

  decomp_comm = MPI_COMM_NULL;

  reference_potentials = 0;
  reference_field = 0;
}

void Configuration::determine_total_duplication()
{
  total_duplications[0] = ((params.periodicity[0])?params.periodic_duplications[0]:1);
  total_duplications[1] = ((params.periodicity[1])?params.periodic_duplications[1]:1);
  total_duplications[2] = ((params.periodicity[2])?params.periodic_duplications[2]:1);

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
  static const size_t READSIZE = 8192;
  ZLIB_IFELSE(gzFile, FILE*) inputfile;
  size_t inputsize, dataread;
  char* inputdata;
  char* datap;
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
      istringstream is(string(attr->value()));
      is.exceptions(istream::failbit | istream::badbit);
      is >> this->error_potential;
    } else if (aname == "error_field") {
      istringstream is(string(attr->value()));
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

    Configuration *config = new Configuration();
    this->configurations.push_back(config);

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

    os.str(string(""));
    os << this->error_potential;
    s = doc.allocate_string(os.str().c_str());
    attr = doc.allocate_attribute("error_potential", s);
    system_node->append_attribute(attr);

    os.str(string(""));
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
