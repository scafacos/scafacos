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

#include "fcs.h"

#include "common.hpp"
#include "common/fcs-common/FCSCommon.h"

using namespace std;


bool parse_value(string s, fcs_int &r) {
  return parse_value<fcs_int>(s, r);
}

bool parse_value(string s, fcs_int &r, char &c)
{
  return parse_value<fcs_int>(s, r, &c);
}

bool parse_value(string s, fcs_float &r) {
  return parse_value<fcs_float>(s, r);
}

bool parse_value(string s, fcs_float &r, char &c) {
  return parse_value<fcs_float>(s, r, &c);
}

fcs_int parse_sequence(string s, fcs_int nmax, fcs_int *rv, char *cv) {
  return parse_sequence<fcs_int>(s, nmax, rv, cv);
}

fcs_int parse_sequence(string s, fcs_int nmax, fcs_float *rv, char *cv) {
  return parse_sequence<fcs_float>(s, nmax, rv, cv);
}


void make_equal_counts_and_displs(fcs_int total_count, fcs_int ncounts, int *counts, int *displs, int *counts3, int *displs3) {
  fcs_int i;

  displs[0] = 0;
  if (displs3) displs3[0] = 0;
  for (i = 0; i < ncounts; ++i) {
    counts[i] = (int) ((fcs_float) total_count / (fcs_float) (ncounts - i) + 0.5);
    if (displs3) counts3[i] = counts[i] * 3;
    if (i > 0) {
      displs[i] = displs[i - 1] + counts[i - 1];
      if (displs3) displs3[i] = 3 * displs[i];
    }
    total_count -= counts[i];
  }
}

void compute_errors(errors_t *e, fcs_int nparticles, 
		    fcs_float *positions, fcs_float *charges, 
		    fcs_float *reference_potentials, fcs_float *reference_field, 
		    fcs_float *result_potentials, fcs_float *result_field, 
		    MPI_Comm comm)
{
  // Compare resulting potential and field values with reference values

  const int idx_potential_error_sqr = 0;
  const int idx_field_error_sqr = 1;
  const int idx_potential_sqr = 2;
  const int idx_field_sqr = 3;
  const int idx_energy_sum = 4;
  const int idx_ref_energy_sum = 5;

  fcs_float local_sum[6], global_sum[6];

  struct {
    fcs_float val;
    int pid;
  } local_max[2], global_max[2];

  fcs_int local_nparticles[3], global_nparticles[3];

  DEBUG(cout << comm_rank << ": compute local error of " << nparticles << " particles" << endl);
  
  e->have_potential_errors = (reference_potentials != NULL && result_potentials != NULL);
  e->have_field_errors = (reference_field != NULL && result_field != NULL);

  // local sums of the squared potential/field
  local_sum[idx_potential_sqr] = 0.0;
  local_sum[idx_field_sqr] = 0.0;

  // local sums of the squared potential/field error
  local_sum[idx_potential_error_sqr] = 0.0;
  local_sum[idx_field_error_sqr] = 0.0;

  // local sums of the energy
  local_sum[idx_energy_sum] = 0.0;
  local_sum[idx_ref_energy_sum] = 0.0;

  // maximum of the squared potential/field error and its particle id
  local_max[idx_potential_error_sqr].val = 0.0;
  local_max[idx_potential_error_sqr].pid = 0;
  local_max[idx_field_error_sqr].val = 0.0;
  local_max[idx_field_error_sqr].pid = 0;

  fcs_int pid_offset = 0;

  MPI_Exscan(&nparticles, &pid_offset, 1, FCS_MPI_INT, MPI_SUM, comm);

  local_nparticles[0] = nparticles;
  local_nparticles[1] = local_nparticles[2] = 0;

  if (e->have_potential_errors)
  for (fcs_int pid = 0; pid < nparticles; pid++) {
    const fcs_float &res_potential = result_potentials[pid];
    const fcs_float &ref_potential = reference_potentials[pid];

    if (isnan(res_potential) || isnan(ref_potential)) continue;

    // compute potential sum
    local_sum[idx_potential_sqr] += res_potential*res_potential;

    // compute potential error sum
    fcs_float potential_error = res_potential - ref_potential;
    fcs_float potential_error_sqr = potential_error*potential_error;
    local_sum[idx_potential_error_sqr] += potential_error_sqr;

    // compute energy sum
    local_sum[idx_energy_sum] += 0.5 * charges[pid] * res_potential;
    local_sum[idx_ref_energy_sum] += 0.5 * charges[pid] * ref_potential;

    // compute potential error max
    if (potential_error_sqr > local_max[idx_potential_error_sqr].val) {
      local_max[idx_potential_error_sqr].val = potential_error_sqr;
      local_max[idx_potential_error_sqr].pid = pid_offset + pid;
    }

    // count valid poential values
    ++local_nparticles[1];
  }

  if (e->have_field_errors)
  for (fcs_int pid = 0; pid < nparticles; pid++) {
    const fcs_float *res_field = &result_field[3*pid];
    const fcs_float *ref_field = &reference_field[3*pid];

    if (isnan(res_field[0]) || isnan(res_field[1]) || isnan(res_field[2]) || isnan(ref_field[0]) || isnan(ref_field[1]) || isnan(ref_field[2])) continue;

    // compute field sum
    local_sum[idx_field_sqr] += res_field[0]*res_field[0] + res_field[1]*res_field[1] + res_field[2]*res_field[2];

    // compute field error sum
    fcs_float dx = res_field[0]-ref_field[0];
    fcs_float dy = res_field[1]-ref_field[1];
    fcs_float dz = res_field[2]-ref_field[2];
    fcs_float field_error_sqr = dx*dx + dy*dy + dz*dz;
    local_sum[idx_field_error_sqr] += field_error_sqr;

    // compute field error max
    if (field_error_sqr > local_max[idx_field_error_sqr].val) {
      local_max[idx_field_error_sqr].val = field_error_sqr;
      local_max[idx_field_error_sqr].pid = pid_offset + pid;
    }

    // count valid field values
    ++local_nparticles[2];
  }

#ifdef FCS_ENABLE_DEBUG
#define NUM_PARTICLES_TO_PRINT 20
  if (comm_rank==MASTER_RANK && (e->have_potential_errors || e->have_field_errors)) {
    for (fcs_int pid = 0; 
	 pid < (nparticles>NUM_PARTICLES_TO_PRINT ? NUM_PARTICLES_TO_PRINT : nparticles); 
	 pid++) {
      const fcs_float *pos = &positions[3*pid];
      const fcs_float &q = charges[pid];
      
      cout << pid << ": "
	   << setw(15) << "position ="
	   << setw(15) << pos[0] << setw(15) << pos[1] << setw(15) << pos[2] << endl;
      cout << pid << ": " 
	   << setw(15) << "q ="
	   << setw(15) << q  << endl;
      
      if (e->have_potential_errors)
	{
	  const fcs_float &res_potential = result_potentials[pid];
	  const fcs_float &ref_potential = reference_potentials[pid];
	  cout << pid << ": " 
	       << setw(15) << "potential ="
	       << setw(15) << res_potential << endl;
	  cout << pid << ": " 
	       << setw(15) << "ref_pot ="
	       << setw(15) << ref_potential << endl;
	}
      
      if (e->have_field_errors)
	{
	  const fcs_float *res_field = &result_field[3*pid];
	  const fcs_float *ref_field = &reference_field[3*pid];
	  cout << pid << ": " 
	       << setw(15) << "field ="
	       << setw(15) << res_field[0] 
	       << setw(15) << res_field[1] 
	       << setw(15) << res_field[2] << endl;
	  cout << pid << ": " 
	       << setw(15) << "ref_field ="
	       << setw(15) << ref_field[0] 
	       << setw(15) << ref_field[1] 
	       << setw(15) << ref_field[2] << endl;
	}
    }
    if (nparticles > NUM_PARTICLES_TO_PRINT)
      cout << "  (" << (nparticles-NUM_PARTICLES_TO_PRINT) << " more on #" 
	   << comm_rank << ")" << endl;
  }
#endif

  MPI_Allreduce(local_nparticles, global_nparticles, 3, FCS_MPI_INT, MPI_SUM, comm);

  MPI_Allreduce(local_sum, global_sum, 6, FCS_MPI_FLOAT, MPI_SUM, comm);
  MPI_Allreduce(local_max, global_max, 2,
#if defined(FCS_FLOAT_IS_FLOAT)
    MPI_FLOAT_INT,
#elif defined(FCS_FLOAT_IS_LONG_DOUBLE)
    MPI_LONG_DOUBLE_INT,
#else
    MPI_DOUBLE_INT,
#endif
    MPI_MAXLOC, comm);

  // sums of the squared potentials/field error
  e->sum_potential_error_sqr = global_sum[idx_potential_error_sqr];
  e->sum_field_error_sqr = global_sum[idx_field_error_sqr];

  // maximum of the squared potentials/field error
  e->max_potential_error_sqr = global_max[idx_potential_error_sqr].val;
  e->max_field_error_sqr = global_max[idx_field_error_sqr].val;

  // pid with the maximal potentials/field error
  e->max_potential_error_pid = global_max[idx_potential_error_sqr].pid;
  e->max_field_error_pid = global_max[idx_field_error_sqr].pid;

  // rms potential/field
  fcs_float rms_potential = (global_nparticles[1])?sqrt(global_sum[idx_potential_sqr] / (fcs_float) global_nparticles[1]):NAN;
  fcs_float rms_field = (global_nparticles[2])?sqrt(global_sum[idx_field_sqr] / (fcs_float) global_nparticles[2]):NAN;

  // absolute rms errors
  e->abs_rms_potential_error = (global_nparticles[1])?sqrt(e->sum_potential_error_sqr / (fcs_float) global_nparticles[1]):NAN;
  e->abs_rms_field_error = (global_nparticles[2])?sqrt(e->sum_field_error_sqr / (fcs_float) global_nparticles[2]):NAN;

  // relative rms errors
  e->rel_rms_potential_error = e->abs_rms_potential_error / rms_potential;
  e->rel_rms_field_error = e->abs_rms_field_error / rms_field;

  // absolute maximum errors
  e->abs_max_potential_error = sqrt(e->max_potential_error_sqr);
  e->abs_max_field_error = sqrt(e->max_field_error_sqr);

  // relative maximum errors
  e->rel_max_potential_error = e->abs_max_potential_error / rms_potential;
  e->rel_max_field_error = e->abs_max_field_error / rms_field;

  // total energy
  e->total_energy = global_sum[idx_energy_sum];
  e->total_energy_ref = global_sum[idx_ref_energy_sum];

  // absolute total energy error
  e->abs_energy_error = fabs(e->total_energy - e->total_energy_ref);
  if (fcs_float_is_zero(e->total_energy_ref))
    e->rel_energy_error = e->abs_energy_error;
  else
    e->rel_energy_error = fabs(e->abs_energy_error / e->total_energy_ref);

  e->total_nparticles = global_nparticles[0];
  e->valid_potentials_errors = global_nparticles[1];
  e->valid_field_errors = global_nparticles[2];
}

void print_errors(errors_t *e, const char *prefix)
{
  // OUTPUT
  if (!e->have_potential_errors && !e->have_field_errors) {
    cout << prefix << "NO ERRORS: No reference or result data available!" << endl;
  } else {
    cout << prefix << "ABSOLUTE ERRORS (from " << e->valid_potentials_errors << " of " << e->total_nparticles << " particles)" << endl;
    if (e->have_field_errors) {
      cout << prefix << setw(30) << "abs_rms_field_error = " 
	   << setw(15) << scientific << e->abs_rms_field_error << endl;
      cout << prefix << setw(30) << "abs_max_field_error = " 
	   << setw(15) << scientific << e->abs_max_field_error << endl;
    }

    if (e->have_potential_errors) {
      cout << prefix << setw(30) << "abs_rms_potential_error = " 
	   << setw(15) << scientific << e->abs_rms_potential_error << endl;
      cout << prefix << setw(30) << "abs_max_potential_error = " 
	   << setw(15) << scientific << e->abs_max_potential_error << endl;
      cout << prefix << setw(30) << "abs_energy_error = "
	   << setw(15) << scientific << e->abs_energy_error << endl;   
    }
  
    cout << endl;
    cout << prefix << "RELATIVE ERRORS (from " << e->valid_potentials_errors << " of " << e->total_nparticles << " particles)" << endl;
    if (e->have_field_errors) {
      cout << prefix << setw(30) << "rel_rms_field_error = " 
	   << setw(15) << scientific << e->rel_rms_field_error << endl;
      cout << prefix << setw(30) << "rel_max_field_error = " 
	   << setw(15) << scientific << e->rel_max_field_error << endl;
    }

    if (e->have_potential_errors) {
      cout << prefix << setw(30) << "rel_rms_potential_error = " 
	   << setw(15) << scientific << e->rel_rms_potential_error << endl;
      cout << prefix << setw(30) << "rel_max_potential_error = " 
	   << setw(15) << scientific << e->rel_max_potential_error << endl;
      cout << prefix << setw(30) << "rel_energy_error = "
	   << setw(15) << scientific << e->rel_energy_error << endl;
    }

    cout << endl;
  }

  if (e->have_potential_errors)
    cout << prefix << setw(30) << "total_energy_ref = " 
	 << setw(15) << scientific << e->total_energy_ref << endl;
  cout << prefix << setw(30) << "total_energy = " 
       << setw(15) << scientific << e->total_energy << endl;
}
