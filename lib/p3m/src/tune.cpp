/*
  Copyright (C) 2011,2012,2013 Olaf Lenz

  This file is part of ScaFaCoS.

  ScaFaCoS is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  ScaFaCoS is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>. 
 */

#include "utils.hpp"
#include "Communication.hpp"
#include "Solver.hpp"
#include <cstdlib>
#include <cstdio>
#include <stdexcept>
#include <list>

#ifdef P3M_ENABLE_DEBUG
#undef P3M_DEBUG
#define P3M_DEBUG(cmd) if (d->comm.onMaster()) { cmd; }
#endif

namespace P3M {
/***************************************************/
/* TYPES AND CONSTANTS */
/***************************************************/
/** Good mesh sizes for fftw
 */
static const p3m_int good_gridsize[] =
{0, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32,
		36, 40, 42, 44, 48, 50, 52, 54, 56, 60, 64,
		66, 70, 72, 78, 80, 84, 88, 90, 96, 98, 100, 104, 108, 110, 112, 120, 126, 128,
		130, 132, 140, 144, 150, 154, 156, 160, 162, 168, 176, 180, 182, 192,
		196, 198, 200, 208, 210, 216, 220, 224, 234, 240, 242, 250, 252, 256,
		260, 264, 270, 280, 288, 294,
		300, 360, 400, 480, 512,
		576, 648, 729, 768, 864, 972, 1024,
		1152, 1296, 1458, 1536,
		1728, 1944, 2048,
		2187, 2304, 2592,
		2916, 3072,
		3456,
		3888};

/** Steps to do when trying to determine smallest possible grid size. */
static const p3m_int step_good_gridsize[] =
{ 0, 15, 26, 44, 58, 72, 78, 83, 90, 94, 98, 101, 103, 104 } ;
static const p3m_int num_steps_good_gridsize = 14;

// fixme
struct tune_params {
	/** charge assignment order ([0,P3M_MAX_CAO]). */
	p3m_int cao;
	/** number of grid points per coordinate direction (>0). */
	p3m_int grid[3];
	/** Ewald splitting parameter */
	p3m_float alpha;

	/** Errors */
	p3m_float rs_error, ks_error, error;
	/** Timings */
	p3m_float timing, timing_near, timing_far;
};
typedef std::list<tune_params*> tune_params_l;

/***************************************************/
/* FORWARD DECLARATIONS OF INTERNAL FUNCTIONS */
/***************************************************/
tune_params*
tune_far(Solver *d,
		p3m_int num_particles, p3m_float *positions, p3m_float *charges,
		p3m_float r_cut);

void
tune_far(Solver *d,
		p3m_int num_particles, p3m_float *positions, p3m_float *charges);

tune_params*
tune_alpha_cao_grid(Solver *d,
		p3m_int num_particles, p3m_float *positions, p3m_float *charges,
		p3m_float r_cut);

tune_params*
tune_cao_grid(Solver *d,
		p3m_int num_particles, p3m_float *positions, p3m_float *charges,
		p3m_float r_cut, p3m_float alpha);

tune_params*
tune_grid(Solver *d, p3m_int num_particles,
		p3m_float *positions, p3m_float *charges,
		tune_params_l &params_to_try);

tune_params*
time_params(Solver *d,
		p3m_int num_particles, p3m_float *positions, p3m_float *charges,
		tune_params_l &params_to_try);

void
count_charges(Solver *d,
		p3m_int num_particles, p3m_float *charges);

/***************************************************/
/* IMPLEMENTATION */
/***************************************************/
void
tune(Solver *d, p3m_int num_particles,
		p3m_float *positions, p3m_float *charges) {

	/* Prepare the communicator before tuning */
	d->comm.prepare(d->box_l);

	Communication &comm = d->comm;

	/* Count the charges */
	p3m_float sum_q2_before = d->sum_q2;
	count_charges(d, num_particles, charges);

	if (!d->comm.onMaster()) {
		int howoften;
		P3M_DEBUG_LOCAL(printf("  %d: How often to run tune_far?\n", d->comm.rank));
		MPI_Bcast(&howoften, 1, MPI_INT, Communication::MPI_MASTER, d->comm.mpicomm);
		P3M_DEBUG_LOCAL(printf("  %d: Running tune_far %d times.\n", d->comm.rank, howoften));
		for (; howoften > 0; howoften--)
			tune_far(d, num_particles, positions, charges);

		d->needs_retune = false;
		return;
	}

	/* Retune if the number of charges has changed */
	if (!(float_is_equal(sum_q2_before, d->sum_q2))) {
		P3M_INFO(printf( "  Number of charges changed, retuning is needed.\n"));
		d->needs_retune = true;
	}

	/* Do not retune if there are no charges */
	if (float_is_zero(d->sum_q2)) {
		P3M_INFO(printf( "  No charges in the system.\n"));
		d->needs_retune = false;
	}

	/* Exit if retuning is unnecessary */
	if (!d->needs_retune) {
		P3M_INFO(printf( "  Retuning is not required.\n"));
		int howoften = 0;
		MPI_Bcast(&howoften, 1, MPI_INT, Communication::MPI_MASTER, d->comm.mpicomm);
		P3M_INFO(printf( "tune() finished.\n"));
		return;
	}

	P3M_INFO(printf( "  Retuning is required.\n"));

	if (!d->tune_r_cut) {
		P3M_INFO(printf( "    r_cut=" FFLOAT " (fixed)\n", d->r_cut));
		int howoften = 1;
		MPI_Bcast(&howoften, 1, MPI_INT, Communication::MPI_MASTER, d->comm.mpicomm);
		tune_far(d, num_particles, positions, charges, d->r_cut);
	} else {
		int howoften = 2;
		MPI_Bcast(&howoften, 1, MPI_INT, Communication::MPI_MASTER, d->comm.mpicomm);

		/* compute the average distance between two charges  */
		p3m_float avg_dist =
				pow((d->box_l[0]*d->box_l[1]*d->box_l[2])
						/ d->sum_qpart, 0.33333);

		/* FIRST ESTIMATE */
		/* set the initial r_cut to 3 times the average distance between
         charges */
		p3m_float r_cut = 3.0 * avg_dist;

		/* tune r_cut to half the box length */
		if (0.5*d->box_l[1]-d->skin < r_cut)
			r_cut = 0.5*d->box_l[0] - d->skin;
		if (0.5*d->box_l[1]-d->skin < r_cut)
			r_cut = 0.5*d->box_l[1] - d->skin;
		if (0.5*d->box_l[2]-d->skin < r_cut)
			r_cut = 0.5*d->box_l[2] - d->skin;

		P3M_INFO(printf( "    r_cut=" FFLOAT " (first estimate)\n", r_cut));

		// @todo get near timing
		tune_params *p =
				tune_far(d, num_particles, positions, charges, r_cut);

		/* SECOND ESTIMATE */
		/* use the fact that we know that timing_near scales like r_cut**3
         and timing_far like r_cut**(-3) to get the second estimate */

		p3m_float rel_timing_diff =
				fabs(p->timing_near - p->timing_far) /
				(p->timing_near + p->timing_far);
		P3M_INFO(printf( "    rel_timing_diff=" FFLOAT "\n", rel_timing_diff));

		p3m_float rcut3 = pow(r_cut, 3);
		p3m_float c_near = p->timing_near/rcut3;
		p3m_float c_far = p->timing_far*rcut3;
		p3m_float rcut_new = pow(c_far/c_near, 1./6.);

		r_cut = rcut_new;
		P3M_INFO(printf( "    r_cut=" FFLOAT " (second estimate)\n", r_cut));

		// @todo get near timing
		// second far tuning
		p = tune_far(d, num_particles, positions, charges, r_cut);

		rel_timing_diff =
				fabs(p->timing_near - p->timing_far) /
				(p->timing_near + p->timing_far);
		P3M_INFO(printf("    rel_timing_diff=" FFLOAT "\n", rel_timing_diff));
		P3M_INFO(printf("    Finished tuning.\n"));
	}

	/* mark that the method was retuned */
	d->needs_retune = false;
}

tune_params*
tune_far(Solver *d,
		p3m_int num_particles, p3m_float *positions, p3m_float *charges,
		p3m_float r_cut) {
	/* r_cut is ignored on the slaves */
	/* Distinguish between two types of parameters:
       Input params:
	 * box length
	 * #charges
	 * cutoff
	 * skin (only if comm is cartesian)
	 * error (default=1e-4)
	 * component flags

       Performance related params:
	 * caf interpolation

       Tuned params:
	 * alpha
	 * cao
	 * grid
	 */
	Communication &comm = d->comm;
	P3M_INFO(printf( "tune_far() started...\n"));

	tune_params *final_params = NULL;

	if (!d->comm.onMaster()) {
		tune_broadcast_slave(d, num_particles, positions, charges);
	} else {

		/* check whether the input parameters are sane */
		if (r_cut < 0.0)
			throw std::domain_error("r_cut is negative!");

		if (float_is_zero(r_cut))
			throw std::domain_error("r_cut is too small!");

		/* check whether cutoff is larger than half a box length */
		if ((r_cut > 0.5*d->box_l[0]) ||
				(r_cut > 0.5*d->box_l[1]) ||
				(r_cut > 0.5*d->box_l[2]))
			throw std::domain_error("r_cut is larger than half a system box length.");

		/* check whether cutoff is larger than domain size */

		try {
			d->r_cut = r_cut;
			final_params =
					tune_alpha_cao_grid(d, num_particles,
							positions, charges, r_cut);

			// @todo: move to tune function
			/* Set and broadcast the final parameters. */
			d->r_cut = r_cut;
			d->alpha = final_params->alpha;
			d->grid[0] = final_params->grid[0];
			d->grid[1] = final_params->grid[1];
			d->grid[2] = final_params->grid[2];
			d->cao = final_params->cao;
			tune_broadcast_command(d, CMD_FINISHED);
		} catch (...) {
			P3M_INFO(printf( "  Tuning failed.\n"));
			tune_broadcast_command(d, CMD_FAILED);
			P3M_INFO(printf( "tune_far() finished.\n"));
			throw;
		}
	}

	P3M_INFO(printf( "  Tuning was successful.\n"));

	/* At the end of retuning, prepare the method */
	d->prepare();

	P3M_INFO(printf( "tune_far() finished.\n"));

	return final_params;
}

/** Slave variant of tune_far. */
void
tune_far(Solver *d,
		p3m_int num_particles, p3m_float *positions, p3m_float *charges) {
	if (d->comm.onMaster()) throw std::runtime_error("tune_far without r_cut cannot be called on master node.");
	tune_far(d, num_particles, positions, charges, 0.0);
}


/* Tune alpha */
tune_params*
tune_alpha_cao_grid(Solver *d,
		p3m_int num_particles, p3m_float *positions, p3m_float *charges,
		p3m_float r_cut) {
	Communication &comm = d->comm;
	if (d->tune_alpha) {
		// fixme
		Parameters p;
		p.cao = d->cao;
		p.alpha = d->alpha;
		p.r_cut = d->r_cut;
		p.grid[0] = d->grid[0];
		p.grid[1] = d->grid[1];
		p.grid[2] = d->grid[2];
		d->errorEstimate->compute_alpha(d->tolerance_field,
				p, d->sum_qpart, d->sum_q2, d->box_l);
		d->alpha = p.alpha;
		P3M_INFO(printf("    => alpha=" FFLOAT "\n", d->alpha));
	} else {
		P3M_INFO(printf("    alpha=" FFLOAT " (fixed)\n", d->alpha));
	}

	return tune_cao_grid(d, num_particles, positions, charges, r_cut, d->alpha);
}

/* Tune cao */
tune_params*
tune_cao_grid(Solver *d,
		p3m_int num_particles, p3m_float *positions, p3m_float *charges,
		p3m_float r_cut, p3m_float alpha) {
	// fixme
	Communication &comm = d->comm;
#ifdef P3M_AD
	const p3m_int min_cao = 2;
#else
	const p3m_int min_cao = 1;
#endif

	tune_params_l params_to_try;

	if (d->tune_cao) {
		P3M_INFO(printf("    Testing cao={ "));
		for (p3m_int cao = CAF::max_cao; cao >= min_cao; --cao) {
			tune_params *p = new tune_params;
			p->cao = cao;
			p->alpha = alpha;
			params_to_try.push_back(p);
			P3M_INFO(printf(FINT " ", cao));
		}
		P3M_INFO(printf("}\n"));
	} else {
		tune_params *p = new tune_params;
		p->cao = d->cao;
		p->alpha = alpha;
		P3M_INFO(printf( "    cao=" FINT " (fixed)\n", d->cao));
	}

	return tune_grid(d, num_particles, positions, charges, params_to_try);
}

/* params should have decreasing cao */
tune_params*
tune_grid(Solver *d, p3m_int num_particles,
		p3m_float *positions, p3m_float *charges,
		tune_params_l &params_to_try) {
	// fixme
	Communication &comm = d->comm;
	if (d->tune_grid) {
		p3m_int step_ix = 0;
		// store the minimal grid size seen so far
		p3m_int min_grid1d = good_gridsize[step_good_gridsize[num_steps_good_gridsize-1]];

		for (tune_params_l::iterator p = params_to_try.begin();
				p != params_to_try.end(); ++p) {
			// Find smallest possible grid for this set of parameters
			// Test whether accuracy can be achieved with given parameters
			d->cao = (*p)->cao;

			p3m_int upper_ix;
			// test this step
			P3M_INFO(printf("    Trying to find grid for r_cut=" FFLOAT ", " \
					"alpha=" FFLOAT ", "                          \
					"cao=" FINT "\n",                             \
					d->r_cut, d->alpha, d->cao));
			p3m_float error, rs_error, ks_error;
			do {
				step_ix++;
				if (step_ix >= num_steps_good_gridsize) break;
				upper_ix = step_good_gridsize[step_ix];
				p3m_int grid1d = good_gridsize[upper_ix];
				P3M_DEBUG(printf("      rough grid=" FINT "\n", grid1d));
				// fixme
				Parameters p2;
				p2.grid[0] = grid1d;
				p2.grid[1] = grid1d;
				p2.grid[2] = grid1d;
				p2.cao = (*p)->cao;
				p2.alpha = (*p)->alpha;
				p2.r_cut = d->r_cut;
				tune_broadcast_command(d, CMD_COMPUTE_ERROR_ESTIMATE);
				d->errorEstimate->compute_master(p2, d->sum_qpart, d->sum_q2, d->box_l, error, rs_error, ks_error);
				(*p)->grid[0] = grid1d;
				(*p)->grid[1] = grid1d;
				(*p)->grid[2] = grid1d;
				(*p)->error = error;
				(*p)->rs_error = rs_error;
				(*p)->ks_error = ks_error;
				P3M_DEBUG(printf("        => error=" FFLOATE "\n", error));
			} while (error > d->tolerance_field);

			// reached largest possible grid, remove the rest of the parameter sets
			if (step_ix >= num_steps_good_gridsize) {
				P3M_INFO(printf("    Too large grid size, skipping rest of parameter sets.\n"));
				for (tune_params_l::iterator pit = params_to_try.end();
						pit != p; --pit) {
					delete[] *pit;
					params_to_try.erase(pit);
				}
				break;
			}

			// error would be small enough at upper_ix, but not at lower_ix,
			// so bisect to find optimal gridsize
			p3m_int lower_ix = step_good_gridsize[step_ix-1];
			while (lower_ix+1 < upper_ix) {
				p3m_int test_ix = (lower_ix+upper_ix)/2;
				p3m_int grid1d = good_gridsize[test_ix];
				P3M_DEBUG(printf("      fine grid=" FINT "\n", grid1d));
				// fixme
				Parameters p2;
				p2.grid[0] = grid1d;
				p2.grid[1] = grid1d;
				p2.grid[2] = grid1d;
				p2.cao = (*p)->cao;
				p2.alpha = (*p)->alpha;
				p2.r_cut = d->r_cut;
                tune_broadcast_command(d, CMD_COMPUTE_ERROR_ESTIMATE);
				d->errorEstimate->compute_master(p2, d->sum_qpart, d->sum_q2, d->box_l, error, rs_error, ks_error);
				(*p)->grid[0] = grid1d;
				(*p)->grid[1] = grid1d;
				(*p)->grid[2] = grid1d;
				(*p)->error = error;
				(*p)->rs_error = rs_error;
				(*p)->ks_error = ks_error;
				P3M_DEBUG(printf("          => error=" FFLOATE "\n", error));
				if (error < d->tolerance_field) {
					// parameters achieve error
					upper_ix = test_ix;
					(*p)->error = error;
					(*p)->rs_error = rs_error;
					(*p)->ks_error = ks_error;
				} else {
					// parameters do not achieve error
					lower_ix = test_ix;
				}
			}

			// now the right size is at upper_ix
			p3m_int grid1d = good_gridsize[upper_ix];

			// store the new grid size and alpha
			if (min_grid1d > grid1d) min_grid1d = grid1d;

			(*p)->grid[0] = grid1d;
			(*p)->grid[1] = grid1d;
			(*p)->grid[2] = grid1d;
			P3M_INFO(printf( "      => grid=" F3INT ", "                      \
					"error=" FFLOATE "\n",                           \
					(*p)->grid[0], (*p)->grid[1], (*p)->grid[2],     \
					(*p)->error));

			// decrease step_ix so that the same step_ix is tested for the
			// next param set
			step_ix--;

			// compare grid size to previous data set
			// if it is larger than any previous size + P3M_MAX_GRID_DIFF, skip it
			if (min_grid1d + P3M_MAX_GRID_DIFF < grid1d) {
				P3M_INFO(printf("      grid too large => removing data set\n"));
				// remove the rest of the params
				for (tune_params_l::iterator pit = p;
						pit != params_to_try.end(); ++pit) {
					delete[] *pit;
					*pit = NULL;
				}
				break;
			}

			if (p != params_to_try.begin()) {
				tune_params_l::iterator prev = p;
				prev--;
				if ((*prev)->grid[0] >= (*p)->grid[0]) {
					P3M_INFO(printf("      better than previous => removing previous data set.\n"));
					delete[] *prev;
					params_to_try.erase(prev);
				}
			}
		}

		// remove empty param sets
		params_to_try.remove(NULL);


	} else {
		// !tune_grid

		for (tune_params_l::iterator p = params_to_try.begin();
				p != params_to_try.end(); ++p) {
			// test whether accuracy can be achieved with given parameters
			// fixme
			Parameters p2;
			p3m_float error, rs_error, ks_error;
			p2.grid[0] = d->grid[0];
			p2.grid[1] = d->grid[1];
			p2.grid[2] = d->grid[2];
			p2.cao = (*p)->cao;
			p2.alpha = (*p)->alpha;
			p2.r_cut = d->r_cut;
			P3M_INFO(printf("    r_cut=" FFLOAT ", "                   \
					"cao=" FINT ", "                           \
					"grid=" F3INT " (fixed)\n",                \
					p2.r_cut, p2.cao,                          \
					p2.grid[0], p2.grid[1], p2.grid[2]));
            tune_broadcast_command(d, CMD_COMPUTE_ERROR_ESTIMATE);
			d->errorEstimate->compute_master(p2, d->sum_qpart, d->sum_q2, d->box_l, error, rs_error, ks_error);

			if (error < d->tolerance_field) {
				// error is small enough for this parameter set, so keep it
				P3M_INFO(printf("    error (%le) < tolerance (%le), keeping params\n", \
						d->error, d->tolerance_field));
				(*p)->grid[0] = d->grid[0];
				(*p)->grid[1] = d->grid[1];
				(*p)->grid[2] = d->grid[2];
				(*p)->alpha = d->alpha;
				(*p)->error = error;
				(*p)->rs_error = rs_error;
				(*p)->ks_error = ks_error;
			} else {
				// otherwise remove this parameter set
				P3M_INFO(printf("    error (%le) > tolerance (%le), removing params\n", \
						error, d->tolerance_field));
				delete[] *p;
				*p = NULL;
			}
		}

		// really remove the removed param sets
		params_to_try.remove(NULL);
	}

	return time_params(d, num_particles, positions, charges, params_to_try);
}

tune_params*
time_params(Solver *d,
		p3m_int num_particles, p3m_float *positions, p3m_float *charges,
		tune_params_l &params_to_try) {
	// fixme
	Communication &comm = d->comm;

	/* Now time the different parameter sets */
	tune_params *best_params = NULL;
	double best_timing = 1.e100;

#ifdef P3M_ENABLE_INFO
	printf("Timing %ld param sets...\n", params_to_try.size());
#endif

	for (tune_params_l::iterator it = params_to_try.begin();
			it != params_to_try.end();
			it++) {
		tune_params &p = **it;
		/* use the parameters */
		d->alpha = p.alpha;
		d->grid[0] = p.grid[0];
		d->grid[1] = p.grid[1];
		d->grid[2] = p.grid[2];
		d->cao = p.cao;

		timing(d, num_particles, positions, charges);
		p.timing = d->timings[TIMING];
		p.timing_near = d->timings[TIMING_NEAR];
		p.timing_far = d->timings[TIMING_FAR];
		P3M_INFO(printf( "  Timing r_cut=" FFLOAT ", "                  \
				"alpha=" FFLOAT ", "                           \
				"grid=" F3INT ", "                             \
				"cao=" FINT " "                                \
				"=> timing=" FFLOAT " "                        \
				"(" FFLOAT " near, " FFLOAT " far)\n",         \
				d->r_cut, d->alpha,                            \
				d->grid[0], d->grid[1], d->grid[2],            \
				d->cao, p.timing,                \
				p.timing_near, p.timing_far));

		if (p.timing < best_timing) {
			best_timing = p.timing;
			best_params = &p;
		}
	}

	if (best_params == NULL)
		throw std::runtime_error("Internal error: Tuning could not get a best timing.");

	for (tune_params_l::iterator p = params_to_try.begin();
			p != params_to_try.end(); ++p) {
		if (*p != best_params)
			delete[] *p;
	}

	return best_params;
}

/** Calculate number of charged particles, the sum of the squared
      charges and the squared sum of the charges. Called in parallel at
      the beginning of tuning. */
void
count_charges(Solver *d, p3m_int num_particles, p3m_float *charges) {
	p3m_float node_sums[3], tot_sums[3];

	for (int i=0; i<3; i++) {
		node_sums[i] = 0.0;
		tot_sums[i] = 0.0;
	}

	for (int i = 0; i < num_particles; i++) {
		if (!float_is_zero(charges[i])) {
			node_sums[0] += 1.0;
			node_sums[1] += SQR(charges[i]);
			node_sums[2] += charges[i];
		}
	}

	MPI_Allreduce(node_sums, tot_sums, 3, P3M_MPI_FLOAT, MPI_SUM, d->comm.mpicomm);
	d->sum_qpart    = (p3m_int)(tot_sums[0]+0.1);
	d->sum_q2       = tot_sums[1];
	d->square_sum_q = SQR(tot_sums[2]);

	P3M_DEBUG(printf("  count_charges(): "			\
			"num_charged_particles=" FINT ", "                   \
			"sum_squared_charges=" FFLOAT ", "                   \
			"sum_charges=" FFLOAT "\n",                          \
			d->sum_qpart, d->sum_q2, sqrt(d->square_sum_q)));
}
}
