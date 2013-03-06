/*
 *  Copyright (C) 2011, 2012, 2013 Michael Hofmann
 *  
 *  This file is part of ScaFaCoS.
 *  
 *  ScaFaCoS is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *  
 *  ScaFaCoS is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser Public License for more details.
 *  
 *  You should have received a copy of the GNU Lesser Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *  

 *  
 *  SL - Sorting Library, michael <dot> hofmann <at> informatik <dot> tu-chemnitz <dot> de
 */

/* src/base/base.c */
  struct {
int dummy_rank;
  } sl;
#ifdef SL_USE_RTI
rti_t rti;
#endif
  struct {
slint_t ip_threshold;
slint_t db_threshold;
slint_t ma_threshold;
  } sr;
  struct {
slint_t threshold;
  } sri;
/* src/base_mpi/base_mpi.c */
#ifdef SL_USE_MPI
  struct {
MPI_Datatype int_datatype;
MPI_Datatype key_datatype;
MPI_Datatype pkey_datatype;
MPI_Datatype pwkey_datatype;
MPI_Datatype index_datatype;
MPI_Datatype weight_datatype;
MPI_Datatype data_datatype[SL_DATA_NMAX + 1];
int rank;
  } mpi;
#endif
#ifdef SL_USE_MPI
  struct {
void *sendrecv_replace_mem;
slint_t sendrecv_replace_memsize;
slint_t sendrecv_replace_mpi_maxsize;
  } me;
#endif
#ifdef SL_USE_MPI
  struct {
double t[2];
slint_t max_nprocs;
slint_t packed;
slint_t minalloc;
double overalloc;
  } meas;
#endif
#ifdef SL_USE_MPI
  struct {
slint_t packed;
slint_t db_packed;
slint_t ip_packed;
  } mea;
#endif
#ifdef SL_USE_MPI
  struct {
#ifdef MSEG_ROOT
int root;
#endif
#ifdef MSEG_BORDER_UPDATE_REDUCTION
double border_update_count_reduction;
double border_update_weight_reduction;
#endif
#ifdef MSEG_FORWARD_ONLY
slint_t forward_only;
#endif
#ifdef MSEG_INFO
slint_t info_rounds;
slint_t *info_finish_rounds;
double info_finish_rounds_avg;
slweight_t info_total_weights;
#endif
slint_t binnings;
slint_t finalize_mode;
  } mseg;
#endif
#ifdef SL_USE_MPI
  struct {
#ifdef MSS_ROOT
int root;
#endif
  } mss;
#endif
#ifdef SL_USE_MPI
  struct {
double t[4];
slint_t sync;
  } msm;
#endif
#ifdef SL_USE_MPI
  struct {
double t[4];
slint_t sync;
partcond_t *r_pc;
  } msp;
#endif
#ifdef SL_USE_MPI
  struct {
double i_t[3];
double p_t[3];
double b_t[3];
slint_t sync;
slint_t i_sync;
slint_t p_sync;
slint_t b_sync;
slint_t back_packed;
  } mssp;
#endif
