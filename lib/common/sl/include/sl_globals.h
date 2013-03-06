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


#ifndef __SL_GLOBALS_H__
#define __SL_GLOBALS_H__


/* src/base/base.c */
extern sl_context_t sl_default_context;
extern const int default_sl_dummy_rank;
#ifdef SL_USE_RTI
extern const rti_t default_rti;
#endif
extern const slint_t default_sr_ip_threshold;
extern const slint_t default_sr_db_threshold;
extern const slint_t default_sr_ma_threshold;
extern const slint_t default_sri_threshold;

/* src/base_mpi/base_mpi.c */
#ifdef SL_USE_MPI
extern const MPI_Datatype default_mpi_int_datatype;
extern const MPI_Datatype default_mpi_key_datatype;
extern const MPI_Datatype default_mpi_pkey_datatype;
extern const MPI_Datatype default_mpi_pwkey_datatype;
extern const MPI_Datatype default_mpi_index_datatype;
extern const MPI_Datatype default_mpi_weight_datatype;
extern const MPI_Datatype default_mpi_data_datatype[];
extern const int default_mpi_rank;
#endif
extern const void *default_me_sendrecv_replace_mem;
extern const slint_t default_me_sendrecv_replace_memsize;
extern const slint_t default_me_sendrecv_replace_mpi_maxsize;
extern const double default_meas_t[];
extern const slint_t default_meas_max_nprocs;
extern const slint_t default_meas_packed;
extern const slint_t default_meas_minalloc;
extern const double default_meas_overalloc;
extern const slint_t default_mea_packed;
extern const slint_t default_mea_db_packed;
extern const slint_t default_mea_ip_packed;
#ifdef MSEG_ROOT
extern const int default_mseg_root;
#endif
#ifdef MSEG_BORDER_UPDATE_REDUCTION
extern const double default_mseg_border_update_count_reduction;
extern const double default_mseg_border_update_weight_reduction;
#endif
#ifdef MSEG_FORWARD_ONLY
extern const slint_t default_mseg_forward_only;
#endif
#ifdef MSEG_INFO
extern const slint_t default_mseg_info_rounds;
extern const slint_t *default_mseg_info_finish_rounds;
extern const double default_mseg_info_finish_rounds_avg;
extern const slweight_t default_mseg_info_total_weights;
#endif
extern const slint_t default_mseg_binnings;
extern const slint_t default_mseg_finalize_mode;
#ifdef MSS_ROOT
extern const int default_mss_root;
#endif
extern const double default_msm_t[];
extern const slint_t default_msm_sync;
extern const double default_msp_t[];
extern const slint_t default_msp_sync;
extern const partcond_t *default_msp_r_pc;
extern const double default_mssp_i_t[];
extern const double default_mssp_p_t[];
extern const double default_mssp_b_t[];
extern const slint_t default_mssp_sync;
extern const slint_t default_mssp_i_sync;
extern const slint_t default_mssp_p_sync;
extern const slint_t default_mssp_b_sync;
extern const slint_t default_mssp_back_packed;


#endif /* __SL_GLOBALS_H__ */
