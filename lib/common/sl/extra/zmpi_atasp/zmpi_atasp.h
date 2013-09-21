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


#ifndef __ZMPI_ATASP_H__
#define __ZMPI_ATASP_H__


#include "spec_public_conf.h"
#include "spec_public.h"


#ifdef ZMPI_PREFIX
# include "zmpi_atasp_rename.h"
#endif


typedef long ZMPI_Count;

typedef struct _ZMPI_Tproc *ZMPI_Tproc;

typedef int ZMPI_TPROC_FN(void *b, ZMPI_Count x, void *tproc_data);
typedef int ZMPI_TPROC_MOD_FN(void *b, ZMPI_Count x, void *tproc_data, void *mod);
typedef int ZMPI_TPROCS_FN(void *b, ZMPI_Count x, void *tproc_data, int *procs);
typedef int ZMPI_TPROCS_MOD_FN(void *b, ZMPI_Count x, void *tproc_data, int *procs, void *mod);

typedef void ZMPI_TPROC_RESET_FN(void *tproc_data);

#define ZMPI_TPROC_RESET_NULL  NULL

typedef struct _ZMPI_Tproc_exdef {
  int type;

  spec_tproc_count_f *tproc_count_db, *tproc_count_ip;
  spec_tproc_rearrange_db_f *tproc_rearrange_db;
  spec_tproc_rearrange_ip_f *tproc_rearrange_ip;

  spec_tproc_mod_count_f *tproc_mod_count_db, *tproc_mod_count_ip;
  spec_tproc_mod_rearrange_db_f *tproc_mod_rearrange_db;
  spec_tproc_mod_rearrange_ip_f *tproc_mod_rearrange_ip;

  spec_tprocs_count_f *tprocs_count_db, *tprocs_count_ip;
  spec_tprocs_rearrange_db_f *tprocs_rearrange_db;
  spec_tprocs_rearrange_ip_f *tprocs_rearrange_ip;

  spec_tprocs_mod_count_f *tprocs_mod_count_db, *tprocs_mod_count_ip;
  spec_tprocs_mod_rearrange_db_f *tprocs_mod_rearrange_db;
  spec_tprocs_mod_rearrange_ip_f *tprocs_mod_rearrange_ip;

} const *ZMPI_Tproc_exdef;

#define ZMPI_TPROC_EXDEF_NULL  NULL

#define ZMPI_TPROC_EXDEF_DEFINE_TPROC(_name_, _tp_, _s_...) \
  SPEC_DEFINE_TPROC(_name_, _tp_, _s_) \
  _s_ const struct _ZMPI_Tproc_exdef _##_name_ = { 1, SPEC_EXT_PARAM_TPROC(_name_), SPEC_EXT_PARAM_TPROC_MOD_NULL, SPEC_EXT_PARAM_TPROCS_NULL, SPEC_EXT_PARAM_TPROCS_MOD_NULL }, *_name_ = &_##_name_;

#define ZMPI_TPROC_EXDEF_DEFINE_TPROC_MOD(_name_, _tp_, _s_...) \
  SPEC_DEFINE_TPROC_MOD(_name_, _tp_, _s_) \
  _s_ const struct _ZMPI_Tproc_exdef _##_name_ = { 2, SPEC_EXT_PARAM_TPROC_NULL, SPEC_EXT_PARAM_TPROC_MOD(_name_), SPEC_EXT_PARAM_TPROCS_NULL, SPEC_EXT_PARAM_TPROCS_MOD_NULL }, *_name_ = &_##_name_;

#define ZMPI_TPROC_EXDEF_DEFINE_TPROCS(_name_, _tp_, _s_...) \
  SPEC_DEFINE_TPROCS(_name_, _tp_, _s_) \
  _s_ const struct _ZMPI_Tproc_exdef _##_name_ = { 3, SPEC_EXT_PARAM_TPROC_NULL, SPEC_EXT_PARAM_TPROC_MOD_NULL, SPEC_EXT_PARAM_TPROCS(_name_), SPEC_EXT_PARAM_TPROCS_MOD_NULL }, *_name_ = &_##_name_;

#define ZMPI_TPROC_EXDEF_DEFINE_TPROCS_MOD(_name_, _tp_, _s_...) \
  SPEC_DEFINE_TPROCS_MOD(_name_, _tp_, _s_) \
  _s_ const struct _ZMPI_Tproc_exdef _##_name_ = { 4, SPEC_EXT_PARAM_TPROC_NULL, SPEC_EXT_PARAM_TPROC_MOD_NULL, SPEC_EXT_PARAM_TPROCS_NULL, SPEC_EXT_PARAM_TPROCS_MOD(_name_) }, *_name_ = &_##_name_;

int ZMPI_Tproc_create_tproc(ZMPI_Tproc *tproc, ZMPI_TPROC_FN *tfn, ZMPI_TPROC_RESET_FN *rfn, ZMPI_Tproc_exdef exdef);
int ZMPI_Tproc_create_tproc_mod(ZMPI_Tproc *tproc, ZMPI_TPROC_MOD_FN *tfn, ZMPI_TPROC_RESET_FN *rfn, ZMPI_Tproc_exdef exdef);
int ZMPI_Tproc_create_tprocs(ZMPI_Tproc *tproc, ZMPI_TPROCS_FN *tfn, ZMPI_TPROC_RESET_FN *rfn, ZMPI_Tproc_exdef exdef);
int ZMPI_Tproc_create_tprocs_mod(ZMPI_Tproc *tproc, ZMPI_TPROCS_MOD_FN *tfn, ZMPI_TPROC_RESET_FN *rfn, ZMPI_Tproc_exdef exdef);
int ZMPI_Tproc_free(ZMPI_Tproc *tproc);

int ZMPI_Tproc_set_neighbors(ZMPI_Tproc tproc, int nneighbors, int *neighbors, MPI_Comm comm);
int ZMPI_Tproc_set_proclists(ZMPI_Tproc tproc, int ndstprocs, int *dstprocs, int nsrcprocs, int *srcprocs, MPI_Comm comm);

typedef int ZMPI_ALLTOALL_SPECIFIC_FN(void *sbuf, int scount, MPI_Datatype stype, void *rbuf, int rcount, MPI_Datatype rtype, ZMPI_Tproc tproc, void *tproc_data, int *received, MPI_Comm comm);

int ZMPI_Alltoall_specific(void *sbuf, int scount, MPI_Datatype stype, void *rbuf, int rcount, MPI_Datatype rtype, ZMPI_Tproc tproc, void *tproc_data, int *received, MPI_Comm comm);

int ZMPI_Neighbor_alltoall_specific(void *sbuf, int scount, MPI_Datatype stype, void *rbuf, int rcount, MPI_Datatype rtype, ZMPI_Tproc tproc, void *tproc_data, int *received, MPI_Comm comm);


#endif /* __ZMPI_ATASP_H__ */
