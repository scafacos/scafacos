/*
 *  Copyright (C) 2011, 2012, 2013, 2014, 2015 Michael Hofmann
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


#if MPI_VERSION >= 3
# define IF_ELSE_MPI_VERSION_3(_if_, _else_)  _if_
#else
# define IF_ELSE_MPI_VERSION_3(_if_, _else_)  _else_
#endif

typedef struct _spec_tproc_t *ZMPI_Tproc;

#define ZMPI_TPROC_NULL  NULL

#if MPI_VERSION < 3
typedef long ZMPI_Count;
#endif

typedef int ZMPI_TPROC_FN(void *b, IF_ELSE_MPI_VERSION_3(MPI_Count, ZMPI_Count) x, void *tproc_data);
typedef int ZMPI_TPROC_MOD_FN(void *b, IF_ELSE_MPI_VERSION_3(MPI_Count, ZMPI_Count) x, void *tproc_data, void *mod);
typedef void ZMPI_TPROCS_FN(void *b, IF_ELSE_MPI_VERSION_3(MPI_Count, ZMPI_Count) x, void *tproc_data, int *nprocs, int *procs);
typedef void ZMPI_TPROCS_MOD_FN(void *b, IF_ELSE_MPI_VERSION_3(MPI_Count, ZMPI_Count) x, void *tproc_data, int *nprocs, int *procs, void *mod);

typedef void ZMPI_TPROC_RESET_FN(void *tproc_data);

#define ZMPI_TPROC_RESET_NULL  NULL

typedef struct _ZMPI_Tproc_exdef
{
  int type;

  spec_tproc_ext_t tproc_ext;
  spec_tproc_mod_ext_t tproc_mod_ext;
  spec_tprocs_ext_t tprocs_ext;
  spec_tprocs_mod_ext_t tprocs_mod_ext;

} const *ZMPI_Tproc_exdef;

#define ZMPI_TPROC_EXDEF_NULL  NULL

/* default */
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

/* fixtype */
#define ZMPI_TPROC_EXDEF_DEFINE_FIXTYPE_TPROC(_name_, _fxt_, _tp_, _s_...) \
  SPEC_DEFINE_FIXED_TPROC(_name_, ZMPI_FIXTYPE_DECLARE, _fxt_, ZMPI_FIXTYPE_CREATE, ZMPI_FIXTYPE_COPY_AT, ZMPI_FIXTYPE_EXCHANGE_AT, ZMPI_FIXTYPE_DESTROY, _tp_, _s_) \
  _s_ const struct _ZMPI_Tproc_exdef _##_name_ = { 1, SPEC_EXT_PARAM_TPROC(_name_), SPEC_EXT_PARAM_TPROC_MOD_NULL, SPEC_EXT_PARAM_TPROCS_NULL, SPEC_EXT_PARAM_TPROCS_MOD_NULL }, *_name_ = &_##_name_;

#define ZMPI_TPROC_EXDEF_DEFINE_FIXTYPE_TPROC_MOD(_name_, _fxt_, _tp_, _s_...) \
  SPEC_DEFINE_FIXED_TPROC_MOD(_name_, ZMPI_FIXTYPE_DECLARE, _fxt_, ZMPI_FIXTYPE_CREATE, ZMPI_FIXTYPE_COPY_AT, ZMPI_FIXTYPE_EXCHANGE_AT, ZMPI_FIXTYPE_DESTROY, _tp_, _s_) \
  _s_ const struct _ZMPI_Tproc_exdef _##_name_ = { 2, SPEC_EXT_PARAM_TPROC_NULL, SPEC_EXT_PARAM_TPROC_MOD(_name_), SPEC_EXT_PARAM_TPROCS_NULL, SPEC_EXT_PARAM_TPROCS_MOD_NULL }, *_name_ = &_##_name_;

#define ZMPI_TPROC_EXDEF_DEFINE_FIXTYPE_TPROCS(_name_, _fxt_, _tp_, _s_...) \
  SPEC_DEFINE_FIXED_TPROCS(_name_, ZMPI_FIXTYPE_DECLARE, _fxt_, ZMPI_FIXTYPE_CREATE, ZMPI_FIXTYPE_COPY_AT, ZMPI_FIXTYPE_EXCHANGE_AT, ZMPI_FIXTYPE_DESTROY, _tp_, _s_) \
  _s_ const struct _ZMPI_Tproc_exdef _##_name_ = { 3, SPEC_EXT_PARAM_TPROC_NULL, SPEC_EXT_PARAM_TPROC_MOD_NULL, SPEC_EXT_PARAM_TPROCS(_name_), SPEC_EXT_PARAM_TPROCS_MOD_NULL }, *_name_ = &_##_name_;

#define ZMPI_TPROC_EXDEF_DEFINE_FIXTYPE_TPROCS_MOD(_name_, _fxt_, _tp_, _s_...) \
  SPEC_DEFINE_FIXED_TPROCS_MOD(_name_, ZMPI_FIXTYPE_DECLARE, _fxt_, ZMPI_FIXTYPE_CREATE, ZMPI_FIXTYPE_COPY_AT, ZMPI_FIXTYPE_EXCHANGE_AT, ZMPI_FIXTYPE_DESTROY, _tp_, _s_) \
  _s_ const struct _ZMPI_Tproc_exdef _##_name_ = { 4, SPEC_EXT_PARAM_TPROC_NULL, SPEC_EXT_PARAM_TPROC_MOD_NULL, SPEC_EXT_PARAM_TPROCS_NULL, SPEC_EXT_PARAM_TPROCS_MOD(_name_) }, *_name_ = &_##_name_;

/* fixsize */
#define ZMPI_TPROC_EXDEF_DEFINE_FIXSIZE_TPROC(_name_, _fxs_, _tp_, _s_...) \
  _s_ const int _name_##_params = _fxs_; \
  SPEC_DEFINE_FIXED_TPROC(_name_, ZMPI_FIXSIZE_DECLARE, _name_##_params, ZMPI_FIXSIZE_CREATE, ZMPI_FIXSIZE_COPY_AT, ZMPI_FIXSIZE_EXCHANGE_AT, ZMPI_FIXSIZE_DESTROY, _tp_, _s_) \
  _s_ const struct _ZMPI_Tproc_exdef _##_name_ = { 1, SPEC_EXT_PARAM_TPROC(_name_), SPEC_EXT_PARAM_TPROC_MOD_NULL, SPEC_EXT_PARAM_TPROCS_NULL, SPEC_EXT_PARAM_TPROCS_MOD_NULL }, *_name_ = &_##_name_;

int ZMPI_Tproc_create_tproc(ZMPI_Tproc *tproc, ZMPI_TPROC_FN *tfn, ZMPI_TPROC_RESET_FN *rfn, const ZMPI_Tproc_exdef exdef);
int ZMPI_Tproc_create_tproc_mod(ZMPI_Tproc *tproc, ZMPI_TPROC_MOD_FN *tfn, ZMPI_TPROC_RESET_FN *rfn, const ZMPI_Tproc_exdef exdef);
int ZMPI_Tproc_create_tprocs(ZMPI_Tproc *tproc, int max_tprocs, ZMPI_TPROCS_FN *tfn, ZMPI_TPROC_RESET_FN *rfn, const ZMPI_Tproc_exdef exdef);
int ZMPI_Tproc_create_tprocs_mod(ZMPI_Tproc *tproc, int max_tprocs, ZMPI_TPROCS_MOD_FN *tfn, ZMPI_TPROC_RESET_FN *rfn, const ZMPI_Tproc_exdef exdef);
int ZMPI_Tproc_free(ZMPI_Tproc *tproc);

int ZMPI_Tproc_set_neighbors(ZMPI_Tproc tproc, int nneighbors, int *neighbors, MPI_Comm comm);
int ZMPI_Tproc_set_proclists(ZMPI_Tproc tproc, int ndstprocs, int *dstprocs, int nsrcprocs, int *srcprocs, MPI_Comm comm);

#define ZMPI_ALLTOALL_SPECIFIC_TYPE_ALLTOALLV    0
#define ZMPI_ALLTOALL_SPECIFIC_TYPE_ALLTOALLW    1
#define ZMPI_ALLTOALL_SPECIFIC_TYPE_PUT          2
#define ZMPI_ALLTOALL_SPECIFIC_TYPE_PUT_2PHASES  3
#define ZMPI_ALLTOALL_SPECIFIC_TYPE_SENDRECV     4

#define ZMPI_ALLTOALL_SPECIFIC_TYPE_DEFAULT  ZMPI_ALLTOALL_SPECIFIC_TYPE_ALLTOALLV

extern int ZMPI_Alltoall_specific_type;

#define ZMPI_NEIGHBOR_ALLTOALL_SPECIFIC_TYPE_DEFAULT  ZMPI_ALLTOALL_SPECIFIC_TYPE_ALLTOALLV

extern int ZMPI_Neighbor_alltoall_specific_type;

#if MPI_VERSION < 3
typedef int ZMPI_Status;
int ZMPI_Get_elements(const ZMPI_Status *status, MPI_Datatype datatype, int *count);
# define ZMPI_STATUS_IGNORE  NULL
#endif

typedef int ZMPI_ALLTOALL_SPECIFIC_FN(void *sbuf, int scount, MPI_Datatype stype, void *rbuf, int rcount, MPI_Datatype rtype, ZMPI_Tproc tproc, void *tproc_data, MPI_Comm comm, IF_ELSE_MPI_VERSION_3(MPI_Status, ZMPI_Status) *status);
int ZMPI_Alltoall_specific(void *sbuf, int scount, MPI_Datatype stype, void *rbuf, int rcount, MPI_Datatype rtype, ZMPI_Tproc tproc, void *tproc_data, MPI_Comm comm, IF_ELSE_MPI_VERSION_3(MPI_Status, ZMPI_Status) *status);
int ZMPI_Neighbor_alltoall_specific(void *sbuf, int scount, MPI_Datatype stype, void *rbuf, int rcount, MPI_Datatype rtype, ZMPI_Tproc tproc, void *tproc_data, MPI_Comm comm, IF_ELSE_MPI_VERSION_3(MPI_Status, ZMPI_Status) *status);

#undef IF_ELSE_MPI_VERSION_3


#endif /* __ZMPI_ATASP_H__ */
