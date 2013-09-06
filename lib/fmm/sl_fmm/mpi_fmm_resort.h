
#ifndef __MPI_FMM_RESORT_H__
#define __MPI_FMM_RESORT_H__


#include "../../common/resort/resort.h"


typedef struct _fcs_fmm_resort_t
{
  MPI_Comm comm;

  fcs_resort_t resort;

} *fcs_fmm_resort_t;

#define FCS_FMM_RESORT_NULL  NULL


void fmm_resort_create(fcs_fmm_resort_t *fr, fcs_int nparticles, MPI_Comm comm);
void fmm_resort_destroy(fcs_fmm_resort_t *fr);


#endif /* __MPI_FMM_RESORT_H__ */
