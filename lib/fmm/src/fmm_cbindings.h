
#ifndef __FMM_CBINDINGS_H__
#define __FMM_CBINDINGS_H__


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


void fmm_cinit(void*);
void fmm_ctune(long long,fcs_float*,fcs_float*,long long,long long,fcs_float,long long,long long*,fcs_float,long long, long long, long long, void*,long long*,long long*);
void fmm_ctunehomogen(void*, long long*, long long*);
void fmm_ccomputewigner(void*,void*,long long);
void fmm_crun(long long,fcs_float*,fcs_float*,fcs_float*,fcs_float*,fcs_float*,long long,long long,fcs_float,long
long,long long*,fcs_float,long long,long long, long long, long long, void*,long long*);
void fmm_cfinalize(void*,long long);
void fmm_cinitload(void *, fcs_float *, fcs_int);
void fmm_csetload(void *, fcs_float);

void fmm_csetpresorted(void *, long long);
void fmm_cinitresort(void *, fcs_fmm_resort_t);
void fmm_csetresort(void *, long long);


#endif /* __FMM_CBINDINGS_H__ */
