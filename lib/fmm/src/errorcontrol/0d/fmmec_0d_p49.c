#include "fmm_c.h"

#if FMM_MAXNMULTIPOLES > 48
const fmm_float_t fmmec_0d_p49[] = {
#include "include/fmmec_0d_p49.h"
};

const fmm_float_t *get_fmmec_0d_p49()
{
  return fmmec_0d_p49;
}
#endif

