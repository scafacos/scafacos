#include "fmm_c.h"

#if FMM_MAXNMULTIPOLES > 30
const fmm_float_t fmmec_0d_p31[] = {
#include "include/fmmec_0d_p31.h"
};

const fmm_float_t *get_fmmec_0d_p31()
{
  return fmmec_0d_p31;
}
#endif

