#include "fmm_c.h"

#if FMM_MAXNMULTIPOLES > 44
const fmm_float_t fmmec_0d_p45[] = {
#include "include/fmmec_0d_p45.h"
};

const fmm_float_t *get_fmmec_0d_p45()
{
  return fmmec_0d_p45;
}
#endif

