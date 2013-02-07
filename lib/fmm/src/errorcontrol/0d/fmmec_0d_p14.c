#include "fmm_c.h"

#if FMM_MAXNMULTIPOLES > 13
const fmm_float_t fmmec_0d_p14[] = {
#include "include/fmmec_0d_p14.h"
};

const fmm_float_t *get_fmmec_0d_p14()
{
  return fmmec_0d_p14;
}
#endif

