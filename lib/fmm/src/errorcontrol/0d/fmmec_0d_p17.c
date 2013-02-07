#include "fmm_c.h"

#if FMM_MAXNMULTIPOLES > 16
const fmm_float_t fmmec_0d_p17[] = {
#include "include/fmmec_0d_p17.h"
};

const fmm_float_t *get_fmmec_0d_p17()
{
  return fmmec_0d_p17;
}
#endif

