#include "fmm_c.h"

#if FMM_MAXNMULTIPOLES > 32
const fmm_float_t fmmec_0d_p33[] = {
#include "include/fmmec_0d_p33.h"
};

const fmm_float_t *get_fmmec_0d_p33()
{
  return fmmec_0d_p33;
}
#endif

