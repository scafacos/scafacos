#include "fmm_c.h"

#if FMM_MAXNMULTIPOLES > 49
const fmm_float_t fmmec_0d_p50[] = {
#include "include/fmmec_0d_p50.h"
};

const fmm_float_t *get_fmmec_0d_p50()
{
  return fmmec_0d_p50;
}
#endif

