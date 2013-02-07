#include "fmm_c.h"

#if FMM_MAXNMULTIPOLES > 37
const fmm_float_t fmmec_0d_p38[] = {
#include "include/fmmec_0d_p38.h"
};

const fmm_float_t *get_fmmec_0d_p38()
{
  return fmmec_0d_p38;
}
#endif

