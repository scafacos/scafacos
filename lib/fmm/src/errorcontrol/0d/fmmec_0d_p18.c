#include "fmm_c.h"

#if FMM_MAXNMULTIPOLES > 17
const fmm_float_t fmmec_0d_p18[] = {
#include "include/fmmec_0d_p18.h"
};

const fmm_float_t *get_fmmec_0d_p18()
{
  return fmmec_0d_p18;
}
#endif

