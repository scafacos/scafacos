#include "fmm_c.h"

#if FMM_MAXNMULTIPOLES > 21
const fmm_float_t fmmec_0d_p22[] = {
#include "include/fmmec_0d_p22.h"
};

const fmm_float_t *get_fmmec_0d_p22()
{
  return fmmec_0d_p22;
}
#endif

