#include "fmm_c.h"

#if FMM_MAXNMULTIPOLES > 34
const fmm_float_t fmmec_0d_p35[] = {
#include "include/fmmec_0d_p35.h"
};

const fmm_float_t *get_fmmec_0d_p35()
{
  return fmmec_0d_p35;
}
#endif

