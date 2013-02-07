#include "fmm_c.h"

#if FMM_MAXNMULTIPOLES > 43
const fmm_float_t fmmec_0d_p44[] = {
#include "include/fmmec_0d_p44.h"
};

const fmm_float_t *get_fmmec_0d_p44()
{
  return fmmec_0d_p44;
}
#endif

