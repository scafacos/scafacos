#include "fmm_c.h"

#if FMM_MAXNMULTIPOLES > 47
const fmm_float_t fmmec_0d_p48[] = {
#include "include/fmmec_0d_p48.h"
};

const fmm_float_t *get_fmmec_0d_p48()
{
  return fmmec_0d_p48;
}
#endif

