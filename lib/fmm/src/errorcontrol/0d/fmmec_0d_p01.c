#include "fmm_c.h"

#if FMM_MAXNMULTIPOLES > 0
const fmm_float_t fmmec_0d_p01[] = {
#include "include/fmmec_0d_p01.h"
};

const fmm_float_t *get_fmmec_0d_p01()
{
  return fmmec_0d_p01;
}
#endif
