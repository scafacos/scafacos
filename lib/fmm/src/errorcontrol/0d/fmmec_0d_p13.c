#include "fmm_c.h"

#if FMM_MAXNMULTIPOLES > 12
const fmm_float_t fmmec_0d_p13[] = {
#include "include/fmmec_0d_p13.h"
};

const fmm_float_t *get_fmmec_0d_p13()
{
  return fmmec_0d_p13;
}
#endif

