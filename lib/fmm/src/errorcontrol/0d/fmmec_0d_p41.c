#include "fmm_c.h"

#if FMM_MAXNMULTIPOLES > 40
const fmm_float_t fmmec_0d_p41[] = {
#include "include/fmmec_0d_p41.h"
};

const fmm_float_t *get_fmmec_0d_p41()
{
  return fmmec_0d_p41;
}
#endif

