#include "fmm_c.h"

#if FMM_MAXNMULTIPOLES > 10
const fmm_float_t fmmec_0d_p11[] = {
#include "include/fmmec_0d_p11.h"
};

const fmm_float_t *get_fmmec_0d_p11()
{
  return fmmec_0d_p11;
}
#endif

