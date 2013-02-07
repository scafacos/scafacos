#include "fmm_c.h"

#if FMM_MAXNMULTIPOLES > 33
const fmm_float_t fmmec_0d_p34[] = {
#include "include/fmmec_0d_p34.h"
};

const fmm_float_t *get_fmmec_0d_p34()
{
  return fmmec_0d_p34;
}
#endif

