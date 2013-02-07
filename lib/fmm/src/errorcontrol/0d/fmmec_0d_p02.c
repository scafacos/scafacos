#include "fmm_c.h"

#if FMM_MAXNMULTIPOLES > 1
const fmm_float_t fmmec_0d_p02[] = {
#include "include/fmmec_0d_p02.h"
};

const fmm_float_t *get_fmmec_0d_p02()
{
  return fmmec_0d_p02;
}
#endif

