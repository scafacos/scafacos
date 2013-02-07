#include "fmm_c.h"

#if FMM_MAXNMULTIPOLES > 4
const fmm_float_t fmmec_0d_p05[] = {
#include "include/fmmec_0d_p05.h"
};

const fmm_float_t *get_fmmec_0d_p05()
{
  return fmmec_0d_p05;
}
#endif

