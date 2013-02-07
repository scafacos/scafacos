#include "fmm_c.h"

#if FMM_MAXNMULTIPOLES > 14
const fmm_float_t fmmec_0d_p15[] = {
#include "include/fmmec_0d_p15.h"
};

const fmm_float_t *get_fmmec_0d_p15()
{
  return fmmec_0d_p15;
}
#endif

