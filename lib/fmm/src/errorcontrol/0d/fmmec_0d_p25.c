#include "fmm_c.h"

#if FMM_MAXNMULTIPOLES > 24
const fmm_float_t fmmec_0d_p25[] = {
#include "include/fmmec_0d_p25.h"
};

const fmm_float_t *get_fmmec_0d_p25()
{
  return fmmec_0d_p25;
}
#endif

