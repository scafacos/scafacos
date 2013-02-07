#include "fmm_c.h"

#if FMM_MAXNMULTIPOLES > 36
const fmm_float_t fmmec_0d_p37[] = {
#include "include/fmmec_0d_p37.h"
};

const fmm_float_t *get_fmmec_0d_p37()
{
  return fmmec_0d_p37;
}
#endif

