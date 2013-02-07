#include "fmm_c.h"

#if FMM_MAXNMULTIPOLES > -1
const fmm_float_t fmmec_0d_p00[] = {
#include "include/fmmec_0d_p00.h"
};

const fmm_float_t *get_fmmec_0d_p00()
{
  return fmmec_0d_p00;
}
#endif

