#include "fmm_c.h"

#if FMM_MAXNMULTIPOLES > 19
const fmm_float_t fmmec_0d_p20[] = {
#include "include/fmmec_0d_p20.h"
};

const fmm_float_t *get_fmmec_0d_p20()
{
  return fmmec_0d_p20;
}
#endif

