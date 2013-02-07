#include "fmm_c.h"

#if FMM_MAXNMULTIPOLES > 35
const fmm_float_t fmmec_0d_p36[] = {
#include "include/fmmec_0d_p36.h"
};

const fmm_float_t *get_fmmec_0d_p36()
{
  return fmmec_0d_p36;
}
#endif

