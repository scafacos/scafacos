#include "fmm_c.h"

#if FMM_MAXNMULTIPOLES > 28
const fmm_float_t fmmec_0d_p29[] = {
#include "include/fmmec_0d_p29.h"
};

const fmm_float_t *get_fmmec_0d_p29()
{
  return fmmec_0d_p29;
}
#endif

