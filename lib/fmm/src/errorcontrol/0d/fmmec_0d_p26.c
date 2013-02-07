#include "fmm_c.h"

#if FMM_MAXNMULTIPOLES > 25
const fmm_float_t fmmec_0d_p26[] = {
#include "include/fmmec_0d_p26.h"
};

const fmm_float_t *get_fmmec_0d_p26()
{
  return fmmec_0d_p26;
}
#endif

