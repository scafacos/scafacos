#include "fmm_c.h"

#if FMM_MAXNMULTIPOLES > 39
const fmm_float_t fmmec_0d_p40[] = {
#include "include/fmmec_0d_p40.h"
};

const fmm_float_t *get_fmmec_0d_p40()
{
  return fmmec_0d_p40;
}
#endif

