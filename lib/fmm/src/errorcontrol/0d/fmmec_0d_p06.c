#include "fmm_c.h"

#if FMM_MAXNMULTIPOLES > 5
const fmm_float_t fmmec_0d_p06[] = {
#include "include/fmmec_0d_p06.h"
};

const fmm_float_t *get_fmmec_0d_p06()
{
  return fmmec_0d_p06;
}
#endif

