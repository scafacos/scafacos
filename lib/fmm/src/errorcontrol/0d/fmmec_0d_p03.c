#include "fmm_c.h"

#if FMM_MAXNMULTIPOLES > 2
const fmm_float_t fmmec_0d_p03[] = {
#include "include/fmmec_0d_p03.h"
};

const fmm_float_t *get_fmmec_0d_p03()
{
  return fmmec_0d_p03;
}
#endif

