#include "fmm_c.h"

#if FMM_MAXNMULTIPOLES > 22
const fmm_float_t fmmec_0d_p23[] = {
#include "include/fmmec_0d_p23.h"
};

const fmm_float_t *get_fmmec_0d_p23()
{
  return fmmec_0d_p23;
}
#endif

