#include "fmm_c.h"

#if FMM_MAXNMULTIPOLES > 31
const fmm_float_t fmmec_0d_p32[] = {
#include "include/fmmec_0d_p32.h"
};

const fmm_float_t *get_fmmec_0d_p32()
{
  return fmmec_0d_p32;
}
#endif

