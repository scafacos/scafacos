#include "fmm_c.h"

#if FMM_MAXNMULTIPOLES > 26
const fmm_float_t fmmec_0d_p27[] = {
#include "include/fmmec_0d_p27.h"
};

const fmm_float_t *get_fmmec_0d_p27()
{
  return fmmec_0d_p27;
}
#endif

