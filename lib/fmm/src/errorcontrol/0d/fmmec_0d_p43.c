#include "fmm_c.h"

#if FMM_MAXNMULTIPOLES > 42
const fmm_float_t fmmec_0d_p43[] = {
#include "include/fmmec_0d_p43.h"
};

const fmm_float_t *get_fmmec_0d_p43()
{
  return fmmec_0d_p43;
}
#endif

