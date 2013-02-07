#include "fmm_c.h"

#if FMM_MAXNMULTIPOLES > 6
const fmm_float_t fmmec_0d_p07[] = {
#include "include/fmmec_0d_p07.h"
};

const fmm_float_t *get_fmmec_0d_p07()
{
  return fmmec_0d_p07;
}
#endif

