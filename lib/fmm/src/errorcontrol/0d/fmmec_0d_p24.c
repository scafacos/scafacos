#include "fmm_c.h"

#if FMM_MAXNMULTIPOLES > 23
const fmm_float_t fmmec_0d_p24[] = {
#include "include/fmmec_0d_p24.h"
};

const fmm_float_t *get_fmmec_0d_p24()
{
  return fmmec_0d_p24;
}
#endif

