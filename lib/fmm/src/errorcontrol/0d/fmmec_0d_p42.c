#include "fmm_c.h"

#if FMM_MAXNMULTIPOLES > 41
const fmm_float_t fmmec_0d_p42[] = {
#include "include/fmmec_0d_p42.h"
};

const fmm_float_t *get_fmmec_0d_p42()
{
  return fmmec_0d_p42;
}
#endif

