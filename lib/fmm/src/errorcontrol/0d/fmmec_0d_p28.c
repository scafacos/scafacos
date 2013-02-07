#include "fmm_c.h"

#if FMM_MAXNMULTIPOLES > 27
const fmm_float_t fmmec_0d_p28[] = {
#include "include/fmmec_0d_p28.h"
};

const fmm_float_t *get_fmmec_0d_p28()
{
  return fmmec_0d_p28;
}
#endif

