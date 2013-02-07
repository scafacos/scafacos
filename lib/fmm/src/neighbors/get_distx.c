#include "../fmm_c.h"

#ifdef FMM_GETDIST
#ifdef FMM_C_CONSTANTS
const fmm_int_t fmmdist_distx[] = {
#include "include/fmmdist_distx.h"
};

const fmm_int_t *get_distx()
{
  return fmmdist_distx;
}
#endif
#endif

