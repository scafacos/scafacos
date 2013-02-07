#include "../fmm_c.h"

#ifdef FMM_GETDISTMS
#ifdef FMM_C_CONSTANTS
const fmm_int_t fmmdist_maxdist[] = {
#include "include/fmmdist_maxdist.h"
};

const fmm_int_t *get_maxdist()
{
  return fmmdist_maxdist;
}
#endif
#endif

