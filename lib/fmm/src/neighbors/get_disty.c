#include "../fmm_c.h"

#ifdef FMM_GETDIST
#ifdef FMM_C_CONSTANTS
const fmm_int_t fmmdist_disty[] = {
#include "include/fmmdist_disty.h"
};

const fmm_int_t *get_disty()
{
  return fmmdist_disty;
}
#endif
#endif

