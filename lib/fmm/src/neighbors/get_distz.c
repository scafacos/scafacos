#include "../fmm_c.h"

#ifdef FMM_GETDIST
#ifdef FMM_C_CONSTANTS
const fmm_int_t fmmdist_distz[] = {
#include "include/fmmdist_distz.h"
};

const fmm_int_t *get_distz()
{
  return fmmdist_distz;
}
#endif
#endif

