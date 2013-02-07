#include "fmm_c.h"

#if FMM_ERRORCONTROL_STAGE > 0
#  if FMM_MAXNMULTIPOLES > 13
   const fmm_float_t fmmec_1d_p14[] = {
#  include "include/fmmec_1d_p14.h"
   };

   const fmm_float_t *get_fmmec_1d_p14()
   {
     return fmmec_1d_p14;
   }
#  endif
#endif
