#include "fmm_c.h"

#if FMM_ERRORCONTROL_STAGE > 0
#  if FMM_MAXNMULTIPOLES > 14
   const fmm_float_t fmmec_1d_p15[] = {
#  include "include/fmmec_1d_p15.h"
   };

   const fmm_float_t *get_fmmec_1d_p15()
   {
     return fmmec_1d_p15;
   }
#  endif
#endif
