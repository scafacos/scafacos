#include "fmm_c.h"

#if FMM_ERRORCONTROL_STAGE > 0
#  if FMM_MAXNMULTIPOLES > 30
   const fmm_float_t fmmec_2d_p31[] = {
#  include "include/fmmec_2d_p31.h"
   };

   const fmm_float_t *get_fmmec_2d_p31()
   {
     return fmmec_2d_p31;
   }
#  endif
#endif
