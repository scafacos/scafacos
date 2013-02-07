#include "fmm_c.h"

#if FMM_ERRORCONTROL_STAGE > 0
#  if FMM_MAXNMULTIPOLES > 0
   const fmm_float_t fmmec_2d_p01[] = {
#  include "include/fmmec_2d_p01.h"
   };

   const fmm_float_t *get_fmmec_2d_p01()
   {
     return fmmec_2d_p01;
   }
#  endif
#endif
