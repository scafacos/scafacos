#include "fmm_c.h"

#if FMM_ERRORCONTROL_STAGE > 0
#  if FMM_MAXNMULTIPOLES > 47
   const fmm_float_t fmmec_2d_p48[] = {
#  include "include/fmmec_2d_p48.h"
   };

   const fmm_float_t *get_fmmec_2d_p48()
   {
     return fmmec_2d_p48;
   }
#  endif
#endif
