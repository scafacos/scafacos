#include "fmm_c.h"

#if FMM_ERRORCONTROL_STAGE > 0
#  if FMM_MAXNMULTIPOLES > 27
   const fmm_float_t fmmec_2d_p28[] = {
#  include "include/fmmec_2d_p28.h"
   };

   const fmm_float_t *get_fmmec_2d_p28()
   {
     return fmmec_2d_p28;
   }
#  endif
#endif
