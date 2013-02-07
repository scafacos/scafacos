#include "../fmm_c.h"

#ifdef FMM_GETDIST
   const fmm_int_t distx_array[] = {
#include "include/distx.h"
   };

   const fmm_int_t *get_distx()
   {
     return distx_array;
   }
#endif
