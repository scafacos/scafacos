#include "../fmm_c.h"

#ifdef FMM_GETDIST
   const fmm_int_t distz_array[] = {
#include "include/disty.h"
   };

   const fmm_int_t *get_distz()
   {
     return distz_array;
   }
#endif
