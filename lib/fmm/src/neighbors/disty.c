#include "../fmm_c.h"

#ifdef FMM_GETDIST
   const fmm_int_t disty_array[] = {
#include "include/disty.h"
   };

   const fmm_int_t *get_disty()
   {
     return disty_array;
   }
#endif
