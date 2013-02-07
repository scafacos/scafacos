#include "../fmm_c.h"

#ifdef FMM_GETDISTMS
   const fmm_int_t dist3dsquared_array[] = {
#include "include/dist3dsquared.h"
   };

   const fmm_int_t *get_dist3dsquared()
   {
     return dist3dsquared_array;
   }
#endif
