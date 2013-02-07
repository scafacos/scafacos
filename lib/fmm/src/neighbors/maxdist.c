#include "../fmm_c.h"

#ifdef FMM_GETDISTMS
   const fmm_int_t maxdist_array[] = {
#include "include/maxdist.h"
   };

   const fmm_int_t *get_maxdist()
   {
     return maxdist_array;
   }
#endif
