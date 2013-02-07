#include "../fmm_c.h"

const fmm_int_t casejump_array[] = {
#include "include/casejump.h"
};

const fmm_int_t *get_cjp512()
{
  return casejump_array;
}
