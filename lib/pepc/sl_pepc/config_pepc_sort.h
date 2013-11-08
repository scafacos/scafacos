
#ifndef __CONFIG_PEPC_SORT_H__
#define __CONFIG_PEPC_SORT_H__


#include "fortran2c_types.h"


#if defined(FCS_ENABLE_DEBUG_SL_PEPC) || 0
# define DO_DEBUG
#endif
#define DEBUG_PRINT_PREFIX  "SL_PEPC_DEBUG: "

#if defined(FCS_ENABLE_INFO_SL_PEPC) || 0
# define DO_INFO
#endif
#define INFO_PRINT_PREFIX  "SL_PEPC_INFO: "

#if defined(FCS_ENABLE_TIMING_SL_PEPC) || 0
# define DO_TIMING
#endif
#define TIMING_PRINT_PREFIX  "SL_PEPC_TIMING: "


/* FCS mode? */
#if defined(fcs_int) && defined(fcs_float)
# define SL_PEPC_PREFIX  fcs_
#endif


#define SORT_RHIGH   -1
#define SORT_RLOW    -1
#define SORT_RWIDTH  -1

#define PART_RHIGH   62
#define PART_RLOW    -1
#define PART_RWIDTH  3


#endif /* __CONFIG_PEPC_SORT_H__ */
