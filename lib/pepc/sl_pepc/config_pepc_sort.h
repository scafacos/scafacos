
#ifndef __CONFIG_PEPC_SORT_H__
#define __CONFIG_PEPC_SORT_H__


#include "fortran2c_types.h"


#if defined(FCS_ENABLE_DEBUG) || 0
# define DO_DEBUG
# define DEBUG_CMD(_cmd_)  Z_MOP(_cmd_)
#else
# define DEBUG_CMD(_cmd_)  Z_NOP()
#endif
#define DEBUG_PRINT_PREFIX  "SL_PEPC_DEBUG: "

#if defined(FCS_ENABLE_INFO) || 0
# define DO_INFO
# define INFO_CMD(_cmd_)  Z_MOP(_cmd_)
#else
# define INFO_CMD(_cmd_)  Z_NOP()
#endif
#define INFO_PRINT_PREFIX  "SL_PEPC_INFO: "

#if defined(FCS_ENABLE_TIMING) || 0
# define DO_TIMING
# define TIMING_CMD(_cmd_)  Z_MOP(_cmd_)
#else
# define TIMING_CMD(_cmd_)  Z_NOP()
#endif
#define TIMING_PRINT_PREFIX  "SL_PEPC_TIMING: "


#endif /* __CONFIG_PEPC_SORT_H__ */
