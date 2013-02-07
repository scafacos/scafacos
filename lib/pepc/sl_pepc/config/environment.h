
#ifdef FCS_ENABLE_TIMING
# define SL_TIMING
# define SL_TIMING_PRINT_PREFIX  "SL_PEPC_TIMING: "
# define SL_TIMING_PRINT(_i_, _s_, _n_, _v_, _r_)  Z_MOP(if ((_r_) == 0) SL_TIMING_PRINT_DEFAULT(_i_, _s_, _n_, _v_, _r_);)
#endif

#ifdef FCS_ENABLE_DEBUG
# define SLDEBUG_OUTPUT  5
#endif
