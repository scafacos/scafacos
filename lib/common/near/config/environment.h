
#ifdef FCS_ENABLE_TIMING_NEAR
# define SL_TIMING
# define SL_TIMING_PRINT_PREFIX  "NEAR_TIMING: "
# define SL_TIMING_PRINT(_i_, _s_, _n_, _v_, _r_)  Z_MOP(if ((_r_) == 0) SL_TIMING_PRINT_DEFAULT(_i_, _s_, _n_, _v_, _r_);)
#endif

#ifdef FCS_ENABLE_DEBUG_NEAR
# define SLDEBUG_OUTPUT  5
#endif
