/* declarations of common subroutines, etc. for use with FFTW
   self-test/benchmark program (see bench.c). */

#include "bench-user.h"
#include "fftw3.h"

#define CONCAT(prefix, name) prefix ## name
#if defined(BENCHFFT_SINGLE)
#define FFTW(x) CONCAT(fcs_fftwf_, x)
#elif defined(BENCHFFT_LDOUBLE)
#define FFTW(x) CONCAT(fcs_fftwl_, x)
#elif defined(BENCHFFT_QUAD)
#define FFTW(x) CONCAT(fcs_fftwq_, x)
#else
#define FFTW(x) CONCAT(fcs_fftw_, x)
#endif

#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */

extern FFTW(plan) mkplan(bench_problem *p, unsigned flags);
extern void initial_cleanup(void);
extern void final_cleanup(void);
extern int import_wisdom(FILE *f);
extern void export_wisdom(FILE *f);

#if defined(HAVE_THREADS) || defined(HAVE_OPENMP)
#  define HAVE_SMP
   extern int threads_ok;
#endif

#ifdef __cplusplus
}  /* extern "C" */
#endif /* __cplusplus */

