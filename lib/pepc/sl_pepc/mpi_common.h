
#ifndef __MPI_COMMON_H__
#define __MPI_COMMON_H__


typedef FINT_TYPE_C finteger_t;
#define finteger_mpi  FINT_TYPE_MPI
#define finteger_fmt  FINT_TYPE_FMT

typedef pepckeys_slint_t slint_t;
#define slint_fmt pepckeys_slint_fmt


#define Z_MOP(_mop_)  do { _mop_ } while (0)
#define Z_NOP()       Z_MOP()


#ifdef DO_DEBUG
# define DEBUG_CMD(_cmd_)  Z_MOP(_cmd_)
#else
# define DEBUG_CMD(_cmd_)  Z_NOP()
#endif

#ifdef DO_INFO
# define INFO_CMD(_cmd_)  Z_MOP(_cmd_)
#else
# define INFO_CMD(_cmd_)  Z_NOP()
#endif

#ifdef DO_TIMING
# define TIMING_DECL(_decl_)       _decl_
# define TIMING_CMD(_cmd_)         Z_MOP(_cmd_)
#else
# define TIMING_DECL(_decl_)
# define TIMING_CMD(_cmd_)         Z_NOP()
#endif
#ifdef DO_TIMING_SYNC
# define TIMING_SYNC(_c_)          TIMING_CMD(MPI_Barrier(_c_);)
#else
# define TIMING_SYNC(_c_)          Z_NOP()
#endif
#define TIMING_START(_t_)          TIMING_CMD(((_t_) = MPI_Wtime());)
#define TIMING_STOP(_t_)           TIMING_CMD(((_t_) = MPI_Wtime() - (_t_));)
#define TIMING_STOP_ADD(_t_, _r_)  TIMING_CMD(((_r_) += MPI_Wtime() - (_t_));)

#ifdef DO_VERBOSE
# define VERBOSE_CMD(_cmd_)       Z_MOP(_cmd_)
#else
# define VERBOSE_CMD(_cmd_)       Z_NOP()
#endif


void sl_pepc_receive_stats(finteger_t nmax, int *scounts, int *sdispls, int *rcounts, int *rdispls, pepckeys_sldata0_t *work, int size, int rank, MPI_Comm comm);
void sl_pepc_border_stats(slint_t nkeys, pepckeys_slkey_t *keys, int size, int rank, MPI_Comm comm);


#endif /* __MPI_COMMON_H__ */
