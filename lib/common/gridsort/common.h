/*
  Copyright (C) 2011, 2012, 2013 Michael Hofmann
  
  This file is part of ScaFaCoS.
  
  ScaFaCoS is free software: you can redistribute it and/or modify
  it under the terms of the GNU Lesser Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  ScaFaCoS is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU Lesser Public License for more details.
  
  You should have received a copy of the GNU Lesser Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/


#if defined(FCS_ENABLE_DEBUG_GRIDSORT)
# define DO_DEBUG
# define DEBUG_CMD(_cmd_)  Z_MOP(_cmd_)
#else
# define DEBUG_CMD(_cmd_)  Z_NOP()
#endif
#define DEBUG_PRINT_PREFIX  "GRIDSORT_DEBUG: "

#if defined(FCS_ENABLE_INFO_GRIDSORT)
# define DO_INFO
# define INFO_CMD(_cmd_)  Z_MOP(_cmd_)
#else
# define INFO_CMD(_cmd_)  Z_NOP()
#endif
#define INFO_PRINT_PREFIX  "GRIDSORT_INFO: "

#if defined(FCS_ENABLE_TIMING_GRIDSORT)
# define DO_TIMING
# define TIMING_CMD(_cmd_)  Z_MOP(_cmd_)
#else
# define TIMING_CMD(_cmd_)  Z_NOP()
#endif
#define TIMING_PRINT_PREFIX  "GRIDSORT_TIMING: "

#define GRIDSORT_INDEX_IS_VALID(_x_)       ((_x_) >= 0)
#define GRIDSORT_INDEX_VAL_PROC(_proc_)    (((fcs_gridsort_index_t) (_proc_)) << 32)
#define GRIDSORT_INDEX_VAL_POS(_pos_)      ((fcs_gridsort_index_t) (_pos_))
#define GRIDSORT_INDEX_VAL(_proc_, _pos_)  (GRIDSORT_INDEX_VAL_PROC(_proc_) + GRIDSORT_INDEX_VAL_POS(_pos_))
#define GRIDSORT_INDEX_GET_PROC(_x_)       ((_x_) >> 32)
#define GRIDSORT_INDEX_GET_POS(_x_)        ((_x_) & 0x00000000FFFFFFFFLL)

#define GRIDSORT_INDEX_STR         "(%lld,%lld)"
#define GRIDSORT_INDEX_PARAM(_x_)  (GRIDSORT_INDEX_IS_VALID(_x_)?GRIDSORT_INDEX_GET_PROC(_x_):-1), (GRIDSORT_INDEX_IS_VALID(_x_)?GRIDSORT_INDEX_GET_POS(_x_):-1)

#define DO_TIMING_SYNC

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

#if 0
# define GRIDSORT_FRONT_TPROC_RANK_CACHE
#endif

#if 1
# define GRIDSORT_PROCLIST
# define GRIDSORT_FRONT_PROCLIST
# define GRIDSORT_BACK_PROCLIST
# define GRIDSORT_RESORT_PROCLIST
# undef GRIDSORT_PROCLIST_VERIFY
#endif

#if 0
# define GRIDSORT_FRONT_TPROC_EXDEF
# define GRIDSORT_BACK_TPROC_EXDEF
#endif
