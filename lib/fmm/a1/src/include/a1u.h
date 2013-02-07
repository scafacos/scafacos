/* -*- Mode: C; c-basic-offset:4 ; -*- */
/*
 *  Copyright (C) 2010 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#if !defined A1U_H_INCLUDED
#define A1U_H_INCLUDED

#include "a1conf.h"

#if defined HAVE_STDIO_H
#include <stdio.h>
#endif /* HAVE_STDIO_H */

#if defined HAVE_STDLIB_H
#include <stdlib.h>
#endif /* HAVE_STDLIB_H */

#if defined HAVE_STDINT_H
#include <stdint.h>
#endif /* HAVE_STDINT_H */

#if defined HAVE_STRING_H
#include <string.h>
#endif /* HAVE_STRING_H */

#if defined HAVE_STRINGS_H
#include <strings.h>
#endif /* HAVE_STRINGS_H */

#if defined HAVE_UNISTD_H
#include <unistd.h>
#endif /* HAVE_UNISTD_H */

#if defined HAVE_STDARG_H
#include <stdarg.h>
#endif /* HAVE_STDARG_H */

#if defined HAVE_SYS_TYPES_H
#include <sys/types.h>
#endif /* HAVE_SYS_TYPES_H */

#if defined HAVE_SYS_STAT_H
#include <sys/stat.h>
#endif /* HAVE_SYS_STAT_H */

#if defined HAVE_TIME_H
#include <time.h>
#endif /* HAVE_TIME_H */

#if defined HAVE_SYS_TIME_H
#include <sys/time.h>
#endif /* HAVE_SYS_TIME_H */

#if defined HAVE_ERRNO_H
#include <errno.h>
#endif /* HAVE_ERRNO_H */

#if !defined HAVE_MACRO_VA_ARGS
#error "VA_ARGS support is required"
#endif /* HAVE_MACRO_VA_ARGS */

/* FIXME: FUNC_ENTER/EXIT can be used for profiling in the future */

#define A1U_FUNC_ENTER(...)
#define A1U_FUNC_EXIT(...)

#if defined HAVE__FUNC__
#define A1U_FUNC __func__
#elif defined HAVE_CAP__FUNC__
#define A1U_FUNC __FUNC__
#elif defined HAVE__FUNCTION__
#define A1U_FUNC __FUNCTION__
#endif

#if defined __FILE__ && defined A1U_FUNC
    #define A1U_error_printf(...)                                          \
        {                                                                   \
            fprintf(stderr, "%s (%s:%d): ", A1U_FUNC, __FILE__, __LINE__);  \
            fprintf(stderr, __VA_ARGS__);                                   \
            fflush(stderr);                                                 \
        }
#elif defined __FILE__
    #define A1U_error_printf(...)                               \
        {                                                        \
            fprintf(stderr, "%s (%d): ", __FILE__, __LINE__);    \
            fprintf(stderr, __VA_ARGS__);                        \
            fflush(stderr);                                      \
        }
#else
    #define A1U_error_printf(...)                                          \
        {                                                                   \
            fprintf(stderr, __VA_ARGS__);                                   \
            fflush(stderr);                                                 \
        }
#endif

#define A1U_output_printf(...)                                         \
    {                                                                   \
        fprintf(stdout, __VA_ARGS__);                                   \
        fflush(stdout);                                                 \
    }

#define A1U_ASSERT_ABORT(x, ...)                                        \
    {                                                                   \
        if (!(x)) {                                                     \
            A1U_error_printf(__VA_ARGS__);                              \
            assert(0);                                                  \
        }                                                               \
    }

#define A1U_ASSERT(x, status)                                           \
    {                                                                   \
        if (!(x)) {                                                     \
            A1U_ERR_SETANDJUMP(status, A1_ERROR,                        \
                               "assert (%s) failed\n", #x);             \
        }                                                               \
    }

#define A1U_WARNING(status, ...)                                      \
    {                                                                   \
        if (status) {                                                   \
            A1U_error_printf(__VA_ARGS__);                              \
        }                                                               \
    }

#define A1U_ERR_ABORT(status, ...)                                      \
    {                                                                   \
        if (status) {                                                   \
            A1U_error_printf(__VA_ARGS__);                              \
            assert(0);                                                  \
        }                                                               \
    }

#define A1U_ERR_POP(status, ...)                                        \
    {                                                                   \
        if (status) {                                                   \
            A1U_error_printf(__VA_ARGS__);                              \
            goto fn_fail;                                               \
        }                                                               \
    }

#define A1U_ERR_SETANDJUMP(status, error, ...)                          \
    {                                                                   \
        status = error;                                                 \
        A1U_ERR_POP(status, __VA_ARGS__);                               \
    }

#define A1U_ERR_CHKANDJUMP(status, chk, error, ...)                     \
    {                                                                   \
        if ((chk))                                                      \
            A1U_ERR_SETANDJUMP(status, error, __VA_ARGS__);             \
    }

#define A1U_ERR_POPANDSTMT(status, stmt, ... )                          \
    {                                                                   \
        if (status) {                                                   \
            A1U_error_printf(__VA_ARGS__);                              \
            stmt;                                                       \
        }                                                               \
    }

#ifdef A1_DEBUG
#define A1U_DEBUG_PRINT(args...)                                  \
do {                                                              \
    int __my_rank;                                                \
    __my_rank = A1_Process_rank(A1_GROUP_WORLD);                  \
    fprintf(stderr, "Debug Message from [%d] :", __my_rank);      \
    fprintf(stderr, args);                                        \
    fflush(stderr);                                               \
} while (0)
#else
#define A1U_DEBUG_PRINT(args...)
#endif

int A1U_Put_memcpy(void* src,
                   void* dst,
                   int bytes);

int A1U_PutS_memcpy(int stride_level,
                    int *block_sizes,
                    void* source_ptr,
                    int *src_stride_ar,
                    void* target_ptr,
                    int *trg_stride_ar);

int A1U_PutV_memcpy(A1_iov_t *iov_ar,
                    int ar_len);

int A1U_Get_memcpy(void* src,
                   void* dst,
                   int bytes);

int A1U_GetS_memcpy(int stride_level,
                    int *block_sizes,
                    void* source_ptr,
                    int *src_stride_ar,
                    void* target_ptr,
                    int *trg_stride_ar);

int A1U_GetV_memcpy(A1_iov_t *iov_ar, int ar_len);

int A1U_Acc_memcpy(void* source_ptr,
                   void* target_ptr,
                   int bytes,
                   A1_datatype_t a1_type,
                   void* scaling);

int A1U_AccS_memcpy(int stride_level,
                    int *block_sizes,
                    void* source_ptr,
                    int *src_stride_ar,
                    void* target_ptr,
                    int *trg_stride_ar,
                    A1_datatype_t a1_type,
                    void* scaling);

int A1U_AccV_memcpy(A1_iov_t *iov_ar,
                    int ar_len,
                    A1_datatype_t a1_type,
                    void* scaling);

#endif /* A1U_H_INCLUDED */
