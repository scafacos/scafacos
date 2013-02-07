
#ifndef __FORTRAN2C_TYPES_H__
#define __FORTRAN2C_TYPES_H__


#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#ifdef DEFAULT_FORTRAN_INTEGER_KIND
# define FINT_DEFAULT  DEFAULT_FORTRAN_INTEGER_KIND
#endif

#ifdef DEFAULT_FORTRAN_REAL_KIND
# define FREAL_DEFAULT  DEFAULT_FORTRAN_REAL_KIND
#endif


#if defined(HAVE_INTTYPES_H) && (!defined(HAVE_MPI) || (((MPI_VERSION >= 2) && (MPI_SUBVERSION >= 2)) || defined(HAVE_MPI_V22_TYPES)))

# include <inttypes.h>

# define FINT4_TYPE_C      int32_t
# define FINT4_TYPE_MPI    MPI_INT32_T
# define FINT4_TYPE_FMT    PRId32

# define FINT8_TYPE_C      int64_t
# define FINT8_TYPE_MPI    MPI_INT64_T
# define FINT8_TYPE_FMT    PRId64

#else

# define FINT4_TYPE_C      int
# define FINT4_TYPE_MPI    MPI_INT
# define FINT4_TYPE_FMT    "d"

# define FINT8_TYPE_C      long long
# define FINT8_TYPE_MPI    MPI_LONG_LONG
# define FINT8_TYPE_FMT    "lld"

#endif

#define FREAL4_TYPE_C      float
#define FREAL4_TYPE_MPI    MPI_FLOAT
#define FREAL4_TYPE_FMT    "f"

#define FREAL8_TYPE_C      double
#define FREAL8_TYPE_MPI    MPI_DOUBLE
#define FREAL8_TYPE_FMT    "f"

#define FREAL16_TYPE_C     long double
#define FREAL16_TYPE_MPI   MPI_LONG_DOUBLE
#define FREAL16_TYPE_FMT   "Lf"


#ifndef FINT_DEFAULT
# define FINT_DEFAULT  8
#endif

#if FINT_DEFAULT == 4
# define FINT_TYPE_C       FINT4_TYPE_C
# define FINT_TYPE_MPI     FINT4_TYPE_MPI
# define FINT_TYPE_FMT     FINT4_TYPE_FMT
#elif FINT_DEFAULT == 8
# define FINT_TYPE_C       FINT8_TYPE_C
# define FINT_TYPE_MPI     FINT8_TYPE_MPI
# define FINT_TYPE_FMT     FINT8_TYPE_FMT
#else
# error "FINT_DEFAULT value not supported!"
#endif

#ifndef FREAL_DEFAULT
# define FREAL_DEFAULT  8
#endif

#if FREAL_DEFAULT == 4
# define FREAL_TYPE_C       FREAL4_TYPE_C
# define FREAL_TYPE_MPI     FREAL4_TYPE_MPI
# define FREAL_TYPE_FMT     FREAL4_TYPE_FMT
#elif FREAL_DEFAULT == 8
# define FREAL_TYPE_C       FREAL8_TYPE_C
# define FREAL_TYPE_MPI     FREAL8_TYPE_MPI
# define FREAL_TYPE_FMT     FREAL8_TYPE_FMT
#elif FREAL_DEFAULT == 16
# define FREAL_TYPE_C       FREAL16_TYPE_C
# define FREAL_TYPE_MPI     FREAL16_TYPE_MPI
# define FREAL_TYPE_FMT     FREAL16_TYPE_FMT
#else
# error "FREAL_DEFAULT value not supported!"
#endif


#endif /* __FORTRAN2C_TYPES_H__ */
