/*
 *  Copyright (C) 2011, 2012, 2013, 2014, 2015, 2016 Michael Hofmann
 *  
 *  This file is part of ScaFaCoS.
 *  
 *  ScaFaCoS is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *  
 *  ScaFaCoS is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser Public License for more details.
 *  
 *  You should have received a copy of the GNU Lesser Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *  

 *  
 *  SL - Sorting Library, michael <dot> hofmann <at> informatik <dot> tu-chemnitz <dot> de
 */


/* _XOPEN_SOURCE >= 500 required for nftw, lstat */
/* _XOPEN_SOURCE >= 600 required for vsscanf */
#define _XOPEN_SOURCE 600

#include "z_pack.h"


#ifdef Z_PACK_MPI

# ifdef HAVE_SPI_KERNEL_INTERFACE_H
#  include <spi/kernel_interface.h>
# endif
# ifdef HAVE_COMMON_BGP_PERSONALITY_H
#  include <common/bgp_personality.h>
# endif
# ifdef HAVE_COMMON_BGP_PERSONALITY_INLINES_H
#  include <common/bgp_personality_inlines.h>
# endif
# ifdef HAVE_MPIX_H
#  include <mpix.h>
# endif


void z_mpi_remap_cart_topology(int from_ndims, int *from_dims, int *from_torus, int *from_pos, int to_ndims, int *to_dims, int *to_torus, int *to_pos) /* z_proto, z_func z_mpi_remap_cart_topology */
{
  int i;
  int torus = 0;
  const int min_ndims = (from_ndims < to_ndims)?from_ndims:to_ndims;

  
  for (i = 0; i < min_ndims; ++i)
  {
    to_dims[i] = from_dims[i];
    to_torus[i] = from_torus[i];
    to_pos[i] = from_pos[i];
    
    if (to_torus[i]) torus = 1;
  }

  for (; i < to_ndims; ++i)
  {
    to_dims[i] = 1;
    to_torus[i] = torus;
    to_pos[i] = 0;
  }

  for (; i < from_ndims; ++i)
  {
    to_dims[to_ndims - 1] *= from_dims[i];
    if (from_torus[i]) to_torus[to_ndims - 1] = 1;
    to_pos[to_ndims - 1] *= from_dims[i]; to_pos[to_ndims - 1] += from_pos[i];
  }
}


#define MAX_CART_NDIMS  6


void z_mpi_get_cart_topology(int *ndims, int *dims, int *torus, int *pos) /* z_proto, z_func z_mpi_get_cart_topology */
{
  int _ndims, _dims[MAX_CART_NDIMS], _torus[MAX_CART_NDIMS], _pos[MAX_CART_NDIMS];

#if defined(HAVE__BGP_PERSONALITY_T)

  /* BlueGene/P */

  _BGP_Personality_t personality;


  _ndims = 4;

  if (dims == NULL || torus == NULL || pos == NULL) goto exit_ndims_only;

  Kernel_GetPersonality(&personality, sizeof(personality));

  _dims[0] = personality.Network_Config.Xnodes;
  _dims[1] = personality.Network_Config.Ynodes;
  _dims[2] = personality.Network_Config.Znodes;

  _torus[0] = BGP_Personality_isTorusX(&personality);
  _torus[1] = BGP_Personality_isTorusY(&personality);
  _torus[2] = BGP_Personality_isTorusZ(&personality);

  _pos[0] = personality.Network_Config.Xcoord;
  _pos[1] = personality.Network_Config.Ycoord;
  _pos[2] = personality.Network_Config.Zcoord;

  switch (personality.Kernel_Config.ProcessConfig)
  {
    case _BGP_PERS_PROCESSCONFIG_SMP:
      _dims[3] = 1;
      break;
    case _BGP_PERS_PROCESSCONFIG_VNM:
      _dims[3] = 4;
      break;
    case _BGP_PERS_PROCESSCONFIG_2x2:
      _dims[3] = 2;
      break;
    default:
      _dims[3] = 1;
      break;
  }

  _torus[3] = (_torus[0] || _torus[1] || _torus[2]);

  _pos[3] = Kernel_PhysicalProcessorID();

#elif defined(HAVE_MPIX_HARDWARE_T)

  /* BlueGene/Q */

  int i;
  MPIX_Hardware_t hw;


  _ndims = MPIX_TORUS_MAX_DIMS + 1;

  if (dims == NULL || torus == NULL || pos == NULL) goto exit_ndims_only;

  MPIX_Hardware(&hw);

  _torus[MPIX_TORUS_MAX_DIMS] = 0;

  for (i = 0; i < MPIX_TORUS_MAX_DIMS; ++i)
  {
    _dims[i] = hw.Size[i];
    _torus[i] = hw.isTorus[i]?1:0;
    _pos[0] = hw.Coords[i];

    if (_torus[i]) _torus[MPIX_TORUS_MAX_DIMS] = 1;
  }

  _dims[MPIX_TORUS_MAX_DIMS] = hw.ppn;

  _pos[MPIX_TORUS_MAX_DIMS] = hw.coreID;

#else

  /* MPI */

  MPI_Comm comm;
  int size, rank;


  _ndims = 3;

  if (dims == NULL || torus == NULL || pos == NULL) goto exit_ndims_only;

  comm = MPI_COMM_WORLD;
  MPI_Comm_size(comm, &size);
  MPI_Comm_rank(comm, &rank);

  _dims[0] = 0;
  _dims[1] = 0;
  _dims[2] = 0;

  MPI_Dims_create(size, 3, _dims);

  _pos[2] = (rank / (1))                       % _dims[2];
  _pos[1] = (rank / (1 * _dims[2]))            % _dims[1];
  _pos[0] = (rank / (1 * _dims[2] * _dims[1])) % _dims[0];

  _torus[0] = 0;
  _torus[1] = 0;
  _torus[2] = 0;
#endif

  if (*ndims <= 0) *ndims = _ndims;

  z_mpi_remap_cart_topology(_ndims, _dims, _torus, _pos, *ndims, dims, torus, pos);

exit_ndims_only:

  *ndims = _ndims;
}


void z_mpi_get_grid4d(int *dims, int *pos) /* z_proto, z_func z_mpi_get_grid4d */
{
  int torus[MAX_CART_NDIMS], ndims = 4;

  z_mpi_get_cart_topology(&ndims, dims, torus, pos);
}

#endif /* Z_PACK_MPI */


#ifdef Z_PACK_ALLOC

#endif /* Z_PACK_ALLOC */


#ifdef Z_PACK_DEBUG

/* z_var z_notice_fstream z_error_fstream z_debug_fstream */
FILE *z_notice_fstream = NULL;
FILE *z_error_fstream = NULL;
FILE *z_debug_fstream = NULL;

#endif /* Z_PACK_DEBUG */


#if 0
#ifdef Z_PACK_TIME

#ifndef Z_PACK_MPI
double z_time_wtime() /* z_proto, z_func z_time_wtime */
{
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return (double) tv.tv_sec + (tv.tv_usec / 1000000.0);
}
#endif

#endif /* Z_PACK_TIME */
#endif


#ifdef Z_PACK_RANDOM

#ifdef Z_RANDOM_REQUIRED

/* "Even Quicker Generator" from Numerical Recipes in C, §7.1 */

unsigned long z_srandom_seed = 0; /* z_var z_srandom_seed */

void z_srandom(unsigned long seed) /* z_proto, z_func z_srandom */
{
  z_srandom_seed = seed;
}

unsigned long z_random() /* z_proto, z_func z_random */
{
  z_srandom_seed = 1664525L * z_srandom_seed + 1013904223L;
  return (z_srandom_seed &= 0xFFFFFFFFL);
}

#endif


void z_srandom64(unsigned long seed) /* z_proto, z_func z_srandom64 */
{
  z_srand(seed);
}


long long z_random64() /* z_proto, z_func z_random64 */
{
  long long x;


  x =  ((long long) z_rand()) << 32;
  x |= ((long long) z_rand()) << 1;

  x |= ((long long) z_rand()) & 0x1;
  
  return x;
}


long long z_random64_minmax(long long min, long long max) /* z_proto, z_func z_random64_minmax */
{
  long double r;
  long long x;


  r = (long double) z_rand() / (long double) Z_RAND_MAX;
  
  x = (long long) (((long double) (max - min) * r) + 0.5);
  
  x += min;

  if (x > max) x = max;
  
  return x;
}


unsigned long long z_random64u() /* z_proto, z_func z_random64u */
{
  unsigned long long x;
  

  x =  ((unsigned long long) z_rand()) << 32;
  x |= ((unsigned long long) z_rand()) << 1;

  x <<= 1;
  x |= ((unsigned long long) z_rand()) & 0x3;
  
  return x;
}


unsigned long long z_random64u_minmax(unsigned long long min, unsigned long long max) /* z_proto, z_func z_random64u_minmax */
{
  long double r;
  unsigned long long x;


  r = (long double) z_rand() / (long double) Z_RAND_MAX;
  
  x = (long long) (((long double) (max - min) * r) + 0.5);
  
  x += min;

  if (x > max) x = max;
  
  return x;
}



#if defined(HAVE_ESSL_H)
# include <essl.h>
#elif defined(HAVE_T4C_H)
# include <trng4c.h>
#endif


#if defined(HAVE_ESSL_H)
double z_nrandom_seed_essl = 1.0; /* z_var z_nrandom_seed_essl */
#endif


void z_nrandom_seed(unsigned long s) /* z_proto, z_func z_nrandom_seed */
{
#if defined(HAVE_ESSL_H)
  z_nrandom_seed_essl = s + 1;
#elif defined(HAVE_T4C_H)
  t4c_dnrand_seed(s);
#endif
}


double z_nrandom() /* z_proto, z_func z_nrandom */
{
  double r[2] = { 0, 0 };

#if defined(HAVE_ESSL_H)
/*  printf("dnrand\n");*/
  dnrand(&z_nrandom_seed_essl, 2, r, NULL, 0);
#elif defined(HAVE_T4C_H)
/*  printf("t4c_dnrand\n");*/
  t4c_dnrand(r);
#endif

/*  printf("z_nrandom: %f\n", r[0]);*/

  return r[0];
}


#if defined(HAVE_ESSL_H)
double z_urandom_seed_essl = 1.0;  /* z_var z_urandom_seed_essl */
#endif


void z_urandom_seed(unsigned long s) /* z_proto, z_func z_urandom_seed */
{
#if defined(HAVE_ESSL_H)
  z_urandom_seed_essl = s + 1;
#elif defined(HAVE_T4C_H)
  t4c_durand_seed(s);
#else
  z_srand(s);
#endif
}


double z_urandom() /* z_proto, z_func z_urandom */
{
  double r[1] = { 0 };

#if defined(HAVE_ESSL_H)
/*  printf("durand\n");*/
  durand(&z_urandom_seed_essl, 1, r);
#elif defined(HAVE_T4C_H)
/*  printf("t4c_durand\n");*/
  t4c_durand(r);
#else
/*  printf("random() / RAND_MAX\n");*/
  r[0] = (double) z_rand() / (double) Z_RAND_MAX;
#endif

/*  printf("z_urandom: %f\n", r[0]);*/

  return r[0];
}

#endif /* Z_PACK_RANDOM */


#ifdef Z_PACK_DIGEST

#ifdef HAVE_GCRYPT_H
# include <gcrypt.h>
#endif


z_int_t z_digest_sum_buffer(const void *buffer, z_int_t length, void *sum) /* z_proto, z_func z_digest_sum_buffer */
{
#if defined(HAVE_GCRYPT_H)

  if (!buffer || !sum) return gcry_md_get_algo_dlen(GCRY_MD_CRC32);

  gcry_md_hash_buffer(GCRY_MD_CRC32, sum, buffer, length);

#elif defined (Z_PACK_CRC32)

  if (!buffer || !sum) return sizeof(z_crc32_t);

  *((z_crc32_t *) sum) = z_crc32_buffer(buffer, length);
  
#endif

/*  printf("z: %.8X\n", *((z_crc32_t *) sum));*/

  return 0;
}


#ifdef HAVE_GCRYPT_H
int z_digest_hash_gcrypt_algo = GCRY_MD_MD5; /* z_var, z_digest_hash_gcrypt_algo */
#endif


void z_digest_hash_open(void **hdl) /* z_proto, z_func z_digest_hash_open */
{
#ifdef HAVE_GCRYPT_H
  gcry_md_hd_t *gcry_hdl;
  
  gcry_hdl = z_alloc(1, sizeof(gcry_md_hd_t));

  gcry_md_open(gcry_hdl, z_digest_hash_gcrypt_algo, 0);

  *hdl = gcry_hdl;
#endif
}


void z_digest_hash_close(void *hdl) /* z_proto, z_func z_digest_hash_close */
{
#ifdef HAVE_GCRYPT_H
  gcry_md_hd_t *gcry_hdl = hdl;

  gcry_md_close(*gcry_hdl);

  z_free(gcry_hdl);
#endif
}


void z_digest_hash_write(void *hdl, const void *buffer, z_int_t length) /* z_proto, z_func z_digest_hash_write */
{
#ifdef HAVE_GCRYPT_H
  gcry_md_hd_t *gcry_hdl = hdl;

  gcry_md_write(*gcry_hdl, buffer, length);
#endif
}


z_int_t z_digest_hash_read(void *hdl, void *hash) /* z_proto, z_func z_digest_hash_read */
{
#ifdef HAVE_GCRYPT_H
  gcry_md_hd_t *gcry_hdl = hdl;

  z_int_t dlen = gcry_md_get_algo_dlen(z_digest_hash_gcrypt_algo);

  if (!hdl || !hash) return dlen;

  memcpy(hash, gcry_md_read(*gcry_hdl, 0), dlen);
#endif

  return 0;
}

#endif /* Z_PACK_DIGEST */


#if defined(Z_PACK_GMP) && defined(HAVE_GMP_H)

#include <gmp.h>


void z_gmp_mpz_set_ull(mpz_t z, unsigned long long v) /* z_proto, z_func z_gmp_mpz_set_ull */
{
/*  printf("z_gmp_mpz_set_ull\n");*/


  if (sizeof(unsigned long long) != sizeof (unsigned long))
  {
    mpz_set_ui(z, (unsigned long) (v >> 32));
    mpz_mul_2exp(z, z, 32);
    mpz_add_ui(z, z, (unsigned long) (v & 0x00000000FFFFFFFFULL));

  } else mpz_set_ui(z, (unsigned long) v);
}


void z_gmp_mpz_set_sll(mpz_t z, long long v) /* z_proto, z_func z_gmp_mpz_set_sll */
{
  unsigned long long uv = (v < 0)?-v:v;


/*  printf("z_gmp_mpz_set_sll\n");*/

  if (sizeof(long long) != sizeof (long))
  {
    mpz_set_ui(z, (unsigned long) (uv >> 32));
    mpz_mul_2exp(z, z, 32);
    mpz_add_ui(z, z, (long) (uv & 0x00000000FFFFFFFFLL));
    
    if (v < 0) mpz_neg(z, z);

  } else mpz_set_si(z, (long) v);
}


unsigned long long z_gmp_mpz_get_ull(mpz_t z) /* z_proto, z_func z_gmp_mpz_get_ull */
{
  mpz_t t;
  unsigned long long r = 0;


/*  printf("z_gmp_mpz_get_ull\n");*/

  if (sizeof(unsigned long long) != sizeof (unsigned long))
  {
    mpz_init(t);

    mpz_tdiv_q_2exp(t, z, 32);
    r = mpz_get_ui(t);
    mpz_tdiv_r_2exp(t, z, 32);
    r = (r << 32) + mpz_get_ui(t);

    mpz_clear(t);

  } else r = mpz_get_ui(z);

  return r;
}


long long z_gmp_mpz_get_sll(mpz_t z) /* z_proto, z_func z_gmp_mpz_get_sll */
{
  long long r = 0;
  mpz_t t;
  int s;


/*  printf("z_gmp_mpz_get_sll\n");*/

  if (sizeof(long long) != sizeof (long))
  {
    s = mpz_sgn(z);

    if (s < 0) mpz_neg(z, z);

    mpz_init(t);

    mpz_tdiv_q_2exp(t, z, 32);
    r = mpz_get_si(t);
    mpz_tdiv_r_2exp(t, z, 32);
    r = (r << 32) + mpz_get_ui(t);

    mpz_clear(t);

    if (s < 0) mpz_neg(z, z);

    r *= s;

  } else r = mpz_get_si(z);

  return r;
}

#endif /* Z_PACK_GMP && HAVE_GMP_H */


#ifdef Z_PACK_FS

#if HAVE_SYS_STAT_H && HAVE_SYS_TYPES_H

#include <sys/stat.h>
#include <sys/types.h>


z_int_t z_fs_exists(const char *pathname) /* z_proto, z_func z_fs_exists */
{
  struct stat s;


  return (stat(pathname, &s) == 0);
}


z_int_t z_fs_is_directory(const char *pathname) /* z_proto, z_func z_fs_is_directory */
{
  struct stat s;


  if (stat(pathname, &s) != 0) return 0;

  return S_ISDIR(s.st_mode)?1:0;
}


z_int_t z_fs_is_file(const char *pathname) /* z_proto, z_func z_fs_is_file */
{
  struct stat s;


  if (stat(pathname, &s) != 0) return 0;

  return S_ISREG(s.st_mode)?1:0;
}


z_int_t z_fs_is_link(const char *pathname) /* z_proto, z_func z_fs_is_link */
{
  struct stat s;


  if (lstat(pathname, &s) != 0) return 0;

  return S_ISLNK(s.st_mode)?1:0;
}


z_int_t z_fs_get_file_size(const char *pathname) /* z_proto, z_func z_fs_get_file_size */
{
  struct stat s;


  if (stat(pathname, &s) != 0) return -1;

  return (z_int_t) s.st_size;
}


#if HAVE_UNISTD_H

#include <unistd.h>


z_int_t z_fs_get_link_target(const char *pathname, char *target, z_int_t size) /* z_proto, z_func z_fs_get_link_target */
{
  ssize_t n = readlink(pathname, target, size);

  if (n < 0) return -1;

  if (n < size) target[n] = '\0';

  return n;
}

#endif /* HAVE_UNISTD_H */


z_int_t z_fs_stat(const char *pathname, z_fs_stat_t *stat)
{
  struct stat s;


  if (lstat(pathname, &s) != 0) return 0;

  stat->is_directory = S_ISDIR(s.st_mode)?1:0;
  stat->is_file = S_ISREG(s.st_mode)?1:0;
  stat->is_link = S_ISLNK(s.st_mode)?1:0;
  stat->file_size = (z_int_t) s.st_size;

  return 1;
}


z_int_t z_fs_mkdir(const char *pathname) /* z_proto, z_func z_fs_mkdir */
{
  const mode_t mode = 0777;


  if (z_fs_is_directory(pathname)) return 1;

  return (mkdir(pathname, mode) == 0);
}

#endif /* HAVE_SYS_STAT_H && HAVE_SYS_TYPES_H */


#if (HAVE_STDLIB_H && HAVE_STRING_H) || STDC_HEADERS

#include <stdlib.h>
#include <string.h>


z_int_t z_fs_mkdir_p(const char *pathname) /* z_proto, z_func z_fs_mkdir_p */
{
  z_int_t ret = 1;
  char *p, *pb, *pe;


  pb = p = strdup(pathname);

  while (ret && pb != NULL)
  {
    pe = strchr(pb, '/');

    if (pe != pb)
    {
      if (pe) *pe = '\0';

      ret = z_fs_mkdir(p);

      if (pe) *pe = '/';
    }

    pb = pe;

    if (pb) ++pb;
  }

  free(p);

  return ret;
}

#endif /* HAVE_STDLIB_H && HAVE_STRING_H */


#if HAVE_STDIO_H || STDC_HEADERS

#include <stdio.h>


z_int_t z_fs_rm(const char *pathname) /* z_proto, z_func z_fs_rm */
{
  return (remove(pathname) == 0);
}


#if HAVE_FTW_H

#include <ftw.h>

/*
static int display_info_fn(const char *fpath, const struct stat *sb, int tflag, struct FTW *ftwbuf)
{
  printf("%-3s %2d %7d   %-40s %d %s\n",
    (tflag == FTW_D)?"d":
    (tflag == FTW_DNR)?"dnr":
    (tflag == FTW_DP)?"dp":
    (tflag == FTW_F)?(
      S_ISBLK(sb->st_mode)?"f b":
      S_ISCHR(sb->st_mode)?"f c":
      S_ISFIFO(sb->st_mode)?"f p":
      S_ISREG(sb->st_mode)?"f r":
      S_ISSOCK(sb->st_mode)?"f s":"f ?"):
    (tflag == FTW_NS)?"ns":
    (tflag == FTW_SL)?"sl":
    (tflag == FTW_SLN)?"sln":"?",
    ftwbuf->level,
    (int) sb->st_size,
    fpath,
    ftwbuf->base,
    fpath + ftwbuf->base);
  return 0;
}
*/


static int rm_r_fn(const char *fpath, const struct stat *sb, int tflag, struct FTW *ftwbuf)
{
/*  printf("remove: %s\n", fpath);*/
  remove(fpath);

  return 0;
}


z_int_t z_fs_rm_r(const char *pathname) /* z_proto, z_func z_fs_rm_r */
{
  return (nftw(pathname, rm_r_fn, 20, FTW_PHYS|FTW_DEPTH) == 0);
}

#endif /* HAVE_FTW_H */

#endif /* HAVE_STDIO_H */

#endif /* Z_PACK_FS */


#ifdef Z_PACK_STDIO

#if ((HAVE_STDIO_H && HAVE_STRING_H) || STDC_HEADERS) && HAVE_STDARG_H

#include <stdio.h>
#include <string.h>
#include <stdarg.h>


int z_snscanf(const char *str, size_t size, const char *format, ...) /* z_proto, z_func z_snscanf */
{
  int n;
  char str2[size + 1];


  strncpy(str2, str, size);
  str2[size] = '\0';

  va_list argp;
  va_start(argp, format);

  n = vsscanf(str2, format, argp);

  va_end(argp);

  return n;
}

#endif /* HAVE_STDIO_H && HAVE_STRING_H && HAVE_STDARG_H */

#endif /* Z_PACK_STDIO */
