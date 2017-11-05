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


#ifndef __SL_TYPES_H__
#define __SL_TYPES_H__


#ifdef SL_USE_MPI
# include <mpi.h>
#endif


/* sl_type slint_t slint */
typedef sl_int_type_c slint_t, slint;  /* deprecated 'slint' */

#define slint_fmt   sl_int_type_fmt    /* sl_macro */

/* sl_type slindex_t */
typedef sl_index_type_c slindex_t;

#define sindex_fmt  sl_index_type_fmt  /* sl_macro */

/* sl_type slkey_t */
typedef sl_key_type_c slkey_t;

/* sl_type slkey_pure_t slpkey_t */
typedef sl_key_pure_type_c slkey_pure_t, slpkey_t;

/* DATAX_TEMPLATE_BEGIN */
/* sl_type sldata0_t */
#ifdef sl_data0_type_c
typedef sl_data0_type_c sldata0_t;
#endif
/* sl_type sldata1_t */
#ifdef sl_data1_type_c
typedef sl_data1_type_c sldata1_t;
#endif
/* sl_type sldata2_t */
#ifdef sl_data2_type_c
typedef sl_data2_type_c sldata2_t;
#endif
/* sl_type sldata3_t */
#ifdef sl_data3_type_c
typedef sl_data3_type_c sldata3_t;
#endif
/* sl_type sldata4_t */
#ifdef sl_data4_type_c
typedef sl_data4_type_c sldata4_t;
#endif
/* sl_type sldata5_t */
#ifdef sl_data5_type_c
typedef sl_data5_type_c sldata5_t;
#endif
/* sl_type sldata6_t */
#ifdef sl_data6_type_c
typedef sl_data6_type_c sldata6_t;
#endif
/* sl_type sldata7_t */
#ifdef sl_data7_type_c
typedef sl_data7_type_c sldata7_t;
#endif
/* sl_type sldata8_t */
#ifdef sl_data8_type_c
typedef sl_data8_type_c sldata8_t;
#endif
/* sl_type sldata9_t */
#ifdef sl_data9_type_c
typedef sl_data9_type_c sldata9_t;
#endif
/* sl_type sldata10_t */
#ifdef sl_data10_type_c
typedef sl_data10_type_c sldata10_t;
#endif
/* sl_type sldata11_t */
#ifdef sl_data11_type_c
typedef sl_data11_type_c sldata11_t;
#endif
/* sl_type sldata12_t */
#ifdef sl_data12_type_c
typedef sl_data12_type_c sldata12_t;
#endif
/* sl_type sldata13_t */
#ifdef sl_data13_type_c
typedef sl_data13_type_c sldata13_t;
#endif
/* sl_type sldata14_t */
#ifdef sl_data14_type_c
typedef sl_data14_type_c sldata14_t;
#endif
/* sl_type sldata15_t */
#ifdef sl_data15_type_c
typedef sl_data15_type_c sldata15_t;
#endif
/* sl_type sldata16_t */
#ifdef sl_data16_type_c
typedef sl_data16_type_c sldata16_t;
#endif
/* sl_type sldata17_t */
#ifdef sl_data17_type_c
typedef sl_data17_type_c sldata17_t;
#endif
/* sl_type sldata18_t */
#ifdef sl_data18_type_c
typedef sl_data18_type_c sldata18_t;
#endif
/* sl_type sldata19_t */
#ifdef sl_data19_type_c
typedef sl_data19_type_c sldata19_t;
#endif
/* DATAX_TEMPLATE_END */

#define SL_DATA_NMAX (0 \
 + 1 \
 + 1 \
 + 1 \
 + 1 \
 + 1 \
 + 1 \
 + 1 \
 + 1 \
 + 1 \
 + 1 \
 + 1 \
 + 1 \
 + 1 \
 + 1 \
 + 1 \
 + 1 \
 + 1 \
 + 1 \
 + 1 \
 + 1 \
)

/* sl_type slweight_t */
typedef sl_weight_type_c slweight_t;

#define slweight_fmt  sl_weight_type_fmt  /* sl_macro */

#if defined(sl_elem_weight) && defined(sl_weight_intequiv)
typedef sl_weight_type_c slcount_t;       /* sl_type slcount_t */
# define slcount_fmt  sl_weight_type_fmt  /* sl_macro */
#else
typedef sl_int_type_c slcount_t;
# define slcount_fmt  sl_int_type_fmt
#endif


/* sl_type _slpwkey_t slpwkey_t */
typedef struct _slpwkey_t
{
  slpkey_t pkey;
  slweight_t weight;

} slpwkey_t;


/* sl_type _elements_t elements_t */
typedef struct _elements_t
{
  slint_t size, max_size;
  slkey_t *keys;

#ifdef SL_INDEX
  slindex_t *indices;
#endif

/* DATAX_TEMPLATE_BEGIN */
#ifdef SL_DATA0
  sldata0_t *data0;
#endif
#ifdef SL_DATA1
  sldata1_t *data1;
#endif
#ifdef SL_DATA2
  sldata2_t *data2;
#endif
#ifdef SL_DATA3
  sldata3_t *data3;
#endif
#ifdef SL_DATA4
  sldata4_t *data4;
#endif
#ifdef SL_DATA5
  sldata5_t *data5;
#endif
#ifdef SL_DATA6
  sldata6_t *data6;
#endif
#ifdef SL_DATA7
  sldata7_t *data7;
#endif
#ifdef SL_DATA8
  sldata8_t *data8;
#endif
#ifdef SL_DATA9
  sldata9_t *data9;
#endif
#ifdef SL_DATA10
  sldata10_t *data10;
#endif
#ifdef SL_DATA11
  sldata11_t *data11;
#endif
#ifdef SL_DATA12
  sldata12_t *data12;
#endif
#ifdef SL_DATA13
  sldata13_t *data13;
#endif
#ifdef SL_DATA14
  sldata14_t *data14;
#endif
#ifdef SL_DATA15
  sldata15_t *data15;
#endif
#ifdef SL_DATA16
  sldata16_t *data16;
#endif
#ifdef SL_DATA17
  sldata17_t *data17;
#endif
#ifdef SL_DATA18
  sldata18_t *data18;
#endif
#ifdef SL_DATA19
  sldata19_t *data19;
#endif
/* DATAX_TEMPLATE_END */

} elements_t;


/* sl_type _packed_element_t packed_element_t */
typedef struct _packed_element_t
{
  slkey_t key;

#ifdef SL_PACKED_INDEX
  slindex_t index;
#endif

/* DATAX_TEMPLATE_BEGIN */
#ifdef SL_DATA0
# ifdef sl_data0_flex
  sldata0_t data0[];
# else
  sldata0_t data0[sl_data0_size_c];
# endif
#endif
#ifdef SL_DATA1
# ifdef sl_data1_flex
  sldata1_t data1[];
# else
  sldata1_t data1[sl_data1_size_c];
# endif
#endif
#ifdef SL_DATA2
# ifdef sl_data2_flex
  sldata2_t data2[];
# else
  sldata2_t data2[sl_data2_size_c];
# endif
#endif
#ifdef SL_DATA3
# ifdef sl_data3_flex
  sldata3_t data3[];
# else
  sldata3_t data3[sl_data3_size_c];
# endif
#endif
#ifdef SL_DATA4
# ifdef sl_data4_flex
  sldata4_t data4[];
# else
  sldata4_t data4[sl_data4_size_c];
# endif
#endif
#ifdef SL_DATA5
# ifdef sl_data5_flex
  sldata5_t data5[];
# else
  sldata5_t data5[sl_data5_size_c];
# endif
#endif
#ifdef SL_DATA6
# ifdef sl_data6_flex
  sldata6_t data6[];
# else
  sldata6_t data6[sl_data6_size_c];
# endif
#endif
#ifdef SL_DATA7
# ifdef sl_data7_flex
  sldata7_t data7[];
# else
  sldata7_t data7[sl_data7_size_c];
# endif
#endif
#ifdef SL_DATA8
# ifdef sl_data8_flex
  sldata8_t data8[];
# else
  sldata8_t data8[sl_data8_size_c];
# endif
#endif
#ifdef SL_DATA9
# ifdef sl_data9_flex
  sldata9_t data9[];
# else
  sldata9_t data9[sl_data9_size_c];
# endif
#endif
#ifdef SL_DATA10
# ifdef sl_data10_flex
  sldata10_t data10[];
# else
  sldata10_t data10[sl_data10_size_c];
# endif
#endif
#ifdef SL_DATA11
# ifdef sl_data11_flex
  sldata11_t data11[];
# else
  sldata11_t data11[sl_data11_size_c];
# endif
#endif
#ifdef SL_DATA12
# ifdef sl_data12_flex
  sldata12_t data12[];
# else
  sldata12_t data12[sl_data12_size_c];
# endif
#endif
#ifdef SL_DATA13
# ifdef sl_data13_flex
  sldata13_t data13[];
# else
  sldata13_t data13[sl_data13_size_c];
# endif
#endif
#ifdef SL_DATA14
# ifdef sl_data14_flex
  sldata14_t data14[];
# else
  sldata14_t data14[sl_data14_size_c];
# endif
#endif
#ifdef SL_DATA15
# ifdef sl_data15_flex
  sldata15_t data15[];
# else
  sldata15_t data15[sl_data15_size_c];
# endif
#endif
#ifdef SL_DATA16
# ifdef sl_data16_flex
  sldata16_t data16[];
# else
  sldata16_t data16[sl_data16_size_c];
# endif
#endif
#ifdef SL_DATA17
# ifdef sl_data17_flex
  sldata17_t data17[];
# else
  sldata17_t data17[sl_data17_size_c];
# endif
#endif
#ifdef SL_DATA18
# ifdef sl_data18_flex
  sldata18_t data18[];
# else
  sldata18_t data18[sl_data18_size_c];
# endif
#endif
#ifdef SL_DATA19
# ifdef sl_data19_flex
  sldata19_t data19[];
# else
  sldata19_t data19[sl_data19_size_c];
# endif
#endif
/* DATAX_TEMPLATE_END */

} packed_element_t;


/* sl_type _packed_elements_t packed_elements_t */
typedef struct _packed_elements_t
{
  slint_t size, max_size;
  
  packed_element_t *elements;
  
} packed_elements_t;


#ifndef SLCINT_T
#define SLCINT_T
typedef long long int slcint_t;
#define slcint_fmt  "ll"
/*#define slcint_sfx  LL*/
#endif


#define SLCM_KEYS     (((slcint_t) 1) << 0)
#define SLCM_INDICES  (((slcint_t) 1) << 1)
#define SLCM_WEIGHTS  (((slcint_t) 1) << 2)

/* DATAX_TEMPLATE_BEGIN */
#define SLCM_DATA0    (((slcint_t) 1) << (3+0))
#define SLCM_DATA1    (((slcint_t) 1) << (3+1))
#define SLCM_DATA2    (((slcint_t) 1) << (3+2))
#define SLCM_DATA3    (((slcint_t) 1) << (3+3))
#define SLCM_DATA4    (((slcint_t) 1) << (3+4))
#define SLCM_DATA5    (((slcint_t) 1) << (3+5))
#define SLCM_DATA6    (((slcint_t) 1) << (3+6))
#define SLCM_DATA7    (((slcint_t) 1) << (3+7))
#define SLCM_DATA8    (((slcint_t) 1) << (3+8))
#define SLCM_DATA9    (((slcint_t) 1) << (3+9))
#define SLCM_DATA10    (((slcint_t) 1) << (3+10))
#define SLCM_DATA11    (((slcint_t) 1) << (3+11))
#define SLCM_DATA12    (((slcint_t) 1) << (3+12))
#define SLCM_DATA13    (((slcint_t) 1) << (3+13))
#define SLCM_DATA14    (((slcint_t) 1) << (3+14))
#define SLCM_DATA15    (((slcint_t) 1) << (3+15))
#define SLCM_DATA16    (((slcint_t) 1) << (3+16))
#define SLCM_DATA17    (((slcint_t) 1) << (3+17))
#define SLCM_DATA18    (((slcint_t) 1) << (3+18))
#define SLCM_DATA19    (((slcint_t) 1) << (3+19))
/* DATAX_TEMPLATE_END */

#define SLCM_DATA     (((slcint_t) 0) \
  |SLCM_DATA0 \
  |SLCM_DATA1 \
  |SLCM_DATA2 \
  |SLCM_DATA3 \
  |SLCM_DATA4 \
  |SLCM_DATA5 \
  |SLCM_DATA6 \
  |SLCM_DATA7 \
  |SLCM_DATA8 \
  |SLCM_DATA9 \
  |SLCM_DATA10 \
  |SLCM_DATA11 \
  |SLCM_DATA12 \
  |SLCM_DATA13 \
  |SLCM_DATA14 \
  |SLCM_DATA15 \
  |SLCM_DATA16 \
  |SLCM_DATA17 \
  |SLCM_DATA18 \
  |SLCM_DATA19 \
  )

#define SLCM_ALL      (~((slcint_t) 0))


/* sl_type _classification_info_t classification_info_t classification_info */
typedef struct _classification_info_t
{
  slint_t nclasses;
  slkey_pure_t *keys;
  slint_t *counts;
  slint_t *masks;

  /* */
  slint_t *all_local_sizes;
  slint_t *local_lt_eq_counts;
  slint_t *all_local_lt_eq_counts;

} classification_info_t, classification_info;  /* deprecated 'classification_info' */


/* key2class, sl_type key2class_f */
typedef slint_t (*key2class_f)(slkey_t *, slint, void *);

/* pivot-element, sl_type pivot_f */
typedef slint_t (*pivot_f)(elements_t *);

/* sorting-network, sl_type sortnet_f sortnet_data_t */
typedef void *sortnet_data_t;
typedef slint_t (*sortnet_f)(slint_t size, slint_t rank, slint_t stage, sortnet_data_t snd, slint_t *up);

/* merge2, sl_type merge2x_f merge2X_f */
typedef slint_t (*merge2x_f)(elements_t *s0, elements_t *s1, elements_t *sx);
typedef slint_t (*merge2X_f)(elements_t *s0, elements_t *s1, elements_t *sx, elements_t *t);

/* sl_type _permute_generic_t permute_generic_t */
typedef struct _permute_generic_t
{
  int type;

  spec_tloc_f *tloc;
  spec_tloc_rearrange_db_f *tloc_rearrange_db;
  spec_tloc_rearrange_ip_f *tloc_rearrange_ip;

  spec_tloc_mod_f *tloc_mod;
  spec_tloc_mod_rearrange_db_f *tloc_mod_rearrange_db;
  spec_tloc_mod_rearrange_ip_f *tloc_mod_rearrange_ip;

} permute_generic_t;

/* sl_macro PERMUTE_GENERIC_DEFINE_TLOC PERMUTE_GENERIC_INIT_TLOC PERMUTE_GENERIC_INIT_EXT_TLOC */
#define PERMUTE_GENERIC_DEFINE_TLOC(_tl_, _s_...)      SPEC_DEFINE_TLOC(_tl_, _tl_, _s_)
#define PERMUTE_GENERIC_INIT_TLOC(_tl_)                { 1, _tl_, SPEC_EXT_PARAM_TLOC_NULL,  NULL, SPEC_EXT_PARAM_TLOC_MOD_NULL }
#define PERMUTE_GENERIC_INIT_EXT_TLOC(_tl_)            { 1, _tl_, SPEC_EXT_PARAM_TLOC(_tl_), NULL, SPEC_EXT_PARAM_TLOC_MOD_NULL }

/* sl_macro PERMUTE_GENERIC_DEFINE_TLOC_MOD PERMUTE_GENERIC_INIT_TLOC_MOD PERMUTE_GENERIC_INIT_EXT_TLOC_MOD */
#define PERMUTE_GENERIC_DEFINE_TLOC_MOD(_tl_, _s_...)  SPEC_DEFINE_TLOC_MOD(_tl_, _tl_, _s_)
#define PERMUTE_GENERIC_INIT_TLOC_MOD(_tl_)            { 2, NULL, SPEC_EXT_PARAM_TLOC_MOD_NULL, _tl_, SPEC_EXT_PARAM_TLOC_MOD_NULL }
#define PERMUTE_GENERIC_INIT_EXT_TLOC_MOD(_tl_)        { 2, NULL, SPEC_EXT_PARAM_TLOC_MOD_NULL, _tl_, SPEC_EXT_PARAM_TLOC_MOD(_tl_) }

/* sl_type _split_generic_t split_generic_t */
typedef struct _split_generic_t
{
  int type;

  slint_t max_tprocs;

  spec_tproc_f *tproc;
  spec_tproc_ext_t tproc_ext;

  spec_tproc_mod_f *tproc_mod;
  spec_tproc_mod_ext_t tproc_mod_ext;

  spec_tprocs_f *tprocs;
  spec_tprocs_ext_t tprocs_ext;

  spec_tprocs_mod_f *tprocs_mod;
  spec_tprocs_mod_ext_t tprocs_mod_ext;

  spec_tproc_reset_f *reset;

} split_generic_t;

/* sl_macro SPLIT_GENERIC_DEFINE_TPROC SPLIT_GENERIC_INIT_TPROC SPLIT_GENERIC_INIT_EXT_TPROC */
#define SPLIT_GENERIC_DEFINE_TPROC(_tp_, _s_...)                SPEC_DEFINE_TPROC(_tp_, _tp_, _s_)
#define SPLIT_GENERIC_INIT_TPROC(_tp_, _r_...)                  { 1, 0, _tp_, SPEC_EXT_PARAM_TPROC_NULL,  NULL, SPEC_EXT_PARAM_TPROC_MOD_NULL, NULL, SPEC_EXT_PARAM_TPROCS_NULL, NULL, SPEC_EXT_PARAM_TPROCS_MOD_NULL, _r_ }
#define SPLIT_GENERIC_INIT_EXT_TPROC(_tp_, _r_...)              { 1, 0, _tp_, SPEC_EXT_PARAM_TPROC(_tp_), NULL, SPEC_EXT_PARAM_TPROC_MOD_NULL, NULL, SPEC_EXT_PARAM_TPROCS_NULL, NULL, SPEC_EXT_PARAM_TPROCS_MOD_NULL, _r_ }

/* sl_macro SPLIT_GENERIC_DEFINE_TPROC_MOD SPLIT_GENERIC_INIT_TPROC_MOD SPLIT_GENERIC_INIT_EXT_TPROC_MOD */
#define SPLIT_GENERIC_DEFINE_TPROC_MOD(_tp_, _s_...)            SPEC_DEFINE_TPROC_MOD(_tp_, _tp_, _s_)
#define SPLIT_GENERIC_INIT_TPROC_MOD(_tp_, _r_...)              { 2, 0, NULL, SPEC_EXT_PARAM_TPROC_NULL, _tp_, SPEC_EXT_PARAM_TPROC_MOD_NULL,  NULL, SPEC_EXT_PARAM_TPROCS_NULL, NULL, SPEC_EXT_PARAM_TPROCS_MOD_NULL, _r_ }
#define SPLIT_GENERIC_INIT_EXT_TPROC_MOD(_tp_, _r_...)          { 2, 0, NULL, SPEC_EXT_PARAM_TPROC_NULL, _tp_, SPEC_EXT_PARAM_TPROC_MOD(_tp_), NULL, SPEC_EXT_PARAM_TPROCS_NULL, NULL, SPEC_EXT_PARAM_TPROCS_MOD_NULL, _r_ }

/* sl_macro SPLIT_GENERIC_DEFINE_TPROCS SPLIT_GENERIC_INIT_TPROCS SPLIT_GENERIC_INIT_EXT_TPROCS */
#define SPLIT_GENERIC_DEFINE_TPROCS(_tp_, _s_...)               SPEC_DEFINE_TPROCS(_tp_, _tp_, _s_)
#define SPLIT_GENERIC_INIT_TPROCS(_tp_, _xtp_, _r_...)          { 3, (_xtp_), NULL, SPEC_EXT_PARAM_TPROC_NULL, NULL, SPEC_EXT_PARAM_TPROC_MOD_NULL, _tp_, SPEC_EXT_PARAM_TPROCS_NULL,  NULL, SPEC_EXT_PARAM_TPROCS_MOD_NULL, _r_ }
#define SPLIT_GENERIC_INIT_EXT_TPROCS(_tp_, _xtp_, _r_...)      { 3, (_xtp_), NULL, SPEC_EXT_PARAM_TPROC_NULL, NULL, SPEC_EXT_PARAM_TPROC_MOD_NULL, _tp_, SPEC_EXT_PARAM_TPROCS(_tp_), NULL, SPEC_EXT_PARAM_TPROCS_MOD_NULL, _r_ }

/* sl_macro SPLIT_GENERIC_DEFINE_TPROCS_MOD SPLIT_GENERIC_INIT_TPROCS_MOD SPLIT_GENERIC_INIT_EXT_TPROCS_MOD */
#define SPLIT_GENERIC_DEFINE_TPROCS_MOD(_tp_, _s_...)           SPEC_DEFINE_TPROCS_MOD(_tp_, _tp_, _s_)
#define SPLIT_GENERIC_INIT_TPROCS_MOD(_tp_, _xtp_, _r_...)      { 4, (_xtp_), NULL, SPEC_EXT_PARAM_TPROC_NULL, NULL, SPEC_EXT_PARAM_TPROC_MOD_NULL, NULL, SPEC_EXT_PARAM_TPROCS_NULL, _tp_, SPEC_EXT_PARAM_TPROCS_MOD_NULL,  _r_ }
#define SPLIT_GENERIC_INIT_EXT_TPROCS_MOD(_tp_, _xtp_, _r_...)  { 4, (_xtp_), NULL, SPEC_EXT_PARAM_TPROC_NULL, NULL, SPEC_EXT_PARAM_TPROC_MOD_NULL, NULL, SPEC_EXT_PARAM_TPROCS_NULL, _tp_, SPEC_EXT_PARAM_TPROCS_MOD(_tp_), _r_ }

/* sl_type tloc_f tloc_mod_f */
typedef slint_t tloc_f(elements_t *b, slint_t x, void *tloc_data);
typedef slint_t tloc_mod_f(elements_t *b, slint_t x, void *tloc_data, elements_t *mod);

/* sl_type tproc_f tproc_mod_f tprocs_f tprocs_mod_f */
typedef int tproc_f(elements_t *b, slint_t x, void *tproc_data);
typedef int tproc_mod_f(elements_t *b, slint_t x, void *tproc_data, elements_t *mod);
typedef void tprocs_f(elements_t *b, slint_t x, void *tproc_data, slint_t *nprocs, int *procs);
typedef void tprocs_mod_f(elements_t *b, slint_t x, void *tproc_data, slint_t *nprocs, int *procs, elements_t *mod);

/* sl_type tproc_reset_f */
typedef void tproc_reset_f(void *tproc_data);

/* sl_macro TPROC_RESET_NULL */
#define TPROC_RESET_NULL  NULL

/* sl_type tproc_t */
typedef struct _spec_tproc_t *tproc_t;

/* sl_type _tproc_exdef tproc_exdef */
typedef struct _tproc_exdef
{
  int type;

  spec_tproc_ext_t tproc_ext;
  spec_tproc_mod_ext_t tproc_mod_ext;
  spec_tprocs_ext_t tprocs_ext;
  spec_tprocs_mod_ext_t tprocs_mod_ext;

} const *tproc_exdef;

/* sl_macro TPROC_EXDEF_NULL */
#define TPROC_EXDEF_NULL  NULL

/* sl_macro TPROC_EXDEF_DEFINE_TPROC TPROC_EXDEF_DEFINE_TPROC_MOD TPROC_EXDEF_DEFINE_TPROCS TPROC_EXDEF_DEFINE_TPROCS_MOD */
#define TPROC_EXDEF_DEFINE_TPROC(_name_, _tp_, _s_...) \
  SPEC_DEFINE_TPROC(_name_, _tp_, _s_) \
  _s_ const struct _tproc_exdef _##_name_ = { 1, SPEC_EXT_PARAM_TPROC(_name_), SPEC_EXT_PARAM_TPROC_MOD_NULL, SPEC_EXT_PARAM_TPROCS_NULL, SPEC_EXT_PARAM_TPROCS_MOD_NULL }, *_name_ = &_##_name_;

#define TPROC_EXDEF_DEFINE_TPROC_MOD(_name_, _tp_, _s_...) \
  SPEC_DEFINE_TPROC_MOD(_name_, _tp_, _s_) \
  _s_ const struct _tproc_exdef _##_name_ = { 2, SPEC_EXT_PARAM_TPROC_NULL, SPEC_EXT_PARAM_TPROC_MOD(_name_), SPEC_EXT_PARAM_TPROCS_NULL, SPEC_EXT_PARAM_TPROCS_MOD_NULL }, *_name_ = &_##_name_;

#define TPROC_EXDEF_DEFINE_TPROCS(_name_, _tp_, _s_...) \
  SPEC_DEFINE_TPROCS(_name_, _tp_, _s_) \
  _s_ const struct _tproc_exdef _##_name_ = { 3, SPEC_EXT_PARAM_TPROC_NULL, SPEC_EXT_PARAM_TPROC_MOD_NULL, SPEC_EXT_PARAM_TPROCS(_name_), SPEC_EXT_PARAM_TPROCS_MOD_NULL }, *_name_ = &_##_name_;

#define TPROC_EXDEF_DEFINE_TPROCS_MOD(_name_, _tp_, _s_...) \
  SPEC_DEFINE_TPROCS_MOD(_name_, _tp_, _s_) \
  _s_ const struct _tproc_exdef _##_name_ = { 4, SPEC_EXT_PARAM_TPROC_NULL, SPEC_EXT_PARAM_TPROC_MOD_NULL, SPEC_EXT_PARAM_TPROCS_NULL, SPEC_EXT_PARAM_TPROCS_MOD(_name_) }, *_name_ = &_##_name_;


/* mpi_elements_alltoall_specific */
#ifndef SL_MEAS_TYPE_ALLTOALLV
# define SL_MEAS_TYPE_ALLTOALLV  0
#endif

#ifndef SL_MEAS_TYPE_SENDRECV
# define SL_MEAS_TYPE_SENDRECV   1
#endif


/* deprecated, sl_type k2c_func pivot_func sn_func m2x_func m2X_func */
typedef key2class_f k2c_func;
typedef pivot_f pivot_func;
typedef sortnet_f sn_func;
typedef merge2x_f m2x_func;
typedef merge2X_f m2X_func;


/* sl_type _mergek_t mergek_t */
typedef struct _mergek_t
{
  sortnet_f sn;
  sortnet_data_t snd;

  merge2x_f m2x;
  elements_t *sx;

} mergek_t;


/* sl_type keys_init_type_t keys_init_data_t */
typedef slint_t keys_init_type_t;
typedef void *keys_init_data_t;

/* sl_type key_set_data_t key_set_f */
typedef void *key_set_data_t;
typedef void (*key_set_f)(slkey_pure_t *k, key_set_data_t d);


#undef SL_EKIT_SET
#define SL_EKIT_SET         1
#undef SL_EKIT_SET_FUNC
#define SL_EKIT_SET_FUNC    2
#undef SL_EKIT_RAND
#define SL_EKIT_RAND        3
#undef SL_EKIT_RAND_QUAD
#define SL_EKIT_RAND_QUAD   4
#undef SL_EKIT_RAND_AND
#define SL_EKIT_RAND_AND    5
#undef SL_EKIT_URAND
#define SL_EKIT_URAND       6
#undef SL_EKIT_NRAND
#define SL_EKIT_NRAND       7


#ifndef SL_EIK_OFFSET
# define SL_EIK_OFFSET     65536LL
#endif

#ifndef SL_EIK_SET
# define SL_EIK_SET        SL_EIK_OFFSET*1
#endif

#ifndef SL_EIK_RAND
# define SL_EIK_RAND       SL_EIK_OFFSET*2
#endif

#ifndef SL_EIK_RAND_QUAD
# define SL_EIK_RAND_QUAD  SL_EIK_OFFSET*3
#endif

#ifndef SL_EIK_RAND_AND
# define SL_EIK_RAND_AND   SL_EIK_OFFSET*4
#endif

#ifndef SL_EIK_RAND_NORM
# define SL_EIK_RAND_NORM  SL_EIK_OFFSET*5
#endif


/* elements_keys_stats */
#ifndef SL_EKS_MIN
# define SL_EKS_MIN   0
#endif

#ifndef SL_EKS_MAX
# define SL_EKS_MAX   1
#endif

#ifndef SL_EKS_SUM
# define SL_EKS_SUM   2
#endif

#ifndef SL_EKS_AVG
# define SL_EKS_AVG   3
#endif

#ifndef SL_EKS_STD
# define SL_EKS_STD   4
#endif

#ifndef SL_EKS_SIZE
# define SL_EKS_SIZE  5
#endif


#ifndef SL_SORTED_IN
# define SL_SORTED_IN   0x1LL
#endif

#ifndef SL_SORTED_OUT
# define SL_SORTED_OUT  0x2LL
#endif


#ifndef SL_MSEG_FM_EXACT
# define SL_MSEG_FM_EXACT         0
#endif
#ifndef SL_MSEG_FM_ALLORNOTHING
# define SL_MSEG_FM_ALLORNOTHING  1
#endif
#ifndef SL_MSEG_FM_MIDDLE
# define SL_MSEG_FM_MIDDLE        2
#endif


/* partition conditions, sl_type _partcond2_t partcond2_t */
typedef struct _partcond2_t
{
  int weighted;
  double min_count, max_count;
  double min_weight, max_weight;
  double min_cpart, max_cpart;
  double min_wpart, max_wpart;

} partcond2_t;


#ifndef SLPC_COUNTS_MM
# define SLPC_COUNTS_MM   0x1
#endif
#ifndef SLPC_COUNTS_LH
# define SLPC_COUNTS_LH   0x2
#endif
#ifndef SLPC_WEIGHTS_MM
# define SLPC_WEIGHTS_MM  0x4
#endif
#ifndef SLPC_WEIGHTS_LH
# define SLPC_WEIGHTS_LH  0x8
#endif

/* partition conditions, sl_type _partcond_t partcond_t partcond_p */
typedef struct _partcond_t
{
  slint_t pcm;
  double count_min, count_max;
  double count_low, count_high;
  double weight_min, weight_max;
  double weight_low, weight_high;

} partcond_t, *partcond_p;


/* internal partition conditions, sl_type _partcond_intern_t partcond_intern_t partcond_intern_p */
typedef struct _partcond_intern_t
{
  slint_t pcm;
  slint_t count_min, count_max;
  slint_t count_low, count_high;
#ifdef elem_weight
  slweight_t weight_min, weight_max;
  slweight_t weight_low, weight_high;
#endif

} partcond_intern_t, *partcond_intern_p;


/* sl_type _parttype_t parttype_t parttype_p */
typedef struct _parttype_t
{
  slint_t type;

} parttype_t, *parttype_p;


/* generic binning method */

/* sl_type _bin_t bin_t */
typedef struct _bin_t
{
  elements_t s;

#ifdef elem_weight
  slweight_t weight;
#endif

} bin_t;


/* sl_type _splitter_t splitter_t */
typedef struct _splitter_t
{
  slint_t n;

  int *displs;
  slkey_pure_t *s;
  slint_t *sn;

} splitter_t;


struct _binning_t;

/* sl_type binning_pre_f binning_exec_f binning_refine_f binning_hit_f binning_finalize_f binning_post_f */
typedef slint_t (*binning_pre_f)(struct _binning_t *bm);
typedef slint_t (*binning_exec_f)(struct _binning_t *bm, bin_t *bin, slcount_t *counts, slweight_t *weights);
typedef slint_t (*binning_refine_f)(struct _binning_t *bm, bin_t *bin, slint_t k, slcount_t *counts, slweight_t *weights, splitter_t *sp, slint_t s, bin_t *new_bin);
typedef slint_t (*binning_hit_f)(struct _binning_t *bm, bin_t *bin, slint_t k, slcount_t *counts, splitter_t *sp, slint_t s);
typedef slint_t (*binning_finalize_f)(struct _binning_t *bm, bin_t *bin, slint_t dc, slweight_t dw, slint_t lc_min, slint_t lc_max, slcount_t *lcs, slweight_t *lws, splitter_t *sp, slint_t s);
typedef slint_t (*binning_post_f)(struct _binning_t *bm);


/* sl_type _binning_data_t binning_data_t */
typedef union _binning_data_t
{
  struct
  {
    slint_t rhigh, rlow, rwidth;
    slint_t rcurrent;
    slkey_pure_t bit_mask;

    elements_t sx;

  } radix;

} binning_data_t;


/* sl_type _binning_t binning_t */
typedef struct _binning_t
{
  slint_t nbins, max_nbins;
  
  binning_pre_f pre;
  binning_exec_f exec;
  binning_refine_f refine;
  binning_hit_f hit;
  binning_finalize_f finalize;
  binning_post_f post;

  slint_t sorted;

  slint_t docounts;
#ifdef elem_weight
  slint_t doweights;
#endif

  binning_data_t bd;

} binning_t;


/* sl_type _local_bins_t local_bins_t */
typedef struct _local_bins_t
{
  binning_t *bm;

  slint_t nbins, max_nbins;
  slint_t nelements;

  slint_t docounts;
#ifdef elem_weight
  slint_t doweights;
#endif

  slint_t nbinnings, max_nbinnings;

  slint_t nbins_new, last_new_b, last_new_k;
  bin_t *bins, *bins_new;
  bin_t *bins0, *bins1;

  slint_t *bcws;

#if defined(elem_weight) && defined(sl_weight_intequiv)
  slint_t cw_factor, w_index, bin_cw_factor;
  slweight_t *cws, *bin_cws;
  slweight_t *prefix_cws;
#else
  slint_t *cs, *bin_cs;
  slint_t *prefix_cs;
# ifdef elem_weight
  slweight_t *ws, *bin_ws;
  slweight_t *prefix_ws;
# endif
#endif

  slint_t last_exec_b;

} local_bins_t;


/* sl_type _global_bins_t global_bins_t */
typedef struct _global_bins_t
{
  binning_t *bm;
  
  local_bins_t lb;

  slint_t *bcws;

#if defined(elem_weight) && defined(sl_weight_intequiv)
  slweight_t *cws;
  slweight_t *prefix_cws;
#else
  slint_t *cs;
  slint_t *prefix_cs;
# ifdef elem_weight
  slweight_t *ws;
  slweight_t *prefix_ws;
# endif
#endif

} global_bins_t;


/* sl_type rti_cmc_t */
typedef struct
{
  slint_t cmp, movek, moved;

} rti_cmc_t;

#ifndef my_rti_ccmp
# define my_rti_ccmp(m)    m.cmc.cmp
# define my_rti_cmovek(m)  m.cmc.movek
# define my_rti_cmoved(m)  m.cmc.moved
#endif


/* sl_type rti_tim_t */
typedef struct
{
  double start, stop;
  double last, cumu;

  slint_t num;

} rti_tim_t[rti_tids];

#ifndef my_rti_tlast
# define my_rti_tlast(m, t)  m.tim[t].last
# define my_rti_tcumu(m, t)  m.tim[t].cumu
# define my_rti_tnum(m, t)   m.tim[t].num
#endif


/* sl_type rti_mem_t */
typedef struct
{
  slint_t nalloc, nfree;
  slint_t max, cur, cur_max;

} rti_mem_t;


/* sl_type rti_t */
typedef struct
{
  /* compare-move-counter */
  rti_cmc_t cmc;
  /* timer */
  rti_tim_t tim;
  /* memory */
  rti_mem_t mem;

} rti_t;

#ifndef my_rti_reset
# define my_rti_reset(m)  memset((void *) &m, 0, sizeof(m))
#endif


#ifdef SL_USE_MPI
/* sl_type _sl_mpi_context_t sl_mpi_context_t */
typedef struct _sl_mpi_context_t
{
  int size, rank;
  MPI_Comm comm;

} sl_mpi_context_t;
#endif


#ifdef SL_USE_OMP
/* sl_type _sl_omp_context_t sl_omp_context_t */
typedef struct _sl_omp_context_t
{
  int thread_num, num_threads, *coop_thread_nums;

} sl_omp_context_t;
#endif


/* sl_type _sl_context_t sl_context_t */
typedef struct _sl_context_t
{
#include "sl_context_struct.h"
} sl_context_t;


#endif /* __SL_TYPES_H__ */
