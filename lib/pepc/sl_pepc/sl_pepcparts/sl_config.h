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


#ifndef __SL_CONFIG_H__
#define __SL_CONFIG_H__


#include "fortran2c_types.h"


/* enable runtime_informations */
/*#define SL_USE_RTI
#define SL_USE_RTI_TIM*/


/* standard (SL) integer data type */
#define sl_int_type_c          long long
#define sl_int_type_mpi        MPI_LONG_LONG
#define sl_int_size_mpi        1
#define sl_int_type_fmt        "lld"


/* index data type */
#define sl_index_type_c        FINT_TYPE_C
#define sl_index_type_mpi      FINT_TYPE_MPI
#define sl_index_size_mpi      1
#define sl_index_type_fmt      FINT_TYPE_FMT

/* use indices */
#define SL_INDEX


/* keys */
#define sl_key_type_c          FINT8_TYPE_C
#define sl_key_type_mpi        FINT8_TYPE_MPI
#define sl_key_size_mpi        1
#define sl_key_type_fmt        FINT8_TYPE_FMT
#define sl_key_integer

/* data0: x */
#define SL_DATA0
#define sl_data0_type_c        FREAL8_TYPE_C
#define sl_data0_size_c        1
#define sl_data0_type_mpi      FREAL8_TYPE_MPI
#define sl_data0_size_mpi      1

/* data1: y */
#define SL_DATA1
#define sl_data1_type_c        FREAL8_TYPE_C
#define sl_data1_size_c        1
#define sl_data1_type_mpi      FREAL8_TYPE_MPI
#define sl_data1_size_mpi      1

/* data2: z */
#define SL_DATA2
#define sl_data2_type_c        FREAL8_TYPE_C
#define sl_data2_size_c        1
#define sl_data2_type_mpi      FREAL8_TYPE_MPI
#define sl_data2_size_mpi      1

/* data3: ux */
#define SL_DATA3
#define sl_data3_type_c        FREAL8_TYPE_C
#define sl_data3_size_c        1
#define sl_data3_type_mpi      FREAL8_TYPE_MPI
#define sl_data3_size_mpi      1

/* data4: uy */
#define SL_DATA4
#define sl_data4_type_c        FREAL8_TYPE_C
#define sl_data4_size_c        1
#define sl_data4_type_mpi      FREAL8_TYPE_MPI
#define sl_data4_size_mpi      1

/* data5: uz */
#define SL_DATA5
#define sl_data5_type_c        FREAL8_TYPE_C
#define sl_data5_size_c        1
#define sl_data5_type_mpi      FREAL8_TYPE_MPI
#define sl_data5_size_mpi      1

/* data6: q */
#define SL_DATA6
#define sl_data6_type_c        FREAL8_TYPE_C
#define sl_data6_size_c        1
#define sl_data6_type_mpi      FREAL8_TYPE_MPI
#define sl_data6_size_mpi      1

/* data7: m */
#define SL_DATA7
#define sl_data7_type_c        FREAL8_TYPE_C
#define sl_data7_size_c        1
#define sl_data7_type_mpi      FREAL8_TYPE_MPI
#define sl_data7_size_mpi      1

/* data8: work */
#define SL_DATA8
#define sl_data8_type_c        FREAL8_TYPE_C
#define sl_data8_size_c        1
#define sl_data8_type_mpi      FREAL8_TYPE_MPI
#define sl_data8_size_mpi      1

/* data9: ex */
#define SL_DATA9
#define sl_data9_type_c        FREAL8_TYPE_C
#define sl_data9_size_c        1
#define sl_data9_type_mpi      FREAL8_TYPE_MPI
#define sl_data9_size_mpi      1

/* data10: ey */
#define SL_DATA10
#define sl_data10_type_c       FREAL8_TYPE_C
#define sl_data10_size_c       1
#define sl_data10_type_mpi     FREAL8_TYPE_MPI
#define sl_data10_size_mpi     1

/* data11: ez */
#define SL_DATA11
#define sl_data11_type_c       FREAL8_TYPE_C
#define sl_data11_size_c       1
#define sl_data11_type_mpi     FREAL8_TYPE_MPI
#define sl_data11_size_mpi     1

/* data12: pelabel */
#define SL_DATA12
#define sl_data12_type_c       FINT_TYPE_C
#define sl_data12_size_c       1
#define sl_data12_type_mpi     FINT_TYPE_MPI
#define sl_data12_size_mpi     1


/* weighted elements */
#define sl_elem_weight(e, at)  ((e)->data8[at])

#define sl_data8_weight

#endif /* __SL_CONFIG_H__ */
