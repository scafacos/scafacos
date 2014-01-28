!
! Copyright (C) 2011, 2012, 2013 Rene Halver, Michael Hofmann
!
! This file is part of ScaFaCoS.
!
! ScaFaCoS is free software: you can redistribute it and/or modify
! it under the terms of the GNU Lesser Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! ScaFaCoS is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser Public License for more details.
!
! You should have received a copy of the GNU Lesser Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>. 
!


#ifndef FCS_FCS4FORTRAN_DEFINITIONS_INCLUDED
#define FCS_FCS4FORTRAN_DEFINITIONS_INCLUDED


!
! @file FCSInterface.h
! @brief public definitions of the ScaFaCoS library interface
! @author Rene Halver, Michael Hofmann
!


!
! @brief definition of bool data type
!
#define fcs_boolean_kind_isoc  fcs_integer_kind_isoc
#define FCS4FORTRAN_TRUE   1
#define FCS4FORTRAN_FALSE  0


!
! @brief definitions of return codes
!
#define FCS4FORTRAN_SUCCESS                    0
#define FCS4FORTRAN_ERROR_NULL_ARGUMENT        1
#define FCS4FORTRAN_ERROR_ALLOC_FAILED         2
#define FCS4FORTRAN_ERROR_WRONG_ARGUMENT       3
#define FCS4FORTRAN_ERROR_MISSING_ELEMENT      4
#define FCS4FORTRAN_ERROR_LOGICAL_ERROR        5
#define FCS4FORTRAN_ERROR_INCOMPATIBLE_METHOD  6
#define FCS4FORTRAN_ERROR_NOT_IMPLEMENTED      7
#define FCS4FORTRAN_ERROR_FORTRAN_CALL         8
#define FCS4FORTRAN_ERROR_RESULT_CREATE        9


!
! @brief definitions of numerical method identifiers
!
#define FCS4FORTRAN_METHOD_NONE   40
#define FCS4FORTRAN_METHOD_FMM    32
#define FCS4FORTRAN_METHOD_P2NFFT 33
#define FCS4FORTRAN_METHOD_PEPC   34
#define FCS4FORTRAN_METHOD_P3M    35
#define FCS4FORTRAN_METHOD_PP3MG  36
#define FCS4FORTRAN_METHOD_VMG    37
#define FCS4FORTRAN_METHOD_DIRECT 38
#define FCS4FORTRAN_METHOD_MEMD   39
#define FCS4FORTRAN_METHOD_MMM1D  41
#define FCS4FORTRAN_METHOD_EWALD  42
#define FCS4FORTRAN_METHOD_MMM2D  43
#define FCS4FORTRAN_METHOD_WOLF   44


!
! @brief definitions of tolerance types
!
#define FCS4FORTRAN_TOLERANCE_TYPE_UNDEFINED      0
#define FCS4FORTRAN_TOLERANCE_TYPE_ENERGY         1
#define FCS4FORTRAN_TOLERANCE_TYPE_ENERGY_REL     2
#define FCS4FORTRAN_TOLERANCE_TYPE_POTENTIAL      3
#define FCS4FORTRAN_TOLERANCE_TYPE_POTENTIAL_REL  4
#define FCS4FORTRAN_TOLERANCE_TYPE_FIELD          5
#define FCS4FORTRAN_TOLERANCE_TYPE_FIELD_REL      6


#endif
