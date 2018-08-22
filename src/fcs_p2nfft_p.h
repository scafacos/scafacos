/*
  Copyright (C) 2011-2013 Michael Pippig
  Copyright (C) 2011-2012 Rene Halver
  Copyright (C) 2011 Sebastian Banert

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



#ifndef FCS_P2NFFT_P_INCLUDED
#define FCS_P2NFFT_P_INCLUDED

#include "fcs_definitions.h"
#include "fcs_result_p.h"
#include "fcs_interface_p.h"

/**
 * @file fcs_p2nfft_p.h
 * @brief file containing all p2nfft specific functions (public version)
 * @author Rene Halver, Michael Pippig
 */

typedef struct fcs_p2nfft_parameters_t *fcs_p2nfft_parameters;

/* Enable definition of multiple wrappers which point to the same internal function.
 * We use this trick to get an P3M-compliant interface but also remain have our own nomenclature. */
#undef FCS_P2NFFT_INTERFACE_WITH_REDIRECTIONS
#define FCS_P2NFFT_INTERFACE_WITH_REDIRECTIONS 1

/* clear macro definitions of headers (redefine as wrappers later) */
#undef FCS_P2NFFT_INTERFACE_WRAPPER_0
#undef FCS_P2NFFT_INTERFACE_WRAPPER_1
#undef FCS_P2NFFT_INTERFACE_WRAPPER_2
#undef FCS_P2NFFT_INTERFACE_WRAPPER_3

/* define macros for definition of wrapper headers */
#define FCS_P2NFFT_INTERFACE_WRAPPER_0(NAME, INAME) \
  FCSResult fcs_p2nfft_ ## NAME(FCS handle);
#define FCS_P2NFFT_INTERFACE_WRAPPER_1(NAME, INAME, TYPE1, ARG1) \
  FCSResult fcs_p2nfft_ ## NAME(FCS handle, TYPE1 ARG1 );
#define FCS_P2NFFT_INTERFACE_WRAPPER_2(NAME, INAME, TYPE1, ARG1, TYPE2, ARG2) \
  FCSResult fcs_p2nfft_ ## NAME(FCS handle, TYPE1 ARG1, TYPE2 ARG2);
#define FCS_P2NFFT_INTERFACE_WRAPPER_3(NAME, INAME, TYPE1, ARG1, TYPE2, ARG2, TYPE3, ARG3) \
  FCSResult fcs_p2nfft_ ## NAME(FCS handle, TYPE1 ARG1, TYPE2 ARG2, TYPE3 ARG3);

/************************************************************
 *     Getter and Setter for P2NFFT, PNFFT, PFFT Parameters 
 ************************************************************/
/* call macros that create the wrapper prototypes */
#include "fcs_p2nfft_wrappers.h"
  
#endif
