/*
 * Copyright (C) 2011-2013 Michael Pippig
 * Copyright (C) 2011 Sebastian Banert
 *
 * This file is part of ScaFaCoS.
 * 
 * ScaFaCoS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ScaFaCoS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *	
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _P2NFFT_PARAMETERS_H
#define _P2NFFT_PARAMETERS_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "FCSCommon.h"
#include "fcs_result.h"

/* clear macro definitions of headers (redefine as wrappers later) */
#undef FCS_P2NFFT_INTERFACE_WRAPPER_0
#undef FCS_P2NFFT_INTERFACE_WRAPPER_1
#undef FCS_P2NFFT_INTERFACE_WRAPPER_2
#undef FCS_P2NFFT_INTERFACE_WRAPPER_3

/* Enable definition of multiple wrappers which point to the same internal function.
 * We use this trick to get an P3M-compliant interface but also remain have our own nomenclature. 
 * However, at this point we must avoid duplicated definition of the interal functions. */
#undef FCS_P2NFFT_INTERFACE_WITH_REDIRECTIONS
#define FCS_P2NFFT_INTERFACE_WITH_REDIRECTIONS 0

/* define macros for definition of wrapper headers */
#define FCS_P2NFFT_INTERFACE_WRAPPER_0(NAME, INAME) \
  FCSResult ifcs_p2nfft_ ## INAME(void *rd, const char* fnc_name);
#define FCS_P2NFFT_INTERFACE_WRAPPER_1(NAME, INAME, TYPE1, ARG1) \
  FCSResult ifcs_p2nfft_ ## INAME(void *rd, const char* fnc_name, TYPE1 ARG1);
#define FCS_P2NFFT_INTERFACE_WRAPPER_2(NAME, INAME, TYPE1, ARG1, TYPE2, ARG2) \
  FCSResult ifcs_p2nfft_ ## INAME(void *rd, const char* fnc_name, TYPE1 ARG1, TYPE2 ARG2);
#define FCS_P2NFFT_INTERFACE_WRAPPER_3(NAME, INAME, TYPE1, ARG1, TYPE2, ARG2, TYPE3, ARG3) \
  FCSResult ifcs_p2nfft_ ## INAME(void *rd, const char* fnc_name, TYPE1 ARG1, TYPE2 ARG2, TYPE3 ARG3);

/************************************************************
 *     Getter and Setter for P2NFFT, PNFFT, PFFT Parameters 
 ************************************************************/
/* call macros that create the wrapper prototypes */
#include "fcs_p2nfft_wrappers.h"

#endif
