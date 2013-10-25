/*
  Copyright (C) 2011-2012 Rene Halver

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

/*! \mainpage ScaFaCoS - Scalable Fast Coulomb Solvers
 *
 *  The project's main web page is located at [www.scafacos.de](http://www.scafacos.de).
 *
 *
 * \section interface_doc Interface documentation
 *
 * Documentation on the frontend interfaces is found in the following files:
 *
 *   Language      | File
 *   ------------- | -------------
 *   C             | FCSInterface_p.h
 *   Fortran       | fcs4fortran.f90 
 *
 */

#ifndef FCS_INCLUDED
#define FCS_INCLUDED

#include "fcs_config.h"

#include "FCSDefinitions.h"
#include "FCSResult_p.h"
#include "FCSInterface_p.h"

#ifdef __cplusplus
extern "C" {
#endif

#ifdef FCS_ENABLE_DIRECT
#include "fcs_direct_p.h"
#endif
#ifdef FCS_ENABLE_EWALD
#include "fcs_ewald_p.h"
#endif
#ifdef FCS_ENABLE_FMM
#include "fcs_fmm_p.h"
#endif
#ifdef FCS_ENABLE_MEMD
#include "fcs_memd_p.h"
#endif
#ifdef FCS_ENABLE_MMM1D
#include "fcs_mmm1d_p.h"
#endif
#ifdef FCS_ENABLE_MMM2D
#include "fcs_mmm2d_p.h"
#endif
#ifdef FCS_ENABLE_P2NFFT
#include "fcs_p2nfft_p.h"
#endif
#ifdef FCS_ENABLE_P3M
#include "fcs_p3m_p.h"
#endif
#ifdef FCS_ENABLE_PEPC
#include "fcs_pepc_p.h"
#endif
#ifdef FCS_ENABLE_PP3MG
#include "fcs_pp3mg_p.h"
#endif
#ifdef FCS_ENABLE_VMG
#include "fcs_vmg_p.h"
#endif
#ifdef FCS_ENABLE_WOLF
#include "fcs_wolf_p.h"
#endif

#ifdef __cplusplus
}
#endif

#endif
