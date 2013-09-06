/*
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

#include "FCSDefinitions.h"

#define FCS_P2NFFT_SQRT2     FCS_CONST(1.4142135623730950488016887242096981)  /* sqrt(2) */
#define FCS_P2NFFT_PI        FCS_CONST(3.1415926535897932384626433832795029)  /* pi */
#define FCS_P2NFFT_1_SQRTPI  FCS_CONST(0.5641895835477562869480794515607726)  /* 1/sqrt(pi) */
#define FCS_P2NFFT_PISQR     FCS_CONST(9.8696044010893586188344909998761511)  /* pi^2 */
#define FCS_P2NFFT_EULER     FCS_CONST(0.5772156649015328606065120900824024)  /* Euler-Mascheroni constant */
#define FCS_P2NFFT_BRILLOUIN 0
#define FCS_P2NFFT_MINGRID 4
#define FCS_P2NFFT_MAXGRID 2048
#define FCS_P2NFFT_MAXCAO 7
#define FCS_P2NFFT_FULL_ESTIMATE_ALPHA_H_THRESHOLD 0.5;

