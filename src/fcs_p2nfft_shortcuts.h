/*
  Copyright (C) 2012-2013 Michael Pippig

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

/* clear macro definitions first since we have to define the shortcuts multiple times */
#undef FCS_P2NFFT_SET_GET_WRAPPER_1
#undef FCS_P2NFFT_SET_GET_TUNE_WRAPPER_1
#undef FCS_P2NFFT_SET_GET_WRAPPER_2
#undef FCS_P2NFFT_SET_GET_TUNE_WRAPPER_2
#undef FCS_P2NFFT_SET_GET_WRAPPER_3
#undef FCS_P2NFFT_SET_GET_TUNE_WRAPPER_3

/* define some shortcuts in order to define setter, getter and set_tune at once */
#define FCS_P2NFFT_SET_GET_WRAPPER_1(NAME, INAME, TYPE, ARG) \
  FCS_P2NFFT_INTERFACE_WRAPPER_1(set_ ## NAME, set_ ## INAME, TYPE,  ARG) \
  FCS_P2NFFT_INTERFACE_WRAPPER_1(get_ ## NAME, get_ ## INAME, TYPE*, ARG)
#define FCS_P2NFFT_SET_GET_TUNE_WRAPPER_1(NAME, INAME, TYPE, ARG) \
  FCS_P2NFFT_SET_GET_WRAPPER_1(NAME, INAME, TYPE, ARG) \
  FCS_P2NFFT_INTERFACE_WRAPPER_0(set_ ## NAME ## _tune, set_ ## INAME ## _tune)

#define FCS_P2NFFT_SET_GET_WRAPPER_2(NAME, INAME, TYPE1, ARG1, TYPE2, ARG2) \
  FCS_P2NFFT_INTERFACE_WRAPPER_2(set_ ## NAME, set_ ## INAME, TYPE1,  ARG1, TYPE2,  ARG2) \
  FCS_P2NFFT_INTERFACE_WRAPPER_2(get_ ## NAME, get_ ## INAME, TYPE1*, ARG1, TYPE2*, ARG2)
#define FCS_P2NFFT_SET_GET_TUNE_WRAPPER_2(NAME, INAME, TYPE1, ARG1, TYPE2, ARG2) \
  FCS_P2NFFT_SET_GET_WRAPPER_2(NAME, INAME, TYPE1, ARG1, TYPE2, ARG2) \
  FCS_P2NFFT_INTERFACE_WRAPPER_0(set_ ## NAME ## _tune, set_ ## INAME ## _tune)

#define FCS_P2NFFT_SET_GET_WRAPPER_3(NAME, INAME, TYPE1, ARG1, TYPE2, ARG2, TYPE3, ARG3) \
  FCS_P2NFFT_INTERFACE_WRAPPER_3(set_ ## NAME, set_ ## INAME, TYPE1,  ARG1, TYPE2,  ARG2, TYPE3,  ARG3) \
  FCS_P2NFFT_INTERFACE_WRAPPER_3(get_ ## NAME, get_ ## INAME, TYPE1*, ARG1, TYPE2*, ARG2, TYPE3*, ARG3)
#define FCS_P2NFFT_SET_GET_TUNE_WRAPPER_3(NAME, INAME, TYPE1, ARG1, TYPE2, ARG2, TYPE3, ARG3) \
  FCS_P2NFFT_SET_GET_WRAPPER_3(NAME, INAME, TYPE1, ARG1, TYPE2, ARG2, TYPE3, ARG3) \
  FCS_P2NFFT_INTERFACE_WRAPPER_0(set_ ## NAME ## _tune, set_ ## INAME ## _tune)

