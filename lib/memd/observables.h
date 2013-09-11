/*
 Copyright (C) 2010/2011/2012 Florian Fahrenberger
 
 This file is part of ScaFaCoS.
 
 ScaFaCoS is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 ScaFaCoS is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>. 
 */

#ifndef _MEMD_OBSERVABLES_H
#define _MEMD_OBSERVABLES_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/** integrates 0.5*D*E over the whole system
 public function!
 @return returns electric energy
 */
fcs_float fcs_memd_electric_energy(void* rawdata);

/** Public funxtion.
 Integrates the B-field over the whole system to get the
 energy of the magnetic field.
 @return returns magnetic energy
 */
fcs_float fcs_memd_magnetic_energy(void* rawdata);

/** print out current setup of maggs method
 @return 0 if successful
 @param interp TCL interpreter handle
 */
//fcs_int tclprint_to_result_Maggs(Tcl_Interp *interp);

#endif