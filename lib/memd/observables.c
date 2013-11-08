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

#include "observables.h"
#include "data_types.h"
#include "helper_functions.h"

/********************************************/
/****** get energy and print out stuff ******/
/********************************************/

/** integrates 0.5*D*E over the whole system
 public function!
 @return returns electric energy
 */
fcs_float fcs_memd_electric_energy(void* rawdata)
{
    memd_struct* memd = (memd_struct*) rawdata;
    
    fcs_int x, y, z, i, k;
    fcs_float localresult = 0.;
    fcs_float globalresult = 0.;
	
    FORALL_INNER_SITES(x, y, z) {
        i = ifcs_memd_get_linear_index(x, y, z, memd->lparams.dim);	  
        FOR3D(k){
            localresult += SQR(memd->Dfield[i*3+k]) / memd->lattice[i].permittivity[k];
        }
    }
    localresult *= 0.5*memd->parameters.a;
    MPI_Allreduce(&localresult,&globalresult,1,FCS_MPI_FLOAT,MPI_SUM,memd->mpiparams.communicator);  
    return globalresult;
}


/** Public funxtion.
 Integrates the B-field over the whole system to get the
 energy of the magnetic field.
 @return returns magnetic energy
 */
fcs_float fcs_memd_magnetic_energy(void* rawdata)
{
    memd_struct* memd = (memd_struct*) rawdata;

    fcs_int x, y, z, i;
    fcs_float result = 0.;
    /*  fcs_float invmass = 1./memd->parameters.f_mass; we have B^~=B*c !!!! */
	
    FORALL_INNER_SITES(x, y, z) {
        i = ifcs_memd_get_linear_index(x, y, z, memd->lparams.dim);	  
        result += SQR(memd->Bfield[i*3]) + SQR(memd->Bfield[i*3+1]) + SQR(memd->Bfield[i*3+2]);
    }
    /* B is rescaled !!! ATTENTION!!! */
    result *= 0.5*memd->parameters.a;
    return result;
}

/** print out current setup of maggs method
 @return 0 if successful
 @param interp TCL interpreter handle
 */
/*
fcs_int tclprint_to_result_Maggs(Tcl_Interp *interp)
{
    char buffer[TCL_DOUBLE_SPACE];
	
    Tcl_PrintDouble(interp, memd->parameters.f_mass, buffer);
    Tcl_AppendResult(interp, "maggs ", buffer, " ", (char *) NULL);
    sprintf(buffer,"%d",memd->parameters.mesh);
    Tcl_AppendResult(interp, buffer, " ", (char *) NULL); 
	
    return TCL_OK;
}
 */

