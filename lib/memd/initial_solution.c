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

#include "initial_solution.h"
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "helper_functions.h"
#include "interpolate.h"
#include <math.h>
#include <stdio.h>
#include "communication.h"
#include "FCSCommon.h"

/***************************************/
/****** initialization procedure: ******/
/***************************************/

/** Calculate the self energy coefficients for the system, if
 corrected with Lattice Green's function. */
void ifcs_memd_calc_self_energy_coeffs(memd_struct* memd)
{
    fcs_float factor, prefac;
    fcs_int px = 0;
    fcs_int py = 0;
    fcs_int pz = 0;
    fcs_int i, j, k, l, m, index;
    fcs_float sx = 0.;
    fcs_float sy = 0.;
    fcs_float sz = 0.;
    fcs_float sxy = 0.;
    fcs_float sxyz = 0.;
    fcs_float nomx, nomy, nomz;
    fcs_float inva = memd->parameters.inva;
    fcs_float invasq = inva*inva;
    fcs_float n[8][SPACE_DIM];// = {{0.,0.,0.}, {1.,0.,0.}, {0.,1.,0.}, {1.,1.,0.}, 
    // {0.,0.,1.}, {1.,0.,1.}, {0.,1.,1.}, {1.,1.,1.}};
    
    m=0;
    l=0;
    index = 0;
    for(i=0;i<2;i++)
        for(j=0;j<2;j++)
            for(k=0;k<2;k++) {
                n[index][2] = m; 
                n[index][1] = l; 
                n[index][0] = k; 
                index++;
            }
    
    factor = M_PI / memd->parameters.mesh;
    prefac = 1. / memd->parameters.mesh;
    prefac = 0.5 * prefac * prefac * prefac * SQR(memd->parameters.prefactor);
    
    for(i=0;i<8;i++) 
    {
        for(j=0;j<8;j++) 
        {
            memd->parameters.alpha[i][j] = 0.;
            for(px = 0; px < memd->parameters.mesh; ++px)
            {
                sx = sin( factor * px );
                sx = sx * sx;
                nomx = 2.*factor*px*(n[i][0] - n[j][0]);
                for(py = 0; py < memd->parameters.mesh; ++py)
                {
                    sy = sin( factor * py );
                    sy = sy * sy; 
                    nomy = 2.*factor*py*(n[i][1] - n[j][1]);
                    sxy = sx + sy;
                    for(pz = 0; pz < memd->parameters.mesh; ++pz)
                    {
                        sz = sin( factor * pz );
                        sz = sz * sz;
                        nomz = 2.*factor*pz*(n[i][2] - n[j][2]);
                        sxyz = sxy + sz;
                        sxyz *= 4.;
                        if(sxyz > 0)
                        {
                            memd->parameters.alpha[i][j] += cos(nomx + nomy + nomz) / sxyz;
                        }
                    }
                }
            }
            /* invasq is needed for the calculation of forces */
            memd->parameters.alpha[i][j] = invasq * prefac * memd->parameters.alpha[i][j];
        }
    }
    
}

/** For energy minimization:
 @return maximum of curl(E) from all sites.
 May also be used to verify constraint surface condition. */
fcs_float ifcs_memd_check_curl_E(memd_struct* memd)
{
    fcs_int i, ix, iy, iz;
    fcs_float curl, maxcurl, gmaxcurl;
    fcs_int* anchor_neighb;
	
    maxcurl = 0.;
	
    FORALL_INNER_SITES(ix, iy, iz) {
        i = ifcs_memd_get_linear_index(ix, iy, iz, memd->lparams.dim); 
        i = ifcs_memd_get_linear_index(ix, iy, iz, memd->lparams.dim); 
        anchor_neighb = memd->neighbor[i];
        curl = memd->Dfield[3*i] + memd->Dfield[3*anchor_neighb[0]+1] 
        - memd->Dfield[3*anchor_neighb[1]] - memd->Dfield[3*i+1];
        curl *= memd->parameters.inva;
        if(fabs(curl)>maxcurl) maxcurl = fabs(curl);
        curl = memd->Dfield[3*i+2] + memd->Dfield[3*anchor_neighb[2]] 
        - memd->Dfield[3*anchor_neighb[0]+2] - memd->Dfield[3*i];
        curl *= memd->parameters.inva;
        if(fabs(curl)>maxcurl) maxcurl = fabs(curl);
        curl = memd->Dfield[3*i+1] + memd->Dfield[3*anchor_neighb[1]+2] 
        - memd->Dfield[3*anchor_neighb[2]+1] - memd->Dfield[3*i+2];
        curl *= memd->parameters.inva;
        if(fabs(curl)>maxcurl) maxcurl = fabs(curl);
    }
    MPI_Allreduce(&maxcurl,&gmaxcurl,1,FCS_MPI_FLOAT,MPI_MAX,memd->mpiparams.communicator);  
    return gmaxcurl;
}

/** If curl(E) is not zero, update plaquette with corrections.
 @param i index of the current lattice site
 @param n coordinate, is the normal direction to the plaquette
 */
void ifcs_memd_perform_rot_move_inplane(memd_struct* memd, fcs_int i, fcs_int n)
{
    fcs_int mue, nue;
    fcs_int * anchor_neighb;
    fcs_float delta;
    fcs_float ROUND_ERR = 0.01*ROUND_ERROR_PREC;
	
    mue = 0; nue = 0;
	
    switch(n) {
        case 0 :
            mue = 1;
            nue = 2;
            break;
			
        case 1 :
            mue = 2;
            nue = 0;
            break;
        case 2 :
            mue = 0;
            nue = 1;
            break;
    }
	
    anchor_neighb = &memd->neighbor[i][0];
	
    delta = memd->Dfield[3*i+mue] + memd->Dfield[3*anchor_neighb[mue]+nue] 
    - memd->Dfield[3*anchor_neighb[nue]+mue] - memd->Dfield[3*i+nue];
    if(fabs(delta)>=ROUND_ERR) {
        delta = -delta/4.; 
        ifcs_memd_update_plaquette(mue, nue, anchor_neighb, i, delta);
    }
}


/** For energy minimization:
 Relax B-fields in all directions.*/
void ifcs_memd_minimize_transverse_field(memd_struct* memd)
{
    fcs_int k, l, m;
    fcs_int i, d;
    fcs_int ind_i, ind_j;
    fcs_int size[2]={0,0};
    fcs_int index = -1; /* force allocation error */
	
    FOR3D(d) {
        switch(d) {
            case 0 :
                size[0] = memd->lparams.size[2];
                size[1] = memd->lparams.size[1];
                break;
            case 1 :
                size[0] = memd->lparams.size[2];
                size[1] = memd->lparams.size[0];
                break;
            case 2 :
                size[0] = memd->lparams.size[1];
                size[1] = memd->lparams.size[0];
                break;
        }
        for(i=0;i<2;i++) {
            /* at first even sites (i==0) then odd */
            for(k=1;k<=memd->lparams.size[d];k++) {
                /* update every plane in direction d */
                ind_i=0;
                for(l=0; l<=size[1]; l++){
                    ind_j=0;
                    for(m=0;m<=size[0]; m++) {
                        switch(d) {
                            case 0 :
                                index=ifcs_memd_get_linear_index(k,l,m,memd->lparams.dim);
                                break;
                            case 1 :
                                index=ifcs_memd_get_linear_index(l,k,m,memd->lparams.dim); 
                                break;
                            case 2 :
                                index=ifcs_memd_get_linear_index(l,m,k,memd->lparams.dim); 
                                break;
                        }
                        if((ind_i+ind_j)%2==i)
                            ifcs_memd_perform_rot_move_inplane(memd, index, d);
                        ind_j++;
                    }
                    ind_i++;
                }   
            }
            /* update boundaries - update halo regions */
            fcs_memd_exchange_surface_patch(memd, memd->Dfield, 3, 0);
        }
    }
}


/** calculates initial electric field configuration.
 currently uses simple and slow method of plaquettes and links.
 energy minimization takes up lots of time. */
void ifcs_memd_calc_init_e_field(memd_struct* memd)
{
    fcs_int xsizeplus, ysizeplus;
    fcs_float localqy, localqz;
    fcs_float qplane, qline;
    fcs_int    i, k, ix, iy, iz;
    fcs_int index = 0;
    fcs_float invasq, tmp_field=0.0;
    fcs_float sqrE, gsqrE, goldE, gavgEx, gavgEy, gavgEz;
    fcs_float qz, qy, qx, avgEx, avgEy, avgEz;
    fcs_float Eall[SPACE_DIM], gEall[SPACE_DIM];
    fcs_float maxcurl;
    MPI_Status status;
    MPI_Comm zplane, yline;
    fcs_int color, rank, dim;

    invasq = memd->parameters.inva*memd->parameters.inva;
	
    xsizeplus = memd->lparams.dim[0];
    ysizeplus = memd->lparams.dim[1];
    
    ifcs_memd_distribute_particle_charges(memd);
	
    dim = memd->mpiparams.node_grid[1]*memd->mpiparams.node_grid[0];
    color = memd->mpiparams.node_pos[2];
    rank  = memd->mpiparams.this_node%dim;
    MPI_Comm_split(memd->mpiparams.communicator, color, rank, &zplane);
    color = memd->mpiparams.node_pos[1];
    rank  = rank%memd->mpiparams.node_grid[0];
    MPI_Comm_split(zplane, color, rank, &yline);


    /* calculate initial solution of Poisson equation */
	
    /* CAUTION: the indexing of the neighbor nodes in Espresso
     *  starts from x left neighbor node
     */
	
    /* get process coordinates */
    if(memd->mpiparams.node_pos[2]!= 0) {
        MPI_Recv(&tmp_field, 1, FCS_MPI_FLOAT, memd->mpiparams.node_neighbors[4], REQ_MAGGS_EQUIL, memd->mpiparams.communicator, &status);
        for(iy=memd->lparams.inner_left_down[1];iy<memd->lparams.inner_up_right[1];iy++) {
            for(ix=memd->lparams.inner_left_down[0];ix<memd->lparams.inner_up_right[0];ix++) {  
                index = ifcs_memd_get_linear_index(ix, iy, memd->lparams.inner_left_down[2], memd->lparams.dim);
                memd->Dfield[3*memd->neighbor[index][ZMINUS]+ZPLUS] = tmp_field;
            }
        }
    }
    
    localqz = 0.;
    for(iz=memd->lparams.inner_left_down[2];iz<memd->lparams.inner_up_right[2];iz++) { /* loop over z-planes */
        localqz = 0.;
        for(iy=memd->lparams.inner_left_down[1];iy<memd->lparams.inner_up_right[1];iy++) {
            for(ix=memd->lparams.inner_left_down[0];ix<memd->lparams.inner_up_right[0];ix++) {  
                index = ifcs_memd_get_linear_index(ix, iy, iz, memd->lparams.dim);
                /* Sum over the charge of all sides in z-plane */
                localqz += memd->lattice[index].charge / memd->lattice[index].permittivity[2];
            }
        } 
		
        MPI_Allreduce(&localqz, &qz, 1, FCS_MPI_FLOAT, MPI_SUM, zplane);
        qz = qz/(memd->parameters.mesh*memd->parameters.mesh);
        qplane = qz*memd->parameters.prefactor*invasq;
        /*    if(fabs(qplane) >= 0.01*ROUND_ERROR_PREC) { */
        for(iy=memd->lparams.inner_left_down[1];iy<memd->lparams.inner_up_right[1];iy++) {
            for(ix=memd->lparams.inner_left_down[0];ix<memd->lparams.inner_up_right[0];ix++) {  
                index = ifcs_memd_get_linear_index(ix, iy, iz, memd->lparams.dim);
                memd->Dfield[3*index+ZPLUS]  = memd->Dfield[3*memd->neighbor[index][ZMINUS]+ZPLUS] + qplane;
                /*	    + qz*memd->parameters.prefactor*invasq; */
            }
        }
        /*    } */
        if(iz>=memd->lparams.inner_up_right[2]-1) {
            if (memd->mpiparams.node_pos[2]<memd->mpiparams.node_grid[2]-1) {
                if(memd->mpiparams.node_grid[2]>1) {
                    MPI_Send(&memd->Dfield[3*index+ZPLUS], 1, FCS_MPI_FLOAT, memd->mpiparams.node_neighbors[5], REQ_MAGGS_EQUIL, memd->mpiparams.communicator); 
                }
            }
            else 
                if (fabs(memd->Dfield[3*index+ZPLUS]) > 100.*ROUND_ERROR_PREC) {
                    fprintf(stderr, "%d: Error in the calculation of Ez(%d,%d,%d)=%f!!\n", 
                            memd->mpiparams.this_node,memd->lattice[index].r[0], memd->lattice[index].r[1], memd->lattice[index].r[2],
                            memd->Dfield[3*index+ZPLUS]);
                    fflush(stderr);
                }
        }
		
        if(memd->mpiparams.node_pos[1]!= 0) {
            MPI_Recv(&tmp_field, 1, FCS_MPI_FLOAT, memd->mpiparams.node_neighbors[2], REQ_MAGGS_EQUIL, memd->mpiparams.communicator, &status);
            for(ix=memd->lparams.inner_left_down[0];ix<memd->lparams.inner_up_right[0];ix++) {  
                index = ifcs_memd_get_linear_index(ix, memd->lparams.inner_left_down[1], iz, memd->lparams.dim);
                memd->Dfield[3*memd->neighbor[index][YMINUS]+YPLUS] = tmp_field;
            }
        }
		
        for(iy=memd->lparams.inner_left_down[1];iy<memd->lparams.inner_up_right[1];iy++) {
            localqy = 0.;
            for(ix=memd->lparams.inner_left_down[0];ix<memd->lparams.inner_up_right[0];ix++) {  
                index = ifcs_memd_get_linear_index(ix, iy, iz, memd->lparams.dim);
                localqy += memd->lattice[index].charge / memd->lattice[index].permittivity[1];
            }
			
            MPI_Allreduce(&localqy, &qy, 1, FCS_MPI_FLOAT, MPI_SUM, yline);
			
            qy = qy/memd->parameters.mesh;
            qline = (qy-qz)*memd->parameters.prefactor*invasq;
            /*      if(fabs(qy-qz)>=ROUND_ERROR_PREC) { */
            for(ix=memd->lparams.inner_left_down[0];ix<memd->lparams.inner_up_right[0];ix++) {  
                index = ifcs_memd_get_linear_index(ix, iy, iz, memd->lparams.dim);
                memd->Dfield[3*index+YPLUS]  = memd->Dfield[3*memd->neighbor[index][YMINUS]+YPLUS] + qline;
                /*	    (qy-qz)*memd->parameters.prefactor*invasq; */
            }
            /*      } */
			
            if(iy>=memd->lparams.inner_up_right[1]-1) {
                if(memd->mpiparams.node_pos[1] < memd->mpiparams.node_grid[1]-1) {
                    if (memd->mpiparams.node_grid[1]>1)
                        MPI_Send(&memd->Dfield[3*index+YPLUS], 1, FCS_MPI_FLOAT, memd->mpiparams.node_neighbors[3], REQ_MAGGS_EQUIL, memd->mpiparams.communicator); 
                }
                else
                    if (fabs(memd->Dfield[3*index+YPLUS]) > 100.*ROUND_ERROR_PREC)
                        fprintf(stderr, "%d: Error in the calculation of Ey(%d,%d,%d)=%f!!\n",
                                memd->mpiparams.this_node, memd->lattice[index].r[0], memd->lattice[index].r[1], memd->lattice[index].r[2],
                                memd->Dfield[3*index+YPLUS]);	  
            }
			
            if(memd->mpiparams.node_pos[0]!= 0) {
                MPI_Recv(&tmp_field, 1, FCS_MPI_FLOAT, memd->mpiparams.node_neighbors[0], REQ_MAGGS_EQUIL, memd->mpiparams.communicator, &status);
                index = ifcs_memd_get_linear_index(memd->lparams.inner_left_down[0], iy, iz, memd->lparams.dim);
                memd->Dfield[3*memd->neighbor[index][XMINUS]+XPLUS] = tmp_field;
            }
			
            for(ix=memd->lparams.inner_left_down[0];ix<memd->lparams.inner_up_right[0];ix++) {  
                index = ifcs_memd_get_linear_index(ix, iy, iz, memd->lparams.dim);
                qx = memd->lattice[index].charge / memd->lattice[index].permittivity[0]; 
                memd->Dfield[3*index+XPLUS] = memd->Dfield[3*memd->neighbor[index][XMINUS]+XPLUS] + 
                (qx-qy)*memd->parameters.prefactor*invasq;
            }
			
            if(ix>=memd->lparams.inner_up_right[0]-1) {
                if(memd->mpiparams.node_pos[0] < memd->mpiparams.node_grid[0]-1) {
                    if(memd->mpiparams.node_grid[0]>1)
                        MPI_Send(&memd->Dfield[3*index+XPLUS], 1, FCS_MPI_FLOAT, memd->mpiparams.node_neighbors[1], REQ_MAGGS_EQUIL, memd->mpiparams.communicator); 
                }
                else
                    if (fabs(memd->Dfield[3*index+XPLUS]) > 100.*ROUND_ERROR_PREC)
                        fprintf(stderr, "%d: Error in the calculation of Ex(%d,%d,%d)=%f!!\n",
                                memd->mpiparams.this_node, memd->lattice[index].r[0], memd->lattice[index].r[1], memd->lattice[index].r[2],
                                memd->Dfield[3*index+XPLUS]);	  
            }
        }    /*** loop over iy */
    }
	
    /* exchange halo-surfaces */
    fcs_memd_exchange_surface_patch(memd, memd->Dfield, 3, 1); 
	
    avgEz = 0.;
    for(iz=memd->lparams.inner_left_down[2];iz<memd->lparams.inner_up_right[2];iz++) {
        index = ifcs_memd_get_linear_index(memd->lparams.inner_left_down[0], memd->lparams.inner_left_down[1], iz, memd->lparams.dim);
        avgEz += memd->Dfield[3*index+ZPLUS];
    }
	
	
    /*  MPI_Barrier(memd->mpiparams.communicator); */
    MPI_Allreduce(&avgEz,&gavgEz,1,FCS_MPI_FLOAT,MPI_SUM,memd->mpiparams.communicator);
    gavgEz = gavgEz/(memd->parameters.mesh*memd->mpiparams.node_grid[0]*memd->mpiparams.node_grid[1]);
	
    FORALL_INNER_SITES(ix, iy,iz) {
        index = ifcs_memd_get_linear_index(ix, iy, iz, memd->lparams.dim);
        memd->Dfield[3*index+ZPLUS] -= gavgEz;
    }
	
    for(iz = memd->lparams.inner_left_down[2];iz<memd->lparams.inner_up_right[2];iz++) {
        avgEy = 0.;  
        for(iy = memd->lparams.inner_left_down[1];iy<memd->lparams.inner_up_right[1];iy++) {
            index = ifcs_memd_get_linear_index(memd->lparams.inner_left_down[0], iy, iz, memd->lparams.dim);
            avgEy += memd->Dfield[3*index+YPLUS];
        }    
		
        MPI_Allreduce(&avgEy, &gavgEy, 1, FCS_MPI_FLOAT, MPI_SUM, zplane);
        gavgEy = gavgEy/(memd->parameters.mesh*memd->mpiparams.node_grid[0]);
		
        for(iy=memd->lparams.inner_left_down[1];iy<memd->lparams.inner_up_right[1];iy++) {
            for(ix=memd->lparams.inner_left_down[0];ix<memd->lparams.inner_up_right[0];ix++)  
                memd->Dfield[3*ifcs_memd_get_linear_index(ix, iy, iz, memd->lparams.dim)+YPLUS] -= gavgEy;
        }
    }
	
    for(iz=memd->lparams.inner_left_down[2];iz<memd->lparams.inner_up_right[2];iz++) {
        for(iy=memd->lparams.inner_left_down[1];iy<memd->lparams.inner_up_right[1];iy++) {
            avgEx = 0.;
            for(ix=memd->lparams.inner_left_down[0];ix<memd->lparams.inner_up_right[0];ix++) {
                avgEx += memd->Dfield[3*ifcs_memd_get_linear_index(ix, iy, iz, memd->lparams.dim)+XPLUS];
            }
			
            MPI_Allreduce(&avgEx, &gavgEx, 1, FCS_MPI_FLOAT, MPI_SUM, yline);
            gavgEx = gavgEx/memd->parameters.mesh;
			
            for(ix=memd->lparams.inner_left_down[0];ix<memd->lparams.inner_up_right[0];ix++)  
                memd->Dfield[3*ifcs_memd_get_linear_index(ix, iy, iz, memd->lparams.dim)+XPLUS] -= gavgEx;
        }
    }
	
	
    /* exchange halo-surfaces */
    fcs_memd_exchange_surface_patch(memd, memd->Dfield, 3, 1);
	
    
	
    MPI_Comm_free(&zplane);
    MPI_Comm_free(&yline);
	
	
	
	
    /* iterative procedure of energy minimization */
    
    sqrE = 0.;
    FORALL_INNER_SITES(ix, iy, iz) {
        i = ifcs_memd_get_linear_index(ix, iy, iz, memd->lparams.dim);
        FOR3D(k) sqrE += SQR(memd->Dfield[3*i+k]);
    }
    
    MPI_Allreduce(&sqrE,&gsqrE,1,FCS_MPI_FLOAT,MPI_SUM,memd->mpiparams.communicator); 
    gsqrE = gsqrE/(SPACE_DIM*memd->parameters.mesh*memd->parameters.mesh*memd->parameters.mesh);  
    
    do {
        goldE = gsqrE;
        sqrE = 0.;
        ifcs_memd_minimize_transverse_field(memd);
		
        FORALL_INNER_SITES(ix, iy, iz) {
            i = ifcs_memd_get_linear_index(ix, iy, iz, memd->lparams.dim);
            FOR3D(k) sqrE += SQR(memd->Dfield[3*i+k]);
        }
        MPI_Allreduce(&sqrE,&gsqrE,1,FCS_MPI_FLOAT,MPI_SUM,memd->mpiparams.communicator); 
        gsqrE = gsqrE/(SPACE_DIM*memd->parameters.mesh*memd->parameters.mesh*memd->parameters.mesh);  
        maxcurl = ifcs_memd_check_curl_E(memd);
		
		
    } while(fabs(maxcurl)>1000000.*ROUND_ERROR_PREC);
	
    
	
	
	
	
	
	
	
    /* exchange halo-surfaces */
	
    FOR3D(k) Eall[k] = 0.;
    FORALL_INNER_SITES(ix, iy, iz) {
        i = ifcs_memd_get_linear_index(ix, iy, iz, memd->lparams.dim);
        FOR3D(k) {
            Eall[k] += memd->Dfield[3*i+k];
        }
    }
	
    MPI_Allreduce(Eall,gEall,3,FCS_MPI_FLOAT,MPI_SUM,memd->mpiparams.communicator);
	
    FOR3D(k) gEall[k] /= (memd->mpiparams.node_grid[0]*memd->mpiparams.node_grid[1]*memd->mpiparams.node_grid[2]*memd->lparams.size[0]*memd->lparams.size[1]*memd->lparams.size[2]);
	
    FORALL_INNER_SITES(ix, iy, iz) {
        i = ifcs_memd_get_linear_index(ix, iy, iz, memd->lparams.dim);
        FOR3D(k) memd->Dfield[3*i+k] -= gEall[k];
    }
    /* exchange whole glue-patch region */
    fcs_memd_exchange_surface_patch(memd, memd->Dfield, 3, 0);
/*
    if(!memd->mpiparams.this_node)
        MAGGS_TRACE(fprintf(stderr, "Ex = %16.12e, Ey = %15.12e, Ez = %15.12e\n", gEall[0], gEall[1], gEall[2]));
*/
}



