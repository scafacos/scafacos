/*
 Copyright (C) 2010,2011,2012,2013 Florian Fahrenberger
 
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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "FCSCommon.h"
#include "communication.h"
#include "helper_functions.h"

void fcs_memd_setup_communicator(memd_struct* memd, MPI_Comm communicator)
{
    /* store given communicator */
    memd->mpiparams.original_comm = communicator;
    MPI_Comm_size(communicator, &memd->mpiparams.size);

    /* Test whether the communicator is cartesian and correct dimensionality */
    int comm_is_cart = 0;
    int status;

    MPI_Topo_test(communicator, &status);
    if (status == MPI_CART) {
        /* Communicator is cartesian, so test dimensionality */
        int ndims;
        MPI_Cartdim_get(communicator, &ndims);
        if (ndims == 3) {
            /* Correct dimensionality, so get grid and test periodicity */
            int periodicity[3];
            MPI_Cart_get(communicator, 3, memd->mpiparams.node_grid, periodicity, memd->mpiparams.node_pos);
            if (periodicity[0] && periodicity[1] && periodicity[2]) {
                /* If periodicity is correct, we can just use this communicator */
                memd->mpiparams.communicator = communicator;
                /* get the rank */
                MPI_Comm_rank(communicator, &memd->mpiparams.this_node);
                comm_is_cart = 1;
            }
        }
    }
    
    /* otherwise, we have to set up the cartesian communicator */
    if (!comm_is_cart) {        
        memd->mpiparams.node_grid[0] = 0.0;
        memd->mpiparams.node_grid[1] = 0.0;
        memd->mpiparams.node_grid[2] = 0.0;
        
        /* compute node grid */
        MPI_Dims_create(memd->mpiparams.size, 3, memd->mpiparams.node_grid);
        /* swap first and last dimension, as MEMD currently wants to have them increasing */
        fcs_int tmp = memd->mpiparams.node_grid[2];
        memd->mpiparams.node_grid[2] = memd->mpiparams.node_grid[0];
        memd->mpiparams.node_grid[0] = tmp;
        
        /* create communicator */
        int periodicity[3] = {1, 1, 1};
        MPI_Cart_create(memd->mpiparams.original_comm, 3, memd->mpiparams.node_grid, periodicity, 1, &memd->mpiparams.communicator);
        
        /* get the rank */
        MPI_Comm_rank(communicator, &memd->mpiparams.this_node);
        /* get node pos */
        MPI_Cart_coords(memd->mpiparams.communicator, memd->mpiparams.this_node, 3, memd->mpiparams.node_pos);
    }
    
    /* fetch neighborhood info */
    for (int dir = 0; dir<3; dir++) {
        MPI_Cart_shift(memd->mpiparams.communicator, dir, 1, 
                       &memd->mpiparams.node_neighbors[2*dir], 
                       &memd->mpiparams.node_neighbors[2*dir+1]);
    }
    
    /* init local points */
    for (int i=0; i< 3; i++) {
//        memd->mpiparams.local_box_l[i] = 0.0;
        memd->mpiparams.my_left[i] = 0.0;
        memd->mpiparams.my_right[i] = 0.0;
    }
    
    /* compute box limits */
    fcs_float local_box_length = 0.0;
    for(fcs_int i = 0; i < 3; i++) {
        local_box_length = memd->parameters.box_length[i] / (fcs_float)memd->mpiparams.node_grid[i];
        memd->mpiparams.my_left[i]   = memd->mpiparams.node_pos[i] * local_box_length;
        memd->mpiparams.my_right[i]  = (memd->mpiparams.node_pos[i]+1) * local_box_length;
    }

}


/** sets up nearest neighbors for each site in linear index */
void ifcs_memd__setup_neighbors(memd_struct* memd)
{
    fcs_int ix = 0;
    fcs_int iy = 0;
    fcs_int iz = 0;
	
    fcs_int xsize = memd->lparams.dim[0];
    fcs_int ysize = memd->lparams.dim[1];
    fcs_int zsize = memd->lparams.dim[2];
	
    fcs_int ixplus = 0;
    fcs_int ixminus = 0;
    fcs_int iyplus = 0;
    fcs_int iyminus = 0;
    fcs_int izplus = 0;
    fcs_int izminus = 0;
	
    fcs_int kount = 0;
	
    fcs_int kountxplus = 0;
    fcs_int kountxminus = 0;
    fcs_int kountyplus = 0;
    fcs_int kountyminus = 0;
    fcs_int kountzplus = 0;
    fcs_int kountzminus = 0;
	
    for (ix = 0; ix < xsize; ix++) 
    {
        ixplus  = ix + 1;
        ixminus = ix - 1;
        for(iy = 0; iy < ysize; iy ++)
        {
            iyplus  = iy + 1;
            iyminus = iy - 1;
            for(iz = 0; iz < zsize; iz ++)
            {
                izplus  = iz + 1;
                izminus = iz - 1;
				
                kount         = ifcs_memd_get_linear_index(ix,      iy,      iz,      memd->lparams.dim);
                kountzplus    = ifcs_memd_get_linear_index(ix,      iy,      izplus,  memd->lparams.dim);
                kountzminus   = ifcs_memd_get_linear_index(ix,      iy,      izminus, memd->lparams.dim);
                kountyplus    = ifcs_memd_get_linear_index(ix,      iyplus,  iz,      memd->lparams.dim);
                kountyminus   = ifcs_memd_get_linear_index(ix,      iyminus, iz,      memd->lparams.dim);
                kountxplus    = ifcs_memd_get_linear_index(ixplus,  iy,      iz,      memd->lparams.dim);
                kountxminus   = ifcs_memd_get_linear_index(ixminus, iy,      iz,      memd->lparams.dim);
				
                if(ixminus < 0)     memd->neighbor[kount][XMINUS] = -1;
                else                memd->neighbor[kount][XMINUS] = kountxminus;
                if(ixplus >= xsize) memd->neighbor[kount][XPLUS]  = -1;
                else                memd->neighbor[kount][XPLUS]  = kountxplus;
				
                if(iyminus < 0)     memd->neighbor[kount][YMINUS] = -1;
                else                memd->neighbor[kount][YMINUS] = kountyminus;
                if(iyplus >= ysize) memd->neighbor[kount][YPLUS]  = -1;
                else                memd->neighbor[kount][YPLUS]  = kountyplus;
				
                if(izminus < 0)     memd->neighbor[kount][ZMINUS] = -1;
                else                memd->neighbor[kount][ZMINUS] = kountzminus;
                if(izplus >= zsize) memd->neighbor[kount][ZPLUS]  = -1;
                else                memd->neighbor[kount][ZPLUS]  = kountzplus;
            }
        }
    }
    return;
}


/** Set up lattice, calculate dimensions and lattice parameters
 Allocate memory for lattice sites and fields */
void fcs_memd_setup_local_lattice(memd_struct* memd)
{
    fcs_int i;
    fcs_int ix = 0;
    fcs_int iy = 0;
    fcs_int iz = 0;
    fcs_int linearindex = 0;
    fcs_int xyzcube;	
        
    xyzcube = 1;
    FOR3D(i) {
        /** inner left down grid point (global index) */
        memd->lparams.inner_left_down[i] = (fcs_int)ceil(memd->mpiparams.my_left[i]*memd->parameters.inva); 
        /** inner up right grid point (global index) */
        memd->lparams.inner_up_right[i] = (fcs_int)floor(memd->mpiparams.my_right[i]*memd->parameters.inva); 
        /** correct roundof errors at boundary */
        if(memd->mpiparams.my_right[i]*memd->parameters.inva-memd->lparams.inner_up_right[i]<ROUND_ERROR_PREC) memd->lparams.inner_up_right[i]--;
        if(1.0+memd->mpiparams.my_left[i]*memd->parameters.inva-memd->lparams.inner_left_down[i]<ROUND_ERROR_PREC) memd->lparams.inner_left_down[i]--;
        /** inner grid dimensions */
        memd->lparams.size[i] = memd->lparams.inner_up_right[i] - memd->lparams.inner_left_down[i] + 1;
        /** spacial position of left down grid point */
        memd->lparams.left_down_position[i] = memd->mpiparams.my_left[i] - memd->parameters.a;  
        /** spacial position of upper right grid point */
        memd->lparams.upper_right_position[i] = memd->mpiparams.my_right[i] + memd->parameters.a;  
        /** left down margin */
        memd->lparams.margin[i*2] = 1;
        /** up right margin */
        memd->lparams.margin[(i*2)+1] = 1;
		
        memd->lparams.dim[i] = memd->lparams.size[i] + memd->lparams.margin[i*2] + memd->lparams.margin[i*2+1];
        xyzcube *= memd->lparams.dim[i];
        /** reduce inner grid indices from global to local */
        memd->lparams.inner_left_down[i] = memd->lparams.margin[i*2];
        memd->lparams.inner_up_right[i] = memd->lparams.margin[i*2]+memd->lparams.size[i];
        memd->lparams.halo_left_down[i] = 0;
        memd->lparams.halo_upper_right[i] = memd->lparams.inner_up_right[i];      
    }
	

    memd->lparams.volume    = xyzcube;
    fcs_int ghost_cube = (memd->lparams.dim[0] + 1) *
        (memd->lparams.dim[1] + 1) *
        (memd->lparams.dim[2] + 1) -
        xyzcube;
    
    /** allocate memory for sites and neighbors */
    memd->lattice  = (t_site*) malloc(xyzcube*sizeof(t_site));
    memd->neighbor = (t_dirs*) malloc(xyzcube*sizeof(t_dirs));
    memd->local_cells = *(memd_cell_list*) calloc(1, sizeof(memd_cell_list));
    memd->ghost_cells = *(memd_cell_list*) calloc(1, sizeof(memd_cell_list));
    memd->local_cells.cell = (memd_cell**) calloc(xyzcube, sizeof(memd_cell));
    memd->ghost_cells.cell = (memd_cell**) calloc(ghost_cube, sizeof(memd_cell));
    
    /** allocate memory for cell contents */
    for (int cid=0; cid<xyzcube;cid++) {
        memd->local_cells.cell[cid] = (memd_cell*) calloc(1, sizeof(memd_cell));
        memd->local_cells.cell[cid]->part = (memd_particle*) calloc(1, sizeof(memd_particle));
    }
    for (int cid=0; cid<ghost_cube;cid++) {
        memd->ghost_cells.cell[cid] = (memd_cell*) calloc(1, sizeof(memd_cell));
        memd->ghost_cells.cell[cid]->part = (memd_particle*) calloc(1, sizeof(memd_particle));
    }
	
//    printf("Setting up lattice %d\n", xyzcube); fflush(stdout);
    
    memd->Bfield   = (fcs_float*) malloc(3*xyzcube*sizeof(fcs_float));
    memd->Dfield   = (fcs_float*) malloc(3*xyzcube*sizeof(fcs_float));

                                             
                                             
    memd->local_cells.n = xyzcube;
    memd->ghost_cells.n = ghost_cube;

    /** set up lattice sites */
    FORALL_SITES(ix, iy, iz) {
        linearindex = ifcs_memd_get_linear_index(ix, iy, iz, memd->lparams.dim);
		
        memd->lattice[linearindex].r[0] = ix;
        memd->lattice[linearindex].r[1] = iy;
        memd->lattice[linearindex].r[2] = iz;
        
        FOR3D(i) {
            memd->Bfield[3*linearindex+i]  = 0.;
            memd->Dfield[3*linearindex+i]  = 0.;
        }
        
        memd->lattice[linearindex].charge = 0.;
		
        /* Here, we need a function to set PERMITTIVITY!!!!! */
        FOR3D(i) memd->lattice[linearindex].permittivity[i] = 1.;
    }
	
    ifcs_memd__setup_neighbors(memd);
}


/** sets up surface patches for all domains.
 @param surface_patch the local surface patch
 */
void ifcs_memd__calc_surface_patches(memd_struct* memd, t_surf_patch* surface_patch)
{
    /* x=memd->lparams.size[0] plane */
    surface_patch[0].offset   = memd->lparams.dim[2]*memd->lparams.dim[1]*memd->lparams.size[0];    /*(size[0],0,0) point */
    surface_patch[0].doffset  = 0;                                             /*(0,0,0) point */
    surface_patch[0].stride   = memd->lparams.dim[2]*memd->lparams.dim[1];
    surface_patch[0].skip     = 0;
    surface_patch[0].nblocks  = 1;
    surface_patch[0].coord[0] = 2;
    surface_patch[0].coord[1] = 1;
    surface_patch[0].volume   = memd->lparams.dim[2]*memd->lparams.dim[1];
	
    /* x=1 plane */
    surface_patch[1].offset   = memd->lparams.dim[2]*memd->lparams.dim[1];                    /*(1,0,0) point */
    surface_patch[1].doffset  = memd->lparams.dim[2]*memd->lparams.dim[1]*memd->lparams.inner_up_right[0];    /*(halo[0],0,0) point */
    surface_patch[1].stride   = memd->lparams.dim[2]*memd->lparams.dim[1];
    surface_patch[1].skip     = 0;
    surface_patch[1].nblocks  = 1;
    surface_patch[1].coord[0] = 2;
    surface_patch[1].coord[1] = 1;
    surface_patch[1].volume   = memd->lparams.dim[2]*memd->lparams.dim[1]; 
	
    /* y=memd->lparams.size[1] plane */
    surface_patch[2].offset   = memd->lparams.dim[2]*memd->lparams.size[1];               /*(0,size[1],0) point */
    surface_patch[2].doffset  = 0;                                          /*(0,0,0) point */
    surface_patch[2].stride   = memd->lparams.dim[2];
    surface_patch[2].skip     = memd->lparams.dim[2]*memd->lparams.dim[1];
    surface_patch[2].nblocks  = memd->lparams.dim[0];  
    surface_patch[2].coord[0] = 2;
    surface_patch[2].coord[1] = 0;
    surface_patch[2].volume   = memd->lparams.dim[2]*memd->lparams.dim[0];
	
    /* y=1 plane */
    surface_patch[3].offset   = memd->lparams.dim[2];                             /*(0,1,0) point */
    surface_patch[3].doffset  = memd->lparams.dim[2]*memd->lparams.inner_up_right[1];             /*(0,inner_up_right[1],0) point */
    surface_patch[3].stride   = memd->lparams.dim[2];
    surface_patch[3].skip     = memd->lparams.dim[2]*memd->lparams.dim[1];
    surface_patch[3].nblocks  = memd->lparams.dim[0];
    surface_patch[3].coord[0] = 2;
    surface_patch[3].coord[1] = 0;
    surface_patch[3].volume   = memd->lparams.dim[2]*memd->lparams.dim[0];
	
    /* z=memd->lparams.size[2] plane */
    surface_patch[4].offset   = memd->lparams.size[2];    /*(0,0,size[2]) point */
    surface_patch[4].doffset  = 0;                 /*(0,0,0) point */
    surface_patch[4].stride   = 1;
    surface_patch[4].skip     = memd->lparams.dim[2];
    surface_patch[4].nblocks  = memd->lparams.dim[0]*memd->lparams.dim[1];
    surface_patch[4].coord[0] = 1;
    surface_patch[4].coord[1] = 0;
    surface_patch[4].volume   = memd->lparams.dim[0]*memd->lparams.dim[1];
	
    /* z=1 plane for z it must be higher*/
    surface_patch[5].offset   = 1;                   /*(0,0,1) point */
    surface_patch[5].doffset  = memd->lparams.inner_up_right[2];     /*(0,0,inner_up_right[2]) point */
    surface_patch[5].stride   = 1;
    surface_patch[5].skip     = memd->lparams.dim[2];
    surface_patch[5].nblocks  = memd->lparams.dim[0]*memd->lparams.dim[1];
    surface_patch[5].coord[0] = 1;
    surface_patch[5].coord[1] = 0;
    surface_patch[5].volume   = memd->lparams.dim[0]*memd->lparams.dim[1];
}


/** sets up MPI communications for domain surfaces */
void ifcs_memd__prepare_surface_planes(fcs_int dim, MPI_Datatype *xy, MPI_Datatype *xz, MPI_Datatype *yz, 
                                  t_surf_patch *surface_patch)
{
    MPI_Type_contiguous(dim*surface_patch[0].stride*sizeof(fcs_float),MPI_BYTE,yz);  
    MPI_Type_commit(yz);
    MPI_Type_vector(surface_patch[4].nblocks, dim*surface_patch[4].stride,
                    dim*surface_patch[4].skip, FCS_MPI_FLOAT,xy);
    MPI_Type_commit(xy); 
    MPI_Type_vector(surface_patch[2].nblocks, dim*surface_patch[2].stride,
                    dim*surface_patch[2].skip, FCS_MPI_FLOAT,xz);
    MPI_Type_commit(xz);
}







/*****************************************/
/****** Surface patch communication ******/
/*****************************************/

/** MPI communication of surface region.
 works for D- and B-fields.
 @param field   Field to communicate. Can be B- or D-field.
 @param dim     Dimension in which to communicate
 @param e_equil Flag if field is already equilibated
 */
void fcs_memd_exchange_surface_patch(memd_struct* memd, fcs_float *field, fcs_int dim, fcs_int e_equil)
{
    static fcs_int init = 1;
    static MPI_Datatype xyPlane,xzPlane,yzPlane; 
    static MPI_Datatype xzPlane2D, xyPlane2D, yzPlane2D;
    /*  fcs_int coord[2]; */
    fcs_int l, s_dir, r_dir;
    /*  fcs_int pos=0; */
    MPI_Status status[2];
    MPI_Request request[]={MPI_REQUEST_NULL, MPI_REQUEST_NULL};
    fcs_int offset, doffset, skip, stride, nblocks;
    /** surface_patch */
    static t_surf_patch  surface_patch[6];
	
    if(init) {
        MPI_Datatype xz_plaq, oneslice;
		
        ifcs_memd__calc_surface_patches(memd, surface_patch);
        ifcs_memd__prepare_surface_planes(dim, &xyPlane, &xzPlane, &yzPlane, surface_patch);
		
        MPI_Type_vector(surface_patch[0].stride, 2, 3, FCS_MPI_FLOAT,&yzPlane2D);    
        MPI_Type_commit(&yzPlane2D);
		
        /* create data type for xz plaquette */
        MPI_Type_create_hvector(2,1*sizeof(fcs_float),2*sizeof(fcs_float), MPI_BYTE, &xz_plaq);
        /* create data type for a 1D section */
        MPI_Type_contiguous(surface_patch[2].stride, xz_plaq, &oneslice); 
        /* create data type for a 2D xz plane */
        MPI_Type_create_hvector(surface_patch[2].nblocks, 1, dim*surface_patch[2].skip*sizeof(fcs_float), oneslice, &xzPlane2D);
        MPI_Type_commit(&xzPlane2D);    
        /* create data type for a 2D xy plane */
        MPI_Type_vector(surface_patch[4].nblocks, 2, dim*surface_patch[4].skip, FCS_MPI_FLOAT, &xyPlane2D);
        MPI_Type_commit(&xyPlane2D); 
		
        init = 0;
    }
	
    
    /** direction loop */
    for(s_dir=0; s_dir < 6; s_dir++) { 
        offset = dim * surface_patch[s_dir].offset;
        doffset= dim * surface_patch[s_dir].doffset;
		
        if(s_dir%2==0) r_dir = s_dir+1;
        else           r_dir = s_dir-1;
        /** pack send halo-plane data */
        if(memd->mpiparams.node_neighbors[s_dir] != memd->mpiparams.this_node) {
            /** communication */
            switch(s_dir) {
                case 0 :
                case 1 :
                    if(e_equil || dim == 1) {
                        MPI_Irecv (&field[doffset],1,yzPlane,memd->mpiparams.node_neighbors[s_dir],REQ_MAGGS_SPREAD,memd->mpiparams.communicator,&request[0]);
                        MPI_Isend(&field[offset],1,yzPlane,memd->mpiparams.node_neighbors[r_dir],REQ_MAGGS_SPREAD,memd->mpiparams.communicator,&request[1]);
                    }
                    else {
                        MPI_Irecv (&field[doffset+1],1,yzPlane2D,memd->mpiparams.node_neighbors[s_dir],REQ_MAGGS_SPREAD,memd->mpiparams.communicator,&request[0]);
                        MPI_Isend(&field[offset+1],1,yzPlane2D,memd->mpiparams.node_neighbors[r_dir],REQ_MAGGS_SPREAD,memd->mpiparams.communicator,&request[1]);
                    }	  
					
                    MPI_Waitall(2,request,status);
                    break;
                case 2 :
                case 3 :
                    if(e_equil || dim == 1) {
                        MPI_Irecv (&field[doffset],1,xzPlane,memd->mpiparams.node_neighbors[s_dir],REQ_MAGGS_SPREAD,memd->mpiparams.communicator,&request[0]);
                        MPI_Isend(&field[offset],1,xzPlane,memd->mpiparams.node_neighbors[r_dir],REQ_MAGGS_SPREAD,memd->mpiparams.communicator,&request[1]);
                    }
                    else {
                        MPI_Irecv (&field[doffset],1,xzPlane2D,memd->mpiparams.node_neighbors[s_dir],REQ_MAGGS_SPREAD,memd->mpiparams.communicator,&request[0]);
                        MPI_Isend(&field[offset],1,xzPlane2D,memd->mpiparams.node_neighbors[r_dir],REQ_MAGGS_SPREAD,memd->mpiparams.communicator,&request[1]);
                    }	  
                    MPI_Waitall(2,request,status);
                    break;
                case 4 :
                case 5 : 
                    if(e_equil || dim == 1) {
                        MPI_Irecv (&field[doffset],1,xyPlane,memd->mpiparams.node_neighbors[s_dir],REQ_MAGGS_SPREAD,memd->mpiparams.communicator,&request[0]);
                        MPI_Isend(&field[offset],1,xyPlane,memd->mpiparams.node_neighbors[r_dir],REQ_MAGGS_SPREAD,memd->mpiparams.communicator,&request[1]);
                    }
                    else {
                        MPI_Irecv (&field[doffset],1,xyPlane2D,memd->mpiparams.node_neighbors[s_dir],REQ_MAGGS_SPREAD,memd->mpiparams.communicator,&request[0]);
                        MPI_Isend(&field[offset],1,xyPlane2D,memd->mpiparams.node_neighbors[r_dir],REQ_MAGGS_SPREAD,memd->mpiparams.communicator,&request[1]);
                    }
                    MPI_Waitall(2,request,status);
                    break;
            }
        }
        
        else {
            /** copy locally */
            skip    = dim * surface_patch[s_dir].skip;
            stride  = dim * surface_patch[s_dir].stride * sizeof(fcs_float);
            nblocks = surface_patch[s_dir].nblocks;
			
            for(l=0; l<nblocks; l++){
                memcpy(&(field[doffset]), &(field[offset]), stride);
                offset  += skip;
                doffset += skip;
            }
			
        }
    }
}




