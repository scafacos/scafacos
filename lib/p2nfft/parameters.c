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

#include "parameters.h"
#include "types.h"
#include <stdlib.h>
#include <string.h>

static fcs_int set_flag(
    unsigned flag, unsigned *flags
    )
{
  if(flag & *flags)
    return 0;

  *flags |= flag;
  return 1;
}

static fcs_int unset_flag(
    unsigned flag, unsigned *flags
    )
{
  if(flag & ~*flags)
    return 0;

  *flags &= ~flag;
  return 1;
}

static fcs_int correct_flag(
    unsigned new_flag, fcs_int set_it,
    unsigned *flags
    )
{
  if(set_it)
    return set_flag(new_flag, flags);
  else
    return unset_flag(new_flag, flags);
}

#define IFCS_P2NFFT_SET_GET_FLAG(METHOD, NAME, FLAG)                            \
FCSResult ifcs_p2nfft_set_ ## METHOD ## NAME(                                   \
    void *rd, const char* fnc_name, fcs_int set_ ## NAME                        \
    ) 	                                                                        \
{                                                                               \
  ifcs_p2nfft_data_struct *d = (ifcs_p2nfft_data_struct*)rd;                    \
  if( rd==NULL )                                                                \
    return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, "Got NULL Pointer."); \
  d->needs_retune |= correct_flag(FLAG, set_ ## NAME,                           \
      &d->METHOD ## flags);                                                     \
  return NULL;                                                                  \
}                                                                               \
FCSResult ifcs_p2nfft_get_ ## METHOD ## NAME(                                   \
    void *rd, const char* fnc_name, fcs_int *get_ ## NAME                       \
    )                                                                           \
{                                                                               \
  ifcs_p2nfft_data_struct *d = (ifcs_p2nfft_data_struct*)rd;                    \
  if( rd==NULL )                                                                \
    return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, "Got NULL Pointer."); \
  *get_ ## NAME =  (d->METHOD ## flags & FLAG) ? 1 : 0;                         \
  return NULL;                                                                  \
}



/*********************************************
 *  Setters and getters for P2NFFT parameters 
 ********************************************/

/* setters/getters for tolerance */
FCSResult ifcs_p2nfft_set_tolerance(
    void *rd, const char* fnc_name, fcs_int tolerance_type, fcs_float tolerance
    )
{
  ifcs_p2nfft_data_struct *d = (ifcs_p2nfft_data_struct*)rd;
  if( rd==NULL )
    return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, "Got NULL Pointer.");

  if(tolerance_type != d->tolerance_type)
    d->needs_retune = 1;
  if (!fcs_float_is_equal(tolerance, d->tolerance))
    d->needs_retune = 1;
  d->tolerance_type = tolerance_type;
  d->tolerance      = tolerance;
  return NULL;
}

FCSResult ifcs_p2nfft_set_tolerance_tune(
    void *rd, const char* fnc_name
    )
{
  ifcs_p2nfft_data_struct *d = (ifcs_p2nfft_data_struct*)rd;
  if( rd==NULL )
    return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, "Got NULL Pointer.");

  d->needs_retune = 1;
  d->tolerance_type = FCS_TOLERANCE_TYPE_UNDEFINED;
  d->tolerance = -1.0;
  return NULL;
}  


FCSResult ifcs_p2nfft_get_tolerance(
    void *rd, const char* fnc_name, fcs_int *tolerance_type, fcs_float *tolerance
    )
{
  ifcs_p2nfft_data_struct *d = (ifcs_p2nfft_data_struct*)rd;
  if( rd==NULL )
    return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, "Got NULL Pointer.");

  *tolerance_type = d->tolerance_type;
  *tolerance      = d->tolerance;
  return NULL;
}

/* Getters and Setters for near field cutoff radius */
FCSResult ifcs_p2nfft_set_r_cut(
    void *rd, const char* fnc_name, fcs_float r_cut
    )
{
  ifcs_p2nfft_data_struct *d = (ifcs_p2nfft_data_struct*)rd;
  if( rd==NULL )
    return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, "Got NULL Pointer.");

  if (!fcs_float_is_equal(r_cut, d->r_cut))
    d->needs_retune = 1;
  d->r_cut = r_cut;
  d->one_over_r_cut = 1.0/r_cut;
  d->tune_r_cut = 0;
  return NULL;
}

FCSResult ifcs_p2nfft_set_r_cut_tune(
    void *rd, const char* fnc_name
    )
{
  ifcs_p2nfft_data_struct *d = (ifcs_p2nfft_data_struct*)rd;
  if( rd==NULL )
    return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, "Got NULL Pointer.");

  d->needs_retune = 1;
  d->tune_r_cut = 1;
  d->r_cut = -1.0;
  return NULL;
}

FCSResult ifcs_p2nfft_get_r_cut(
    void *rd, const char* fnc_name, fcs_float *r_cut
    )
{
  ifcs_p2nfft_data_struct *d = (ifcs_p2nfft_data_struct*)rd;
  if( rd==NULL )
    return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, "Got NULL Pointer.");

  *r_cut = d->r_cut;
  return NULL;
}

/* Getters and Setters for scaled near field cutoff radius */
FCSResult ifcs_p2nfft_set_epsI(
    void *rd, const char* fnc_name, fcs_float epsI
    )
{
  ifcs_p2nfft_data_struct *d = (ifcs_p2nfft_data_struct*)rd;
  if( rd==NULL )
    return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, "Got NULL Pointer.");

  if( epsI >= 0.5 )
    return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, "Near field cutoff 'epsI' does not hold epsI < 0.5");

  if (!fcs_float_is_equal(epsI, d->epsI))
    d->needs_retune = 1;
  d->epsI = epsI;
  d->tune_epsI = 0;
  return NULL;
}

FCSResult ifcs_p2nfft_set_epsI_tune(
    void *rd, const char* fnc_name
    )
{
  ifcs_p2nfft_data_struct *d = (ifcs_p2nfft_data_struct*)rd;
  if( rd==NULL )
    return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, "Got NULL Pointer.");

  d->needs_retune = 1;
  d->tune_epsI = 1;
  d->epsI = -1.0;
  return NULL;
}

FCSResult ifcs_p2nfft_get_epsI(
    void *rd, const char* fnc_name, fcs_float *epsI
    )
{
  ifcs_p2nfft_data_struct *d = (ifcs_p2nfft_data_struct*)rd;
  if( rd==NULL )
    return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, "Got NULL Pointer.");

  *epsI = d->epsI;
  return NULL;
}


/* Getters and Setters for scaled far field regularization border */
FCSResult ifcs_p2nfft_set_epsB(
    void *rd, const char* fnc_name, fcs_float epsB
    )
{
  ifcs_p2nfft_data_struct *d = (ifcs_p2nfft_data_struct*)rd;
  if( rd==NULL )
    return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, "Got NULL Pointer.");

  if( epsB >= 0.5 )
    return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, "Far field border 'epsB' does not hold epsB < 0.5");

  if (!fcs_float_is_equal(epsB, d->epsB))
    d->needs_retune = 1;
  d->epsB = epsB;
  d->tune_epsB = 0;
  return NULL;
}

FCSResult ifcs_p2nfft_set_epsB_tune(
    void *rd, const char* fnc_name
    )
{
  ifcs_p2nfft_data_struct *d = (ifcs_p2nfft_data_struct*)rd;
  if( rd==NULL )
    return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, "Got NULL Pointer.");

  d->needs_retune = 1;
  d->tune_epsB = 1;
  d->epsB = -1.0;
  return NULL;
}

FCSResult ifcs_p2nfft_get_epsB(
    void *rd, const char* fnc_name, fcs_float *epsB
    )
{
  ifcs_p2nfft_data_struct *d = (ifcs_p2nfft_data_struct*)rd;
  if( rd==NULL )
    return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, "Got NULL Pointer.");

  *epsB = d->epsB;
  return NULL;
}


/* Getter and Setter for far field continuation value c used by taylor2p */
FCSResult ifcs_p2nfft_set_c(
    void *rd, const char* fnc_name, fcs_float c
    )
{
  ifcs_p2nfft_data_struct *d = (ifcs_p2nfft_data_struct*)rd;
  if( rd==NULL )
    return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, "Got NULL Pointer.");

  if (!fcs_float_is_equal(c, d->c))
    d->needs_retune = 1;
  d->c = c;
  d->tune_c = 0;
  return NULL;
}

FCSResult ifcs_p2nfft_set_c_tune(
    void *rd, const char* fnc_name
    )
{
  ifcs_p2nfft_data_struct *d = (ifcs_p2nfft_data_struct*)rd;
  if( rd==NULL )
    return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, "Got NULL Pointer.");

  d->needs_retune = 1;
  d->tune_c = 1;
  d->c = 0.0;
  return NULL;
}

FCSResult ifcs_p2nfft_get_c(
    void *rd, const char* fnc_name, fcs_float *c
    )
{
  ifcs_p2nfft_data_struct *d = (ifcs_p2nfft_data_struct*)rd;
  if( rd==NULL )
    return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, "Got NULL Pointer.");

  *c = d->c;
  return NULL;
}


/* Getters and Setters for Ewald splitting parameter */
FCSResult ifcs_p2nfft_set_alpha(
    void *rd, const char* fnc_name, fcs_float alpha
    )
{
  ifcs_p2nfft_data_struct *d = (ifcs_p2nfft_data_struct*)rd;
  if( rd==NULL )
    return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, "Got NULL Pointer.");

  if (!fcs_float_is_equal(alpha, d->alpha))
    d->needs_retune = 1;
  d->alpha = alpha;
  d->tune_alpha = 0;
  return NULL;
}

FCSResult ifcs_p2nfft_set_alpha_tune(
    void *rd, const char* fnc_name
    )
{
  ifcs_p2nfft_data_struct *d = (ifcs_p2nfft_data_struct*)rd;
  if( rd==NULL )
    return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, "Got NULL Pointer.");

  d->needs_retune = 1;
  d->tune_alpha = 1;
  return NULL;
}

FCSResult ifcs_p2nfft_get_alpha(
    void *rd, const char* fnc_name, fcs_float *alpha
    )
{
  ifcs_p2nfft_data_struct *d = (ifcs_p2nfft_data_struct*)rd;
  if( rd==NULL )
    return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, "Got NULL Pointer.");

  *alpha = d->alpha;
  return NULL;
}

/* Getters and Setters for charge assignment order */
FCSResult ifcs_p2nfft_set_cao(
    void *rd, const char* fnc_name, fcs_int cao
    )
{
  return ifcs_p2nfft_set_pnfft_m(rd, fnc_name, (cao+1)/2);
}

FCSResult ifcs_p2nfft_set_cao_tune(
    void *rd, const char* fnc_name
    )
{
  return ifcs_p2nfft_set_pnfft_m_tune(rd, fnc_name);
}

FCSResult ifcs_p2nfft_get_cao(
    void *rd, const char* fnc_name, fcs_int *cao
    )
{
  fcs_int m=0;
  FCSResult res;

  res = ifcs_p2nfft_get_pnfft_m(rd, fnc_name, &m);
  if(res != NULL)
    return res;

  *cao = 2*m;
  return res;
}

/* Getters and Setters for order of P2NFFT near and far field interpolation */
FCSResult ifcs_p2nfft_set_interpolation_order(
    void *rd, const char* fnc_name, fcs_int intpol_order
    )
{
  ifcs_p2nfft_data_struct *d = (ifcs_p2nfft_data_struct*)rd;
  if( rd==NULL )
    return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, "Got NULL Pointer.");

  if( intpol_order > 3 )
    return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, "Interpolation order not yet supported.");

  /* any negative value turns off interpolation */
  if( intpol_order < 0 ){
    if( d->interpolation_order >=0 )
      d->needs_retune = 1;
  } else {
    if(intpol_order != d->interpolation_order)
      d->needs_retune = 1;
  }
  d->interpolation_order = intpol_order;

  return NULL;
}

FCSResult ifcs_p2nfft_get_interpolation_order(
    void *rd, const char* fnc_name, fcs_int *intpol_order
    )
{
  ifcs_p2nfft_data_struct *d = (ifcs_p2nfft_data_struct*)rd;
  if( rd==NULL )
    return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, "Got NULL Pointer.");

  *intpol_order = d->interpolation_order;
  return NULL;
}

/* Getters and Setters for P2NFFT near field regularization flag */
FCSResult ifcs_p2nfft_set_reg_near(
    void *rd, const char* fnc_name, fcs_int reg
    )
{
  ifcs_p2nfft_data_struct *d = (ifcs_p2nfft_data_struct*)rd;
  if( rd==NULL )
    return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, "Got NULL Pointer.");

  if( (reg != FCS_P2NFFT_REG_NEAR_DEFAULT) && (reg != FCS_P2NFFT_REG_NEAR_CG) && (reg != FCS_P2NFFT_REG_NEAR_T2P) )
    return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, "Unknown regularization.");
  
  if(reg != d->reg_near)
    d->needs_retune = 1;
  d->reg_near = reg;

  return NULL;
}

FCSResult ifcs_p2nfft_set_reg_near_by_name(
    void *rd, const char* fnc_name, char *reg_name
    )
{
  unsigned reg_flag;

  if (strcmp(reg_name,"default") == 0)
    reg_flag = FCS_P2NFFT_REG_NEAR_DEFAULT;
  else if (strcmp(reg_name,"cg") == 0)
    reg_flag = FCS_P2NFFT_REG_NEAR_CG;
  else if (strcmp(reg_name,"t2p") == 0)
    reg_flag = FCS_P2NFFT_REG_NEAR_T2P;
  else /* unknown regularization */
    return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, "Unknown regularization.");

  return ifcs_p2nfft_set_reg_near(rd, fnc_name, reg_flag);
}

FCSResult ifcs_p2nfft_get_reg_near(
    void *rd, const char* fnc_name, fcs_int *reg
    )
{
  ifcs_p2nfft_data_struct *d = (ifcs_p2nfft_data_struct*)rd;
  if( rd==NULL )
    return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, "Got NULL Pointer.");

  *reg = d->reg_near;
  return NULL;
}

/* Getters and Setters for P2NFFT far field regularization flag */
FCSResult ifcs_p2nfft_set_reg_far(
    void *rd, const char* fnc_name, fcs_int reg
    )
{
  ifcs_p2nfft_data_struct *d = (ifcs_p2nfft_data_struct*)rd;
  if( rd==NULL )
    return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, "Got NULL Pointer.");

  if( (reg != FCS_P2NFFT_REG_FAR_DEFAULT)
      && (reg != FCS_P2NFFT_REG_FAR_RAD_CG)
      && (reg != FCS_P2NFFT_REG_FAR_RAD_T2P_SYM)
      && (reg != FCS_P2NFFT_REG_FAR_RAD_T2P_EC)
      && (reg != FCS_P2NFFT_REG_FAR_RAD_T2P_IC)
      && (reg != FCS_P2NFFT_REG_FAR_REC_T2P_SYM)
      && (reg != FCS_P2NFFT_REG_FAR_REC_T2P_EC)
      && (reg != FCS_P2NFFT_REG_FAR_REC_T2P_IC)   
    )
    return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, "Unknown far field regularization.");
  
  if(reg != d->reg_far)
    d->needs_retune = 1;
  d->reg_far = reg;

  return NULL;
}

FCSResult ifcs_p2nfft_set_reg_far_by_name(
    void *rd, const char* fnc_name, char *reg_name
    )
{
  unsigned reg_flag;

  if (strcmp(reg_name,"default") == 0)
    reg_flag = FCS_P2NFFT_REG_FAR_DEFAULT;
  else if (strcmp(reg_name,"rad_cg") == 0)
    reg_flag = FCS_P2NFFT_REG_FAR_RAD_CG;
  else if (strcmp(reg_name,"rad_t2p_sym") == 0)
    reg_flag = FCS_P2NFFT_REG_FAR_RAD_T2P_SYM;
  else if (strcmp(reg_name,"rad_t2p_ec") == 0)
    reg_flag = FCS_P2NFFT_REG_FAR_RAD_T2P_EC;
  else if (strcmp(reg_name,"rad_t2p_ic") == 0)
    reg_flag = FCS_P2NFFT_REG_FAR_RAD_T2P_IC;
  else if (strcmp(reg_name,"rec_t2p_sym") == 0)
    reg_flag = FCS_P2NFFT_REG_FAR_REC_T2P_SYM;
  else if (strcmp(reg_name,"rec_t2p_ec") == 0)
    reg_flag = FCS_P2NFFT_REG_FAR_REC_T2P_EC;
  else if (strcmp(reg_name,"rec_t2p_ic") == 0)
    reg_flag = FCS_P2NFFT_REG_FAR_REC_T2P_IC;
  else /* unknown regularization */
    return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, "Unknown regularization.");

  return ifcs_p2nfft_set_reg_far(rd, fnc_name, reg_flag);
}

FCSResult ifcs_p2nfft_get_reg_far(
    void *rd, const char* fnc_name, fcs_int *reg
    )
{
  ifcs_p2nfft_data_struct *d = (ifcs_p2nfft_data_struct*)rd;
  if( rd==NULL )
    return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, "Got NULL Pointer.");

  *reg = d->reg_far;
  return NULL;
}

/* Getters and Setters for polynomial degree of near field regularization */
FCSResult ifcs_p2nfft_set_p(
    void *rd, const char* fnc_name, fcs_int p
    )
{
  ifcs_p2nfft_data_struct *d = (ifcs_p2nfft_data_struct*)rd;
  if( rd==NULL )
    return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, "Got NULL Pointer.");

  if (!fcs_float_is_equal(p, d->p))
    d->needs_retune = 1;
  d->p = p;
  d->tune_p = 0;
  return NULL;
}

FCSResult ifcs_p2nfft_set_p_tune(
    void *rd, const char* fnc_name
    )
{
  ifcs_p2nfft_data_struct *d = (ifcs_p2nfft_data_struct*)rd;
  if( rd==NULL )
    return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, "Got NULL Pointer.");

  d->needs_retune = 1;
  d->tune_p = 1;
  return NULL;
}

FCSResult ifcs_p2nfft_get_p(
    void *rd, const char* fnc_name, fcs_int *p
    )
{
  ifcs_p2nfft_data_struct *d = (ifcs_p2nfft_data_struct*)rd;
  if( rd==NULL )
    return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, "Got NULL Pointer.");

  *p = d->p;
  return NULL;
}

/* setters/getters for virial */
FCSResult ifcs_p2nfft_require_virial(
    void *rd, const char* fnc_name, fcs_int flag
    )
{
  ifcs_p2nfft_data_struct *d = (ifcs_p2nfft_data_struct*)rd;
  if( rd==NULL )
    return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, "Got NULL Pointer.");

  if(flag){
    if(d->virial == NULL)
      d->virial = (fcs_float*) malloc(sizeof(fcs_float) * 9);
  } else {
    if(d->virial != NULL){
      free(d->virial);
      d->virial = NULL;
    }
  }
  return NULL;
}

FCSResult ifcs_p2nfft_get_virial(
    void *rd, const char* fnc_name, fcs_float *virial
    )
{
  ifcs_p2nfft_data_struct *d = (ifcs_p2nfft_data_struct*)rd;
  if( rd==NULL )
    return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, "Got NULL Pointer.");

  if(!virial)
    return fcsResult_create(FCS_NULL_ARGUMENT,fnc_name,"Supplied virial pointer must not be a null pointer."); 
  if(!d->virial)
    return fcsResult_create(FCS_NULL_ARGUMENT,fnc_name,"Virial computation is not activated. Use fcs_require_virial to do so."); 

  if(d->virial != NULL){
    for (fcs_int t=0; t < 9; t++)
      virial[t] = d->virial[t];
  }
  return NULL;
}

FCSResult ifcs_p2nfft_virial_is_active(
    void *rd, const char* fnc_name, fcs_int *is_active
    )
{
  ifcs_p2nfft_data_struct *d = (ifcs_p2nfft_data_struct*)rd;
  if( rd==NULL )
    return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, "Got NULL Pointer.");

  *is_active = d->virial != NULL;
  return NULL;
}


/* Getters and Setters for ignore tolerance flag */
IFCS_P2NFFT_SET_GET_FLAG(, ignore_tolerance,      FCS_P2NFFT_IGNORE_TOLERANCE)


/********************************************
 *  Setters and getter for PNFFT parameters 
 *******************************************/

/* Getters and Setters for FFT grid size */
FCSResult ifcs_p2nfft_set_pnfft_N(
    void *rd, const char* fnc_name, fcs_int N0, fcs_int N1, fcs_int N2
    )
{
  ifcs_p2nfft_data_struct *d = (ifcs_p2nfft_data_struct*)rd;
  if( rd==NULL )
    return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, "Got NULL Pointer.");

  d->N[0] = N0; d->N[1] = N1; d->N[2] = N2;
  d->needs_retune = 1; 
  d->tune_N = 0;
  return NULL;
}

FCSResult ifcs_p2nfft_set_pnfft_N_tune(
    void *rd, const char* fnc_name
    )
{
  ifcs_p2nfft_data_struct *d = (ifcs_p2nfft_data_struct*)rd;
  if( rd==NULL )
    return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, "Got NULL Pointer.");

  d->needs_retune = 1; 
  d->tune_N = 1;
  return NULL;
}

FCSResult ifcs_p2nfft_get_pnfft_N(
    void *rd, const char* fnc_name,
    fcs_int *N0, fcs_int *N1, fcs_int *N2
    )
{
  ifcs_p2nfft_data_struct *d = (ifcs_p2nfft_data_struct*)rd;
  if( rd==NULL )
    return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, "Got NULL Pointer.");

  *N0 = d->N[0]; *N1 = d->N[1]; *N2 = d->N[2]; 
  return NULL;
}

/* Getters and Setters for oversampled FFT grid size */
FCSResult ifcs_p2nfft_set_pnfft_n(
    void *rd, const char* fnc_name, fcs_int n0, fcs_int n1, fcs_int n2
    )
{
  ifcs_p2nfft_data_struct *d = (ifcs_p2nfft_data_struct*)rd;
  if( rd==NULL )
    return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, "Got NULL Pointer.");

  d->n[0] = n0; d->n[1] = n1; d->n[2] = n2;
  d->needs_retune = 1; 
  d->tune_n = 0;
  return NULL;
}

FCSResult ifcs_p2nfft_set_pnfft_n_tune(
    void *rd, const char* fnc_name
    )
{
  ifcs_p2nfft_data_struct *d = (ifcs_p2nfft_data_struct*)rd;
  if( rd==NULL )
    return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, "Got NULL Pointer.");

  d->needs_retune = 1; 
  d->tune_n = 1;
  return NULL;
}

FCSResult ifcs_p2nfft_get_pnfft_n(
    void *rd, const char* fnc_name,
    fcs_int *n0, fcs_int *n1, fcs_int *n2
    )
{
  ifcs_p2nfft_data_struct *d = (ifcs_p2nfft_data_struct*)rd;
  if( rd==NULL )
    return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, "Got NULL Pointer.");

  *n0 = d->n[0]; *n1 = d->n[1]; *n2 = d->n[2]; 
  return NULL;
}

/* Getters and Setters for PNFFT window function */
FCSResult ifcs_p2nfft_set_pnfft_window(
    void *rd, const char* fnc_name, fcs_int window
    )
{
  ifcs_p2nfft_data_struct *d = (ifcs_p2nfft_data_struct*)rd;
  if( rd==NULL )
    return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, "Got NULL Pointer.");

  if((window < 0) || (4 < window) )
    return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, "Unknown window function.");

  if ( window != d->pnfft_window )
    d->needs_retune = 1;

  d->pnfft_window = window;
  return NULL;
}

FCSResult ifcs_p2nfft_set_pnfft_window_by_name(
    void *rd, const char* fnc_name, char *window_name
    )
{
  fcs_int window;

  if (strcmp(window_name,"gaussian") == 0)
    window = 0;
  else if (strcmp(window_name,"bspline") == 0)
    window = 1;
  else if (strcmp(window_name,"sinc") == 0)
    window = 2;
  else if (strcmp(window_name,"kaiser") == 0)
    window = 3;
  else if (strcmp(window_name,"bessel_i0") == 0)
    window = 4;
  else /* unknown window */
    return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, "Unknown window function.");

  return ifcs_p2nfft_set_pnfft_window(rd, fnc_name, window);
}

FCSResult ifcs_p2nfft_get_pnfft_window(
    void *rd, const char* fnc_name, fcs_int *window
    )
{
  ifcs_p2nfft_data_struct *d = (ifcs_p2nfft_data_struct*)rd;
  if( rd==NULL )
    return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, "Got NULL Pointer.");

  *window = d->pnfft_window;
  return NULL;
}

/* Getters and Setters for charge assignment order */
FCSResult ifcs_p2nfft_set_pnfft_m(
    void *rd, const char* fnc_name, fcs_int m
    )
{
  ifcs_p2nfft_data_struct *d = (ifcs_p2nfft_data_struct*)rd;
  if( rd==NULL )
    return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, "Got NULL Pointer.");

  if (!fcs_float_is_equal(m, d->m))
    d->needs_retune = 1;
  d->m = m;
  d->tune_m = 0;
  return NULL;
}

FCSResult ifcs_p2nfft_set_pnfft_m_tune(
    void *rd, const char* fnc_name
    )
{
  ifcs_p2nfft_data_struct *d = (ifcs_p2nfft_data_struct*)rd;
  if( rd==NULL )
    return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, "Got NULL Pointer.");

  d->needs_retune = 1;
  d->tune_m = 1;
  return NULL;
}

FCSResult ifcs_p2nfft_get_pnfft_m(
    void *rd, const char* fnc_name, fcs_int *m
    )
{
  ifcs_p2nfft_data_struct *d = (ifcs_p2nfft_data_struct*)rd;
  if( rd==NULL )
    return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, "Got NULL Pointer.");

  *m = d->m;
  return NULL;
}

/* Getters and Setters for PNFFT window interpolation order */
FCSResult ifcs_p2nfft_set_pnfft_interpolation_order(
    void *rd, const char* fnc_name, fcs_int intpol_order
    )
{
  ifcs_p2nfft_data_struct *d = (ifcs_p2nfft_data_struct*)rd;
  if( rd==NULL )
    return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, "Got NULL Pointer.");

  if( intpol_order > 3 )
    return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, "Interpolation order not yet supported.");

  /* any negative value turns off interpolation */
  if( intpol_order < 0 ){
    if( d->pnfft_interpolation_order >=0 )
      d->needs_retune = 1;
  } else {
    if(intpol_order != d->pnfft_interpolation_order)
      d->needs_retune = 1;
  }
  d->pnfft_interpolation_order = intpol_order;

  return NULL;
}

FCSResult ifcs_p2nfft_get_pnfft_interpolation_order(
    void *rd, const char* fnc_name, fcs_int *intpol_order
    )
{
  ifcs_p2nfft_data_struct *d = (ifcs_p2nfft_data_struct*)rd;
  if( rd==NULL )
    return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, "Got NULL Pointer.");

  *intpol_order = d->pnfft_interpolation_order;
  return NULL;
}

/* Getters and Setters for PNFFT flags */
IFCS_P2NFFT_SET_GET_FLAG(pnfft_, pre_phi_hat,      PNFFT_PRE_PHI_HAT)
IFCS_P2NFFT_SET_GET_FLAG(pnfft_, fg_psi,           PNFFT_FG_PSI)
IFCS_P2NFFT_SET_GET_FLAG(pnfft_, fft_in_place,     PNFFT_FFT_IN_PLACE)
IFCS_P2NFFT_SET_GET_FLAG(pnfft_, sort_nodes,       PNFFT_SORT_NODES)
IFCS_P2NFFT_SET_GET_FLAG(pnfft_, interlaced,       PNFFT_INTERLACED)
IFCS_P2NFFT_SET_GET_FLAG(pnfft_, grad_ik,          PNFFT_GRAD_IK)
IFCS_P2NFFT_SET_GET_FLAG(pnfft_, pre_psi,          PNFFT_PRE_PSI)
IFCS_P2NFFT_SET_GET_FLAG(pnfft_, pre_fg_psi,       PNFFT_PRE_FG_PSI)
IFCS_P2NFFT_SET_GET_FLAG(pnfft_, pre_full_psi,     PNFFT_PRE_FULL_PSI)
IFCS_P2NFFT_SET_GET_FLAG(pnfft_, pre_full_fg_psi,  PNFFT_PRE_FULL_FG_PSI)
IFCS_P2NFFT_SET_GET_FLAG(pnfft_, real_f,           PNFFT_REAL_F)


/********************************************
 *  Setters and getter for PFFT parameters 
 *******************************************/

/* Getters and Setters for PFFT patience */
FCSResult ifcs_p2nfft_set_pfft_patience(
    void *rd, const char* fnc_name, fcs_int patience
    )
{
  ifcs_p2nfft_data_struct *d = (ifcs_p2nfft_data_struct*)rd;
  if( rd==NULL )
    return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, "Got NULL Pointer.");

  if( (patience < 0) || (3 < patience) )
    return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, "Patience not supported.");

  if(d->pfft_patience != patience)
    d->needs_retune = 1;

  d->pfft_patience = patience;
  return NULL;
}

FCSResult ifcs_p2nfft_get_pfft_patience(
    void *rd, const char* fnc_name, fcs_int *patience
    )
{
  ifcs_p2nfft_data_struct *d = (ifcs_p2nfft_data_struct*)rd;
  if( rd==NULL )
    return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, "Got NULL Pointer.");

  *patience = d->pfft_patience;
  return NULL;
}

FCSResult ifcs_p2nfft_set_pfft_patience_by_name(
    void *rd, const char* fnc_name, char *patience_name
    )
{
  fcs_int patience;

  if (strcmp(patience_name,"estimate") == 0)
    patience = 0;
  else if (strcmp(patience_name,"measure") == 0)
    patience = 1;
  else if (strcmp(patience_name,"patient") == 0)
    patience = 2;
  else if (strcmp(patience_name,"exhaustive") == 0)
    patience = 3;
  else /* unknown patience */
    return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, "Unknown patience function.");

  return ifcs_p2nfft_set_pfft_patience(rd, fnc_name, patience);
}

/* Getters and Setters for PFFT flags */
IFCS_P2NFFT_SET_GET_FLAG(pfft_,  preserve_input,   PFFT_PRESERVE_INPUT)
IFCS_P2NFFT_SET_GET_FLAG(pfft_,  tune,             PFFT_TUNE)


