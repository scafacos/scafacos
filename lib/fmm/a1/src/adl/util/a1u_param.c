/* -*- Mode: C; c-basic-offset:4 ; -*- */
/*
 *  Copyright (C) 2010 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#include "a1.h"
#include "a1d.h"
#include "a1u.h"

/* FIXME: move this to a more appropriate place once all
 *         this A1U/A1D settings crap is sorted out */

#define A1C_NETWORK_BYPASS 1
/* units of bytes.  Very hardware specific but need to reorganize
 * code before moving into device layer. */
#define A1C_NETWORK_BYPASS_UPPER_LIMIT_1D 32768
#define A1C_NETWORK_BYPASS_UPPER_LIMIT_ND 32768

#define A1C_ARMCI_STRICT_ORDERING 0

A1U_Settings_t a1u_settings;

int A1U_Read_parameters(void)
{
    int result = A1_SUCCESS;
    char* value = NULL;

    A1U_FUNC_ENTER();

    a1u_settings.network_bypass = A1C_NETWORK_BYPASS;
    if ((value = getenv("A1_NETWORK_BYPASS")) != NULL)
    {
        a1u_settings.network_bypass = atoi(value);
    }

    /* The threshold BELOW which we do NIC-bypass.  We do this because
     * some architectures (BG/P) have a DMA that beats CPU-based
     * intranode transfers for large buffers.
     */
    a1u_settings.network_bypass_upper_limit_1d = A1C_NETWORK_BYPASS_UPPER_LIMIT_1D;
    if ((value = getenv("A1_NETWORK_BYPASS_UPPER_LIMIT_1D")) != NULL)
    {
        a1u_settings.network_bypass_upper_limit_1d = atoi(value);
    }
    /* For strided, the threshold is much higher. */
    a1u_settings.network_bypass_upper_limit_Nd = A1C_NETWORK_BYPASS_UPPER_LIMIT_ND;
    if ((value = getenv("A1_NETWORK_BYPASS_UPPER_LIMIT_ND")) != NULL)
    {
        a1u_settings.network_bypass_upper_limit_Nd = atoi(value);
    }
    /* If bypass is off, just set upper limit to zero so we always
     * use the NIC.  We do not query network_bypass in contiguous ops. */
    if (a1u_settings.network_bypass == 0)
    {
        a1u_settings.network_bypass_upper_limit_1d = 0;
        a1u_settings.network_bypass_upper_limit_Nd = 0;
    }

    a1u_settings.armci_strict_ordering = A1C_ARMCI_STRICT_ORDERING;
    if ((value = getenv("A1_ARMCI_STRICT_ORDERING")) != NULL)
    {
        a1u_settings.armci_strict_ordering = atoi(value);
    }
    fn_exit: A1U_FUNC_EXIT();
    return result;

    fn_fail: goto fn_exit;
}

int A1U_Print_parameters(void)
{
    int result = A1_SUCCESS;

    A1U_FUNC_ENTER();

    if ( 0 == A1D_Process_id(A1_GROUP_WORLD) )
    {
        A1U_output_printf("=============== A1U Parameters ================\n");
        A1U_output_printf("These are device-independent settings.\n");

        if ( 1==a1u_settings.armci_strict_ordering )
        {
            A1U_output_printf("ARMCI strict ordering        = %s\n","ON");
        }
        else if ( 0==a1u_settings.armci_strict_ordering)
        {
            A1U_output_printf("ARMCI strict ordering        = %s\n","OFF");
        }
        else
        {
            A1U_output_printf("ARMCI strict ordering        = %s\n","WTF");
        }

        if ( 1==a1u_settings.network_bypass )
        {
            A1U_output_printf("NIC bypass                   = %s\n","ON");
            A1U_output_printf("NIC bypass upper limit (1D)  = %u\n",a1u_settings.network_bypass_upper_limit_1d);
            A1U_output_printf("NIC bypass upper limit (ND)  = %u\n",a1u_settings.network_bypass_upper_limit_Nd);

        }
        else if ( 0==a1u_settings.network_bypass)
        {
            A1U_output_printf("Network bypass               = %s\n","OFF");
        }
        else
        {
            A1U_output_printf("Network bypass               = %s\n","WTF");
        }
        A1U_output_printf("===============================================\n\n\n");
        fflush(stdout);
    }

  fn_exit:
    A1U_FUNC_EXIT();
    return result;

  fn_fail:
    goto fn_exit;
}
