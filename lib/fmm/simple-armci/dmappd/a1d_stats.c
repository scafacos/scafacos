/********************************************************************
 * The following is a notice of limited availability of the code, and disclaimer
 * which must be included in the prologue of the code and in all source listings
 * of the code.
 *
 * Copyright (c) 2010 Argonne Leadership Computing Facility, Argonne National Laboratory
 *
 * Permission is hereby granted to use, reproduce, prepare derivative works, and
 * to redistribute to others.
 *
 *                 LICENSE
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * - Redistributions of source code must retain the above copyright
 *   notice, this list of conditions and the following disclaimer.
 *
 *  - Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer listed
 *    in this license in the documentation and/or other materials
 *    provided with the distribution.
 *
 *  - Neither the name of the copyright holders nor the names of its
 *    contributors may be used to endorse or promote products derived from
 *    this software without specific prior written permission.
 *
 * The copyright holders provide no reassurances that the source code
 * provided does not infringe any patent, copyright, or any other
 * intellectual property rights of third parties.  The copyright holders
 * disclaim any liability to any recipient for claims brought against
 * recipient by any third party for infringement of that parties
 * intellectual property rights.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *********************************************************************/

#include "a1d_headers.h"
#include "a1d_globals.h"

#ifdef __CRAYXE
void A1D_Parse_error(dmapp_return_t rc)
{
    switch (rc)
    {
        case DMAPP_RC_SUCCESS:
            //fprintf(stderr,"DMAPP_RC_SUCCESS\n");
            break;

        case DMAPP_RC_INVALID_PARAM:
            fprintf(stderr,"DMAPP_RC_INVALID_PARAM\n");
            break;

        case DMAPP_RC_ALIGNMENT_ERROR:
            fprintf(stderr,"DMAPP_RC_ALIGNMENT_ERROR\n");
            break;

        case DMAPP_RC_TRANSACTION_ERROR:
            fprintf(stderr,"DMAPP_RC_TRANSACTION_ERROR\n");
            break;

        case DMAPP_RC_RESOURCE_ERROR:
            fprintf(stderr,"DMAPP_RC_RESOURCE_ERROR\n");
            break;

        case DMAPP_RC_PERMISSION_ERROR:
            fprintf(stderr,"DMAPP_RC_PERMISSION_ERROR\n");
            break;

        case DMAPP_RC_NO_SPACE:
            fprintf(stderr,"DMAPP_RC_NO_SPACE\n");
            break;

        case DMAPP_RC_NOT_DONE:
            fprintf(stderr,"DMAPP_RC_NOT_DONE\n");
            break;

        case DMAPP_RC_NOT_SUPPORTED:
            fprintf(stderr,"DMAPP_RC_NOT_SUPPORTED\n");
            break;

        case DMAPP_RC_NOT_FOUND:
            fprintf(stderr,"DMAPP_RC_NOT_FOUND\n");
            break;

        case DMAPP_RC_BUSY:
            fprintf(stderr,"DMAPP_RC_BUSY\n");
            break;

        case DMAPP_RC_NOT_USED:
            fprintf(stderr,"DMAPP_RC_NOT_USED\n");
            break;

        default:
            fprintf(stderr,"Unknown DMAPP error code. \n");
            break;
    }
    fflush(stderr);
    return;
}
#endif

void A1D_Print_stats()
{
    return;
}







