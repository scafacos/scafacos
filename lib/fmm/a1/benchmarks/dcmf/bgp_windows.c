/*

 Jeff's test of BGP ProcessWindows

 HOW TO COMPILE (cannot use mpicc):
 mpixlc_r windows.c -I/bgsys/drivers/ppcfloor/arch/include

 HOW TO RUN:
 qsub -n 1 -t 10 --env BG_PROCESSWINDOWS=1 ./a.out
 qsub -n 1 -t 10 --env BG_PROCESSWINDOWS=2 ./a.out
 qsub -n 1 -t 10 --env BG_PROCESSWINDOWS=3 ./a.out
 qsub -n 1 -t 10 --env BG_PROCESSWINDOWS=4 ./a.out

 See: http://dcmf.anl-external.org/docs/ccmi:dcmf-bgp/kernel__interface_8h.html

 __INLINE__ int Kernel_SetProcessWindow(int tlbslot, uint64_t window_paddr, size_t window_reqsize, uint32_t window_permissions, uint32_t *window_actualvaddr, uint64_t *window_actualpaddr, size_t *window_actualsize)
 Sets a virtual memory window for the process based on a user supplied physical address and tlb slot.

 __INLINE__ int Kernel_GetProcessWindow(int tlbslot, uint32_t *window_actualvaddr, uint64_t *window_actualpaddr, size_t *window_actualsize)
 Returns size of the process memory window that was set by the _SetProcessWindow.

 __INLINE__ int Kernel_GetProcessWindowSlotRange(int *minslot, int *maxslot)
 Returns the range of available TLB slots for use by Kernel_SetProcessWindow.

 */

#include <stdio.h>
#include <stdlib.h>
#include "../../../arch/include/spi/kernel_interface.h"

int main(int argc, void* argv[])
{
    int minslot, maxslot;
    Kernel_GetProcessWindowSlotRange(&minslot, &maxslot);
    printf("minslot = %d\n", minslot);
    printf("maxslot = %d\n", maxslot);

    return (0);
}
;
