#include <armci.h>

int parmci_calls = 0;

int ARMCI_Init(void) {
  parmci_calls++;
  return PARMCI_Init();
}
   
int ARMCI_Finalize(void) {
  parmci_calls++;
  return PARMCI_Finalize();
}

void ARMCI_Barrier(void) {
  parmci_calls++;
  /* C (ISO 9899:1999) 6.8.6.4/1:
   * "A return statement with an expression shall not appear in a function whose return type is void." */
  PARMCI_Barrier();
  return;
}

int ARMCI_Get(void *src, void *dst, int size, int target) {
    parmci_calls++;
    return PARMCI_Get(src, dst, size, target);
}

int ARMCI_Put(void *src, void *dst, int size, int target) {
    parmci_calls++;
    return PARMCI_Put(src, dst, size, target);
}
