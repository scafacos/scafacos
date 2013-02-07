#include <stdlib.h>
/* TODO Scafacos configre should get alignmenta automatically*/
#define FMM_ALIGNMENT 16

void* dummy_malloc(int bsize, int *ierr)
{
   void *ptr = malloc(bsize);

   if (ptr) {
     *ierr = 0;
   } else {
     *ierr = 1;
   }

   return ptr;
}

void* dummy_malloc_aligned(int bsize, int *ierr)
{
   void ** ptr;
   int err = posix_memalign((void**)&ptr, FMM_ALIGNMENT, (size_t) bsize);

   if (!err) {
     *ierr = 0;
   } else {
     *ierr = 1;
   }

   return ptr;
}

void dummy_free(void *ptr)
{
free(ptr);
}
