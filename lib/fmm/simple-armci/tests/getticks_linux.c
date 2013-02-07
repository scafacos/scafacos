#include <time.h>

unsigned long long getticks(void)
{
    return (unsigned long long) clock();
}

