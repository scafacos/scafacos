#include <stdio.h>
#include <stdlib.h>
#include <math.h>


int main(int argc, char* argv[])
{
    int i, dim=1024;
    double* a;
    posix_memalign((void**)&a, 128*sizeof(double), dim*sizeof(double));

    printf("%d\n",(int)a); 
    printf("%d\n",(int)&a[0]); 
    printf("%d\n",(int)&a[1]); 
    printf("%d\n",((int)a)%16); 
    printf("%d\n",((int)&a[0])%16); 
    printf("%d\n",((int)&a[1])%16); 

    fprintf(stderr,"ALL DONE\n");

    return(0);
}


