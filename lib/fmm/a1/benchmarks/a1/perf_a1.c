/* 
 *    Author: Jialin Ju, PNNL
 */

/* $Id: perf.c,v 1.21 2006-09-12 23:21:21 andriy Exp $ */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <unistd.h>
#include <math.h>
#include "a1.h"

#define SIZE 550
#define MAXPROC 8
#define CHUNK_NUM 28

#ifndef ABS
#define ABS(a) ((a)>0? (a): -(a))
#endif

int CHECK_RESULT=0;

int chunk[CHUNK_NUM] = {1,3,4,6,9,12,16,20,24,30,40,48,52,64,78,91,104,
                        128,142,171,210,256,300,353,400,440,476,512};

char check_type[15];
int nproc, me;
int warn_accuracy=0;

void fill_array(double *arr, int count, int which);
void check_result(double *src_buf, double *dst_buf, int *stride, int *count,
                  int stride_levels);
void acc_array(double scale, double *array1, double *array2, int *stride,
               int *count, int stride_levels);

#define TIMER A1_Time_seconds

double time_get(double *src_buf, double *dst_buf, int chunk, int loop,
                int proc, int levels)
{
    int i, bal = 0;
    
    int stride[2];
    int count[2];
    int stride_levels = levels;
    double *tmp_buf, *tmp_buf_ptr;
    
    double start_time, stop_time, total_time = 0;

    stride[0] = SIZE * sizeof(double);
    count[0] = chunk * sizeof(double); count[1] = chunk;

    if(CHECK_RESULT) {
        tmp_buf = (double *)malloc(SIZE * SIZE * sizeof(double));
        assert(tmp_buf != NULL);

        fill_array(tmp_buf, SIZE*SIZE, proc);
        tmp_buf_ptr = tmp_buf;
    }
    
    start_time = TIMER();
    for(i=0; i<loop; i++) {
         
        if(levels)
           A1_GetS(proc, stride_levels, count, src_buf, stride, dst_buf, stride);
        else
           A1_Get(proc, src_buf, dst_buf, count[0]);

        if(CHECK_RESULT) {
            sprintf(check_type, "A1_GetS:");
            check_result(tmp_buf_ptr, dst_buf, stride, count, stride_levels);
        }
        
        /* prepare next src and dst ptrs: avoid cache locality */
        if(bal == 0) {
            src_buf += 128;
            dst_buf += 128;
            if(CHECK_RESULT) tmp_buf_ptr += 128;
            bal = 1;
        } else {
            src_buf -= 128;
            dst_buf -= 128;
            if(CHECK_RESULT) tmp_buf_ptr -= 128;
            bal = 0;
        }
    }
    stop_time = TIMER();
    total_time = (stop_time - start_time);

    if(CHECK_RESULT) free(tmp_buf);

    if(total_time == 0.0){
       total_time=0.000001; /* workaround for inaccurate timers */
       warn_accuracy++;
    }
    return(total_time/loop);
}

double time_put(double *src_buf, double *dst_buf, int chunk, int loop,
                int proc, int levels)
{
    int i, bal = 0;

    int stride[2];
    int count[2];
    int stride_levels = levels;
    double *tmp_buf;

    double start_time, stop_time, total_time = 0;

    stride[0] = SIZE * sizeof(double);
    count[0] = chunk * sizeof(double); count[1] = chunk;

    if(CHECK_RESULT) {
        tmp_buf = (double *)malloc(SIZE * SIZE * sizeof(double));
        assert(tmp_buf != NULL);
    }
    
    start_time = TIMER();
    for(i=0; i<loop; i++) {
        if(levels)
           A1_PutS(proc, stride_levels, count, src_buf, stride, dst_buf, stride);
        else
           A1_Put(proc, src_buf, dst_buf, count[0]);

        if(CHECK_RESULT) {
            A1_Flush_group(A1_GROUP_WORLD); 
            A1_GetS(proc, stride_levels, count, dst_buf, stride, tmp_buf, stride);

            sprintf(check_type, "A1_PutS:");
            check_result(tmp_buf, src_buf, stride, count, stride_levels);
        }
        
        /* prepare next src and dst ptrs: avoid cache locality */
        if(bal == 0) {
            src_buf += 128;
            dst_buf += 128;
            bal = 1;
        } else {
            src_buf -= 128;
            dst_buf -= 128;
            bal = 0;
        }
    }
    stop_time = TIMER();
    total_time = (stop_time - start_time);

    if(CHECK_RESULT) free(tmp_buf);
    
    if(total_time == 0.0){ 
       total_time=0.000001; /* workaround for inaccurate timers */
       warn_accuracy++;
    }
    return(total_time/loop);
}

double time_put_remote(double *src_buf, double *dst_buf, int chunk, int loop,
                       int proc, int levels)
{
    int i, bal = 0;

    int stride[2];
    int count[2];
    int stride_levels = levels;
    double *tmp_buf;

    double start_time, stop_time, total_time = 0;

    stride[0] = SIZE * sizeof(double);
    count[0] = chunk * sizeof(double); count[1] = chunk;

    if(CHECK_RESULT) {
        tmp_buf = (double *)malloc(SIZE * SIZE * sizeof(double));
        assert(tmp_buf != NULL);
    }
    
    start_time = TIMER();
    for(i=0; i<loop; i++) {
        if(levels)
           A1_PutS(proc, stride_levels, count, src_buf, stride, dst_buf, stride);
        else
           A1_Put(proc, src_buf, dst_buf, count[0]);

        A1_Flush_group(A1_GROUP_WORLD);

        if(CHECK_RESULT) {
            A1_GetS(proc, stride_levels, count, dst_buf, stride, tmp_buf, stride);

            sprintf(check_type, "A1_PutS:");
            check_result(tmp_buf, src_buf, stride, count, stride_levels);
        }
        
        /* prepare next src and dst ptrs: avoid cache locality */
        if(bal == 0) {
            src_buf += 128;
            dst_buf += 128;
            bal = 1;
        } else {
            src_buf -= 128;
            dst_buf -= 128;
            bal = 0;
        }
    }
    stop_time = TIMER();
    total_time = (stop_time - start_time);

    if(CHECK_RESULT) free(tmp_buf);
    
    if(total_time == 0.0){ 
       total_time=0.000001; /* workaround for inaccurate timers */
       warn_accuracy++;
    }
    return(total_time/loop);
}

double time_acc(double *src_buf, double *dst_buf, int chunk, int loop,
                int proc, int levels)
{
    int i, bal = 0;

    int stride[2];
    int count[2];
    int stride_levels = levels;
    double *before_buf, *after_buf;
    
    double start_time, stop_time, total_time = 0;

    stride[0] = SIZE * sizeof(double);
    count[0] = chunk * sizeof(double); count[1] = chunk;

    if(CHECK_RESULT) {
        before_buf = (double *)malloc(SIZE * SIZE * sizeof(double));
        assert(before_buf != NULL);
        after_buf = (double *)malloc(SIZE * SIZE * sizeof(double));
        assert(after_buf != NULL);
    }
    
    start_time = TIMER();
    for(i=0; i<loop; i++) {
        double scale = (double)i;

        if(CHECK_RESULT) {
            A1_GetS(proc, stride_levels, count, dst_buf, stride, before_buf, stride);

            acc_array(scale, before_buf, src_buf, stride, count,stride_levels);
        }

        A1_PutAccS(proc, stride_levels, count, src_buf, stride, dst_buf, stride,
                   A1_DOUBLE, &scale);

        if(CHECK_RESULT) {
            A1_Flush_group(A1_GROUP_WORLD);
            A1_GetS(proc, stride_levels, count, dst_buf, stride, after_buf, stride);
            
            sprintf(check_type, "A1_AccS:");
            check_result(after_buf, before_buf, stride, count, stride_levels);
        }
        
        /* prepare next src and dst ptrs: avoid cache locality */
        if(bal == 0) {
            src_buf += 128;
            dst_buf += 128;
            bal = 1;
        } else {
            src_buf -= 128;
            dst_buf -= 128;
            bal = 0;
        }
    }
    stop_time = TIMER();
    total_time = (stop_time - start_time);

    if(CHECK_RESULT) { free(before_buf); free(after_buf); }
    
    if(total_time == 0.0){ 
       total_time=0.000001; /* workaround for inaccurate timers */
       warn_accuracy++;
    }
    return(total_time/loop);
}

double time_acc_remote(double *src_buf, double *dst_buf, int chunk, int loop,
                       int proc, int levels)
{
    int i, bal = 0;

    int stride[2];
    int count[2];
    int stride_levels = levels;
    double *before_buf, *after_buf;
    
    double start_time, stop_time, total_time = 0;

    stride[0] = SIZE * sizeof(double);
    count[0] = chunk * sizeof(double); count[1] = chunk;

    if(CHECK_RESULT) {
        before_buf = (double *)malloc(SIZE * SIZE * sizeof(double));
        assert(before_buf != NULL);
        after_buf = (double *)malloc(SIZE * SIZE * sizeof(double));
        assert(after_buf != NULL);
    }
    
    start_time = TIMER();
    for(i=0; i<loop; i++) {
        double scale = (double)i;

        if(CHECK_RESULT) {
            A1_GetS(proc, stride_levels, count, dst_buf, stride, before_buf, stride);

            acc_array(scale, before_buf, src_buf, stride, count,stride_levels);
        }

        A1_PutAccS(proc, stride_levels, count, src_buf, stride, dst_buf, stride,
                   A1_DOUBLE, &scale);

        A1_Flush_group(A1_GROUP_WORLD);

        if(CHECK_RESULT) {
            A1_GetS(proc, stride_levels, count, dst_buf, stride, after_buf, stride);
            
            sprintf(check_type, "A1_AccS:");
            check_result(after_buf, before_buf, stride, count, stride_levels);
        }
        
        /* prepare next src and dst ptrs: avoid cache locality */
        if(bal == 0) {
            src_buf += 128;
            dst_buf += 128;
            bal = 1;
        } else {
            src_buf -= 128;
            dst_buf -= 128;
            bal = 0;
        }
    }
    stop_time = TIMER();
    total_time = (stop_time - start_time);

    if(CHECK_RESULT) { free(before_buf); free(after_buf); }
    
    if(total_time == 0.0){ 
       total_time=0.000001; /* workaround for inaccurate timers */
       warn_accuracy++;
    }
    return(total_time/loop);
}

void test_1D()
{
    int i;
    int src, dst;
    int ierr;
    double *buf;
    void *ptr[MAXPROC], *get_ptr[MAXPROC];

    /* find who I am and the dst process */
    src = me;
    
    /* memory allocation */
    if(me == 0) {
        buf = (double *)malloc(SIZE * SIZE * sizeof(double));
        assert(buf != NULL);
    }
   
    ierr = A1_Alloc_segment(&ptr[me], (SIZE * SIZE * sizeof(double))); 
    assert(ierr == 0); 
    ierr = A1_Exchange_segments(A1_GROUP_WORLD, ptr);    
    assert(ierr == 0); 
    ierr = A1_Alloc_segment(&get_ptr[me], (SIZE * SIZE * sizeof(double)));      
    assert(ierr == 0); 
    ierr = A1_Exchange_segments(A1_GROUP_WORLD, get_ptr);
    assert(ierr == 0); 

    /* A1 - initialize the data window */
    fill_array(ptr[me], SIZE*SIZE, me);
    fill_array(get_ptr[me], SIZE*SIZE, me);
    A1_Barrier_group(A1_GROUP_WORLD);
    
    /* only the proc 0 does the work */
    if(me == 0) {
        if(!CHECK_RESULT){
          printf("  section               get                 put");
          printf("                 acc\n");
          printf("bytes   loop       usec      MB/s       usec      MB/s");
          printf("       usec      MB/s\n");
          printf("------- ------  --------  --------  --------  --------");
          printf("  --------  --------\n");
          fflush(stdout);
        }
        
        for(i=0; i<CHUNK_NUM; i++) {
            int loop;
            int bytes = chunk[i] * chunk[i] * sizeof(double);
            
            double t_get = 0, t_put = 0, t_acc = 0;
            double latency_get, latency_put, latency_acc;
            double bandwidth_get, bandwidth_put, bandwidth_acc;
            
            loop = (SIZE * SIZE) / (chunk[i] * chunk[i]);
            loop = (int)sqrt((double)loop);
            if(loop<2)loop=2;
            
            for(dst=1; dst<nproc; dst++) {
                /* strided get */
                fill_array(buf, SIZE*SIZE, me*10);
                t_get += time_get((double *)(get_ptr[dst]), (double *)buf,
                                  chunk[i]*chunk[i], loop, dst, 0);
                
                /* strided put */
                fill_array(buf, SIZE*SIZE, me*10);
                t_put += time_put((double *)buf, (double *)(ptr[dst]),
                                  chunk[i]*chunk[i], loop, dst, 0);
                
                /* strided acc */
                fill_array(buf, SIZE*SIZE, me*10);
                t_acc += time_acc((double *)buf, (double *)(ptr[dst]),
                                  chunk[i]*chunk[i], loop, dst, 0);
            }
            
            latency_get = t_get/(nproc - 1);
            latency_put = t_put/(nproc - 1);
            latency_acc = t_acc/(nproc - 1);
            
            bandwidth_get = (bytes * (nproc - 1) * 1e-6)/t_get;
            bandwidth_put = (bytes * (nproc - 1) * 1e-6)/t_put;
            bandwidth_acc = (bytes * (nproc - 1) * 1e-6)/t_acc;

            /* print */
            if(!CHECK_RESULT)
                   printf("%d\t%d\t %7.2lf %9.2lf %9.2lf %9.2lf %9.2lf  %9.2lf\n",
                   bytes, loop, latency_get*1000000, bandwidth_get,
                   latency_put*1000000, bandwidth_put, latency_acc*1000000, bandwidth_acc);
        }
    }
    else sleep(60);
    
    A1_Flush_group(A1_GROUP_WORLD);
    A1_Barrier_group(A1_GROUP_WORLD);
    
    /* cleanup */
    A1_Release_segments(A1_GROUP_WORLD, get_ptr[me]);
    A1_Free_segment(get_ptr[me]);
    A1_Release_segments(A1_GROUP_WORLD, ptr[me]);
    A1_Free_segment(ptr[me]);
    
    if(me == 0) free(buf);
}

void test_1D_remote()
{
    int i;
    int src, dst;
    int ierr;
    double *buf;
    void *ptr[MAXPROC], *get_ptr[MAXPROC];

    /* find who I am and the dst process */
    src = me;
    
    /* memory allocation */
    if(me == 0) {
        buf = (double *)malloc(SIZE * SIZE * sizeof(double));
        assert(buf != NULL);
    }
   
    ierr = A1_Alloc_segment(&ptr[me], (SIZE * SIZE * sizeof(double))); 
    assert(ierr == 0); 
    ierr = A1_Exchange_segments(A1_GROUP_WORLD, ptr);    
    assert(ierr == 0); 
    ierr = A1_Alloc_segment(&get_ptr[me], (SIZE * SIZE * sizeof(double)));      
    assert(ierr == 0); 
    ierr = A1_Exchange_segments(A1_GROUP_WORLD, get_ptr);
    assert(ierr == 0); 

    /* A1 - initialize the data window */
    fill_array(ptr[me], SIZE*SIZE, me);
    fill_array(get_ptr[me], SIZE*SIZE, me);
    A1_Barrier_group(A1_GROUP_WORLD);
    
    /* only the proc 0 does the work */
    if(me == 0) {
        if(!CHECK_RESULT){
          printf("  section               get                 put");
          printf("                 acc\n");
          printf("bytes   loop       usec      MB/s       usec      MB/s");
          printf("       usec      MB/s\n");
          printf("------- ------  --------  --------  --------  --------");
          printf("  --------  --------\n");
          fflush(stdout);
        }
        
        for(i=0; i<CHUNK_NUM; i++) {
            int loop;
            int bytes = chunk[i] * chunk[i] * sizeof(double);
            
            double t_get = 0, t_put = 0, t_acc = 0;
            double latency_get, latency_put, latency_acc;
            double bandwidth_get, bandwidth_put, bandwidth_acc;
            
            loop = (SIZE * SIZE) / (chunk[i] * chunk[i]);
            loop = (int)sqrt((double)loop);
            if(loop<2)loop=2;
            
            for(dst=1; dst<nproc; dst++) {
                /* strided get */
                fill_array(buf, SIZE*SIZE, me*10);
                t_get += time_get((double *)(get_ptr[dst]), (double *)buf,
                                  chunk[i]*chunk[i], loop, dst, 0);
                
                /* strided put */
                fill_array(buf, SIZE*SIZE, me*10);
                t_put += time_put_remote((double *)buf, (double *)(ptr[dst]),
                                          chunk[i]*chunk[i], loop, dst, 0);
                
                /* strided acc */
                fill_array(buf, SIZE*SIZE, me*10);
                t_acc += time_acc_remote((double *)buf, (double *)(ptr[dst]),
                                        chunk[i]*chunk[i], loop, dst, 0);
            }
            
            latency_get = t_get/(nproc - 1);
            latency_put = t_put/(nproc - 1);
            latency_acc = t_acc/(nproc - 1);
            
            bandwidth_get = (bytes * (nproc - 1) * 1e-6)/t_get;
            bandwidth_put = (bytes * (nproc - 1) * 1e-6)/t_put;
            bandwidth_acc = (bytes * (nproc - 1) * 1e-6)/t_acc;

            /* print */
            if(!CHECK_RESULT)
                   printf("%d\t%d\t %7.2lf %9.2lf %9.2lf %9.2lf %9.2lf  %9.2lf\n",
                   bytes, loop, latency_get*1000000, bandwidth_get,
                   latency_put*1000000, bandwidth_put, latency_acc*1000000, bandwidth_acc);
        }
    }
    else sleep(60);
    
    A1_Flush_group(A1_GROUP_WORLD);
    A1_Barrier_group(A1_GROUP_WORLD);
    
    /* cleanup */
    A1_Release_segments(A1_GROUP_WORLD, get_ptr[me]);
    A1_Free_segment(get_ptr[me]);
    A1_Release_segments(A1_GROUP_WORLD, ptr[me]);
    A1_Free_segment(ptr[me]);
    
    if(me == 0) free(buf);
}

void test_2D()
{
    int i;
    int src, dst;
    int ierr;
    double *buf;
    void *ptr[MAXPROC], *get_ptr[MAXPROC];

    /* find who I am and the dst process */
    src = me;
    
    if(me == 0) {
        buf = (double *)malloc(SIZE * SIZE * sizeof(double));
        assert(buf != NULL);
    }

    ierr = A1_Alloc_segment(&ptr[me], (SIZE * SIZE * sizeof(double)));
    assert(ierr == 0);
    ierr = A1_Exchange_segments(A1_GROUP_WORLD, ptr);
    assert(ierr == 0);
    ierr = A1_Alloc_segment(&get_ptr[me], (SIZE * SIZE * sizeof(double)));
    assert(ierr == 0);
    ierr = A1_Exchange_segments(A1_GROUP_WORLD, get_ptr);
    assert(ierr == 0);
    
    /* A1 - initialize the data window */
    fill_array(ptr[me], SIZE*SIZE, me);
    fill_array(get_ptr[me], SIZE*SIZE, me);

    A1_Barrier_group(A1_GROUP_WORLD);
    
    /* only the proc 0 doest the work */
    /* print the title */
    if(me == 0) {
        if(!CHECK_RESULT){
           printf("  section               get                 put");
           printf("                 acc\n");
           printf("bytes   loop       usec      MB/s       usec      MB/s");
           printf("       usec      MB/s\n");
           printf("------- ------  --------  --------  --------  --------");
           printf("  --------  --------\n");
           fflush(stdout);
        }
        
        for(i=0; i<CHUNK_NUM; i++) {
            int loop;
            int bytes = chunk[i] * chunk[i] * sizeof(double);

            double t_get = 0, t_put = 0, t_acc = 0;
            double latency_get, latency_put, latency_acc;
            double bandwidth_get, bandwidth_put, bandwidth_acc;
            
            loop = SIZE / chunk[i];
            if(loop<2)loop=2;

            for(dst=1; dst<nproc; dst++) {
                /* strided get */
                fill_array(buf, SIZE*SIZE, me*10);
                t_get += time_get((double *)(get_ptr[dst]), (double *)buf,
                                 chunk[i], loop, dst, 1);
 
                /* strided put */
                fill_array(buf, SIZE*SIZE, me*10);
                t_put += time_put((double *)buf, (double *)(ptr[dst]),
                                 chunk[i], loop, dst, 1);
                
                /* strided acc */
                fill_array(buf, SIZE*SIZE, me*10);
                t_acc += time_acc((double *)buf, (double *)(ptr[dst]),
                                 chunk[i], loop, dst, 1);
            }
            
            latency_get = t_get/(nproc - 1);
            latency_put = t_put/(nproc - 1);
            latency_acc = t_acc/(nproc - 1);
            
            bandwidth_get = (bytes * (nproc - 1) * 1e-6)/t_get;
            bandwidth_put = (bytes * (nproc - 1) * 1e-6)/t_put;
            bandwidth_acc = (bytes * (nproc - 1) * 1e-6)/t_acc;

            /* print */
            if(!CHECK_RESULT)
                   printf("%d\t%d\t %7.2lf %9.2lf %9.2lf %9.2lf %9.2lf  %9.2lf\n",
                   bytes, loop, latency_get*1000000, bandwidth_get,
                   latency_put*1000000, bandwidth_put, latency_acc*1000000, bandwidth_acc);
        }
    }
    else sleep(60);
    
    A1_Flush_group(A1_GROUP_WORLD);
    A1_Barrier_group(A1_GROUP_WORLD);

    /* cleanup */
    A1_Release_segments(A1_GROUP_WORLD, get_ptr[me]);
    A1_Free_segment(get_ptr[me]);
    A1_Release_segments(A1_GROUP_WORLD, ptr[me]);
    A1_Free_segment(ptr[me]);

    if(me == 0) free(buf);

}

void test_2D_remote()
{
    int i;
    int src, dst;
    int ierr;
    double *buf;
    void *ptr[MAXPROC], *get_ptr[MAXPROC];

    /* find who I am and the dst process */
    src = me;
    
    if(me == 0) {
        buf = (double *)malloc(SIZE * SIZE * sizeof(double));
        assert(buf != NULL);
    }

    ierr = A1_Alloc_segment(&ptr[me], (SIZE * SIZE * sizeof(double)));
    assert(ierr == 0);
    ierr = A1_Exchange_segments(A1_GROUP_WORLD, ptr);
    assert(ierr == 0);
    ierr = A1_Alloc_segment(&get_ptr[me], (SIZE * SIZE * sizeof(double)));
    assert(ierr == 0);
    ierr = A1_Exchange_segments(A1_GROUP_WORLD, get_ptr);
    assert(ierr == 0);
    
    /* A1 - initialize the data window */
    fill_array(ptr[me], SIZE*SIZE, me);
    fill_array(get_ptr[me], SIZE*SIZE, me);

    A1_Barrier_group(A1_GROUP_WORLD);
    
    /* only the proc 0 doest the work */
    /* print the title */
    if(me == 0) {
        if(!CHECK_RESULT){
           printf("  section               get                 put");
           printf("                 acc\n");
           printf("bytes   loop       usec      MB/s       usec      MB/s");
           printf("       usec      MB/s\n");
           printf("------- ------  --------  --------  --------  --------");
           printf("  --------  --------\n");
           fflush(stdout);
        }
        
        for(i=0; i<CHUNK_NUM; i++) {
            int loop;
            int bytes = chunk[i] * chunk[i] * sizeof(double);

            double t_get = 0, t_put = 0, t_acc = 0;
            double latency_get, latency_put, latency_acc;
            double bandwidth_get, bandwidth_put, bandwidth_acc;
            
            loop = SIZE / chunk[i];
            if(loop<2)loop=2;

            for(dst=1; dst<nproc; dst++) {
                /* strided get */
                fill_array(buf, SIZE*SIZE, me*10);
                t_get += time_get((double *)(get_ptr[dst]), (double *)buf,
                                 chunk[i], loop, dst, 1);
 
                /* strided put */
                fill_array(buf, SIZE*SIZE, me*10);
                t_put += time_put_remote((double *)buf, (double *)(ptr[dst]),
                                         chunk[i], loop, dst, 1);
                
                /* strided acc */
                fill_array(buf, SIZE*SIZE, me*10);
                t_acc += time_acc_remote((double *)buf, (double *)(ptr[dst]),
                                         chunk[i], loop, dst, 1);
            }
            
            latency_get = t_get/(nproc - 1);
            latency_put = t_put/(nproc - 1);
            latency_acc = t_acc/(nproc - 1);
            
            bandwidth_get = (bytes * (nproc - 1) * 1e-6)/t_get;
            bandwidth_put = (bytes * (nproc - 1) * 1e-6)/t_put;
            bandwidth_acc = (bytes * (nproc - 1) * 1e-6)/t_acc;

            /* print */
            if(!CHECK_RESULT)
                   printf("%d\t%d\t %7.2lf %9.2lf %9.2lf %9.2lf %9.2lf  %9.2lf\n",
                   bytes, loop, latency_get*1000000, bandwidth_get,
                   latency_put*1000000, bandwidth_put, latency_acc*1000000, bandwidth_acc);
        }
    }
    else sleep(60);
    
    A1_Flush_group(A1_GROUP_WORLD);
    A1_Barrier_group(A1_GROUP_WORLD);

    /* cleanup */
    A1_Release_segments(A1_GROUP_WORLD, get_ptr[me]);
    A1_Free_segment(get_ptr[me]);
    A1_Release_segments(A1_GROUP_WORLD, ptr[me]);
    A1_Free_segment(ptr[me]);

    if(me == 0) free(buf);

}
    
int main(int argc, char **argv)
{

    /* initialize A1 */
    A1_Initialize(A1_THREAD_SINGLE);

    me = A1_Process_id(A1_GROUP_WORLD);
    nproc = A1_Process_total(A1_GROUP_WORLD);

    if(nproc < 2 || nproc> MAXPROC) {
        if(me == 0)
            fprintf(stderr,
                    "USAGE: 2 <= processes <= %d - got %d\n", MAXPROC, nproc);
        A1_Barrier_group(A1_GROUP_WORLD);
        exit(0);
    }
    
    if(!me)printf("\n             Performance of Basic Blocking Communication Operations\n");
    A1_Barrier_group(A1_GROUP_WORLD);
    
    CHECK_RESULT=1; test_1D(); CHECK_RESULT=0; /* warmup run */

    /* test 1 dimension array */
    if(!me)printf("\n\t\t\tContiguous Data Transfer\n");
    test_1D();

    /* test 1 dimension array */
    if(!me)printf("\n\t\t\tContiguous Data Transfer - Remote completion\n");
    test_1D_remote();
    
    /* test 2 dimension array */
    if(!me)printf("\n\t\t\tStrided Data Transfer\n");
    test_2D();

    /* test 2 dimension array */
    if(!me)printf("\n\t\t\tStrided Data Transfer - Remote completion\n");
    test_2D_remote();

    A1_Barrier_group(A1_GROUP_WORLD);
    if(me == 0){
       if(warn_accuracy) 
          printf("\nWARNING: Your timer does not have sufficient accuracy for this test (%d)\n",warn_accuracy);
       printf("\n\n------------ Now we test the same data transfer for correctness ----------\n");
       fflush(stdout);
    }

    A1_Barrier_group(A1_GROUP_WORLD);
    CHECK_RESULT=1;
    if(!me)printf("\n\t\t\tContiguous Data Transfer\n");
    test_1D();
    if(me == 0) printf("OK\n");
    A1_Barrier_group(A1_GROUP_WORLD);
    if(!me)printf("\n\t\t\tStrided Data Transfer\n");
    test_2D();
    if(me == 0) printf("OK\n\n\nTests Completed.\n");
    A1_Barrier_group(A1_GROUP_WORLD);

    /* done */
    A1_Finalize();
    return(0);
}    

void fill_array(double *arr, int count, int which)
{
    int i;

    for(i=0; i<count; i++) arr[i] = i * 8.23 + which * 2.89;
}

void check_result(double *src_buf, double *dst_buf, int *stride, int *count,
                  int stride_levels)
{
    int i, j, size;
    long idx;
    int n1dim;  /* number of 1 dim block */
    int bvalue[4], bunit[4];

    /* number of n-element of the first dimension */
    n1dim = 1;
    for(i=1; i<=stride_levels; i++)
        n1dim *= count[i];

    /* calculate the destination indices */
    bvalue[0] = 0; bvalue[1] = 0; bunit[0] = 1; bunit[1] = 1;
    for(i=2; i<=stride_levels; i++) {
        bvalue[i] = 0;
        bunit[i] = bunit[i-1] * count[i-1];
    }

    for(i=0; i<n1dim; i++) {
        idx = 0;
        for(j=1; j<=stride_levels; j++) {
            idx += bvalue[j] * stride[j-1];
            if((i+1) % bunit[j] == 0) bvalue[j]++;
            if(bvalue[j] > (count[j]-1)) bvalue[j] = 0;
        }
        
        size = count[0] / sizeof(double);
        for(j=0; j<size; j++)
            if(ABS(((double *)((char *)src_buf+idx))[j] - 
               ((double *)((char *)dst_buf+idx))[j]) > 0.000001 ){
                fprintf(stdout,"Error:%s comparison failed: (%d) (%f :%f) %d\n",
                        check_type, j, ((double *)((char *)src_buf+idx))[j],
                        ((double *)((char *)dst_buf+idx))[j], count[0]);
                A1_Abort(0, "failed");
            }
    }
}

/* array1 = array1 + array2 * scale */
void acc_array(double scale, double *array1, double *array2, int *stride,
               int *count, int stride_levels)
{
        int i, j, size;
    long idx;
    int n1dim;  /* number of 1 dim block */
    int bvalue[4], bunit[4];

    /* number of n-element of the first dimension */
    n1dim = 1;
    for(i=1; i<=stride_levels; i++)
        n1dim *= count[i];

    /* calculate the destination indices */
    bvalue[0] = 0; bvalue[1] = 0; bunit[0] = 1; bunit[1] = 1;
    for(i=2; i<=stride_levels; i++) {
        bvalue[i] = 0;
        bunit[i] = bunit[i-1] * count[i-1];
    }

    for(i=0; i<n1dim; i++) {
        idx = 0;
        for(j=1; j<=stride_levels; j++) {
            idx += bvalue[j] * stride[j-1];
            if((i+1) % bunit[j] == 0) bvalue[j]++;
            if(bvalue[j] > (count[j]-1)) bvalue[j] = 0;
        }

        size = count[0] / sizeof(double);
        for(j=0; j<size; j++)
            ((double *)((char *)array1+idx))[j] +=
                ((double *)((char *)array2+idx))[j] * scale;

    }
}
