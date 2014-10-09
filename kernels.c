/********************************************************
 * Kernels to be optimized for the CS:APP Performance Lab
 ********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <sys/times.h>
#include <time.h>
#include "defs.h"

/* Misc constants */
#define BSIZE 32     /* cache block size in bytes */     
#define MIN(a,b) (((a) < (b)) ? (a) : (b))
#define DEBUG_PROFILE_TIME 0
#define DEBUG_PROFILE_TIMES 0
#define DEBUG_PRINT DEBUG_PROFILE_TIME || DEBUG_PROFILE_TIMES
#define log(args ...) do {{printf(args); printf("\n");}} while (0)

/* 
 * ECE454 Students: 
 * Please fill in the following team struct 
 */
team_t team = {
    "PrinterOnFire",            /* Team name */

    "Jun Tao Luo",              /* First member full name */
    "jt.luo@mail.utoronto.ca",  /* First member email address */

    "Chang Liu",                /* Second member full name (leave blank if none) */
    "cha.liu@mail.utoronto.ca"  /* Second member email addr (leave blank if none) */
};

/***************
 * ROTATE KERNEL
 ***************/

/******************
 * Helper Functions
 ******************/

unsigned get_seconds() {
    struct tms t;
    times(&t);
    
    return t.tms_utime;
}

// for benchmarking 
void copy_block(int dim, pixel *src, pixel *dst, int B) 
{
    int i, j, ii, jj, ii_limit, jj_limit;

    for (i = 0; i < dim; i+=B)
    {
        ii_limit = MIN(dim, i+B);
        for (j = 0; j < dim; j+=B)
        {
            jj_limit = MIN(dim, j+B);

            int index;
            
            for (ii = i; ii < ii_limit; ii++)
            {
                index = j*dim + ii;
                
                for (jj = j; jj < jj_limit; jj++)
                {
                    dst[index] = src[index];
                    index ++;
                }
            }
        }
    }
}

// two step part 1
void transpose_block(int dim, pixel *src, pixel *dst, int B) 
{
    int i, j, ii, jj, ii_limit, jj_limit, dst_index, src_index;

    for (i = 0; i < dim; i+=B)
    {
        ii_limit = MIN(dim, i+B);
        for (j = 0; j < dim; j+=B)
        {
            jj_limit = MIN(dim, j+B);
            
            for (ii = i; ii < ii_limit; ii++)
            {
                dst_index = j*dim + ii;
                src_index = ii*dim + j;
                
                for (jj = j; jj < jj_limit; jj++)
                {
                    dst[dst_index] = src[src_index];
                    src_index ++;
                    dst_index += dim;
                }
            }
        }
    }
}

// two step part 2: in-place version
void flip_ip(int dim, pixel *src) 
{
    int i, j, midpoint = dim/2;
    pixel *a=src + dim*(dim - 1), *b=src, temp;

    for (i = 0; i < midpoint; i++)
    {
        for (j = 0; j < dim; j+=4) 
        {
            temp = a[0];
            a[0] = b[0];
            b[0] = temp;

            temp = a[1];
            a[1] = b[1];
            b[1] = temp;

            temp = a[2];
            a[2] = b[2];
            b[2] = temp;

            temp = a[3];
            a[3] = b[3];
            b[3] = temp;

            a += 4;
            b += 4;
        }
        a -= dim<<1;
    }
}


/******************************************************
 * Your different versions of the rotate kernel go here
 ******************************************************/

/* 
 * new_rotate - The naive baseline version of rotate with optimizations
 */
char new_rotate_descr[] = "new_rotate: Naive baseline implementation with optimizations";
void new_rotate(int dim, pixel *src, pixel *dst, int iB, int jB) 
{
    int i, j, ii, jj, ii_limit, jj_limit;

    for (i = 0; i < dim; i+=iB)
    {
        ii_limit = MIN(dim, i+iB);
        for (j = 0; j < dim; j+=jB)
        {
            jj_limit = MIN(dim, j+jB);
            for (ii = i; ii < ii_limit; ii++)
            {
                for (jj = j; jj < jj_limit; jj++)
                {
	                dst[RIDX(dim-1-ii, jj, dim)] = src[RIDX(jj, ii, dim)];
                }
                ii++;
                for (jj = j; jj < jj_limit; jj++)
                {
	                dst[RIDX(dim-1-ii, jj, dim)] = src[RIDX(jj, ii, dim)];
                }
                ii++;
                for (jj = j; jj < jj_limit; jj++)
                {
	                dst[RIDX(dim-1-ii, jj, dim)] = src[RIDX(jj, ii, dim)];
                }
                ii++;
                for (jj = j; jj < jj_limit; jj++)
                {
	                dst[RIDX(dim-1-ii, jj, dim)] = src[RIDX(jj, ii, dim)];
                }
            }
        }
    }
}

// does not help since loops are already unrolled
char naive_rotate_unrolled_descr[] = "naive_rotate_unrolled: Naive baseline implementation with unrolling";
void naive_rotate_unrolled(int dim, pixel *src, pixel *dst) 
{
    int i, j;

    for (i = 0; i < dim; i++)
        for (j = 0; j < dim; j+=4)
        {
            dst[RIDX(dim-1-j, i, dim)] = src[0];
            dst[RIDX(dim-1-(j + 1), i, dim)] = src[1];
            dst[RIDX(dim-1-(j + 2), i, dim)] = src[2];
            dst[RIDX(dim-1-(j + 3), i, dim)] = src[3];
//             dst[RIDX(dim-1-(j + 4), i, dim)] = *(src + 4);
//             dst[RIDX(dim-1-(j + 5), i, dim)] = *(src + 5);
//             dst[RIDX(dim-1-(j + 6), i, dim)] = *(src + 6);
//             dst[RIDX(dim-1-(j + 7), i, dim)] = *(src + 7);
//             dst[RIDX(dim-1-(j + 8), i, dim)] = *(src + 8);
//             dst[RIDX(dim-1-(j + 9), i, dim)] = *(src + 9);
//             dst[RIDX(dim-1-(j + 10), i, dim)] = *(src + 10);
//             dst[RIDX(dim-1-(j + 11), i, dim)] = *(src + 11);
//             dst[RIDX(dim-1-(j + 12), i, dim)] = *(src + 12);
//             dst[RIDX(dim-1-(j + 13), i, dim)] = *(src + 13);
//             dst[RIDX(dim-1-(j + 14), i, dim)] = *(src + 14);
//             dst[RIDX(dim-1-(j + 15), i, dim)] = *(src + 15);
            src += 4;
        }
}

char default_naive_rotate_descr[] = "default_naive_rotate: Naive baseline implementation";
char default_naive_rotate_descr[] = "default_naive_rotate: Naive baseline implementation";
void default_naive_rotate(int dim, pixel *src, pixel *dst) 
{
    int i, j;

    for (i = 0; i < dim; i++)
	for (j = 0; j < dim; j++)
	    dst[RIDX(dim-1-j, i, dim)] = src[RIDX(i, j, dim)];
}

void default_naive_rotate(int dim, pixel *src, pixel *dst) 
{
    int i, j;

    for (i = 0; i < dim; i++)
	for (j = 0; j < dim; j++)
	    dst[RIDX(dim-1-j, i, dim)] = src[RIDX(i, j, dim)];
}

// for debugging and profiling
void rotate_core(int dim, pixel *src, pixel *dst, int D) 
{
    if (dim < 256)
    {
#if DEBUG_PRINT
        log("doing improved naive rotate");
#endif
        new_rotate(dim, src, dst, D*2, D/2);
    }
    else
    {
#if DEBUG_PRINT
        log("doing 2 step rotate on input size %d", dim);
#endif

#if DEBUG_PROFILE_TIME
        clock_t test, start, middle, end;
        
        test = clock();
        copy_block(dim, src, dst, D);
        start = clock();
#elif DEBUG_PROFILE_TIMES
        unsigned test, start, middle, end;
        
        test = get_seconds();
        copy_block(dim, src, dst, D);
        start = get_seconds();
#endif

        // transpose the image in src and store
        transpose_block(dim, src, dst, D);

#if DEBUG_PROFILE_TIME
        middle = clock();
#elif DEBUG_PROFILE_TIMES
        middle = get_seconds();
#endif

        // perform row exchange on the transposed image
        flip_ip(dim, dst);

#if DEBUG_PROFILE_TIME
        end = clock();
        log("test took %u transpose took %u and flip took %u", (float)start - (float)test, (float)middle - (float)start, (float)end - (float)middle);
#elif DEBUG_PROFILE_TIMES
        end = get_seconds();
        log("test took %u transpose took %u and flip took %u", start - test, middle - start, end - middle);
#endif
    }
}

/*
 * ECE 454 Students: Write your rotate functions here:
 */ 

/* 
 * rotate - Your current working version of rotate
 * IMPORTANT: This is the version you will be graded on
 */
char rotate_descr[] = "rotate: Current working version";
void rotate(int dim, pixel *src, pixel *dst) 
{
    // negativ/zero dims are invalid
    if (dim <= 0) return;

    // from empirical results
    int tile_i_size = (dim > 724) ? 128 : 256;
    int tile_j_size = (dim > 724) ? 8 : 16;

    // the best tiles are not squares because we want to keep the j tiles to be 8
    // which is the number of pixels that fit in a cache line
    new_rotate(dim, src, dst, tile_i_size, tile_j_size);
}

// searching for optimal cache sizes
char rotate_two_descr[] = "second attempt";
void attempt_two(int dim, pixel *src, pixel *dst) 
{
    rotate_core(dim, src, dst, 32);
}

char rotate_three_descr[] = "third attempt";
void attempt_three(int dim, pixel *src, pixel *dst) 
{
    rotate_core(dim, src, dst, 64);
}

char rotate_four_descr[] = "fourth attempt";
void attempt_four(int dim, pixel *src, pixel *dst) 
{
    rotate_core(dim, src, dst, 96);
}

char rotate_five_descr[] = "fifth attempt";
void attempt_five(int dim, pixel *src, pixel *dst) 
{
    rotate_core(dim, src, dst, 128);
}


/*********************************************************************
 * register_rotate_functions - Register all of your different versions
 *     of the rotate kernel with the driver by calling the
 *     add_rotate_function() for each test function. When you run the
 *     driver program, it will test and report the performance of each
 *     registered test function.  
 *********************************************************************/

void register_rotate_functions() 
{
    // benchmark this should yield a speedup of 1.0 for an unloaded computer
//     add_rotate_function(&default_naive_rotate, default_naive_rotate_descr);   
    add_rotate_function(&rotate, rotate_descr);   

//     add_rotate_function(&attempt_two, rotate_two_descr);   
//     add_rotate_function(&attempt_three, rotate_three_descr);   
//     add_rotate_function(&attempt_four, rotate_four_descr);   
//     add_rotate_function(&attempt_five, rotate_five_descr);   
    //add_rotate_function(&attempt_six, rotate_six_descr);   
    //add_rotate_function(&attempt_seven, rotate_seven_descr);   
    //add_rotate_function(&attempt_eight, rotate_eight_descr);   
    //add_rotate_function(&attempt_nine, rotate_nine_descr);   
    //add_rotate_function(&attempt_ten, rotate_ten_descr);   
    //add_rotate_function(&attempt_eleven, rotate_eleven_descr);   

    /* ... Register additional rotate functions here */
}

