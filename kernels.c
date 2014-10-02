/********************************************************
 * Kernels to be optimized for the CS:APP Performance Lab
 ********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "defs.h"

/* Misc constants */
#define BSIZE 32     /* cache block size in bytes */     
#define MIN(a,b) (((a) < (b)) ? (a) : (b))
/* 
 * ECE454 Students: 
 * Please fill in the following team struct 
 */
team_t team = {
    "PrinterOnFire",            /* Team name */

    "Jun Tao Luo",              /* First member full name */
    "jt.luo@mail.utoronto.ca",  /* First member email address */

    "Chang Liu",                /* Second member full name (leave blank if none) */
    ""                    /* Second member email addr (leave blank if none) */
};

/***************
 * ROTATE KERNEL
 ***************/

/******************
 * Helper Functions
 ******************/
void swap(pixel *a, pixel *b)
{
    pixel temp = *a;
    *a = *b;
    *b = temp;
}

void transpose(int dim, pixel *src, pixel *dst) 
{
    int i, j;
    pixel *src_px=src, *dst_col=dst, *dst_px=dst;

    // top right portion
    for (i = 0; i < dim; i++)
    {
        dst_px = dst_col++;
        for (j = 0; j < dim; j++)
        {
            *dst_px = *(src_px++);
            dst_px += dim;
        }
    }
}

void transpose_block(int dim, pixel *src, pixel *dst, int B) 
{
    int i, j, ii, jj, ii_limit, jj_limit;

    // top right portion
    for (i = 0; i < dim; i+=B)
    {
        ii_limit = MIN(dim, i+B);
        for (j = 0; j < dim; j+=B)
        {
            jj_limit = MIN(dim, j+B);
            for (ii = i; ii < ii_limit; ii++)
            {
                for (jj = j; jj < jj_limit; jj++)
                {
                    dst[RIDX(jj, ii, dim)] = src[RIDX(ii, jj, dim)];
                }
            }
        }
    }
}

void xrows(int dim, pixel *src, pixel *dst) 
{
    int i, j;

    for (i = 0; i < dim; i++)
        for (j = 0; j < dim; j++)
            dst[RIDX(dim-i-1, j, dim)] = src[RIDX(i, j, dim)];
}

// in-place version
void xrows_ip(int dim, pixel *src) 
{
    int i, j, midpoint = dim/2;
    pixel *a=src + dim*(dim - 1), *b=src, temp;

    for (i = 0; i < midpoint; i++)
    {
        for (j = 0; j < dim; j+=4) 
        {
            temp = *a;
            *a = *b;
            *b = temp;

            a++;
            b++;

            temp = *a;
            *a = *b;
            *b = temp;

            a++;
            b++;

            temp = *a;
            *a = *b;
            *b = temp;

            a++;
            b++;

            temp = *a;
            *a = *b;
            *b = temp;

            a++;
            b++;
        }
        a -= dim<<1;
    }
}


/******************************************************
 * Your different versions of the rotate kernel go here
 ******************************************************/

/* 
 * naive_rotate - The naive baseline version of rotate 
 */
char naive_rotate_descr[] = "naive_rotate: Naive baseline implementation";
void naive_rotate(int dim, pixel *src, pixel *dst) 
{
    int i, j;

    for (i = 0; i < dim; i++)
        for (j = 0; j < dim; j+=4)
        {
            dst[RIDX(dim-1-j, i, dim)] = *(src);
            dst[RIDX(dim-1-(j + 1), i, dim)] = *(src + 1);
            dst[RIDX(dim-1-(j + 2), i, dim)] = *(src + 2);
            dst[RIDX(dim-1-(j + 3), i, dim)] = *(src + 3);
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

char naive_rotate_old_descr[] = "naive_rotate: Naive baseline implementation";
void naive_rotate_old(int dim, pixel *src, pixel *dst) 
{
    int i, j;

    for (i = 0; i < dim; i++)
	for (j = 0; j < dim; j++)
	    dst[RIDX(dim-1-j, i, dim)] = src[RIDX(i, j, dim)];
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
    if (dim < 256)
        naive_rotate(dim, src, dst);
    else
    {
        // transpose the image in src and store
        transpose_block(dim, src, dst, 16);

        // perform row exchange on the transposed image
        xrows_ip(dim, dst);
    }
}


char rotate_two_descr[] = "second attempt";
void attempt_two(int dim, pixel *src, pixel *dst) 
{
    if (dim < 256)
        naive_rotate(dim, src, dst);
    else
    {
        // transpose the image in src and store
        transpose_block(dim, src, dst, 24);

        // perform row exchange on the transposed image
        xrows_ip(dim, dst);
    }
}

char rotate_three_descr[] = "third attempt";
void attempt_three(int dim, pixel *src, pixel *dst) 
{
    if (dim < 256)
        naive_rotate(dim, src, dst);
    else
    {
        // transpose the image in src and store
        transpose_block(dim, src, dst, 32);

        // perform row exchange on the transposed image
        xrows_ip(dim, dst);
    }
}

char rotate_four_descr[] = "fourth attempt";
void attempt_four(int dim, pixel *src, pixel *dst) 
{
    if (dim < 256)
        naive_rotate(dim, src, dst);
    else
    {
        // transpose the image in src and store
        transpose_block(dim, src, dst, 40);

        // perform row exchange on the transposed image
        xrows_ip(dim, dst);
    }
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
    //add_rotate_function(&naive_rotate, naive_rotate_descr);   
    add_rotate_function(&rotate, rotate_descr);   

    add_rotate_function(&attempt_two, rotate_two_descr);   
    add_rotate_function(&attempt_three, rotate_three_descr);   
    add_rotate_function(&attempt_four, rotate_four_descr);   
    //add_rotate_function(&attempt_five, rotate_five_descr);   
    //add_rotate_function(&attempt_six, rotate_six_descr);   
    //add_rotate_function(&attempt_seven, rotate_seven_descr);   
    //add_rotate_function(&attempt_eight, rotate_eight_descr);   
    //add_rotate_function(&attempt_nine, rotate_nine_descr);   
    //add_rotate_function(&attempt_ten, rotate_ten_descr);   
    //add_rotate_function(&attempt_eleven, rotate_eleven_descr);   

    /* ... Register additional rotate functions here */
}

