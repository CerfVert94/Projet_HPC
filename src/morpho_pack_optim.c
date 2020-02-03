#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include "nrdef.h"
#include "nrutil.h"
#include <util.h>
#include <img.h>
#include "vnrdef.h"
#include "vnrutil.h"
#include <img_SIMD.h>
#include <morpho.h>

#define scalar_and3(input, col)                  (input[col - 1] & input[col + 0] & input[col + 1])
#define scalar_and5(input, col) (input[col - 2] & input[col - 1] & input[col + 0] & input[col + 1] & input[col + 2])
#define scalar_or3(input, col)                   (input[col - 1] | input[col + 0] | input[col + 1])
#define scalar_or5(input, col)  (input[col - 2] | input[col - 1] | input[col + 0] | input[col + 1] | input[col + 2])

#define scalar_and3x3(input, col) scalar_and3((input)[-1], col) &\
                                  scalar_and3((input)[ 0], col) &\
                                  scalar_and3((input)[ 1], col) 

#define scalar_and5x5(input, col) scalar_and5((input)[-2], col) &\
                                  scalar_and5((input)[-1], col) &\
                                  scalar_and5((input)[ 0], col) &\
                                  scalar_and5((input)[ 1], col) &\
                                  scalar_and5((input)[ 2], col)

#define scalar_or3x3(input, col)  scalar_or3((input)[-1], col) |\
                                  scalar_or3((input)[ 0], col) |\
                                  scalar_or3((input)[ 1], col)

#define scalar_or5x5(input, col)  scalar_or5((input)[-2], col) |\
                                  scalar_or5((input)[-1], col) |\
                                  scalar_or5((input)[ 0], col) |\
                                  scalar_or5((input)[ 1], col) |\
                                  scalar_or5((input)[ 2], col)

#include "nrdef.h"
#include "nrutil.h"
#include <img.h>
#include <stdlib.h>
#include <stdio.h>
//#include <malloc.h>
#include <morpho.h>
// #include <test_morpho.h>
#include <util.h>

#define scalar_and3(input, col)                  (input[col - 1] & input[col + 0] & input[col + 1])
#define scalar_and5(input, col) (input[col - 2] & input[col - 1] & input[col + 0] & input[col + 1] & input[col + 2])
#define scalar_or3(input, col)                   (input[col - 1] | input[col + 0] | input[col + 1])
#define scalar_or5(input, col)  (input[col - 2] | input[col - 1] | input[col + 0] | input[col + 1] | input[col + 2])

#define scalar_and3x3(input, col) scalar_and3((input)[-1], col) &\
                                  scalar_and3((input)[ 0], col) &\
                                  scalar_and3((input)[ 1], col) 

#define scalar_and5x5(input, col) scalar_and5((input)[-2], col) &\
                                  scalar_and5((input)[-1], col) &\
                                  scalar_and5((input)[ 0], col) &\
                                  scalar_and5((input)[ 1], col) &\
                                  scalar_and5((input)[ 2], col)

#define scalar_or3x3(input, col)  scalar_or3((input)[-1], col) |\
                                  scalar_or3((input)[ 0], col) |\
                                  scalar_or3((input)[ 1], col)

#define scalar_or5x5(input, col)  scalar_or5((input)[-2], col) |\
                                  scalar_or5((input)[-1], col) |\
                                  scalar_or5((input)[ 0], col) |\
                                  scalar_or5((input)[ 1], col) |\
                                  scalar_or5((input)[ 2], col)

#define SE_NRL -1
#define SE_NRH  1
#define SE_NCL -1
#define SE_NCH  1


/*******************************************/
/******* Optimisation : Loop Unroll ********/
/*******************************************/

void ui8matrix_dilation_hpacked_divide_row_and_conquer (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **temp_buffer, uint8 **Y)
{
	const long order = 3;
	
	uint8 *temp_row, *temp_row0, *temp_row1, *temp_row2, *out_row, *in_row;
	long row, col, col_prime;
	long packed_ncl, packed_nch;
	int nb_threads;
	uint8 x0, x1, x2;
 
	
	for (row = nrl - 1; row < (nrh + 1) + 1; row ++)
	{
		temp_row = temp_buffer[row];
		in_row   = X[row];

		for(col = ncl; col < nch + 1; col++) {
			x1 = in_row[col + 0];
			x0 = x1 << 1 | (in_row[col + 1] >> 7) & 0x1; 
			x2 = x1 >> 1 | (in_row[col - 1] & 0x1) << 7;
			temp_row[col] = x0 | x1 | x2;
			
		}
	}
	
	for (row = nrl; row < nrh + 1; row ++)
	{
		temp_row0 =     temp_buffer[row - 1];
		temp_row1 =     temp_buffer[row + 0];
		temp_row2 =     temp_buffer[row + 1];
		out_row   =     Y[row];
		for(col = ncl; col < nch + 1; col++) {
			out_row[col] = temp_row0[col]|
						   temp_row1[col]|
						   temp_row2[col];
		}
	}
}
