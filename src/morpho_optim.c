#include "nrdef.h"
#include "nrutil.h"
#include <img.h>
#include <stdlib.h>
#include <stdio.h>
//#include <malloc.h>
#include <morpho.h>
#include <test_morpho.h>
#include <util.h>

// #define scalar_and3(input, col) input[col - 2] & input[col - 1] & input[col + 0] & input[col + 1] & input[col + 2]
// #define scalar_and5(input, col)                  input[col - 1] & input[col + 0] & input[col + 1]
// #define scalar_or3 (input, col) input[col - 2] | input[col - 1] | input[col + 0] | input[col + 1] | input[col + 2]
// #define scalar_or5 (input, col)                  input[col - 1] | input[col + 0] | input[col + 1]

// #define scalar_and3x3(input, col) scalar_and3(input[-1], col) &\
//                                   scalar_and3(input[ 0], col) &\
//                                   scalar_and3(input[ 1], col)
// #define scalar_and5x5(input, col) scalar_and3(input[-2], col) &\
//                                   scalar_and3(input[-1], col) &\
//                                   scalar_and3(input[ 0], col) &\
//                                   scalar_and3(input[ 1], col) &\
//                                   scalar_and3(input[ 2], col)

// #define scalar_or3x3(input, col)  scalar_or3(input[-1], col) |\
//                                   scalar_or3(input[ 0], col) |\
//                                   scalar_or3(input[ 1], col)
// #define scalar_or5x5(input, col)  scalar_or3(input[-2], col) |\
//                                   scalar_or3(input[-1], col) |\
//                                   scalar_or3(input[ 0], col) |\
//                                   scalar_or3(input[ 1], col) |\
//                                   scalar_or3(input[ 2], col)

#define scalar_and3(input)            (input[-1] & input[0] & input[1])
#define scalar_and5(input) (input[-2] & input[-1] & input[0] & input[1] & input[2])

#define scalar_or3(input)              (input[-1] | input[0] | input[1])
#define scalar_or5(input)  (input[-2] | input[-1] | input[0] | input[1] | input[2])

#define scalar_and3x3(input) scalar_and3(input[-1]) &\
                             scalar_and3(input[ 0]) &\
                             scalar_and3(input[ 1])
#define scalar_and5x5(input) scalar_and5(input[-2]) &\
                             scalar_and5(input[-1]) &\
                             scalar_and5(input[ 0]) &\
                             scalar_and5(input[ 1]) &\
                             scalar_and5(input[ 2])

#define scalar_or3x3(input)  scalar_or3(input[-1]) |\
                             scalar_or3(input[ 0]) |\
                             scalar_or3(input[ 1])
#define scalar_or5x5(input)  scalar_or5(input[-2]) |\
                             scalar_or5(input[-1]) |\
                             scalar_or5(input[ 0]) |\
                             scalar_or5(input[ 1]) |\
                             scalar_or5(input[ 2])
void ui8matrix_erosion_3x3_noloop(uint8** ppInput, long nrl, long nrh, long ncl, long nch, p_struct_elem_dim s, uint8 **ppOutput)
{
    const long order = 3;
	long row, col, x, y;
    // Erode
	for (row = nrl; row < nrh + 1; row++)
		for (col = ncl; col < nch + 1; col++) 
            ppOutput[row][col] = scalar_and3x3((&ppInput[row]));
}
void ui8matrix_dilation_3x3_noloop(uint8** ppInput, long nrl, long nrh, long ncl, long nch, p_struct_elem_dim s, uint8 **ppOutput)
{
	long row, col, x, y;
	// dilate
	for (row = nrl; row < nrh + 1; row++)
		for (col = ncl; col < nch + 1; col++)
            ppOutput[row][col] = scalar_or3x3((&ppInput[row]));
}

void ui8matrix_erosion_5x5_noloop(uint8** ppInput, long nrl, long nrh, long ncl, long nch, p_struct_elem_dim s, uint8 **ppOutput)
{
    const long order = 3;
	long row, col, x, y;
    // Erode
	for (row = nrl; row < nrh + 1; row++)
		for (col = ncl; col < nch + 1; col++) 
            ppOutput[row][col] = scalar_and5x5((&ppInput[row]));
}
void ui8matrix_dilation_5x5_noloop(uint8** ppInput, long nrl, long nrh, long ncl, long nch, p_struct_elem_dim s, uint8 **ppOutput)
{
	long row, col, x, y;
	// dilate
	for (row = nrl; row < nrh + 1; row++)
		for (col = ncl; col < nch + 1; col++) 
            ppOutput[row][col] = scalar_or5x5((&ppInput[row]));
}


inline void ui8matrix_sequence_naive_inline(uint8** ppInput, long nrl, long nrh, long ncl, long nch, p_struct_elem_dim s, uint8 **ppOutput)
{
	long row, col, x, y;
	long snrl = s->nrl, snrh = s->nrh, sncl = s->ncl, snch = s->nch;
	uint8 **ppPreOutput0, **ppPreOutput1;

	ppPreOutput0 = ui8matrix(nrl + snrl, nrh + snrh, ncl + sncl, nch + snch);
    ppPreOutput1 = ui8matrix(nrl + snrl, nrh + snrh, ncl + sncl, nch + snch);

	for (row = nrl; row < nrh + 1; row++){
		for (col = ncl; col < nch + 1; col++) {
			// Copy the input to the output
			ppPreOutput0[row][col] = ppInput[row][col];

			// Apply morpho to the output
			for (y = snrl; y < snrh + 1; y++ ) 
				for (x = sncl; x < snch + 1; x++ ) 
					ppPreOutput0[row][col] &= ppInput[row + y][col + x];
		}
	}


	for (row = nrl; row < nrh + 1; row++){
		for (col = ncl; col < nch + 1; col++) {
			// Copy the input to the output
			ppPreOutput1[row][col] = ppPreOutput0[row][col];
			
			// Apply morpho to the output
			for (y = snrl; y < snrh + 1; y++ ) 
				for (x = sncl; x < snch + 1; x++ ) 
					ppPreOutput1[row][col] |= ppPreOutput0[row + y][col + x];
		}
	}
	
	for (row = nrl; row < nrh + 1; row++){
		for (col = ncl; col < nch + 1; col++) {
			// Copy the input to the output
			ppPreOutput0[row][col] = ppPreOutput1[row][col];
			
			// Apply morpho to the output
			for (y = snrl; y < snrh + 1; y++ ) 
				for (x = sncl; x < snch + 1; x++ ) 
					ppPreOutput0[row][col] |= ppPreOutput1[row + y][col + x];
		}
	}

	for (row = nrl; row < nrh + 1; row++){
		for (col = ncl; col < nch + 1; col++) {
			// Copy the input to the output
			ppOutput[row][col] = ppPreOutput0[row][col];

			// Apply morpho to the output
			for (y = snrl; y < snrh + 1; y++ ) 
				for (x = sncl; x < snch + 1; x++ ) 
					ppOutput[row][col] &= ppPreOutput0[row + y][col + x];
		}
	}
	free_ui8matrix(ppPreOutput0, nrl + snrl, nrh + snrh, ncl + sncl, nch + snch);
	free_ui8matrix(ppPreOutput1, nrl + snrl, nrh + snrh, ncl + sncl, nch + snch);

}