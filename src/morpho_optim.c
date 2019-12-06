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

void ui8matrix_erosion_LU3x3(uint8** ppInput, long nrl, long nrh, long ncl, long nch, p_struct_elem_dim s, uint8 **ppOutput)
{
	
	long row, col, x, y;
    // Erode
	for (row = nrl; row < nrh + 1; row++)
		for (col = ncl; col < nch + 1; col++) 
            ppOutput[row][col] = scalar_and3x3(&ppInput[row], col);
}
void ui8matrix_dilation_LU3x3(uint8** ppInput, long nrl, long nrh, long ncl, long nch, p_struct_elem_dim s, uint8 **ppOutput)
{
	long row, col, x, y;
	// dilate
	for (row = nrl; row < nrh + 1; row++)
		for (col = ncl; col < nch + 1; col++)
            ppOutput[row][col] = scalar_or3x3(&ppInput[row], col);
}

void ui8matrix_erosion_LU5x5(uint8** ppInput, long nrl, long nrh, long ncl, long nch, p_struct_elem_dim s, uint8 **ppOutput)
{
	long row, col, x, y;
    // Erode
	for (row = nrl; row < nrh + 1; row++)
		for (col = ncl; col < nch + 1; col++) 
            ppOutput[row][col] = scalar_and5x5(&ppInput[row], col);
}
void ui8matrix_dilation_LU5x5(uint8** ppInput, long nrl, long nrh, long ncl, long nch, p_struct_elem_dim s, uint8 **ppOutput)
{
	long row, col, x, y;
	// dilate
	for (row = nrl; row < nrh + 1; row++)
		for (col = ncl; col < nch + 1; col++) 
            ppOutput[row][col] = scalar_or5x5(&ppInput[row], col);
			
}

void ui8matrix_sequence_LU3x3(uint8** ppInput, long nrl, long nrh, long ncl, long nch, p_struct_elem_dim s, uint8 **ppOutput)
{
	uint8 **ppPreOutput0, **ppPreOutput1;
	ppPreOutput0 = ui8matrix(nrl + s->nrl, nrh + s->nrh, ncl + s->ncl, nch + s->nch);
	ppPreOutput1 = ui8matrix(nrl + s->nrl, nrh + s->nrh, ncl + s->ncl, nch + s->nch);
	ui8matrix_erosion_LU3x3 (ppInput     , nrl, nrh, ncl, nch, s, ppPreOutput0);
	ui8matrix_dilation_LU3x3(ppPreOutput0, nrl, nrh, ncl, nch, s, ppPreOutput1);
	ui8matrix_dilation_LU3x3(ppPreOutput1, nrl, nrh, ncl, nch, s, ppPreOutput0);
	ui8matrix_erosion_LU3x3 (ppPreOutput0, nrl, nrh, ncl, nch, s, ppOutput);
	free_ui8matrix(ppPreOutput0, nrl + s->nrl, nrh + s->nrh, ncl + s->ncl, nch + s->nch);
	free_ui8matrix(ppPreOutput1, nrl + s->nrl, nrh + s->nrh, ncl + s->ncl, nch + s->nch);
}

void ui8matrix_dilation_LU3x3_order3 (uint8** ppInput, long nrl, long nrh, long ncl, long nch, p_struct_elem_dim s, uint8 **ppOutput)
{
	const long order = 3;
	long row, col, x, y;
	long y0, y1, y2, x0, x1, x2, x3, r;
	uint8 *pRow0, *pRow1, *pRow2;

	r = (nrh + 1)  % order;
	// Erode	
	for (row = nrl; row < nrh + 1 - r; row += order) {
		pRow0 = ppInput[row + 0];
		pRow1 = ppInput[row + 1];
		pRow2 = ppInput[row + 2];
		for (col = ncl; col < nch + 1; col++) {
			x0 = scalar_or3(pRow0, col);
			x1 = scalar_or3(pRow1, col);
			x2 = scalar_or3(pRow2, col);

			y0 = scalar_or3(ppInput[row - 1], col) | x0 | x1;
			y1 = 						        x0 | x1 | x2;
			y2 = 						        x1 | x2 | scalar_or3(ppInput[row + 3], col);

			ppOutput[row + 0][col] = y0;
			ppOutput[row + 1][col] = y1;
			ppOutput[row + 2][col] = y2;
		}
	}
	// printf("%ld %ld %ld\n", r, nrh + 1 - r, nrh);
	switch(r) {
		case 2:
			for (col = ncl; col < nch + 1; col++) 
				ppOutput[row + 1][col] = scalar_or3x3(&ppInput[row + 1], col);
		case 1:
			for (col = ncl; col < nch + 1; col++) 
				ppOutput[row][col] = scalar_or3x3(&ppInput[row], col);
				
				break;
		default:
			break;
	}
	// if (nrh > 0){
	// 	for (row = nrl; row < nrh + 1; row++){
	// 		for (col = ncl; col < nch + 1; col++) 
	// 			printf("%u", ppOutput[row][col]);
	// 		putchar('\n');
	// 	}
	// 	// getchar();
	// }
}

void ui8matrix_erosion_LU3x3_order3 (uint8** ppInput, long nrl, long nrh, long ncl, long nch, p_struct_elem_dim s, uint8 **ppOutput)
{
	const long order = 3;
	long row, col, x, y;
	long y0, y1, y2, x0, x1, x2, x3, r;
	uint8 *pRow0, *pRow1, *pRow2;

	r = (nrh + 1) % order;
	// Erode	
	for (row = nrl; row < nrh + 1 - r; row += order) {
		pRow0 = ppInput[row + 0];
		pRow1 = ppInput[row + 1];
		pRow2 = ppInput[row + 2];
		for (col = ncl; col < nch + 1; col++) {
			x0 = scalar_and3(pRow0, col);
			x1 = scalar_and3(pRow1, col);
			x2 = scalar_and3(pRow2, col);

			y0 = scalar_and3(ppInput[row - 1], col) & x0 & x1;
			y1 = 							     x0 & x1 & x2;
			y2 = 							     x1 & x2 & scalar_and3(ppInput[row + 3], col);

			ppOutput[row + 0][col] = y0;
			ppOutput[row + 1][col] = y1;
			ppOutput[row + 2][col] = y2;
		}
	}
	switch(r) {
		case 2:
			for (col = ncl; col < nch + 1; col++) 
				ppOutput[row + 1][col] = scalar_and3x3(&ppInput[row + 1], col);
		case 1:
			for (col = ncl; col < nch + 1; col++) 
				ppOutput[row][col] = scalar_and3x3(&ppInput[row], col);
				break;
		default:
			break;
	}


}
void ui8matrix_dilation_LU3x3_order3_RR (uint8** ppInput, long nrl, long nrh, long ncl, long nch, p_struct_elem_dim s, uint8 **ppOutput)
{
	
	const long order = 3;
	long row, col, x, y;
	long y0, y1, y2, y3, y4, x0, x1, x2, x3, x4, r;
	uint8 *pRow0, *pRow1, *pRow2, *pRow3, *pRow4;

	r = (nrh + 1)  % order;
	// Erode	
	pRow0 = ppInput[nrl - 1];
	for (row = nrl; row < nrh + 1 - r; row+=order) {
		pRow1 = ppInput[row + 0];
		pRow2 = ppInput[row + 1];
		pRow3 = ppInput[row + 2];
		pRow4 = ppInput[row + 3];
		for (col = ncl; col < nch + 1; col++) {
			x0 = scalar_or3(pRow0, col);
			x1 = scalar_or3(pRow1, col);
			x2 = scalar_or3(pRow2, col);
			x3 = scalar_or3(pRow3, col);
			x4 = scalar_or3(pRow4, col);

			y0 = x0 | x1 | x2;
			y1 = x1 | x2 | x3;
			y2 = x2 | x3 | x4;

			ppOutput[row + 0][col] = y0;
			ppOutput[row + 1][col] = y1;
			ppOutput[row + 2][col] = y2;
		}
		pRow0 = pRow4;
		/**
		 * 3, 4, 5
		 * y0 =  x[-1][-1] + x[-1][ 0] x[-1][ 1]
		 * 		 x[ 0][-1] + x[ 0][ 0] x[ 0][ 1]
		 *       x[+1][-1] + x[+1][ 0] x[+1][ 1]
		 * 
		 * y0 =  x[ 0][-1] + x[ 0][ 0] x[ 0][ 1]
		 *       x[+1][-1] + x[+1][ 0] x[+1][ 1]
		*        x[+2][ 0] + x[+2][ 1] x[+2][ 2]
		 * **/
		// if (nrh > 0){
		// 	for (y = nrl; y < nrh + 1; y++){
		// 		for (x = ncl; x < nch + 1; x++) 
		// 			printf("%u", ppOutput[y][col]);
		// 		putchar('\n');
		// 	}
		// 	getchar();
		// }
	}
	// printf("%ld %ld %ld\n", r, nrh + 1 - r, nrh);
	switch(r) {
		case 2:
			for (col = ncl; col < nch + 1; col++) 
				ppOutput[row + 1][col] = scalar_or3x3(&ppInput[row + 1], col);
		case 1:
			for (col = ncl; col < nch + 1; col++) 
				ppOutput[row][col] = scalar_or3x3(&ppInput[row], col);
				
				break;
		default:
			break;
	}	
		// if (nrh > 0){
		// 	for (y = nrl; y < nrh + 1; y++){
		// 		for (x = ncl; x < nch + 1; x++) 
		// 			printf("%u", ppOutput[y][x]);
		// 		putchar('\n');
		// 	}
		// 	getchar();
		// }
}

void ui8matrix_erosion_LU3x3_order3_RR (uint8** ppInput, long nrl, long nrh, long ncl, long nch, p_struct_elem_dim s, uint8 **ppOutput)
{
	const long order = 3;
	long row, col, x, y;
	long y0, y1, y2, x0, x1, x2, x3, r;
	uint8 *pRow0, *pRow1, *pRow2;

	r = (nrh + 1)  % order;
	// Erode	
	for (row = nrl; row < nrh + 1 - r; row += order) {
		pRow0 = ppInput[row + 0];
		pRow1 = ppInput[row + 1];
		pRow2 = ppInput[row + 2];
		for (col = ncl; col < nch + 1; col++) {
			x0 = scalar_and3(pRow0, col);
			x1 = scalar_and3(pRow1, col);
			x2 = scalar_and3(pRow2, col);

			y0 = scalar_and3(ppInput[row - 1], col) & x0 & x1;
			y1 = 							     x0 & x1 & x2;
			y2 = 							     x1 & x2 & scalar_and3(ppInput[row + 3], col);

			ppOutput[row + 0][col] = y0;
			ppOutput[row + 1][col] = y1;
			ppOutput[row + 2][col] = y2;
		}
	}
	switch(r) {
		case 2:
			for (col = ncl; col < nch + 1; col++) 
				ppOutput[row - 2][col] = scalar_and3x3(&ppInput[row - 2], col);
		case 1:
			for (col = ncl; col < nch + 1; col++) 
				ppOutput[row - 1][col] = scalar_and3x3(&ppInput[row - 1], col);
				break;
		default:
			break;
	}

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