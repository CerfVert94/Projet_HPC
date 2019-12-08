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
void ui8matrix_dilation_LU3x3_O1xO1(uint8** ppInput, long nrl, long nrh, long ncl, long nch, p_struct_elem_dim s, uint8 **ppOutput)
{
	long row = nrl, col = ncl, x, y;
	// dilate
	for (row = nrl; row < nrh + 1; row++)
		for (col = ncl; col < nch + 1; col++)
            ppOutput[row][col] = scalar_or3x3(&ppInput[row], col);
}

void ui8matrix_dilation_LU5x5(uint8** ppInput, long nrl, long nrh, long ncl, long nch, p_struct_elem_dim s, uint8 **ppOutput)
{
	long row = nrl, col = ncl, x, y;
	// dilate
	for (row = nrl; row < nrh + 1; row++)
		for (col = ncl; col < nch + 1; col++) 
            ppOutput[row][col] = scalar_or5x5(&ppInput[row], col);
			
}


void ui8matrix_dilation_LU3x3_O3xO1 (uint8** ppInput, long nrl, long nrh, long ncl, long nch, p_struct_elem_dim s, uint8 **ppOutput)
{
	const long order = 3;
	long row = nrl, col = ncl, x, y;
	uint8 y0, y1, y2, x0, x1, x2, x3, r;
	uint8 *pRow0, *pRow1, *pRow2, *pRow3, *pRow4;

	r = (nrh + 1)  % order;
	
	for (row = nrl; row < nrh + 1 - r; row += order) {
		pRow0 = ppInput[row - 1];
		pRow1 = ppInput[row + 0];
		pRow2 = ppInput[row + 1];
		pRow3 = ppInput[row + 2];
		pRow4 = ppInput[row + 3];
		for (col = ncl; col < nch + 1; col++) {
			x0 = scalar_or3(pRow1, col);
			x1 = scalar_or3(pRow2, col);
			x2 = scalar_or3(pRow3, col);

			ppOutput[row + 0][col] = scalar_or3(pRow0, col) | x0 | x1;
			ppOutput[row + 1][col] = 				     x0 | x1 | x2;
			ppOutput[row + 2][col] = 				     x1 | x2 | scalar_or3(pRow4, col);
		}
	}
	
	switch(r) {
		case 2:
			pRow0 = ppInput[row - 1];
			pRow1 = ppInput[row + 0];
			pRow2 = ppInput[row + 1];
			pRow3 = ppInput[row + 2];
			for (col = ncl; col < nch + 1; col++) {
				x0 = scalar_or3(pRow1, col);
				x1 = scalar_or3(pRow2, col);
				ppOutput[row + 0][col] = scalar_or3(pRow0, col) | x0 | x1;
				ppOutput[row + 1][col] = 				     x0 | x1 | scalar_or3(pRow3, col);
			}
		case 1:
			pRow0 = ppInput[row - 1];
			pRow1 = ppInput[row + 0];
			pRow2 = ppInput[row + 1];
			for (col = ncl; col < nch + 1; col++) {
				ppOutput[row + 0][col] = scalar_or3(pRow0, col) | 
										 scalar_or3(pRow1, col) |
										 scalar_or3(pRow2, col);
				
			}
			break;
		default:
			break;
	}
}


void ui8matrix_dilation_LU3x3_O1xO3 (uint8** ppInput, long nrl, long nrh, long ncl, long nch, p_struct_elem_dim s, uint8 **ppOutput)
{
	const long order = 3;
	long row = nrl, col = ncl, x, y, r;
	uint8 y0, y1, y2, x0, x1, x2, x3, x4, x5, x6, x7, x8;
	uint8 *pRow0, *pRow1, *pRow2, *pRow3, *pRow4;

	r = (nrh + 1)  % order;
	
	for (row = nrl; row < nrh + 1; row ++) {
		pRow0 = ppInput[row - 1];
		pRow1 = ppInput[row + 0];
		pRow2 = ppInput[row + 1];
		for (col = ncl; col < nch + 1 - r; col += order) {
			// x0 = scalar_or3(pRow0, col); x3 = scalar_or3(pRow0, col + 1); x6 = scalar_or3(pRow0, col + 2);
			// x1 = scalar_or3(pRow1, col); x4 = scalar_or3(pRow1, col + 1); x7 = scalar_or3(pRow1, col + 2);
			// x2 = scalar_or3(pRow2, col); x5 = scalar_or3(pRow2, col + 1); x8 = scalar_or3(pRow2, col + 2);

			ppOutput[row + 0][col + 0] = scalar_or3(pRow0, col + 0)|
										 scalar_or3(pRow1, col + 0)|
										 scalar_or3(pRow2, col + 0);
			ppOutput[row + 0][col + 1] = scalar_or3(pRow0, col + 1)|
 										 scalar_or3(pRow1, col + 1)|
 										 scalar_or3(pRow2, col + 1);
			ppOutput[row + 0][col + 2] = scalar_or3(pRow0, col + 2)|
										 scalar_or3(pRow1, col + 2)|
										 scalar_or3(pRow2, col + 2);
		}
	}

	switch(r) {
		case 2:
			for (row = nrl; row < nrh + 1; row ++) {
				pRow0 = ppInput[row - 1];
				pRow1 = ppInput[row + 0];
				pRow2 = ppInput[row + 1];
				// x0 = scalar_or3(pRow0, col); x3 = scalar_or3(pRow0, col + 1); 
				// x1 = scalar_or3(pRow1, col); x4 = scalar_or3(pRow1, col + 1); 
				// x2 = scalar_or3(pRow2, col); x5 = scalar_or3(pRow2, col + 1); 

				ppOutput[row + 0][col + 0] = scalar_or3(pRow0, col + 0)|
 											 scalar_or3(pRow1, col + 0)|
 											 scalar_or3(pRow2, col + 0);
				ppOutput[row + 0][col + 1] = scalar_or3(pRow0, col + 1)|
											 scalar_or3(pRow1, col + 1)|
											 scalar_or3(pRow2, col + 1);
			}
			break;
		case 1:
			for (row = nrl; row < nrh + 1; row ++) {
				// pRow0 = ;
				// pRow1 = ;
				// pRow2 = ;
				ppOutput[row + 0][col + 0] = scalar_or3(ppInput[row - 1], col) | 
											 scalar_or3(ppInput[row + 0], col) | 
											 scalar_or3(ppInput[row + 1], col);
			}
			break;
		default:
			break;
	}
}



void ui8matrix_dilation_LU3x3_O1xO3_RR (uint8** ppInput, long nrl, long nrh, long ncl, long nch, p_struct_elem_dim s, uint8 **ppOutput)
{
	const long order = 3;
	long row = nrl, col = ncl, x, y, r;
	uint8 y0, y1, y2, x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14;
	uint8 *pRow0, *pRow1, *pRow2, *pRow3, *pRow4;

	r = (nrh + 1)  % order;
	pRow0 = ppInput[row - 1];
	pRow1 = ppInput[row + 0];
	for (row = nrl; row < nrh + 1; row ++) {
		pRow2 = ppInput[row + 1];

		x0 = pRow1[-1]; x3 = pRow1[ 0]; 
		x1 = pRow2[-1]; x4 = pRow2[ 0]; 
		for (col = ncl; col < nch + 1 - r; col += order) {
			x6 = pRow1[col + 1]; x9  = pRow1[col + 2]; x12 = pRow1[col + 3];
			x7 = pRow2[col + 1]; x10 = pRow2[col + 2]; x13 = pRow2[col + 3];

			ppOutput[row][col + 0] = scalar_or3(pRow0, col + 0) | (x0 | x3 | x6 ) | (x1 | x4  | x7);
			ppOutput[row][col + 1] = scalar_or3(pRow0, col + 1) | (x3 | x6 | x9 ) | (x4 | x7  | x10);
			ppOutput[row][col + 2] = scalar_or3(pRow0, col + 2) | (x6 | x9 | x12) | (x7 | x10 | x13);
			
			x0 = x9 ; x3 = x12;
			x1 = x10; x4 = x13;
		}
		pRow0 = pRow1;
		pRow1 = pRow2;
	}

	
	switch(r) {
	case 2:
		pRow2 = ppInput[row + 1];
		x0 = pRow1[-1]; x3 = pRow1[ 0]; x6 = pRow1[col + 1];
		x1 = pRow2[-1]; x4 = pRow2[ 0]; x7 = pRow2[col + 1];
		for (col; col < nch + 1; col ++) {
			x9  = pRow1[col + 2]; 
			x10 = pRow2[col + 2]; 
			x11 = pRow3[col + 2]; 
			ppOutput[row + 0][col + 0] = scalar_or3(pRow0, col + 0) | (x0 | x3 | x6) | (x1 | x4 | x7);
			ppOutput[row + 1][col + 0] = 		     (x0 | x3 | x6) | (x1 | x4 | x7) | (x2 | x5 | x8);

			ppOutput[row + 0][col + 1] = scalar_or3(pRow0, col + 1) | (x3 | x6 | x9 ) | (x4 | x7 | x10);
			ppOutput[row + 1][col + 1] =   		    (x3 | x6 | x9 ) | (x4 | x7 | x10) | (x5 | x8 | x11);
		
			x0 = x3; x3 = x6; x6 = x9 ;
			x1 = x4; x4 = x7; x7 = x10;
			x2 = x5; x5 = x8; x8 = x11;
		}
		break;
	case 1:
		pRow2 = ppInput[row + 1];
		x0 = pRow1[-1]; x3 = pRow1[ 0]; 
		x1 = pRow2[-1]; x4 = pRow2[ 0]; 
		for (col; col < nch + 1; col ++) {
			x6 = pRow1[col + 1];
			x7 = pRow2[col + 1];
			x8 = pRow3[col + 1];
			ppOutput[row + 0][col + 0] = scalar_or3(pRow0, col + 0) | (x0 | x3 | x6) | (x1 | x4 | x7);
			ppOutput[row + 1][col + 0] = 		     (x0 | x3 | x6) | (x1 | x4 | x7) | (x2 | x5 | x8);
			x0 = x3; x3 = x6; 
			x1 = x4; x4 = x7; 
			x2 = x5; x5 = x8; 
		}
		break;
	default:
		break;
	}
}


void ui8matrix_dilation_LU3x3_O3xO3 (uint8** ppInput, long nrl, long nrh, long ncl, long nch, p_struct_elem_dim s, uint8 **ppOutput)
{
	const long order = 3;
	long row = nrl, col = ncl, x, y;
	uint8 y0, y1, y2, rr, cr;
	uint8 x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14;
	uint8 *pRow0, *pRow1, *pRow2, *pRow3, *pRow4;

	rr = (nrh + 1) % order;
	cr = (nch + 1) % order;
	
	for (row = nrl; row < nrh + 1 - rr; row += order) {
		pRow0 = ppInput[row - 1];
		pRow1 = ppInput[row + 0];
		pRow2 = ppInput[row + 1];
		pRow3 = ppInput[row + 2];
		pRow4 = ppInput[row + 3];
		for (col = ncl; col < nch + 1 - cr; col += order) {
			x0 = scalar_or3(pRow1, col); x3 = scalar_or3(pRow1, col + 1); x6 = scalar_or3(pRow1, col + 2);
			x1 = scalar_or3(pRow2, col); x4 = scalar_or3(pRow2, col + 1); x7 = scalar_or3(pRow2, col + 2);
			x2 = scalar_or3(pRow3, col); x5 = scalar_or3(pRow3, col + 1); x8 = scalar_or3(pRow3, col + 2);

			ppOutput[row + 0][col + 0] = scalar_or3(pRow0, col + 0) | x0 | x1;
			ppOutput[row + 1][col + 0] = 				         x0 | x1 | x2;
			ppOutput[row + 2][col + 0] = 				         x1 | x2 | scalar_or3(pRow4, col + 0);


			ppOutput[row + 0][col + 1] = scalar_or3(pRow0, col + 1) | x3 | x4;
			ppOutput[row + 1][col + 1] = 				         x3 | x4 | x5;
			ppOutput[row + 2][col + 1] = 				         x4 | x5 | scalar_or3(pRow4, col + 1);


			ppOutput[row + 0][col + 2] = scalar_or3(pRow0, col + 2) | x6 | x7;
			ppOutput[row + 1][col + 2] = 				         x6 | x7 | x8;
			ppOutput[row + 2][col + 2] = 				         x7 | x8 | scalar_or3(pRow4, col + 2);

		}
	}
	pRow0 = ppInput[row - 1];
	pRow1 = ppInput[row + 0];
	pRow2 = ppInput[row + 1];
	pRow3 = ppInput[row + 2];
	switch (rr) {
		case 2 :
			switch(cr){
				case 2 :
					for (col; col < nch + 1; col ++) {
						x0 = scalar_or3(pRow1, col); x3 = scalar_or3(pRow1, col + 1); 
						x1 = scalar_or3(pRow2, col); x4 = scalar_or3(pRow2, col + 1); 

						ppOutput[row + 0][col + 0] = scalar_or3(pRow0, col + 0) | x0 | x1;
						ppOutput[row + 1][col + 0] = 				         x0 | x1 | scalar_or3(pRow3, col);


						ppOutput[row + 0][col + 1] = scalar_or3(pRow0, col + 1) | x3 | x4;
						ppOutput[row + 1][col + 1] = 				         x3 | x4 | scalar_or3(pRow3, col + 1);
					}
				break;
				case 1: 
					for (col; col < nch + 1; col ++) {
						x0 = scalar_or3(pRow1, col); 
						x1 = scalar_or3(pRow2, col); 

						ppOutput[row + 0][col + 0] = scalar_or3(pRow0, col + 0) | x0 | x1;
						ppOutput[row + 1][col + 0] = 				         x0 | x1 | scalar_or3(pRow3, col);
					}
				break;
				default :
				break;
			}
		case 1:
			switch(cr) {
				case 2:
					for (col; col < nch + 1; col ++) {
						ppOutput[row + 0][col + 0] = scalar_or3(pRow0, col + 0) | scalar_or3(pRow1, col + 0) | scalar_or3(pRow2, col + 0);
						ppOutput[row + 0][col + 1] = scalar_or3(pRow0, col + 1) | scalar_or3(pRow1, col + 1) | scalar_or3(pRow2, col + 1);
					}
				break;
				case 1:
					for (col; col < nch + 1; col ++) {
						ppOutput[row + 0][col + 0] = scalar_or3(pRow0, col + 0) | 
												 	 scalar_or3(pRow1, col + 0) | 
													 scalar_or3(pRow2, col + 0);
					}
				break;
			}
		break;
		default:
		break;
	}
	// ui8matrix_dilation_LU3x3(ppInput, row, nrh, col, nrh, s, ppOutput);
}

void ui8matrix_dilation_LU3x3_O3xO3_RR(uint8** ppInput, long nrl, long nrh, long ncl, long nch, p_struct_elem_dim s, uint8 **ppOutput)
{
	const long order = 3;
	long row = nrl, col = ncl, x, y;
	uint8 y0, y1, y2, rr, cr;
	uint8 x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14;
	uint8 *pRow0, *pRow1, *pRow2, *pRow3, *pRow4;

	rr = (nrh + 1) % order;
	cr = (nch + 1) % order;
	
	pRow0 = ppInput[row - 1];
	pRow1 = ppInput[row + 0];
	for (row = nrl; row < nrh + 1 - rr; row += order) {
		pRow2 = ppInput[row + 1];
		pRow3 = ppInput[row + 2];
		pRow4 = ppInput[row + 3];
		
		x0 = pRow1[-1]; x3 = pRow1[ 0]; 
		x1 = pRow2[-1]; x4 = pRow2[ 0]; 
		x2 = pRow3[-1]; x5 = pRow3[ 0]; 
		for (col = ncl; col < nch + 1 - cr; col += order) {
			x6 = pRow1[col + 1]; x9  = pRow1[col + 2]; x12 = pRow1[col + 3];
			x7 = pRow2[col + 1]; x10 = pRow2[col + 2]; x13 = pRow2[col + 3];
			x8 = pRow3[col + 1]; x11 = pRow3[col + 2]; x14 = pRow3[col + 3];
			ppOutput[row + 0][col + 0] = scalar_or3(pRow0, col + 0) | (x0 | x3 | x6) | (x1 | x4 | x7);
			ppOutput[row + 1][col + 0] = 		     (x0 | x3 | x6) | (x1 | x4 | x7) | (x2 | x5 | x8);
			ppOutput[row + 2][col + 0] = 		     (x1 | x4 | x7) | (x2 | x5 | x8) | scalar_or3(pRow4, col + 0);

			
			
			
			ppOutput[row + 0][col + 1] = scalar_or3(pRow0, col + 1) | (x3 | x6 | x9 ) | (x4 | x7 | x10);
			ppOutput[row + 1][col + 1] =   		    (x3 | x6 | x9 ) | (x4 | x7 | x10) | (x5 | x8 | x11);
			ppOutput[row + 2][col + 1] =   		    (x4 | x7 | x10) | (x5 | x8 | x11) | scalar_or3(pRow4, col + 1);
			
			
			
			ppOutput[row + 0][col + 2] = scalar_or3(pRow0, col + 2) | (x6 | x9  | x12) | (x7 | x10 | x13);
			ppOutput[row + 1][col + 2] = 		   (x6 | x9  | x12) | (x7 | x10 | x13) | (x8 | x11 | x14);
			ppOutput[row + 2][col + 2] = 		   (x7 | x10 | x13) | (x8 | x11 | x14) | scalar_or3(pRow4, col + 2);
			
			// printf("%u %u %u %u %u\n", x0, x3, x6, x9 , x12);
			// printf("%u %u %u %u %u\n", x1, x4, x7, x10, x13);
			// printf("%u %u %u %u %u\n", x2, x5, x8, x11, x14);	
			// display_ui8matrix(ppOutput, nrl, nrh, ncl, nch, "%u", "debug");
			// getchar();
			x0 = x9 ; x3 = x12;
			x1 = x10; x4 = x13;
			x2 = x11; x5 = x14;
		}
		pRow0 = pRow3;
		pRow1 = pRow4;
	}
	switch(rr) {
	case 2: 
		pRow2 = ppInput[row + 1];
		pRow3 = ppInput[row + 2];
		switch(cr) {
		case 2:
			x0 = pRow1[-1]; x3 = pRow1[ 0]; x6 = pRow1[col + 1];
			x1 = pRow2[-1]; x4 = pRow2[ 0]; x7 = pRow2[col + 1];
			x2 = pRow3[-1]; x5 = pRow3[ 0]; x8 = pRow3[col + 1];
			for (col; col < nch + 1; col ++) {
				x9  = pRow1[col + 2]; 
				x10 = pRow2[col + 2]; 
				x11 = pRow3[col + 2]; 
				ppOutput[row + 0][col + 0] = scalar_or3(pRow0, col + 0) | (x0 | x3 | x6) | (x1 | x4 | x7);
				ppOutput[row + 1][col + 0] = 		     (x0 | x3 | x6) | (x1 | x4 | x7) | (x2 | x5 | x8);

				ppOutput[row + 0][col + 1] = scalar_or3(pRow0, col + 1) | (x3 | x6 | x9 ) | (x4 | x7 | x10);
				ppOutput[row + 1][col + 1] =   		    (x3 | x6 | x9 ) | (x4 | x7 | x10) | (x5 | x8 | x11);
			
				x0 = x3; x3 = x6; x6 = x9 ;
				x1 = x4; x4 = x7; x7 = x10;
				x2 = x5; x5 = x8; x8 = x11;
			}
			break;
		case 1:
			x0 = pRow1[-1]; x3 = pRow1[ 0]; 
			x1 = pRow2[-1]; x4 = pRow2[ 0]; 
			x2 = pRow3[-1]; x5 = pRow3[ 0]; 
			for (col; col < nch + 1; col ++) {
				x6 = pRow1[col + 1];
				x7 = pRow2[col + 1];
				x8 = pRow3[col + 1];
				ppOutput[row + 0][col + 0] = scalar_or3(pRow0, col + 0) | (x0 | x3 | x6) | (x1 | x4 | x7);
				ppOutput[row + 1][col + 0] = 		     (x0 | x3 | x6) | (x1 | x4 | x7) | (x2 | x5 | x8);
				x0 = x3; x3 = x6; 
				x1 = x4; x4 = x7; 
				x2 = x5; x5 = x8; 
			}
			break;
		default:
			break;
		}
		break;
	case 1: 
		pRow2 = ppInput[row + 1];
		switch(cr) {
		case 2: 
			x0 = pRow1[-1]; x3 = pRow1[ 0]; x6 = pRow1[col + 1];
			x1 = pRow2[-1]; x4 = pRow2[ 0]; x7 = pRow2[col + 1];
			for (col; col < nch + 1; col ++) {
				x9  = pRow1[col + 2]; 
				x10 = pRow2[col + 2]; 
				ppOutput[row + 0][col + 0] = scalar_or3(pRow0, col + 0) | (x0 | x3 | x6) | (x1 | x4 | x7);
				ppOutput[row + 0][col + 1] = scalar_or3(pRow0, col + 1) | (x3 | x6 | x9) | (x4 | x7 | x10);
			
				x0 = x3; x3 = x6; x6 = x9 ;
				x1 = x4; x4 = x7; x7 = x10;
				x2 = x5; x5 = x8; x8 = x11;
			}
			break;
		case 1:
			x0 = pRow1[-1]; x3 = pRow1[ 0]; 
			x1 = pRow2[-1]; x4 = pRow2[ 0]; 
			for (col; col < nch + 1; col ++) {
				x6 = pRow1[col + 1];
				x7 = pRow2[col + 1];
				ppOutput[row + 0][col + 0] = scalar_or3(pRow0, col + 0) | (x0 | x3 | x6) | (x1 | x4 | x7);
				x0 = x3; x3 = x6; 
				x1 = x4; x4 = x7; 
			}
			break;
		default:
			break;
		}
		break;
	default:
		break;
	}
		
}


void ui8matrix_dilation_LU3x3_O3xO1_RR(uint8** ppInput, long nrl, long nrh, long ncl, long nch, p_struct_elem_dim s, uint8 **ppOutput)
{
	const long order = 3;
	long row = nrl, col = ncl, x, y;
	uint8 y0, y1, y2, r;
	uint8 x0, x1, x2, x3, x4, x5, x6, x7, x8;
	uint8 *pRow0, *pRow1, *pRow2, *pRow3, *pRow4;

	r = (nrh + 1)  % order;
	
	pRow0 = ppInput[row - 1];
	pRow1 = ppInput[row + 0];
	for (row = nrl; row < nrh + 1 - r; row += order) {
		pRow2 = ppInput[row + 1]; 
		pRow3 = ppInput[row + 2]; 
		pRow4 = ppInput[row + 3];
		x0 = pRow1[-1]; x3 = pRow1[ 0];
		x1 = pRow2[-1]; x4 = pRow2[ 0];
		x2 = pRow3[-1]; x5 = pRow3[ 0];
		for (col = ncl; col < nch + 1; col++) {
			x6 = pRow1[col + 1]; //scalar_or3(pRow1, col);
			x7 = pRow2[col + 1]; //scalar_or3(pRow2, col);
			x8 = pRow3[col + 1]; //scalar_or3(pRow3, col);
			// y0 = 
			// y1 = 
			// y2 = 

			ppOutput[row + 0][col] = scalar_or3(pRow0, col) | (x0 | x3 | x6) | (x1 | x4 | x7);
			ppOutput[row + 1][col] =         (x0 | x3 | x6) | (x1 | x4 | x7) | (x2 | x5 | x8);
			ppOutput[row + 2][col] = 		 (x1 | x4 | x7) | (x2 | x5 | x8) | scalar_or3(pRow4, col);
			x0 = x3; x3 = x6;
			x1 = x4; x4 = x7;
			x2 = x5; x5 = x8;
		}
		pRow0 = pRow3;
		pRow1 = pRow4;
	}
	pRow2 = ppInput[row + 1];
	pRow3 = ppInput[row + 2];
	pRow4 = ppInput[row + 3];
	switch(r) {
		case 2:
			x0 = pRow1[-1]; x3 = pRow1[ 0];
			x1 = pRow2[-1]; x4 = pRow2[ 0];
			x2 = pRow3[-1]; x5 = pRow3[ 0];
			for (col = ncl; col < nch + 1; col++) {
				x6 = pRow1[col + 1]; //scalar_or3(pRow1, col);
				x7 = pRow2[col + 1]; //scalar_or3(pRow2, col);
				x8 = pRow3[col + 1]; //scalar_or3(pRow3, col);
				ppOutput[row + 0][col] = scalar_or3(pRow0, col) | (x0 | x3 | x6) | (x1 | x4 | x7);;
				ppOutput[row + 1][col] =         (x0 | x3 | x6) | (x1 | x4 | x7) | (x2 | x5 | x8);;
				x0 = x3; x3 = x6;
				x1 = x4; x4 = x7;
				x2 = x5; x5 = x8;
			}
		case 1:
			x0 = pRow1[-1]; x3 = pRow1[ 0];
			x1 = pRow2[-1]; x4 = pRow2[ 0];
			for (col = ncl; col < nch + 1; col++) {
				x6 = pRow1[col + 1]; //scalar_or3(pRow1, col);
				x7 = pRow2[col + 1]; //scalar_or3(pRow2, col);
				ppOutput[row + 0][col] = scalar_or3(pRow0, col) | (x0 | x3 | x6) | (x1 | x4 | x7);
				x0 = x3; x3 = x6;
				x1 = x4; x4 = x7;
			}				
			break;
		default:
			break;
	}
}



// void ui8matrix_erosion_LU3x3(uint8** ppInput, long nrl, long nrh, long ncl, long nch, p_struct_elem_dim s, uint8 **ppOutput)
// {
	
// 	long row = nrl, col = ncl, x, y;
//     // Erode
// 	for (row = nrl; row < nrh + 1; row++)
// 		for (col = ncl; col < nch + 1; col++) 
//             ppOutput[row][col] = scalar_and3x3(&ppInput[row], col);
// }
// void ui8matrix_erosion_LU5x5(uint8** ppInput, long nrl, long nrh, long ncl, long nch, p_struct_elem_dim s, uint8 **ppOutput)
// {
// 	long row = nrl, col = ncl, x, y;
//     // Erode
// 	for (row = nrl; row < nrh + 1; row++)
// 		for (col = ncl; col < nch + 1; col++) 
//             ppOutput[row][col] = scalar_and5x5(&ppInput[row], col);
// }
// void ui8matrix_sequence_LU3x3(uint8** ppInput, long nrl, long nrh, long ncl, long nch, p_struct_elem_dim s, uint8 **ppOutput)
// {
// 	uint8 **ppPreOutput0, **ppPreOutput1;
// 	ppPreOutput0 = ui8matrix(nrl + s->nrl, nrh + s->nrh, ncl + s->ncl, nch + s->nch);
// 	ppPreOutput1 = ui8matrix(nrl + s->nrl, nrh + s->nrh, ncl + s->ncl, nch + s->nch);
// 	ui8matrix_erosion_LU3x3 (ppInput     , nrl, nrh, ncl, nch, s, ppPreOutput0);
// 	ui8matrix_dilation_LU3x3_O1xO1(ppPreOutput0, nrl, nrh, ncl, nch, s, ppPreOutput1);
// 	ui8matrix_dilation_LU3x3_O1xO1(ppPreOutput1, nrl, nrh, ncl, nch, s, ppPreOutput0);
// 	ui8matrix_erosion_LU3x3 (ppPreOutput0, nrl, nrh, ncl, nch, s, ppOutput);
// 	free_ui8matrix(ppPreOutput0, nrl + s->nrl, nrh + s->nrh, ncl + s->ncl, nch + s->nch);
// 	free_ui8matrix(ppPreOutput1, nrl + s->nrl, nrh + s->nrh, ncl + s->ncl, nch + s->nch);
// }

