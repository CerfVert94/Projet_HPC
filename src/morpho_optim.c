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
	uint8 *row0, *row1, *row2, *row3, *row4;

	r = (nrh + 1)  % order;
	
	for (row = nrl; row < nrh + 1 - r; row += order) {
		row0 = ppInput[row - 1];
		row1 = ppInput[row + 0];
		row2 = ppInput[row + 1];
		row3 = ppInput[row + 2];
		row4 = ppInput[row + 3];
		for (col = ncl; col < nch + 1; col++) {
			x0 = scalar_or3(row1, col);
			x1 = scalar_or3(row2, col);
			x2 = scalar_or3(row3, col);

			ppOutput[row + 0][col] = scalar_or3(row0, col) | x0 | x1;
			ppOutput[row + 1][col] = 				     x0 | x1 | x2;
			ppOutput[row + 2][col] = 				     x1 | x2 | scalar_or3(row4, col);
		}
	}
	
	switch(r) {
		case 2:
			row0 = ppInput[row - 1];
			row1 = ppInput[row + 0];
			row2 = ppInput[row + 1];
			row3 = ppInput[row + 2];
			for (col = ncl; col < nch + 1; col++) {
				x0 = scalar_or3(row1, col);
				x1 = scalar_or3(row2, col);
				ppOutput[row + 0][col] = scalar_or3(row0, col) | x0 | x1;
				ppOutput[row + 1][col] = 				     x0 | x1 | scalar_or3(row3, col);
			}
		case 1:
			row0 = ppInput[row - 1];
			row1 = ppInput[row + 0];
			row2 = ppInput[row + 1];
			for (col = ncl; col < nch + 1; col++) {
				ppOutput[row + 0][col] = scalar_or3(row0, col) | 
										 scalar_or3(row1, col) |
										 scalar_or3(row2, col);
				
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
	uint8 *row0, *row1, *row2, *row3, *row4;

	r = (nch + 1)  % order;
	
	for (row = nrl; row < nrh + 1; row ++) {
		row0 = ppInput[row - 1];
		row1 = ppInput[row + 0];
		row2 = ppInput[row + 1];
		for (col = ncl; col < nch + 1 - r; col += order) {
			// x0 = scalar_or3(row0, col); x3 = scalar_or3(row0, col + 1); x6 = scalar_or3(row0, col + 2);
			// x1 = scalar_or3(row1, col); x4 = scalar_or3(row1, col + 1); x7 = scalar_or3(row1, col + 2);
			// x2 = scalar_or3(row2, col); x5 = scalar_or3(row2, col + 1); x8 = scalar_or3(row2, col + 2);

			ppOutput[row + 0][col + 0] = scalar_or3(row0, col + 0)|
										 scalar_or3(row1, col + 0)|
										 scalar_or3(row2, col + 0);
			ppOutput[row + 0][col + 1] = scalar_or3(row0, col + 1)|
 										 scalar_or3(row1, col + 1)|
 										 scalar_or3(row2, col + 1);
			ppOutput[row + 0][col + 2] = scalar_or3(row0, col + 2)|
										 scalar_or3(row1, col + 2)|
										 scalar_or3(row2, col + 2);
		}
	}

	switch(r) {
		case 2:
			for (row = nrl; row < nrh + 1; row ++) {
				row0 = ppInput[row - 1];
				row1 = ppInput[row + 0];
				row2 = ppInput[row + 1];

				ppOutput[row + 0][col + 0] = scalar_or3(row0, col + 0)|
 											 scalar_or3(row1, col + 0)|
 											 scalar_or3(row2, col + 0);
				ppOutput[row + 0][col + 1] = scalar_or3(row0, col + 1)|
											 scalar_or3(row1, col + 1)|
											 scalar_or3(row2, col + 1);
			}
			break;
		case 1:
			for (row = nrl; row < nrh + 1; row ++) {
				ppOutput[row + 0][col + 0] = scalar_or3(ppInput[row - 1], col) | 
											 scalar_or3(ppInput[row + 0], col) | 
											 scalar_or3(ppInput[row + 1], col);
			}
			break;
		default:
			break;
	}
}

void ui8matrix_dilation_LU3x3_O3xO3_1 (uint8** ppInput, long nrl, long nrh, long ncl, long nch, p_struct_elem_dim s, uint8 **ppOutput)
{
	const long order = 3;
	long row = nrl, col = ncl, x, y;
	uint8 y0, y1, y2, rr, cr;
	uint8 x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14;
	uint8 *row0, *row1, *row2, *row3, *row4;
	uint8 *out_row0, *out_row1, *out_row2;

	rr = (nrh + 1) % order;
	cr = (nch + 1) % order;
	
	for (row = nrl; row < nrh + 1 - rr; row += order) {
		row0 = ppInput[row - 1];
		row1 = ppInput[row + 0];
		row2 = ppInput[row + 1];
		row3 = ppInput[row + 2];
		row4 = ppInput[row + 3];
		
		for (col = ncl; col < nch + 1 - cr; col += order) {
			x0 = scalar_or3(row1, col); x1 = scalar_or3(row1, col + 1); x2 = scalar_or3(row1, col + 2);
			x3 = scalar_or3(row2, col); x4 = scalar_or3(row2, col + 1); x5 = scalar_or3(row2, col + 2);
			x6 = scalar_or3(row3, col); x7 = scalar_or3(row3, col + 1); x8 = scalar_or3(row3, col + 2);

			ppOutput[row + 0][col + 0] = scalar_or3(row0, col + 0) | x0 | x3;
			ppOutput[row + 1][col + 0] = 				         x0 | x3 | x6;
			ppOutput[row + 2][col + 0] = 				         x3 | x6 | scalar_or3(row4, col + 0);


			ppOutput[row + 0][col + 1] = scalar_or3(row0, col + 1) | x1 | x4;
			ppOutput[row + 1][col + 1] = 				         x1 | x4 | x7;
			ppOutput[row + 2][col + 1] = 				         x4 | x7 | scalar_or3(row4, col + 1);


			ppOutput[row + 0][col + 2] = scalar_or3(row0, col + 2) | x2 | x5;
			ppOutput[row + 1][col + 2] = 				         x2 | x5 | x8;
			ppOutput[row + 2][col + 2] = 				         x5 | x8 | scalar_or3(row4, col + 2);

		}
	}
	switch (rr) {
		case 2 :
			switch(cr){
				case 2 :
					ui8matrix_dilation_LU3x3_O3xO1(ppInput, nrl, nrh, nch - 1, nch, s, ppOutput);
					ui8matrix_dilation_LU3x3_O1xO3(ppInput, nrh - 1, nrh, ncl, nch, s, ppOutput);
				break;
				case 1: 
					ui8matrix_dilation_LU3x3_O3xO1(ppInput, nrl, nrh, nch, nch, s, ppOutput);		
					ui8matrix_dilation_LU3x3_O1xO3(ppInput, nrh - 1, nrh, ncl, nch, s, ppOutput);	
				break;
				default :
					ui8matrix_dilation_LU3x3_O1xO3(ppInput, nrh - 1, nrh, ncl, nch, s, ppOutput);
				break;
			}
		case 1:
			switch(cr){
				case 2 :
					ui8matrix_dilation_LU3x3_O3xO1(ppInput, nrl, nrh, nch - 1, nch, s, ppOutput);
					ui8matrix_dilation_LU3x3_O1xO3(ppInput, nrh, nrh, ncl, nch, s, ppOutput);
				break;
				case 1: 
					ui8matrix_dilation_LU3x3_O3xO1(ppInput, nrl, nrh, nch, nch, s, ppOutput);		
					ui8matrix_dilation_LU3x3_O1xO3(ppInput, nrh, nrh, ncl, nch, s, ppOutput);	
				break;
				default :
					ui8matrix_dilation_LU3x3_O1xO3(ppInput, nrh, nrh, ncl, nch, s, ppOutput);
				break;
			}
		break;
		default:
			switch(cr){
				case 2 :
					ui8matrix_dilation_LU3x3_O3xO1(ppInput, nrl, nrh, nch - 1, nch, s, ppOutput);
				break;
				case 1: 
					ui8matrix_dilation_LU3x3_O3xO1(ppInput, nrl, nrh, nch, nch, s, ppOutput);	
				break;
				default :
				break;
			}
		break;
	}
}
void ui8matrix_dilation_LU3x3_O3xO3 (uint8** ppInput, long nrl, long nrh, long ncl, long nch, p_struct_elem_dim s, uint8 **ppOutput)
{
	const long order = 3;
	long row = nrl, col = ncl, x, y;
	uint8 y0, y1, y2, rr, cr;
	uint8 x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14;
	uint8 *row0, *row1, *row2, *row3, *row4;

	rr = (nrh + 1) % order;
	cr = (nch + 1) % order;
	
	for (row = nrl; row < nrh + 1 - rr; row += order) {
		row0 = ppInput[row - 1];
		row1 = ppInput[row + 0];
		row2 = ppInput[row + 1];
		row3 = ppInput[row + 2];
		row4 = ppInput[row + 3];
		for (col = ncl; col < nch + 1 - cr; col += order) {
			x0 = scalar_or3(row1, col); x1 = scalar_or3(row1, col + 1); x2 = scalar_or3(row1, col + 2);
			x3 = scalar_or3(row2, col); x4 = scalar_or3(row2, col + 1); x5 = scalar_or3(row2, col + 2);
			x6 = scalar_or3(row3, col); x7 = scalar_or3(row3, col + 1); x8 = scalar_or3(row3, col + 2);

			ppOutput[row + 0][col + 0] = scalar_or3(row0, col + 0) | x0 | x3;
			ppOutput[row + 1][col + 0] = 				         x0 | x3 | x6;
			ppOutput[row + 2][col + 0] = 				         x3 | x6 | scalar_or3(row4, col + 0);


			ppOutput[row + 0][col + 1] = scalar_or3(row0, col + 1) | x1 | x4;
			ppOutput[row + 1][col + 1] = 				         x1 | x4 | x7;
			ppOutput[row + 2][col + 1] = 				         x4 | x7 | scalar_or3(row4, col + 1);


			ppOutput[row + 0][col + 2] = scalar_or3(row0, col + 2) | x2 | x5;
			ppOutput[row + 1][col + 2] = 				         x2 | x5 | x8;
			ppOutput[row + 2][col + 2] = 				         x5 | x8 | scalar_or3(row4, col + 2);

		}
	}
	switch (rr) {
		case 2 :
			switch(cr){
				case 2 :
					ui8matrix_dilation_LU3x3_O3xO1(ppInput, nrl, nrh, nch - 1, nch, s, ppOutput);
					ui8matrix_dilation_LU3x3_O1xO3(ppInput, nrh - 1, nrh, ncl, nch, s, ppOutput);
				break;
				case 1: 
					ui8matrix_dilation_LU3x3_O3xO1(ppInput, nrl, nrh, nch, nch, s, ppOutput);		
					ui8matrix_dilation_LU3x3_O1xO3(ppInput, nrh - 1, nrh, ncl, nch, s, ppOutput);	
				break;
				default :
					ui8matrix_dilation_LU3x3_O1xO3(ppInput, nrh - 1, nrh, ncl, nch, s, ppOutput);
				break;
			}
		case 1:
			switch(cr){
				case 2 :
					ui8matrix_dilation_LU3x3_O3xO1(ppInput, nrl, nrh, nch - 1, nch, s, ppOutput);
					ui8matrix_dilation_LU3x3_O1xO3(ppInput, nrh, nrh, ncl, nch, s, ppOutput);
				break;
				case 1: 
					ui8matrix_dilation_LU3x3_O3xO1(ppInput, nrl, nrh, nch, nch, s, ppOutput);		
					ui8matrix_dilation_LU3x3_O1xO3(ppInput, nrh, nrh, ncl, nch, s, ppOutput);	
				break;
				default :
					ui8matrix_dilation_LU3x3_O1xO3(ppInput, nrh, nrh, ncl, nch, s, ppOutput);
				break;
			}
		break;
		default:
			switch(cr){
				case 2 :
					ui8matrix_dilation_LU3x3_O3xO1(ppInput, nrl, nrh, nch - 1, nch, s, ppOutput);
				break;
				case 1: 
					ui8matrix_dilation_LU3x3_O3xO1(ppInput, nrl, nrh, nch, nch, s, ppOutput);	
				break;
				default :
				break;
			}
		break;
	}
	// ui8matrix_dilation_LU3x3(ppInput, row, nrh, col, nrh, s, ppOutput);
}

void ui8matrix_dilation_LU3x3_O3xO3_RR(uint8** ppInput, long nrl, long nrh, long ncl, long nch, p_struct_elem_dim s, uint8 **ppOutput)
{
	const long order = 3;
	long row = nrl, col = ncl, x, y;
	uint8 y0, y1, y2, y3, rr, cr;
	uint8 x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14;
	uint8 *row0, *row1, *row2, *row3, *row4;

	rr = (nrh + 1) % order;
	cr = (nch + 1) % order;
	
	row0 = ppInput[row - 1];
	row1 = ppInput[row + 0];
	for (row = nrl; row < nrh + 1 - rr; row += order) {
		row2 = ppInput[row + 1];
		row3 = ppInput[row + 2];
		row4 = ppInput[row + 3];
		
		x0 = row1[-1]; x3 = row1[ 0]; 
		x1 = row2[-1]; x4 = row2[ 0]; 
		x2 = row3[-1]; x5 = row3[ 0]; 
		for (col = ncl; col < nch + 1 - cr; col += order) {
			x6 = row1[col + 1]; x9  = row1[col + 2]; x12 = row1[col + 3];
			x7 = row2[col + 1]; x10 = row2[col + 2]; x13 = row2[col + 3];
			x8 = row3[col + 1]; x11 = row3[col + 2]; x14 = row3[col + 3];

			y1 = (x0 | x3 | x6);
			y2 = (x1 | x4 | x7);
			y3 = (x2 | x5 | x8);
			ppOutput[row + 0][col + 0] = scalar_or3(row0, col + 0) | y1 | y2;
			ppOutput[row + 1][col + 0] = 						 y1 | y2 | y3;
			ppOutput[row + 2][col + 0] = 						 y2 | y3 | scalar_or3(row4, col + 0);

			y1 = (x3 | x6 | x9);
			y2 = (x4 | x7 | x10);
			y3 = (x5 | x8 | x11);
			ppOutput[row + 0][col + 1] = scalar_or3(row0, col + 1) | y1 | y2;
			ppOutput[row + 1][col + 1] = 						 y1 | y2 | y3;
			ppOutput[row + 2][col + 1] = 						 y2 | y3 | scalar_or3(row4, col + 1);
			
			
			y1 = (x6 | x9  | x12);
			y2 = (x7 | x10 | x13);
			y3 = (x8 | x11 | x14);			
			ppOutput[row + 0][col + 2] = scalar_or3(row0, col + 2) | y1 | y2;
			ppOutput[row + 1][col + 2] = 						 y1 | y2 | y3;
			ppOutput[row + 2][col + 2] = 						 y2 | y3 | scalar_or3(row4, col + 2);
			
			// printf("%u %u %u %u %u\n", x0, x3, x6, x9 , x12);
			// printf("%u %u %u %u %u\n", x1, x4, x7, x10, x13);
			// printf("%u %u %u %u %u\n", x2, x5, x8, x11, x14);	
			// display_ui8matrix(ppOutput, nrl, nrh, ncl, nch, "%u", "debug");
			// getchar();
			x0 = x9 ; x3 = x12;
			x1 = x10; x4 = x13;
			x2 = x11; x5 = x14;
		}
		row0 = row3;
		row1 = row4;
	}
	switch(rr) {
	case 2: 
		row2 = ppInput[row + 1];
		row3 = ppInput[row + 2];
		switch(cr) {
		case 2:
			x0 = row1[-1]; x3 = row1[ 0]; x6 = row1[col + 1];
			x1 = row2[-1]; x4 = row2[ 0]; x7 = row2[col + 1];
			x2 = row3[-1]; x5 = row3[ 0]; x8 = row3[col + 1];
			for (col; col < nch + 1; col ++) {
				x9  = row1[col + 2]; 
				x10 = row2[col + 2]; 
				x11 = row3[col + 2]; 
				ppOutput[row + 0][col + 0] = scalar_or3(row0, col + 0) | (x0 | x3 | x6) | (x1 | x4 | x7);
				ppOutput[row + 1][col + 0] = 		     (x0 | x3 | x6) | (x1 | x4 | x7) | (x2 | x5 | x8);

				ppOutput[row + 0][col + 1] = scalar_or3(row0, col + 1) | (x3 | x6 | x9 ) | (x4 | x7 | x10);
				ppOutput[row + 1][col + 1] =   		    (x3 | x6 | x9 ) | (x4 | x7 | x10) | (x5 | x8 | x11);
			
				x0 = x3; x3 = x6; x6 = x9 ;
				x1 = x4; x4 = x7; x7 = x10;
				x2 = x5; x5 = x8; x8 = x11;
			}
			break;
		case 1:
			x0 = row1[-1]; x3 = row1[ 0]; 
			x1 = row2[-1]; x4 = row2[ 0]; 
			x2 = row3[-1]; x5 = row3[ 0]; 
			for (col; col < nch + 1; col ++) {
				x6 = row1[col + 1];
				x7 = row2[col + 1];
				x8 = row3[col + 1];
				ppOutput[row + 0][col + 0] = scalar_or3(row0, col + 0) | (x0 | x3 | x6) | (x1 | x4 | x7);
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
		row2 = ppInput[row + 1];
		switch(cr) {
		case 2: 
			x0 = row1[-1]; x3 = row1[ 0]; x6 = row1[col + 1];
			x1 = row2[-1]; x4 = row2[ 0]; x7 = row2[col + 1];
			for (col; col < nch + 1; col ++) {
				x9  = row1[col + 2]; 
				x10 = row2[col + 2]; 
				ppOutput[row + 0][col + 0] = scalar_or3(row0, col + 0) | (x0 | x3 | x6) | (x1 | x4 | x7);
				ppOutput[row + 0][col + 1] = scalar_or3(row0, col + 1) | (x3 | x6 | x9) | (x4 | x7 | x10);
			
				x0 = x3; x3 = x6; x6 = x9 ;
				x1 = x4; x4 = x7; x7 = x10;
				x2 = x5; x5 = x8; x8 = x11;
			}
			break;
		case 1:
			x0 = row1[-1]; x3 = row1[ 0]; 
			x1 = row2[-1]; x4 = row2[ 0]; 
			for (col; col < nch + 1; col ++) {
				x6 = row1[col + 1];
				x7 = row2[col + 1];
				ppOutput[row + 0][col + 0] = scalar_or3(row0, col + 0) | (x0 | x3 | x6) | (x1 | x4 | x7);
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
	uint8 *row0, *row1, *row2, *row3, *row4;

	r = (nrh + 1)  % order;
	
	row0 = ppInput[row - 1];
	row1 = ppInput[row + 0];
	for (row = nrl; row < nrh + 1 - r; row += order) {
		row2 = ppInput[row + 1]; 
		row3 = ppInput[row + 2]; 
		row4 = ppInput[row + 3];
		x0 = row1[-1]; x3 = row1[ 0];
		x1 = row2[-1]; x4 = row2[ 0];
		x2 = row3[-1]; x5 = row3[ 0];
		for (col = ncl; col < nch + 1; col++) {
			x6 = row1[col + 1]; //scalar_or3(row1, col);
			x7 = row2[col + 1]; //scalar_or3(row2, col);
			x8 = row3[col + 1]; //scalar_or3(row3, col);
			y0 = (x0 | x3 | x6); // 2 calc
			y1 = (x1 | x4 | x7); // 2 calc
			y2 = (x2 | x5 | x8); // 2 calc

			ppOutput[row + 0][col] = scalar_or3(row0, col) | y0 | y1; 						// 3 + 2    =  5
			ppOutput[row + 1][col] =         			 y0 | y1 | y2;						//     2    =  2
			ppOutput[row + 2][col] = 		 			 y1 | y2 | scalar_or3(row4, col);  //     2 + 3=  5
			// 6 + 12 = 18 calc
			x0 = x3; x3 = x6;
			x1 = x4; x4 = x7;
			x2 = x5; x5 = x8;
		}
		row0 = row3;
		row1 = row4;
	}
	row2 = ppInput[row + 1];
	row3 = ppInput[row + 2];
	row4 = ppInput[row + 3];
	switch(r) {
		case 2:
			x0 = row1[-1]; x3 = row1[ 0];
			x1 = row2[-1]; x4 = row2[ 0];
			x2 = row3[-1]; x5 = row3[ 0];
			for (col = ncl; col < nch + 1; col++) {
				x6 = row1[col + 1]; //scalar_or3(row1, col);
				x7 = row2[col + 1]; //scalar_or3(row2, col);
				x8 = row3[col + 1]; //scalar_or3(row3, col);
				ppOutput[row + 0][col] = scalar_or3(row0, col) | (x0 | x3 | x6) | (x1 | x4 | x7);;
				ppOutput[row + 1][col] =         (x0 | x3 | x6) | (x1 | x4 | x7) | (x2 | x5 | x8);;
				x0 = x3; x3 = x6;
				x1 = x4; x4 = x7;
				x2 = x5; x5 = x8;
			}
		case 1:
			x0 = row1[-1]; x3 = row1[ 0];
			x1 = row2[-1]; x4 = row2[ 0];
			for (col = ncl; col < nch + 1; col++) {
				x6 = row1[col + 1]; //scalar_or3(row1, col);
				x7 = row2[col + 1]; //scalar_or3(row2, col);
				ppOutput[row + 0][col] = scalar_or3(row0, col) | (x0 | x3 | x6) | (x1 | x4 | x7);
				x0 = x3; x3 = x6;
				x1 = x4; x4 = x7;
			}				
			break;
		default:
			break;
	}
}

// Refer to the DEPRECATED VERSIONS OF ui8matrix_dilation_LU3x3_O1xO3_RR to see other versions of column loop unrolling
void ui8matrix_dilation_LU3x3_O1xO3_RR (uint8** ppInput, long nrl, long nrh, long ncl, long nch, p_struct_elem_dim s, uint8 **ppOutput)
{
	const long order = 3;
	long row = nrl, col = ncl, r;
	uint8 y0, y1, y2;
	uint8 x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14;
	uint8 *row0, *row1, *row2, *row3, *row4;

	r = (nch + 1) % order;
	row0 = ppInput[row - 1];
	row1 = ppInput[row + 0];
	for (row = nrl; row < nrh + 1; row++) {
		row2 = ppInput[row + 1];

		x0  = row0[-1]; x1  = row0[ 0]; 
		x5  = row1[-1]; x6  = row1[ 0]; 
		x10 = row2[-1]; x11 = row2[ 0]; 
		for (col = ncl; col < nch + 1 - r; col += order) {
			x2  = row0[col + 1]; x3  = row0[col + 2]; x4  = row0[col + 3];
			x7  = row1[col + 1]; x8  = row1[col + 2]; x9  = row1[col + 3];
			x12 = row2[col + 1]; x13 = row2[col + 2]; x14 = row2[col + 3];
			

			y1 = x1 | x2 | x6 | x7 | x11 | x12; // 5 calcs
			ppOutput[row][col + 0] = y1 | x0 | x5 | x10; // 3 calc	
			ppOutput[row][col + 1] = y1 | x3 | x8 | x13; // 3 calc	
			ppOutput[row][col + 2] =     (x2 | x3 | x4) | (x8 | x7 | x9) | (x12 | x13 | x14); //8 calc 

			// 8 + 6 + 5 = 21 calcs
			x0  = x3;  x1  = x4;
			x5  = x8;  x6  = x9;
			x10 = x13; x11 = x14;
			
		}
		row0 = row1;
		row1 = row2;
	}

		
	row = nrl;
	row0 = ppInput[nrl - 1];
	row1 = ppInput[nrl + 0];
	
	switch(r) {
	case 2:
		for (row; row < nrh + 1; row++) {
			row2 = ppInput[row + 1];

			y1 = x1 | x2 | x6 | x7 | x11 | x12; // 5
			x0  = row0[col - 1]; x1  = row0[col +  0]; x2  = row0[col + 1]; x3  = row0[col + 2];
			x5  = row1[col - 1]; x6  = row1[col +  0]; x7  = row1[col + 1]; x8  = row1[col + 2];
			x10 = row2[col - 1]; x11 = row2[col +  0]; x12 = row2[col + 1]; x13 = row2[col + 2];
			ppOutput[row][col + 0] = y1 | x0 | x5 | x10; // 3 calc	
			ppOutput[row][col + 1] = y1 | x3 | x8 | x13; // 3 calc	
			row0 = row1;
			row1 = row2;
		}
		break;
	case 1:
		for (row; row < nrh + 1; row++) {			
			row2 = ppInput[row + 1];
			x0  = row0[col - 1]; x1  = row0[col +  0]; x2  = row0[col + 1];
			x5  = row1[col - 1]; x6  = row1[col +  0]; x7  = row1[col + 1];
			x10 = row2[col - 1]; x11 = row2[col +  0]; x12 = row2[col + 1];
			ppOutput[row][col + 0] = (x0 | x1 | x2 ) | (x5 | x6 | x7) | (x10 | x11 | x12);
			row0 = row1;
			row1 = row2;
		}
		break;
	default:
		break;
	}
}


void ui8matrix_dilation_pipeline3x3(uint8** ppInput, long nrl, long nrh, long ncl, long nch, p_struct_elem_dim s, uint8 **ppOutput)
{
	long row = nrl, col = ncl, x, y;
	uint8 **temp = ui8matrix(nrl + s->nrl, nrh + s->nrh, ncl + s->ncl, nch + s->nch);
	
	// Prologue	
	for (col = ncl; col < nch + 1; col++) {
    	temp[row - 1][col] = scalar_or3(ppInput[row - 1], col);
		temp[row + 0][col] = scalar_or3(ppInput[row + 0], col);
	}	

	for (row = row + 1; row < nrh + 1; row++){
		for (col = ncl; col < nch + 1; col++) {
			temp[row + 0][col] = scalar_or3(ppInput[row + 0], col);
			ppOutput[row - 1][col] = temp[row - 2][ col] |
									 temp[row - 1][ col] |
									 temp[row - 0][ col];
		}
	}

	for (col = ncl; col < nch + 1; col++) {

    	ppOutput[nrh + 0][col] =    temp[nrh - 1][ col] |
						   	        temp[nrh + 0][ col] |
			          scalar_or3(ppInput[nrh + 1], col);
	}
	free_ui8matrix(temp, nrl + s->nrl, nrh + s->nrh, ncl + s->ncl, nch + s->nch);
}
void ui8matrix_dilation_pipeline_LU3x3_O3xO1 (uint8** ppInput, long nrl, long nrh, long ncl, long nch, p_struct_elem_dim s, uint8 **ppOutput)
{
	const long order = 3;
	long row = nrl, col = ncl, x, y, r = 0;
	uint8 **temp = ui8matrix(nrl + s->nrl, nrh + s->nrh, ncl + s->ncl, nch + s->nch);
	uint8 *row0, *row1, *row2, *row3, x0, x1, x2;
	uint8 *temp_row0, *temp_row1, *temp_row2, *temp_row3, *temp_row4;
	uint8 *out_row0, *out_row1, *out_row2, *out_row3;
	
	r = (nrh + 1) % order;
	// Prologue	
	row0 = ppInput[row - 1];
	row1 = ppInput[row + 0];
	row2 = ppInput[row + 1];
	
	temp_row1 = temp[row - 1];
	temp_row2 = temp[row + 0];
	temp_row3 = temp[row + 1];
	if (nrh - nrl > 0){
		row3 = ppInput[row + 2];
		temp_row4 = temp[row + 2];
		for (col = ncl; col < nch + 1; col++) {
			temp_row1[col] = scalar_or3(row0, col);
			temp_row2[col] = scalar_or3(row1, col);
			temp_row3[col] = scalar_or3(row2, col);
			temp_row4[col] = scalar_or3(row3, col);
		}	
	}
	else {
		for (col = ncl; col < nch + 1; col++) {
			temp_row1[col] = scalar_or3(row0, col);
			temp_row2[col] = scalar_or3(row1, col);
			temp_row3[col] = scalar_or3(row2, col);
		}	
	}

	for (row = row + 1; row < nrh + 1 - r; row += order){
		row0 = ppInput[row - 1];
		row1 = ppInput[row + 0];
		row2 = ppInput[row + 1];
		row3 = ppInput[row + 2];
		temp_row0 = temp[row - 2];
		temp_row1 = temp[row - 1];
		temp_row2 = temp[row + 0];
		temp_row3 = temp[row + 1];
		temp_row4 = temp[row + 2];
		out_row0 = ppOutput[row - 1];
		out_row1 = ppOutput[row + 0];
		out_row2 = ppOutput[row + 1];
		for (col = ncl; col < nch + 1; col++) {
			temp_row4[col] = scalar_or3(row3, col);
			x0 = temp_row1[col];
			x1 = scalar_or3(row1, col);
			x2 = scalar_or3(row2, col);
			out_row0[col] = temp_row0[ col] | x0 | x1;
			out_row1[col] = 			 x0 | x1 | x2;
			out_row2[col] = 			 x1 | x2 | temp_row4[ col];	
			temp_row2[col] = x1;
			temp_row3[col] = x2;
		}
	}
	
	row1 = ppInput[nrh + 0];
	row2 = ppInput[nrh + 1];
	temp_row1 = temp[nrh - 1];
	temp_row2 = temp[nrh + 0];
	out_row1 = ppOutput[nrh + 0];
	switch (r) {
		case 2: 
			out_row0 = ppOutput[row - 1];
			temp_row0 = temp[row - 2];
			for (col = ncl; col < nch + 1; col++) {
				x1 = scalar_or3(row1, col);
				x2 = temp_row2[col];
				out_row0[col] = temp_row0[col] | x1 | x2;
				out_row1[col] = 			x1 | x2 | scalar_or3(row2, col);
				temp[nrh + 0][col] = x1;
			}
		case 1: 
			for (col = ncl; col < nch + 1; col++) {
				out_row1[col] = temp_row1[col] |
								temp_row2[col] |
			  				  	scalar_or3(row2, col);
			}
	}
	free_ui8matrix(temp, nrl + s->nrl, nrh + s->nrh, ncl + s->ncl, nch + s->nch);
}

void ui8matrix_dilation_pipeline_LU3x3_O3xO1_RR (uint8** ppInput, long nrl, long nrh, long ncl, long nch, p_struct_elem_dim s, uint8 **ppOutput)
{
	const long order = 3;
	long row = nrl, col = ncl, x, y, r = 0;
	uint8 **temp = ui8matrix(nrl + s->nrl, nrh + s->nrh, ncl + s->ncl, nch + s->nch);
	uint8 *row0, *row1, *row2, *row3, x0, x1, x2, x3;
	uint8 *temp_row0, *temp_row1, *temp_row2, *temp_row3, *temp_row4;
	uint8 *out_row0, *out_row1, *out_row2, *out_row3;
	
	r = (nrh + 1) % order;
	// Prologue	
	row0 = ppInput[row - 1];
	row1 = ppInput[row + 0];
	row2 = ppInput[row + 1];
	
	temp_row1 = temp[row - 1];
	temp_row2 = temp[row + 0];
	temp_row3 = temp[row + 1];
	if (nrh - nrl > 0){
		row3 = ppInput[row + 2];
		temp_row4 = temp[row + 2];
		for (col = ncl; col < nch + 1; col++) {
			temp_row1[col] = scalar_or3(row0, col);
			temp_row2[col] = scalar_or3(row1, col);
			temp_row3[col] = scalar_or3(row2, col);
			temp_row4[col] = scalar_or3(row3, col);
		}	
	}
	else {
		for (col = ncl; col < nch + 1; col++) {
			temp_row1[col] = scalar_or3(row0, col);
			temp_row2[col] = scalar_or3(row1, col);
			temp_row3[col] = scalar_or3(row2, col);
		}	
	}

	row = row + 1;
	row0 = ppInput[row - 1];
	row1 = ppInput[row + 0];
	temp_row0 = temp[row - 2];
	temp_row1 = temp[row - 1];
	
	for (row; row < nrh + 1 - r; row += order){
		row2 = ppInput[row + 1];
		row3 = ppInput[row + 2];
		temp_row2 = temp[row + 0];
		temp_row3 = temp[row + 1];
		temp_row4 = temp[row + 2];

		out_row0 = ppOutput[row - 1];
		out_row1 = ppOutput[row + 0];
		out_row2 = ppOutput[row + 1];
		for (col = ncl; col < nch + 1; col++) {
			x0 = temp_row1[col];
			x1 = scalar_or3(row1, col);
			x2 = scalar_or3(row2, col);
			x3 = scalar_or3(row3, col);
			out_row0[col] = x3 | x0 | x1;
			out_row1[col] = x0 | x1 | x2;
			out_row2[col] = x1 | x2 | x3;	
			temp_row2[col] = x1;
			temp_row3[col] = x2;
			temp_row4[col] = x3;
		}
		row0 = row2;
		row1 = row3;
		temp_row0 = temp_row3;
		temp_row1 = temp_row4;
	}
	
	row1 = ppInput[nrh + 0];
	row2 = ppInput[nrh + 1];
	temp_row1 = temp[nrh - 1];
	temp_row2 = temp[nrh + 0];
	out_row1 = ppOutput[nrh + 0];

	switch (r) {
		case 2: 
			out_row0 = ppOutput[row - 1];
			temp_row0 = temp[row - 2];
			for (col = ncl; col < nch + 1; col++) {
				x1 = scalar_or3(row1, col);
				x2 = temp_row2[col];
				out_row0[col] = temp_row0[col] | x1 | x2;
				out_row1[col] = 			x1 | x2 | scalar_or3(row2, col);
				temp[nrh + 0][col] = x1;
			}
		case 1: 
			for (col = ncl; col < nch + 1; col++) {
				out_row1[col] = temp_row1[col] |
								temp_row2[col] |
			  				  	scalar_or3(row2, col);
			}
	}
	free_ui8matrix(temp, nrl + s->nrl, nrh + s->nrh, ncl + s->ncl, nch + s->nch);
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




// void ui8matrix_dilation_LU3x3_O1xO3_RR_second_best (uint8** ppInput, long nrl, long nrh, long ncl, long nch, p_struct_elem_dim s, uint8 **ppOutput)
// {
// 	const long order = 3;
// 	long row = nrl, col = ncl, x, y, r;
// 	uint8 y0, y1, y2, y3, y4, y5;
// 	uint8 x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16, x17, x18, x19;
// 	uint8 *row0, *row1, *row2, *row3, *row4;

// 	r = (nch + 1) % order;
// 	row0 = ppInput[row - 1];
// 	row1 = ppInput[row + 0];
// 	for (row = nrl; row < nrh + 1; row++) {
// 		row2 = ppInput[row + 1];

// 		x0  = row0[-1]; x1  = row0[ 0]; 
// 		x5  = row1[-1]; x6  = row1[ 0]; 
// 		x10 = row2[-1]; x11 = row2[ 0]; 
// 		for (col = ncl; col < nch + 1 - r; col += order) {
// 			x2  = row0[col + 1]; x3  = row0[col + 2]; x4  = row0[col + 3];
// 			x7  = row1[col + 1]; x8  = row1[col + 2]; x9  = row1[col + 3];
// 			x12 = row2[col + 1]; x13 = row2[col + 2]; x14 = row2[col + 3];
// 			// y0 = x1  | x2; 
// 			// y1 = x2  | x3;
// 			// y2 = x6  | x7;
// 			// y3 = x7  | x8;
// 			// y4 = x11 | x12; 
// 			// y5 = x12 | x13;
// 			// ppOutput[row][col + 0] = (x0 | x1 | x2 ) | (x5 | x6 | x7) | (x10 | x11 | x12);
// 			// ppOutput[row][col + 1] = (x1 | x2 | x3 ) | (x6 | x7 | x8) | (x11 | x12 | x13);
// 			// ppOutput[row][col + 2] = (x2 | x3 | x4 ) | (x7 | x8 | x9) | (x12 | x13 | x14);
			
// 			ppOutput[row][col + 0] = x0 | x1 | x2 |           x5 | x6 | x7 |           x10 | x11 | x12; 			//8 calc -> 3 calc
// 			ppOutput[row][col + 1] =      x1 | x2 | x3 |           x6 | x7 | x8 |            x11 | x12 | x13; 		//8 calc -> 3 calc
// 			ppOutput[row][col + 2] =           x2 | x3 | x4 |           x7 | x8 | x9 |             x12 | x13 | x14; //8 calc -> 3 calc
// 			x0  = x3;  x1  = x4;
// 			x5  = x8;  x6  = x9;
// 			x10 = x13; x11 = x14;
			
// 		}
// 		row0 = row1;
// 		row1 = row2;
// 	}

		
// 	row = nrl;
// 	row0 = ppInput[nrl - 1];
// 	row1 = ppInput[nrl + 0];
	
// 	switch(r) {
// 	case 2:
// 		for (row; row < nrh + 1; row++) {
// 			row2 = ppInput[row + 1];
// 			x0  = row0[col - 1]; x1  = row0[col +  0]; x2  = row0[col + 1]; x3  = row0[col + 2];
// 			x5  = row1[col - 1]; x6  = row1[col +  0]; x7  = row1[col + 1]; x8  = row1[col + 2];
// 			x10 = row2[col - 1]; x11 = row2[col +  0]; x12 = row2[col + 1]; x13 = row2[col + 2];
// 			ppOutput[row][col + 0] = (x0 | x1 | x2 ) | (x5 | x6 | x7) | (x10 | x11 | x12);
// 			ppOutput[row][col + 1] = (x1 | x2 | x3 ) | (x6 | x7 | x8) | (x11 | x12 | x13);
// 			row0 = row1;
// 			row1 = row2;
// 		}
// 		break;
// 	case 1:
// 		for (row; row < nrh + 1; row++) {			
// 			row2 = ppInput[row + 1];
// 			x0  = row0[col - 1]; x1  = row0[col +  0]; x2  = row0[col + 1];
// 			x5  = row1[col - 1]; x6  = row1[col +  0]; x7  = row1[col + 1];
// 			x10 = row2[col - 1]; x11 = row2[col +  0]; x12 = row2[col + 1];
// 			ppOutput[row][col + 0] = (x0 | x1 | x2 ) | (x5 | x6 | x7) | (x10 | x11 | x12);
// 			row0 = row1;
// 			row1 = row2;
// 		}
// 		break;
// 	default:
// 		break;
// 	}
// }

// DEPRECATED VERSIONS OF ui8matrix_dilation_LU3x3_O1xO3_RR
// void ui8matrix_dilation_LU3x3_O1xO3_RR_worst (uint8** ppInput, long nrl, long nrh, long ncl, long nch, p_struct_elem_dim s, uint8 **ppOutput)
// {
// 	const long order = 3;
// 	long row = nrl, col = ncl, x, y, r;
// 	uint8 y0, y1, y2, y3, y4, y5;
// 	uint8 x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16, x17, x18, x19;
// 	uint8 *row0, *row1, *row2, *row3, *row4;

// 	r = (nch + 1) % order;
// 	row0 = ppInput[row - 1];
// 	row1 = ppInput[row + 0];
// 	for (row = nrl; row < nrh + 1; row++) {
// 		row2 = ppInput[row + 1];

// 		x0  = row0[-1]; x1  = row0[ 0]; 
// 		x5  = row1[-1]; x6  = row1[ 0]; 
// 		x10 = row2[-1]; x11 = row2[ 0]; 
// 		for (col = ncl; col < nch + 1 - r; col += order) {
// 			x2  = row0[col + 1]; x3  = row0[col + 2]; x4  = row0[col + 3];
// 			x7  = row1[col + 1]; x8  = row1[col + 2]; x9  = row1[col + 3];
// 			x12 = row2[col + 1]; x13 = row2[col + 2]; x14 = row2[col + 3];
			
// 			y1 = x1 | x2 | x6 | x7 | x11 | x12; // 5 calcs
// 			y2 = x3 | x8 | x13;					// 2 calcs
// 			ppOutput[row][col + 0] = x0 | x5 | y1 | x10; //3 calcs
// 			ppOutput[row][col + 1] =      y2 | y1; 		 //1 calc
// 			ppOutput[row][col + 2] =      y2 | x2 | x4          | x7 | x9             | x12 | x14; //6 calcs
// 			// 5 + 2 + 3 + 1 + 6 = 17						

// 			x0  = x3;  x1  = x4;
// 			x5  = x8;  x6  = x9;
// 			x10 = x13; x11 = x14;
			
// 		}
// 		row0 = row1;
// 		row1 = row2;
// 	}

		
// 	row = nrl;
// 	row0 = ppInput[nrl - 1];
// 	row1 = ppInput[nrl + 0];
	
// 	switch(r) {
// 	case 2:
// 		for (row; row < nrh + 1; row++) {
// 			row2 = ppInput[row + 1];
// 			x0  = row0[col - 1]; x1  = row0[col +  0]; x2  = row0[col + 1]; x3  = row0[col + 2];
// 			x5  = row1[col - 1]; x6  = row1[col +  0]; x7  = row1[col + 1]; x8  = row1[col + 2];
// 			x10 = row2[col - 1]; x11 = row2[col +  0]; x12 = row2[col + 1]; x13 = row2[col + 2];
// 			y1 = x1 | x2 | x6 | x7 | x11 | x12; // 5 calcs
// 			y2 = x3 | x8 | x13;					// 2 calcs
// 			ppOutput[row][col + 0] = x0 | x5 | y1 | x10; //3 calcs
// 			ppOutput[row][col + 1] =      y2 | y1; 		 //1 calc
// 			row0 = row1;
// 			row1 = row2;
// 		}
// 		break;
// 	case 1:
// 		for (row; row < nrh + 1; row++) {			
// 			row2 = ppInput[row + 1];
// 			x0  = row0[col - 1]; x1  = row0[col +  0]; x2  = row0[col + 1];
// 			x5  = row1[col - 1]; x6  = row1[col +  0]; x7  = row1[col + 1];
// 			x10 = row2[col - 1]; x11 = row2[col +  0]; x12 = row2[col + 1];
// 			ppOutput[row][col + 0] = (x0 | x1 | x2 ) | (x5 | x6 | x7) | (x10 | x11 | x12);
// 			row0 = row1;
// 			row1 = row2;
// 		}
// 		break;
// 	default:
// 		break;
// 	}
// }