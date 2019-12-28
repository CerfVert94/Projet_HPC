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


void ui8matrix_dilation_LU3x3_O1xO3_RR_OMP(uint8** ppInput, long nrl, long nrh, long ncl, long nch, p_struct_elem_dim s, uint8 **ppOutput)
{
	const long order = 3;
	long row = nrl, col = ncl, r;
	uint8 y0, y1, y2;
	uint8 x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14;
	uint8 *row0, *row1, *row2, *row3, *row4;
	uint8 *out_row0;

	r = (nch + 1) % order;
	row0 = ppInput[nrl - 1];
	row1 = ppInput[nrl + 0];
	for (row = nrl; row < nrh + 1; row++) {
		// row0 = ppInput[row - 1];
		// row1 = ppInput[row + 0];
		row2 = ppInput[row + 1];

		out_row0 = ppOutput[row + 0];
		x0  = row0[ncl - 1]; x1  = row0[ncl + 0];
		x5  = row1[ncl - 1]; x6  = row1[ncl + 0];
		x10 = row2[ncl - 1]; x11 = row2[ncl + 0];
		
		for (col = ncl; col < nch + 1 - r; col += order) {
			x0  = row0[col - 1]; x1  = row0[col + 0]; x2  = row0[col + 1]; x3  = row0[col + 2]; x4  = row0[col + 3];
			x5  = row1[col - 1]; x6  = row1[col + 0]; x7  = row1[col + 1]; x8  = row1[col + 2]; x9  = row1[col + 3];
			x10 = row2[col - 1]; x11 = row2[col + 0]; x12 = row2[col + 1]; x13 = row2[col + 2]; x14 = row2[col + 3];
			
			y0 = x1  | x2  ;
			y1 = x6  | x7  ;
			y2 = x11 | x12 ;

			out_row0[col + 0] = (x0  | y0) | 
								(x5  | y1) | 
								(x10 | y2);

			out_row0[col + 1] = (y0 | x3) | 
								(y1 | x8) | 
								(y2 | x13);
								
			out_row0[col + 2] = (x2  | x3  | x4) | 
								(x7  | x8  | x9) | 
								(x12 | x13 | x14); 

			x0  = x3;  x1  = x4;
			x5  = x8;  x6  = x9;
			x10 = x13; x11 = x14;
			
		}
		row0 = row1;
		row1 = row2;
	}

	switch(r) {
		case 2:		
		for (row = nrl; row < nrh + 1; row++) {
			row0 = ppInput[row - 1];
			row1 = ppInput[row + 0];
			row2 = ppInput[row + 1];

			out_row0 = ppOutput[row + 0];
		
			// y1 = x1 | x2 | x6 | x7 | x11 | x12; // 5
			x0  = row0[nch - 2]; x1  = row0[nch - 1]; x2  = row0[nch + 0]; x3  = row0[nch + 1]; 
			x5  = row1[nch - 2]; x6  = row1[nch - 1]; x7  = row1[nch + 0]; x8  = row1[nch + 1]; 
			x10 = row2[nch - 2]; x11 = row2[nch - 1]; x12 = row2[nch + 0]; x13 = row2[nch + 1]; 

			out_row0[nch - 1] = (x0  | x1  | x2) | 
								(x5  | x6  | x7) | 
								(x10 | x11 | x12);

			out_row0[nch + 0] = (x1  | x2  | x3) | 
			  					(x6  | x7  | x8) | 
								(x11 | x12 | x13);
			
		}
		break;
	case 1:
		for (row = nrl; row < nrh + 1; row++) {
			row0 = ppInput[row - 1];
			row1 = ppInput[row + 0];
			row2 = ppInput[row + 1];

			out_row0 = ppOutput[row + 0];
		
			// x1  = row0[nch - 1] | x2  = row0[nch + 0] | x3  = row0[nch + 1]; 
			// x6  = row1[nch - 1] | x7  = row1[nch + 0] | x8  = row1[nch + 1]; 
			// x11 = row2[nch - 1] | x12 = row2[nch + 0] | x13 = row2[nch + 1]; 

			out_row0[nch + 0] = (row0[nch - 1] |  row0[nch + 0] | row0[nch + 1])| 
			  					(row1[nch - 1] |  row1[nch + 0] | row1[nch + 1])| 
								(row2[nch - 1] |  row2[nch + 0] | row2[nch + 1]);
			
		}
		break;
	default:
		break;
	}
}

void ui8matrix_dilation_LU3x3_O3xO3_RR_OMP(uint8** ppInput, long nrl, long nrh, long ncl, long nch, p_struct_elem_dim s, uint8 **ppOutput)
{
	const long order = 3;
	long row = nrl, col = ncl, x, y, rr, cr;
	// uint8 y00, y10, y20, y30, y40, y01, y11, y21, y31, y41, y02, y12, y22, y32, y42;
	uint8 y0, y1, y2, y3, y4, y5, y6, y7, y8, y9, y10, y11, y12, y13, y14;
	uint8 x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16, x17, x18, x19, x20, x21, x22, x23, x24, x25;
	uint8 r0c0, r1c0, r2c0, r0c1,r1c1, r2c1, r0c2, r1c2, r2c2;
	uint8 *row0, *row1, *row2, *row3, *row4;
	uint8 *out_row0, *out_row1, *out_row2;

	rr = (nrh + 1) % order;
	cr = (nch + 1) % order;
	
	row0 = ppInput[row - 1];
	row1 = ppInput[row + 0];

	for (row = nrl; row < nrh + 1 - rr; row += order) {
		row2 = ppInput[row + 1];
		row3 = ppInput[row + 2];
		row4 = ppInput[row + 3];
		out_row0 = ppOutput[row + 0];
		out_row1 = ppOutput[row + 1];
		out_row2 = ppOutput[row + 2];


		x0  = row0[ncl - 1]; x1  = row0[ncl + 0]; x2  = row0[ncl + 1]; x3  = row0[ncl + 2]; x4  = row0[ncl + 3];
		x5  = row1[ncl - 1]; x6  = row1[ncl + 0]; x7  = row1[ncl + 1]; x8  = row1[ncl + 2]; x9  = row1[ncl + 3];
		x10 = row2[ncl - 1]; x11 = row2[ncl + 0]; x12 = row2[ncl + 1]; x13 = row2[ncl + 2]; x14 = row2[ncl + 3];
		x15 = row3[ncl - 1]; x16 = row3[ncl + 0]; x17 = row3[ncl + 1]; x18 = row3[ncl + 2]; x19 = row3[ncl + 3];
		x20 = row4[ncl - 1]; x21 = row4[ncl + 0]; x22 = row4[ncl + 1]; x23 = row4[ncl + 2]; x24 = row4[ncl + 3];
		y0  = x0  | x5  | x10; y1  = x1  | x6  | x11;
		y5  = x5  | x10 | x15; y6  = x6  | x11 | x16;
		y10 = x10 | x15 | x20; y11 = x11 | x16 | x21;

		for (col = ncl; col < nch + 1 - cr; col += order) {
			// x2  = row0[col + 1]; x3  = row0[col + 2]; x4  = row0[col + 3];
			// x7  = row1[col + 1]; x8  = row1[col + 2]; x9  = row1[col + 3];
			// x12 = row2[col + 1]; x13 = row2[col + 2]; x14 = row2[col + 3];
			// x17 = row3[col + 1]; x18 = row3[col + 2]; x19 = row3[col + 3];
			// x22 = row4[col + 1]; x23 = row4[col + 2]; x24 = row4[col + 3];
			x0  = row0[col - 1]; x1  = row0[col + 0]; x2  = row0[col + 1]; x3  = row0[col + 2]; x4  = row0[col + 3];
			x5  = row1[col - 1]; x6  = row1[col + 0]; x7  = row1[col + 1]; x8  = row1[col + 2]; x9  = row1[col + 3];
			x10 = row2[col - 1]; x11 = row2[col + 0]; x12 = row2[col + 1]; x13 = row2[col + 2]; x14 = row2[col + 3];
			x15 = row3[col - 1]; x16 = row3[col + 0]; x17 = row3[col + 1]; x18 = row3[col + 2]; x19 = row3[col + 3];
			x20 = row4[col - 1]; x21 = row4[col + 0]; x22 = row4[col + 1]; x23 = row4[col + 2]; x24 = row4[col + 3];


			y2  = x2  | x7  | x12; y3  = x3  | x8  | x13; y4  = x4  | x9  | x14; 
			y7  = x7  | x12 | x17; y8  = x8  | x13 | x18; y9  = x9  | x14 | x19; 
			y12 = x12 | x17 | x22; y13 = x13 | x18 | x23; y14 = x14 | x19 | x24; 
 
			out_row0[col + 0] = y0 | y1 | y2; 
			out_row0[col + 1] = y1 | y2 | y3; 
			out_row0[col + 2] = y2 | y3 | y4;
			y0 = y3; y1 = y4;

			out_row1[col + 0] = y5 | y6 | y7;
			out_row1[col + 1] = y6 | y7 | y8;
			out_row1[col + 2] = y7 | y8 | y9;
			y5 = y8; y6 = y9;
			
			out_row2[col + 0] = y10 | y11 | y12;
			out_row2[col + 1] = y11 | y12 | y13;
			out_row2[col + 2] = y12 | y13 | y14;
			y10 = y13; y11 = y14;
		}
		row0 = row3;
		row1 = row4;
	}
	switch (rr) {
		case 2 :
			switch(cr){
				case 2 :
					ui8matrix_dilation_LU3x3_O3xO1_RR(ppInput, nrl, nrh, nch - 1, nch, s, ppOutput);
					ui8matrix_dilation_LU3x3_O1xO3_RR(ppInput, nrh - 1, nrh, ncl, nch, s, ppOutput);
				break;
				case 1: 
					ui8matrix_dilation_LU3x3_O3xO1_RR(ppInput, nrl, nrh, nch, nch, s, ppOutput);		
					ui8matrix_dilation_LU3x3_O1xO3_RR(ppInput, nrh - 1, nrh, ncl, nch, s, ppOutput);	
				break;
				default :
					ui8matrix_dilation_LU3x3_O1xO3_RR(ppInput, nrh - 1, nrh, ncl, nch, s, ppOutput);
				break;
			}
			break;
		case 1:
			switch(cr){
				case 2 :
					ui8matrix_dilation_LU3x3_O3xO1_RR(ppInput, nrl, nrh, nch - 1, nch, s, ppOutput);
					ui8matrix_dilation_LU3x3_O1xO3_RR(ppInput, nrh, nrh, ncl, nch, s, ppOutput);
				break;
				case 1: 
					ui8matrix_dilation_LU3x3_O3xO1_RR(ppInput, nrl, nrh, nch, nch, s, ppOutput);		
					ui8matrix_dilation_LU3x3_O1xO3_RR(ppInput, nrh, nrh, ncl, nch, s, ppOutput);	
					// display_ui8matrix(ppOutput, nrl, nrh, ncl, nch, "%u", "rr=1; cr=1");
				break;
				default :
					ui8matrix_dilation_LU3x3_O1xO3_RR(ppInput, nrh, nrh, ncl, nch, s, ppOutput);
				break;
			}
		break;
		default:
			switch(cr){
				case 2 :
					ui8matrix_dilation_LU3x3_O3xO1_RR(ppInput, nrl, nrh, nch - 1, nch, s, ppOutput);
				break;
				case 1: 
					ui8matrix_dilation_LU3x3_O3xO1_RR(ppInput, nrl, nrh, nch, nch, s, ppOutput);	
				break;
				default :
				break;
			}
		break;
	}
		
}

void ui8matrix_dilation_pipeline_LU3x3_O1xO3_RR_OMP (uint8** ppInput, long nrl, long nrh, long ncl, long nch, p_struct_elem_dim s, uint8 **ppOutput)
{
	const long order = 3;
	long row = nrl, col = ncl, x, y, r = 0;
	uint8 **temp = ui8matrix(nrl, nrh, ncl + s->ncl, nch + s->nch);
	uint8 *row0, *row1, *row2, *row3, x0, x1, x2, x3, x4, x5;
	uint8 *temp_row0, *temp_row1, *temp_row2, *temp_row3, *temp_row4;
	uint8 *out_row0, *out_row1, *out_row2, *out_row3;
	
	r = (nch + 1) % order;
	// Prologue	
	row0 = ppInput[nrl - 1];
	row1 = ppInput[nrl + 0];
	for (row = nrl; row < nrh + 1; row++) {
		row2 = ppInput[row + 1];
		temp_row1 = temp[row + 0];

		temp_row1[ncl - 1] = row0[ncl - 1] | row1[ncl - 1] | row2[ncl - 1];
		temp_row1[ncl + 0] = row0[ncl + 0] | row1[ncl + 0] | row2[ncl + 0];
		row0 = row1;
		row1 = row2;
	}

	row0 = ppInput[nrl - 1];
	row1 = ppInput[nrl + 0];
	for (row = nrl; row < nrh + 1; row ++){
		temp_row1 = temp[row + 0];
		row2 = ppInput[row + 1];

		out_row1 = ppOutput[row]; 
		x0 = temp_row1[ncl - 1];		// LOAD => A
		x1 = temp_row1[ncl - 0];		// LOAD => B
		
		for (col = ncl; col < nch + 1 - r; col += order) {
			x2 = row0[col + 1] | row1[col + 1] | row2[col + 1];
			x3 = row0[col + 2] | row1[col + 2] | row2[col + 2];
			x4 = row0[col + 3] | row1[col + 3] | row2[col + 3];
			
			out_row1 [col + 0] = x0 | x1 | x2; // A | B | C
			out_row1 [col + 1] = x1 | x2 | x3; // B | C | D
			out_row1 [col + 2] = x2 | x3 | x4; // C | D | E

			x0 = x3; 						  // STORE => D
			x1 = x4; 						  // STORE => E 
											  // NEXT  :  F 
											  // NEXT  :  G 
											  // NEXT  :  H
		}
		row0 = row1;
		row1 = row2;
	}
	switch (r) {
		case 2: 
			row0 = ppInput[nrl - 1];
			row1 = ppInput[nrl + 0];
			for (row = nrl; row < nrh + 1; row++) {
				row2 = ppInput[row + 1];
				temp_row1 = temp[row + 0];
				out_row1 = ppOutput[row]; 

				x1 = row0[nch - 1] | row1[nch - 1] | row2[nch - 1]|row0[nch + 0] | row1[nch + 0] | row2[nch + 0];
				// x2 = ;

				out_row1 [nch - 1] = row0[nch - 2] | row1[nch - 2] | row2[nch - 2] | x1 ; // A | B | C
				out_row1 [nch + 0] = x1 | row0[nch + 1] | row1[nch + 0] | row2[nch + 0]; // B | C | D
				row0 = row1;
				row1 = row2;
			}
			break;
		case 1: 
			row0 = ppInput[nrl - 1];
			row1 = ppInput[nrl + 0];
			for (row = nrl; row < nrh + 1; row++) {
				row2 = ppInput[row + 1];
				// temp_row1 = temp[row + 0];
				ppOutput[row][nch] = row0[nch - 1] | row1[nch - 1] | row2[nch - 1] | row0[nch - 0] | row1[nch - 0] | row2[nch - 0] | row0[nch + 1] | row1[nch + 1] | row2[nch + 1]; // A | B | C
				// out_row1 [nch + 0] = x1 | x2 | row0[nch + 1] | row1[nch + 0] | row2[nch + 0]; // B | C | D
				row0 = row1;
				row1 = row2;
			}
			break;
		default:
			break;
	}
	// display_ui8matrix(ppOutput, nrl - 1, nrh + 1, ncl - 1, nch + 1, "%u", "TEst");
	free_ui8matrix(temp, nrl, nrh, ncl + s->ncl, nch + s->nch);
	
}

void ui8matrix_dilation_LU3x3_O1xO3_OMP (uint8** ppInput, long nrl, long nrh, long ncl, long nch, p_struct_elem_dim s, uint8 **ppOutput)
{
	const long order = 3;
	long row = nrl, col = ncl, x, y, r;
	uint8 y0, y1, y2, x0, x1, x2, x3, x4, x5, x6, x7, x8;
	uint8 *row0, *row1, *row2, *row3, *row4;
	uint8 *out_row0;
	r = (nch + 1)  % order;
	
	for (row = nrl; row < nrh + 1; row ++) {
		row0 = ppInput[row - 1];
		row1 = ppInput[row + 0];
		row2 = ppInput[row + 1];
		out_row0 = ppOutput[row + 0];
		for (col = ncl; col < nch + 1 - r; col += order) {
			// x0 = scalar_or3(row0, col); x3 = scalar_or3(row0, col + 1); x6 = scalar_or3(row0, col + 2);
			// x1 = scalar_or3(row1, col); x4 = scalar_or3(row1, col + 1); x7 = scalar_or3(row1, col + 2);
			// x2 = scalar_or3(row2, col); x5 = scalar_or3(row2, col + 1); x8 = scalar_or3(row2, col + 2);

			out_row0[col + 0] = scalar_or3(row0, col + 0)|
								scalar_or3(row1, col + 0)|
								scalar_or3(row2, col + 0);
			out_row0[col + 1] = scalar_or3(row0, col + 1)|
 								scalar_or3(row1, col + 1)|
 								scalar_or3(row2, col + 1);
			out_row0[col + 2] = scalar_or3(row0, col + 2)|
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
				out_row0 = ppOutput[row + 0];

				out_row0[nch - 1] = scalar_or3(row0, nch - 1)|
 									scalar_or3(row1, nch - 1)|
 									scalar_or3(row2, nch - 1);
				out_row0[nch + 0] = scalar_or3(row0, nch + 0)|
									scalar_or3(row1, nch + 0)|
									scalar_or3(row2, nch + 0);
			}
			break;
		case 1:
			for (row = nrl; row < nrh + 1; row ++) {
				ppOutput[row + 0][nch + 0] = scalar_or3(ppInput[row - 1], nch) | 
										  	 scalar_or3(ppInput[row + 0], nch) | 
										  	 scalar_or3(ppInput[row + 1], nch);
			}
			// ui8matrix_dilation_LU3x3_O1xO1(ppInput, nrl, nrh, nch, nch, s, ppOutput);
			break;
		default:
			break;
	}
}


void ui8matrix_dilation_LU3x3_O1xO3_RR (uint8** ppInput, long nrl, long nrh, long ncl, long nch, p_struct_elem_dim s, uint8 **ppOutput)
{
	const long order = 3;
	long row = nrl, col = ncl, r;
	uint8 y0, y1, y2;
	uint8 x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14;
	uint8 *row0, *row1, *row2, *row3, *row4;
	uint8 *out_row0;

	r = (nch + 1) % order;
	row0 = ppInput[nrl - 1];
	row1 = ppInput[nrl + 0];
	for (row = nrl; row < nrh + 1; row++) {
		// row0 = ppInput[row - 1];
		// row1 = ppInput[row + 0];
		row2 = ppInput[row + 1];

		out_row0 = ppOutput[row + 0];
		x0  = row0[ncl - 1]; x1  = row0[ncl + 0];
		x5  = row1[ncl - 1]; x6  = row1[ncl + 0];
		x10 = row2[ncl - 1]; x11 = row2[ncl + 0];
		for (col = ncl; col < nch + 1 - r; col += order) {
			x0  = row0[col - 1]; x1  = row0[col + 0]; x2  = row0[col + 1]; x3  = row0[col + 2]; x4  = row0[col + 3];
			x5  = row1[col - 1]; x6  = row1[col + 0]; x7  = row1[col + 1]; x8  = row1[col + 2]; x9  = row1[col + 3];
			x10 = row2[col - 1]; x11 = row2[col + 0]; x12 = row2[col + 1]; x13 = row2[col + 2]; x14 = row2[col + 3];
			
			y0 = x1  | x2  ;
			y1 = x6  | x7  ;
			y2 = x11 | x12 ;

			out_row0[col + 0] = (x0  | y0) | 
								(x5  | y1) | 
								(x10 | y2);

			out_row0[col + 1] = (y0 | x3) | 
								(y1 | x8) | 
								(y2 | x13);
								
			out_row0[col + 2] = (x2  | x3  | x4) | 
								(x7  | x8  | x9) | 
								(x12 | x13 | x14); 

			// x0  = row0[col - 1]; x1  = row0[col + 0]; x2  = row0[col + 1]; x3  = row0[col + 2]; x4  = row0[col + 3];
			// x5  = row1[col - 1]; x6  = row1[col + 0]; x7  = row1[col + 1]; x8  = row1[col + 2]; x9  = row1[col + 3];
			// x10 = row2[col - 1]; x11 = row2[col + 0]; x12 = row2[col + 1]; x13 = row2[col + 2]; x14 = row2[col + 3];
			

			// out_row0[col + 0] = (x0 | x1 | x2) | (x5 | x6 | x7) | (x10 | x11 | x12);
			// out_row0[col + 1] = (x1 | x2 | x3) | (x6 | x7 | x8) | (x11 | x12 | x13);
			// out_row0[col + 2] = (x2 | x3 | x4) | (x7 | x8 | x9) | (x12 | x13 | x14); 
			x0  = x3;  x1  = x4;
			x5  = x8;  x6  = x9;
			x10 = x13; x11 = x14;
			
		}
		row0 = row1;
		row1 = row2;
	}

	switch(r) {
		case 2:

		
		for (row = nrl; row < nrh + 1; row++) {
			row0 = ppInput[row - 1];
			row1 = ppInput[row + 0];
			row2 = ppInput[row + 1];

			out_row0 = ppOutput[row + 0];
		
			// y1 = x1 | x2 | x6 | x7 | x11 | x12; // 5
			x0  = row0[nch - 2]; x1  = row0[nch - 1]; x2  = row0[nch + 0]; x3  = row0[nch + 1]; 
			x5  = row1[nch - 2]; x6  = row1[nch - 1]; x7  = row1[nch + 0]; x8  = row1[nch + 1]; 
			x10 = row2[nch - 2]; x11 = row2[nch - 1]; x12 = row2[nch + 0]; x13 = row2[nch + 1]; 

			out_row0[nch - 1] = (x0  | x1  | x2) | 
								(x5  | x6  | x7) | 
								(x10 | x11 | x12);

			out_row0[nch + 0] = (x1  | x2  | x3) | 
			  					(x6  | x7  | x8) | 
								(x11 | x12 | x13);
			
		}
		break;
	case 1:
		for (row = nrl; row < nrh + 1; row++) {
			row0 = ppInput[row - 1];
			row1 = ppInput[row + 0];
			row2 = ppInput[row + 1];

			out_row0 = ppOutput[row + 0];
		
			// x1  = row0[nch - 1] | x2  = row0[nch + 0] | x3  = row0[nch + 1]; 
			// x6  = row1[nch - 1] | x7  = row1[nch + 0] | x8  = row1[nch + 1]; 
			// x11 = row2[nch - 1] | x12 = row2[nch + 0] | x13 = row2[nch + 1]; 

			out_row0[nch + 0] = (row0[nch - 1] |  row0[nch + 0] | row0[nch + 1])| 
			  					(row1[nch - 1] |  row1[nch + 0] | row1[nch + 1])| 
								(row2[nch - 1] |  row2[nch + 0] | row2[nch + 1]);
			
		}
		break;
	default:
		break;
	}
}

void ui8matrix_dilation_LU3x3_O3xO3_RR(uint8** ppInput, long nrl, long nrh, long ncl, long nch, p_struct_elem_dim s, uint8 **ppOutput)
{
	const long order = 3;
	long row = nrl, col = ncl, x, y, rr, cr;
	// uint8 y00, y10, y20, y30, y40, y01, y11, y21, y31, y41, y02, y12, y22, y32, y42;
	uint8 y0, y1, y2, y3, y4, y5, y6, y7, y8, y9, y10, y11, y12, y13, y14;
	uint8 x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16, x17, x18, x19, x20, x21, x22, x23, x24, x25;
	uint8 r0c0, r1c0, r2c0, r0c1,r1c1, r2c1, r0c2, r1c2, r2c2;
	uint8 *row0, *row1, *row2, *row3, *row4;
	uint8 *out_row0, *out_row1, *out_row2;

	rr = (nrh + 1) % order;
	cr = (nch + 1) % order;
	
	row0 = ppInput[row - 1];
	row1 = ppInput[row + 0];

	for (row = nrl; row < nrh + 1 - rr; row += order) {
		row2 = ppInput[row + 1];
		row3 = ppInput[row + 2];
		row4 = ppInput[row + 3];
		out_row0 = ppOutput[row + 0];
		out_row1 = ppOutput[row + 1];
		out_row2 = ppOutput[row + 2];


		x0  = row0[ncl - 1]; x1  = row0[ncl + 0]; x2  = row0[ncl + 1]; x3  = row0[ncl + 2]; x4  = row0[ncl + 3];
		x5  = row1[ncl - 1]; x6  = row1[ncl + 0]; x7  = row1[ncl + 1]; x8  = row1[ncl + 2]; x9  = row1[ncl + 3];
		x10 = row2[ncl - 1]; x11 = row2[ncl + 0]; x12 = row2[ncl + 1]; x13 = row2[ncl + 2]; x14 = row2[ncl + 3];
		x15 = row3[ncl - 1]; x16 = row3[ncl + 0]; x17 = row3[ncl + 1]; x18 = row3[ncl + 2]; x19 = row3[ncl + 3];
		x20 = row4[ncl - 1]; x21 = row4[ncl + 0]; x22 = row4[ncl + 1]; x23 = row4[ncl + 2]; x24 = row4[ncl + 3];
		y0  = x0  | x5  | x10; y1  = x1  | x6  | x11;
		y5  = x5  | x10 | x15; y6  = x6  | x11 | x16;
		y10 = x10 | x15 | x20; y11 = x11 | x16 | x21;

		for (col = ncl; col < nch + 1 - cr; col += order) {
			// x2  = row0[col + 1]; x3  = row0[col + 2]; x4  = row0[col + 3];
			// x7  = row1[col + 1]; x8  = row1[col + 2]; x9  = row1[col + 3];
			// x12 = row2[col + 1]; x13 = row2[col + 2]; x14 = row2[col + 3];
			// x17 = row3[col + 1]; x18 = row3[col + 2]; x19 = row3[col + 3];
			// x22 = row4[col + 1]; x23 = row4[col + 2]; x24 = row4[col + 3];
			x0  = row0[col - 1]; x1  = row0[col + 0]; x2  = row0[col + 1]; x3  = row0[col + 2]; x4  = row0[col + 3];
			x5  = row1[col - 1]; x6  = row1[col + 0]; x7  = row1[col + 1]; x8  = row1[col + 2]; x9  = row1[col + 3];
			x10 = row2[col - 1]; x11 = row2[col + 0]; x12 = row2[col + 1]; x13 = row2[col + 2]; x14 = row2[col + 3];
			x15 = row3[col - 1]; x16 = row3[col + 0]; x17 = row3[col + 1]; x18 = row3[col + 2]; x19 = row3[col + 3];
			x20 = row4[col - 1]; x21 = row4[col + 0]; x22 = row4[col + 1]; x23 = row4[col + 2]; x24 = row4[col + 3];


			y2  = x2  | x7  | x12; y3  = x3  | x8  | x13; y4  = x4  | x9  | x14; 
			y7  = x7  | x12 | x17; y8  = x8  | x13 | x18; y9  = x9  | x14 | x19; 
			y12 = x12 | x17 | x22; y13 = x13 | x18 | x23; y14 = x14 | x19 | x24; 
 
			out_row0[col + 0] = y0 | y1 | y2; 
			out_row0[col + 1] = y1 | y2 | y3; 
			out_row0[col + 2] = y2 | y3 | y4;
			y0 = y3; y1 = y4;

			out_row1[col + 0] = y5 | y6 | y7;
			out_row1[col + 1] = y6 | y7 | y8;
			out_row1[col + 2] = y7 | y8 | y9;
			y5 = y8; y6 = y9;
			
			out_row2[col + 0] = y10 | y11 | y12;
			out_row2[col + 1] = y11 | y12 | y13;
			out_row2[col + 2] = y12 | y13 | y14;
			y10 = y13; y11 = y14;
		}
		row0 = row3;
		row1 = row4;
	}
	switch (rr) {
		case 2 :
			switch(cr){
				case 2 :
					ui8matrix_dilation_LU3x3_O3xO1_RR(ppInput, nrl, nrh, nch - 1, nch, s, ppOutput);
					ui8matrix_dilation_LU3x3_O1xO3_RR(ppInput, nrh - 1, nrh, ncl, nch, s, ppOutput);
				break;
				case 1: 
					ui8matrix_dilation_LU3x3_O3xO1_RR(ppInput, nrl, nrh, nch, nch, s, ppOutput);		
					ui8matrix_dilation_LU3x3_O1xO3_RR(ppInput, nrh - 1, nrh, ncl, nch, s, ppOutput);	
				break;
				default :
					ui8matrix_dilation_LU3x3_O1xO3_RR(ppInput, nrh - 1, nrh, ncl, nch, s, ppOutput);
				break;
			}
			break;
		case 1:
			switch(cr){
				case 2 :
					ui8matrix_dilation_LU3x3_O3xO1_RR(ppInput, nrl, nrh, nch - 1, nch, s, ppOutput);
					ui8matrix_dilation_LU3x3_O1xO3_RR(ppInput, nrh, nrh, ncl, nch, s, ppOutput);
				break;
				case 1: 
					ui8matrix_dilation_LU3x3_O3xO1_RR(ppInput, nrl, nrh, nch, nch, s, ppOutput);		
					ui8matrix_dilation_LU3x3_O1xO3_RR(ppInput, nrh, nrh, ncl, nch, s, ppOutput);	
					// display_ui8matrix(ppOutput, nrl, nrh, ncl, nch, "%u", "rr=1; cr=1");
				break;
				default :
					ui8matrix_dilation_LU3x3_O1xO3_RR(ppInput, nrh, nrh, ncl, nch, s, ppOutput);
				break;
			}
		break;
		default:
			switch(cr){
				case 2 :
					ui8matrix_dilation_LU3x3_O3xO1_RR(ppInput, nrl, nrh, nch - 1, nch, s, ppOutput);
				break;
				case 1: 
					ui8matrix_dilation_LU3x3_O3xO1_RR(ppInput, nrl, nrh, nch, nch, s, ppOutput);	
				break;
				default :
				break;
			}
		break;
	}
		
}

void ui8matrix_dilation_pipeline_LU3x3_O1xO3_RR (uint8** ppInput, long nrl, long nrh, long ncl, long nch, p_struct_elem_dim s, uint8 **ppOutput)
{
	const long order = 3;
	long row = nrl, col = ncl, x, y, r = 0;
	uint8 **temp = ui8matrix(nrl, nrh, ncl + s->ncl, nch + s->nch);
	uint8 *row0, *row1, *row2, *row3, x0, x1, x2, x3, x4, x5;
	uint8 *temp_row0, *temp_row1, *temp_row2, *temp_row3, *temp_row4;
	uint8 *out_row0, *out_row1, *out_row2, *out_row3;
	
	r = (nch + 1) % order;
	// Prologue	
	row0 = ppInput[nrl - 1];
	row1 = ppInput[nrl + 0];
	for (row = nrl; row < nrh + 1; row++) {
		row2 = ppInput[row + 1];
		temp_row1 = temp[row + 0];

		temp_row1[ncl - 1] = row0[ncl - 1] | row1[ncl - 1] | row2[ncl - 1];
		temp_row1[ncl + 0] = row0[ncl + 0] | row1[ncl + 0] | row2[ncl + 0];
		row0 = row1;
		row1 = row2;
	}

	row0 = ppInput[nrl - 1];
	row1 = ppInput[nrl + 0];
	for (row = nrl; row < nrh + 1; row ++){
		temp_row1 = temp[row + 0];
		row2 = ppInput[row + 1];

		out_row1 = ppOutput[row]; 
		x0 = temp_row1[ncl - 1];		// LOAD => A
		x1 = temp_row1[ncl - 0];		// LOAD => B
		
		for (col = ncl; col < nch + 1 - r; col += order) {
			x2 = row0[col + 1] | row1[col + 1] | row2[col + 1];
			x3 = row0[col + 2] | row1[col + 2] | row2[col + 2];
			x4 = row0[col + 3] | row1[col + 3] | row2[col + 3];
			
			out_row1 [col + 0] = x0 | x1 | x2; // A | B | C
			out_row1 [col + 1] = x1 | x2 | x3; // B | C | D
			out_row1 [col + 2] = x2 | x3 | x4; // C | D | E

			x0 = x3; 						  // STORE => D
			x1 = x4; 						  // STORE => E 
											  // NEXT  :  F 
											  // NEXT  :  G 
											  // NEXT  :  H
		}
		row0 = row1;
		row1 = row2;
	}
	switch (r) {
		case 2: 
			row0 = ppInput[nrl - 1];
			row1 = ppInput[nrl + 0];
			for (row = nrl; row < nrh + 1; row++) {
				row2 = ppInput[row + 1];
				temp_row1 = temp[row + 0];
				out_row1 = ppOutput[row]; 

				x1 = row0[nch - 1] | row1[nch - 1] | row2[nch - 1]|row0[nch + 0] | row1[nch + 0] | row2[nch + 0];
				// x2 = ;

				out_row1 [nch - 1] = row0[nch - 2] | row1[nch - 2] | row2[nch - 2] | x1 ; // A | B | C
				out_row1 [nch + 0] = x1 | row0[nch + 1] | row1[nch + 0] | row2[nch + 0]; // B | C | D
				row0 = row1;
				row1 = row2;
			}
			break;
		case 1: 
			row0 = ppInput[nrl - 1];
			row1 = ppInput[nrl + 0];
			for (row = nrl; row < nrh + 1; row++) {
				row2 = ppInput[row + 1];
				// temp_row1 = temp[row + 0];
				ppOutput[row][nch] = row0[nch - 1] | row1[nch - 1] | row2[nch - 1] | row0[nch - 0] | row1[nch - 0] | row2[nch - 0] | row0[nch + 1] | row1[nch + 1] | row2[nch + 1]; // A | B | C
				// out_row1 [nch + 0] = x1 | x2 | row0[nch + 1] | row1[nch + 0] | row2[nch + 0]; // B | C | D
				row0 = row1;
				row1 = row2;
			}
			break;
		default:
			break;
	}
	// display_ui8matrix(ppOutput, nrl - 1, nrh + 1, ncl - 1, nch + 1, "%u", "TEst");
	free_ui8matrix(temp, nrl, nrh, ncl + s->ncl, nch + s->nch);
	
}

void ui8matrix_dilation_LU3x3_O1xO3 (uint8** ppInput, long nrl, long nrh, long ncl, long nch, p_struct_elem_dim s, uint8 **ppOutput)
{
	const long order = 3;
	long row = nrl, col = ncl, x, y, r;
	uint8 y0, y1, y2, x0, x1, x2, x3, x4, x5, x6, x7, x8;
	uint8 *row0, *row1, *row2, *row3, *row4;
	uint8 *out_row0;
	r = (nch + 1)  % order;
	
	for (row = nrl; row < nrh + 1; row ++) {
		row0 = ppInput[row - 1];
		row1 = ppInput[row + 0];
		row2 = ppInput[row + 1];
		out_row0 = ppOutput[row + 0];
		for (col = ncl; col < nch + 1 - r; col += order) {
			// x0 = scalar_or3(row0, col); x3 = scalar_or3(row0, col + 1); x6 = scalar_or3(row0, col + 2);
			// x1 = scalar_or3(row1, col); x4 = scalar_or3(row1, col + 1); x7 = scalar_or3(row1, col + 2);
			// x2 = scalar_or3(row2, col); x5 = scalar_or3(row2, col + 1); x8 = scalar_or3(row2, col + 2);

			out_row0[col + 0] = scalar_or3(row0, col + 0)|
								scalar_or3(row1, col + 0)|
								scalar_or3(row2, col + 0);
			out_row0[col + 1] = scalar_or3(row0, col + 1)|
 								scalar_or3(row1, col + 1)|
 								scalar_or3(row2, col + 1);
			out_row0[col + 2] = scalar_or3(row0, col + 2)|
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
				out_row0 = ppOutput[row + 0];

				out_row0[nch - 1] = scalar_or3(row0, nch - 1)|
 									scalar_or3(row1, nch - 1)|
 									scalar_or3(row2, nch - 1);
				out_row0[nch + 0] = scalar_or3(row0, nch + 0)|
									scalar_or3(row1, nch + 0)|
									scalar_or3(row2, nch + 0);
			}
			break;
		case 1:
			for (row = nrl; row < nrh + 1; row ++) {
				ppOutput[row + 0][nch + 0] = scalar_or3(ppInput[row - 1], nch) | 
										  	 scalar_or3(ppInput[row + 0], nch) | 
										  	 scalar_or3(ppInput[row + 1], nch);
			}
			// ui8matrix_dilation_LU3x3_O1xO1(ppInput, nrl, nrh, nch, nch, s, ppOutput);
			break;
		default:
			break;
	}
}



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
void ui8matrix_dilation_LU3x3_O3xO1_ver1 (uint8** ppInput, long nrl, long nrh, long ncl, long nch, p_struct_elem_dim s, uint8 **ppOutput)
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
			break;
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
void ui8matrix_dilation_LU3x3_O1xO3_ver1 (uint8** ppInput, long nrl, long nrh, long ncl, long nch, p_struct_elem_dim s, uint8 **ppOutput)
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
void ui8matrix_dilation_LU3x3_O3xO3_ver1 (uint8** ppInput, long nrl, long nrh, long ncl, long nch, p_struct_elem_dim s, uint8 **ppOutput)
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
		break;
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
void ui8matrix_dilation_LU3x3_O3xO1_RR_ver1(uint8** ppInput, long nrl, long nrh, long ncl, long nch, p_struct_elem_dim s, uint8 **ppOutput)
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
			break;
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

void ui8matrix_dilation_LU3x3_O1xO3_RR_ver1 (uint8** ppInput, long nrl, long nrh, long ncl, long nch, p_struct_elem_dim s, uint8 **ppOutput)
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

		x0  = row0[ncl -1]; x1  = row0[ncl + 0]; 
		x5  = row1[ncl -1]; x6  = row1[ncl + 0]; 
		x10 = row2[ncl -1]; x11 = row2[ncl + 0]; 
		for (col = ncl; col < nch + 1 - r; col += order) {
			x2  = row0[col + 1]; x3  = row0[col + 2]; x4  = row0[col + 3];
			x7  = row1[col + 1]; x8  = row1[col + 2]; x9  = row1[col + 3];
			x12 = row2[col + 1]; x13 = row2[col + 2]; x14 = row2[col + 3];
			

			y1 = x1 | x2 | x6 | x7 | x11 | x12; // 5 calcs
			ppOutput[row][col + 0] = (x0 		  ) | (x5          ) | (x10            ) | y1; // 3 calc	
			ppOutput[row][col + 1] = (          x3) | (          x8) | (            x13) | y1; // 3 calc	
			ppOutput[row][col + 2] = (x2 | x3 | x4) | (x7 | x8 | x9) | (x12 | x13 | x14); //8 calc 

			// 8 + 6 + 5 = 21 calcs
			x0  = x3;  x1  = x4;
			x5  = x8;  x6  = x9;
			x10 = x13; x11 = x14;
			
		}
		row0 = row1;
		row1 = row2;
	}

	row0 = ppInput[nrl - 1];
	row1 = ppInput[nrl + 0];
	
	switch(r) {
	case 2:
		// ui8matrix_dilation_LU3x3_O1xO3_RR(ppInput, nrl, nrh, nch - 1, nch, s, ppOutput);
		for (row = nrl; row < nrh + 1; row++) {
			row2 = ppInput[row + 1];

			x0  = row0[nch - 2]; x1  = row0[nch - 1];  x2  = row0[nch + 0]; x3  = row0[nch + 1]; 
			x5  = row1[nch - 2]; x6  = row1[nch - 1];  x7  = row1[nch + 0]; x8  = row1[nch + 1]; 
			x10 = row2[nch - 2]; x11 = row2[nch - 1];  x12 = row2[nch + 0]; x13 = row2[nch + 1];

			y1 = x1 | x2 | x6 | x7 | x11 | x12; // 5 calcs
			ppOutput[row][nch - 1] = (x0 		  ) | (x5          ) | (x10            ) | y1; // 3 calc	
			ppOutput[row][nch + 0] = (          x3) | (          x8) | (            x13) | y1; // 3 calc	
		
			row0 = row1;
			row1 = row2;
		}

		
		break;
	case 1:
		for (row = nrl; row < nrh + 1; row++) {
			row2 = ppInput[row + 1];

			x1  = row0[nch - 1];  x2  = row0[nch + 0]; x3  = row0[nch + 1]; 
			x6  = row1[nch - 1];  x7  = row1[nch + 0]; x8  = row1[nch + 1]; 
			x11 = row2[nch - 1];  x12 = row2[nch + 0]; x13 = row2[nch + 1];

			ppOutput[row][nch + 0] = (x1 | x2 | x3) | (x6 | x7 | x8) | (x11 | x12 | x13); // 3 calc

			row0 = row1;
			row1 = row2;
		}
		break;
	default:
		break;
	}
}

void ui8matrix_dilation_LU3x3_O3xO1 (uint8** ppInput, long nrl, long nrh, long ncl, long nch, p_struct_elem_dim s, uint8 **ppOutput)
{
	const long order = 3;
	long row = nrl, col = ncl, x, y;
	uint8 y0, y1, y2, x0, x1, x2, x3, r;
	uint8 *row0, *row1, *row2, *row3, *row4;
	uint8 *out_row0, *out_row1, *out_row2;

	r = (nrh + 1)  % order;
	
	for (row = nrl; row < nrh + 1 - r; row += order) {
		row0 = ppInput[row - 1];
		row1 = ppInput[row + 0];
		row2 = ppInput[row + 1];
		row3 = ppInput[row + 2];
		row4 = ppInput[row + 3];

		out_row0 = ppOutput[row + 0];
		out_row1 = ppOutput[row + 1];
		out_row2 = ppOutput[row + 2];
		
		for (col = ncl; col < nch + 1; col++) {
			x0 = scalar_or3(row1, col);
			x1 = scalar_or3(row2, col);
			x2 = scalar_or3(row3, col);

			out_row0[col] = scalar_or3(row0, col) | x0 | x1;
			out_row1[col] = 			       x0 | x1 | x2;
			out_row2[col] = 			       x1 | x2 | scalar_or3(row4, col);
		}
	}
	
	switch(r) {
		case 2:
			row0 = ppInput[row - 1];
			row1 = ppInput[row + 0];
			row2 = ppInput[row + 1];
			row3 = ppInput[row + 2];
			out_row0 = ppOutput[row + 0];
			out_row1 = ppOutput[row + 1];
			for (col = ncl; col < nch + 1; col++) {
				x0 = scalar_or3(row1, col);
				x1 = scalar_or3(row2, col);
				out_row0[col] = scalar_or3(row0, col) | x0 | x1;
				out_row1[col] = 				     x0 | x1 | scalar_or3(row3, col);
			}
			break;
		case 1:
			row0 = ppInput[row - 1];
			row1 = ppInput[row + 0];
			row2 = ppInput[row + 1];
			out_row0 = ppOutput[row + 0];
			for (col = ncl; col < nch + 1; col++) {
				out_row0[col] = scalar_or3(row0, col) | 
										 scalar_or3(row1, col) |
										 scalar_or3(row2, col);
				
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
		out_row0 = ppOutput[row + 0];
		out_row1 = ppOutput[row + 1];
		out_row2 = ppOutput[row + 2];
		for (col = ncl; col < nch + 1 - cr; col += order) {
			x0 = scalar_or3(row1, col); x1 = scalar_or3(row1, col + 1); x2 = scalar_or3(row1, col + 2);
			x3 = scalar_or3(row2, col); x4 = scalar_or3(row2, col + 1); x5 = scalar_or3(row2, col + 2);
			x6 = scalar_or3(row3, col); x7 = scalar_or3(row3, col + 1); x8 = scalar_or3(row3, col + 2);

			out_row0[col + 0] = scalar_or3(row0, col + 0) | x0 | x3;
			out_row1[col + 0] = 				         x0 | x3 | x6;
			out_row2[col + 0] = 				         x3 | x6 | scalar_or3(row4, col + 0);


			out_row0[col + 1] = scalar_or3(row0, col + 1) | x1 | x4;
			out_row1[col + 1] = 				         x1 | x4 | x7;
			out_row2[col + 1] = 				         x4 | x7 | scalar_or3(row4, col + 1);


			out_row0[col + 2] = scalar_or3(row0, col + 2) | x2 | x5;
			out_row1[col + 2] = 				         x2 | x5 | x8;
			out_row2[col + 2] = 				         x5 | x8 | scalar_or3(row4, col + 2);

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
			break;
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
void ui8matrix_dilation_LU3x3_O3xO3_RR_ver1(uint8** ppInput, long nrl, long nrh, long ncl, long nch, p_struct_elem_dim s, uint8 **ppOutput)
{
	const long order = 3;
	long row = nrl, col = ncl, x, y, rr, cr;
	// uint8 y00, y10, y20, y30, y40, y01, y11, y21, y31, y41, y02, y12, y22, y32, y42;
	uint8 y0, y1, y2, y3, y4, y5, y6, y7, y8, y9, y10, y11, y12, y13, y14;
	uint8 x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16, x17, x18, x19, x20, x21, x22, x23, x24, x25;
	uint8 r0c0, r1c0, r2c0, r0c1,r1c1, r2c1, r0c2, r1c2, r2c2;
	uint8 *row0, *row1, *row2, *row3, *row4;
	uint8 *out_row0, *out_row1, *out_row2;

	rr = (nrh + 1) % order;
	cr = (nch + 1) % order;

	row0 = ppInput[nrl - 1]; x0  = row0[ncl - 1]; x1  = row0[ncl + 0];
	row1 = ppInput[nrl + 0]; x5  = row1[ncl - 1]; x6  = row1[ncl + 0];
		
	for (col = ncl; col < nch + 1 - cr; col += order) {
		row0 = ppInput[nrl - 1];
		row1 = ppInput[nrl + 0]; 
		
		x2  = row0[col + 1]; x3  = row0[col + 2]; x4  = row0[col + 3];
		x7  = row1[col + 1]; x8  = row1[col + 2]; x9  = row1[col + 3];
		y0  = x0  | x1  | x2 ; y1  = x1  | x2  | x3 ; y2  = x2  | x3  | x4 ; 
		y3  = x5  | x6  | x7 ; y4  = x6  | x7  | x8 ; y5  = x7  | x8  | x9 ;
			
		for (row = nrl; row < nrh + 1 - rr; row += order) {
			row2 = ppInput[row + 1];
			row3 = ppInput[row + 2];
			row4 = ppInput[row + 3];
			out_row0 = ppOutput[row + 0];
			out_row1 = ppOutput[row + 1];
			out_row2 = ppOutput[row + 2];

			x10 = row2[col - 1]; x11 = row2[col + 0]; x12 = row2[col + 1]; x13 = row2[col + 2]; x14 = row2[col + 3];
			x15 = row3[col - 1]; x16 = row3[col + 0]; x17 = row3[col + 1]; x18 = row3[col + 2]; x19 = row3[col + 3];
			x20 = row4[col - 1]; x21 = row4[col + 0]; x22 = row4[col + 1]; x23 = row4[col + 2]; x24 = row4[col + 3];

			y6  = x10 | x11 | x12; y7  = x11 | x12 | x13; y8  = x12 | x13 | x14;
			y9  = x15 | x16 | x17; y10 = x16 | x17 | x18; y11 = x17 | x18 | x19;
			y12 = x20 | x21 | x22; y13 = x21 | x22 | x23; y14 = x22 | x23 | x24;

			out_row0[col + 0] = y0 | y3 | y6;
			out_row0[col + 1] = y1 | y4 | y7;
			out_row0[col + 2] = y2 | y5 | y8;

			out_row1[col + 0] = y3 | y6 | y9 ;
			out_row1[col + 1] = y4 | y7 | y10;
			out_row1[col + 2] = y5 | y8 | y11;

			out_row2[col + 0] = y6 | y9  | y12;
			out_row2[col + 1] = y7 | y10 | y13;
			out_row2[col + 2] = y8 | y11 | y14;
			
			y0  = y9 ; y1  = y10; y2  = y11; 
			y3  = y12; y4  = y13; y5  = y14;

			row0 = row3;
			row1 = row4;
		}
		x0 = x3; x1 = x4;
		x5 = x8; x6 = x9;
	}
	switch (rr) {
		case 2 :
			switch(cr){
				case 2 :
					ui8matrix_dilation_LU3x3_O3xO1_RR(ppInput, nrl, nrh, nch - 1, nch, s, ppOutput);
					ui8matrix_dilation_LU3x3_O1xO3_RR(ppInput, nrh - 1, nrh, ncl, nch, s, ppOutput);
				break;
				case 1: 
					ui8matrix_dilation_LU3x3_O3xO1_RR(ppInput, nrl, nrh, nch, nch, s, ppOutput);		
					ui8matrix_dilation_LU3x3_O1xO3_RR(ppInput, nrh - 1, nrh, ncl, nch, s, ppOutput);	
				break;
				default :
					ui8matrix_dilation_LU3x3_O1xO3_RR(ppInput, nrh - 1, nrh, ncl, nch, s, ppOutput);
				break;
			}
			break;
		case 1:
			switch(cr){
				case 2 :
					ui8matrix_dilation_LU3x3_O3xO1_RR(ppInput, nrl, nrh, nch - 1, nch, s, ppOutput);
					ui8matrix_dilation_LU3x3_O1xO3_RR(ppInput, nrh, nrh, ncl, nch, s, ppOutput);
				break;
				case 1: 
					ui8matrix_dilation_LU3x3_O3xO1_RR(ppInput, nrl, nrh, nch, nch, s, ppOutput);	
					ui8matrix_dilation_LU3x3_O1xO3_RR(ppInput, nrh, nrh, ncl, nch, s, ppOutput);	
				break;
				default :
					ui8matrix_dilation_LU3x3_O1xO3_RR(ppInput, nrh, nrh, ncl, nch, s, ppOutput);
				break;
			}
		break;
		default:
			switch(cr){
				case 2 :
					ui8matrix_dilation_LU3x3_O3xO1_RR(ppInput, nrl, nrh, nch - 1, nch, s, ppOutput);
				break;
				case 1: 
					ui8matrix_dilation_LU3x3_O3xO1_RR(ppInput, nrl, nrh, nch, nch, s, ppOutput);	
				break;
				default :
				break;
			}
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
	uint8 *out_row0, *out_row1, *out_row2;

	r = (nrh + 1)  % order;
	
	row0 = ppInput[row - 1];
	row1 = ppInput[row + 0];
	for (row = nrl; row < nrh + 1 - r; row += order) {
		row2 = ppInput[row + 1]; 
		row3 = ppInput[row + 2]; 
		row4 = ppInput[row + 3];
		x0 = row1[ncl - 1]; x3 = row1[ncl + 0];
		x1 = row2[ncl - 1]; x4 = row2[ncl + 0];
		x2 = row3[ncl - 1]; x5 = row3[ncl + 0];
		out_row0 = ppOutput[row + 0];
		out_row1 = ppOutput[row + 1];
		out_row2 = ppOutput[row + 2];
		for (col = ncl; col < nch + 1; col++) {
			x6 = row1[col + 1]; //scalar_or3(row1, col);
			x7 = row2[col + 1]; //scalar_or3(row2, col);
			x8 = row3[col + 1]; //scalar_or3(row3, col);
			y0 = (x0 | x3 | x6); // 2 calc
			y1 = (x1 | x4 | x7); // 2 calc
			y2 = (x2 | x5 | x8); // 2 calc

			out_row0[col] = scalar_or3(row0, col) | y0 | y1; 						// 3 + 2    =  5
			out_row1[col] =         			 y0 | y1 | y2;						//     2    =  2
			out_row2[col] = 		 			 y1 | y2 | scalar_or3(row4, col);  //     2 + 3=  5
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
	out_row0 = ppOutput[nrh - 1];
	out_row1 = ppOutput[nrh + 0];
	switch(r) {
		case 2:
			x0 = row1[ncl - 1]; x3 = row1[ncl + 0];
			x1 = row2[ncl - 1]; x4 = row2[ncl + 0];
			x2 = row3[ncl - 1]; x5 = row3[ncl + 0];
			for (col = ncl; col < nch + 1; col++) {
				x6 = row1[col + 1]; //scalar_or3(row1, col);
				x7 = row2[col + 1]; //scalar_or3(row2, col);
				x8 = row3[col + 1]; //scalar_or3(row3, col);
				out_row0[col] = scalar_or3(row0, col) | (x0 | x3 | x6) | (x1 | x4 | x7);;
				out_row1[col] =         (x0 | x3 | x6) | (x1 | x4 | x7) | (x2 | x5 | x8);;
				x0 = x3; x3 = x6;
				x1 = x4; x4 = x7;
				x2 = x5; x5 = x8;
			}
			break;
		case 1:
			x0 = row1[ncl - 1]; x3 = row1[ncl + 0];
			x1 = row2[ncl - 1]; x4 = row2[ncl + 0];
			for (col = ncl; col < nch + 1; col++) {
				x6 = row1[col + 1]; //scalar_or3(row1, col);
				x7 = row2[col + 1]; //scalar_or3(row2, col);
				out_row1[col] = scalar_or3(row0, col) | (x0 | x3 | x6) | (x1 | x4 | x7);
				x0 = x3; x3 = x6;
				x1 = x4; x4 = x7;
			}				
			break;
		default:
			break;
	}
}


void ui8matrix_dilation_LU3x3_O1xO3_PJ (uint8** ppInput, long nrl, long nrh, long ncl, long nch, p_struct_elem_dim s, uint8 **ppOutput)
{
	const long order = 3;
	long row = nrl, col = ncl, r;
	uint8 y0, y1, y2;
	uint8 x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14;
	uint8 *row0, *row1, *row2, *row3, *row4;
	uint8 *out_row0;

	r = (nch + 1) % order;
	row0 = ppInput[nrl - 1];
	row1 = ppInput[nrl + 0];
	for (row = nrl; row < nrh + 1; row++) {
		row2 = ppInput[row + 1];

		out_row0 = ppOutput[row + 0];
		x0  = row0[ncl - 1]; x1  = row0[ncl + 0];
		x5  = row1[ncl - 1]; x6  = row1[ncl + 0];
		x10 = row2[ncl - 1]; x11 = row2[ncl + 0];
		for (col = ncl; col < nch + 1 - r; col += order) {
			x2  = row0[col + 1]; x3  = row0[col + 2]; x4  = row0[col + 3];
			x7  = row1[col + 1]; x8  = row1[col + 2]; x9  = row1[col + 3];
			x12 = row2[col + 1]; x13 = row2[col + 2]; x14 = row2[col + 3];
			
			out_row0[col + 0] = (x0 | x1 | x2) | (x5 | x6 | x7) | (x10 | x11 | x12);
			out_row0[col + 1] = (x1 | x2 | x3) | (x6 | x7 | x8) | (x11 | x12 | x13);
			out_row0[col + 2] = (x2 | x3 | x4) | (x7 | x8 | x9) | (x12 | x13 | x14); 

			x0  = x3;  x1  = x4;
			x5  = x8;  x6  = x9;
			x10 = x13; x11 = x14;
			
		}
		row0 = row1;
		row1 = row2;
	}

	switch(r) {
		case 2:

		for (row = nrl; row < nrh + 1; row++) {
			row0 = ppInput[row - 1];
			row1 = ppInput[row + 0];
			row2 = ppInput[row + 1];

			out_row0 = ppOutput[row + 0];
		
			// y1 = x1 | x2 | x6 | x7 | x11 | x12; // 5
			x0  = row0[nch - 2]; x1  = row0[nch - 1]; x2  = row0[nch + 0]; x3  = row0[nch + 1]; 
			x5  = row1[nch - 2]; x6  = row1[nch - 1]; x7  = row1[nch + 0]; x8  = row1[nch + 1]; 
			x10 = row2[nch - 2]; x11 = row2[nch - 1]; x12 = row2[nch + 0]; x13 = row2[nch + 1]; 

			out_row0[nch - 1] = (x0  | x1  | x2) | 
								(x5  | x6  | x7) | 
								(x10 | x11 | x12);

			out_row0[nch + 0] = (x1  | x2  | x3) | 
			  					(x6  | x7  | x8) | 
								(x11 | x12 | x13);
			
		}
		break;
	case 1:
		for (row = nrl; row < nrh + 1; row++) {
			row0 = ppInput[row - 1];
			row1 = ppInput[row + 0];
			row2 = ppInput[row + 1];

			out_row0 = ppOutput[row + 0];
		
			// x1  = row0[nch - 1] | x2  = row0[nch + 0] | x3  = row0[nch + 1]; 
			// x6  = row1[nch - 1] | x7  = row1[nch + 0] | x8  = row1[nch + 1]; 
			// x11 = row2[nch - 1] | x12 = row2[nch + 0] | x13 = row2[nch + 1]; 

			out_row0[nch + 0] = (row0[nch - 1] |  row0[nch + 0] | row0[nch + 1])| 
			  					(row1[nch - 1] |  row1[nch + 0] | row1[nch + 1])| 
								(row2[nch - 1] |  row2[nch + 0] | row2[nch + 1]);
			
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
	uint8 *out_row0;
	// Prologue	
	for (col = ncl; col < nch + 1; col++) {
    	temp[row - 1][col] = scalar_or3(ppInput[row - 1], col);
		temp[row + 0][col] = scalar_or3(ppInput[row + 0], col);
	}	

	for (row = row + 1; row < nrh + 1; row++){
		out_row0 = ppOutput[row - 1];
		for (col = ncl; col < nch + 1; col++) {
			temp[row + 0][col] = scalar_or3(ppInput[row + 0], col);
			out_row0[col] = temp[row - 2][ col] |
									 temp[row - 1][ col] |
									 temp[row - 0][ col];
		}
	}
	out_row0 = ppOutput[nrh + 0];
	for (col = ncl; col < nch + 1; col++) {

    	out_row0[col] =    temp[nrh - 1][ col] |
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
	for (col = ncl; col < nch + 1; col++) {
		temp_row1[col] = scalar_or3(row0, col);
		temp_row2[col] = scalar_or3(row1, col);
	}

	for (row = row + 1; row < nrh + 1 - r; row += order){
		row0 = ppInput[row - 1];   //  0 : b0 b1 b2 ...
		row1 = ppInput[row + 0];   //  1 : c0 c1 c2 ...
		row2 = ppInput[row + 1];   //  2 : d0 d1 d2 ...
		row3 = ppInput[row + 2];   //  3 : e0 e1 e2 ...

		temp_row0 = temp[row - 2]; // -1 : A    <= (a0 a1 a2 ...) Prologued
		temp_row1 = temp[row - 1]; //  0 : B    <= (b0 b1 b2 ...) Prologued
		temp_row2 = temp[row + 0]; //  1 : NULL <= (c0 c1 c2 ...)
		temp_row3 = temp[row + 1]; //  2 : NULL <= (d0 d1 d2 ...)
		temp_row4 = temp[row + 2]; //  3 : NULL <= (e0 e1 e2 ...)
		out_row0 = ppOutput[row - 1]; 
		out_row1 = ppOutput[row + 0];
		out_row2 = ppOutput[row + 1];

		for (col = ncl; col < nch + 1; col++) {
			//  temp_row0[ col]			   LOAD => A										
			x0 = temp_row1[col];		// LOAD => B
			x1 = scalar_or3(row1, col); // OR3  => C
			x2 = scalar_or3(row2, col); // OR3  => D
			x3 = scalar_or3(row3, col); // OR3  => E
			out_row0[col] = temp_row0[ col] | x0 | x1; // A | B | C
			out_row1[col] = 			 x0 | x1 | x2; // B | C | D
			out_row2[col] = 			 x1 | x2 | x3; // C | D | E
			temp_row2[col] = x1;		// STORE => C 
			temp_row3[col] = x2;		// STORE => D
			temp_row4[col] = x3;		// STORE => E 
										// NEXT  :  F 
										// NEXT  :  G 
										// NEXT  :  H
		}
	}
	
	
	switch (r) {
		case 2: 
			row1 = ppInput[nrh - 1];   // -1 : (f0 f1 f2 ...)
			row2 = ppInput[nrh + 0];   //  0 : (g0 g1 g2 ...)
			row3 = ppInput[nrh + 1];   //  1 : (h0 h1 h2 ...)
			temp_row0 = temp[nrh - 2]; //  E    <= (e0 e1 e2 ...)
			out_row0 = ppOutput[nrh - 1];
			out_row1 = ppOutput[nrh + 0];
			for (col = ncl; col < nch + 1; col++) {
				//  temp_row0[ col]		   		 LOAD => E		
				x0 = scalar_or3(row1, col);	  // OR3  => F	
				x1 = scalar_or3(row2, col);	  // OR3  => G
				out_row0[col] = temp_row0[ col] | x0 | x1; // A | B | C
				out_row1[col] = 			 x0 | x1 | scalar_or3(row3, col); // B | C | D	
			}
			break;
		case 1: 
			row2 = ppInput[nrh + 0];   //  0 : (f0 f1 f2 ...)
			row3 = ppInput[nrh + 1];   //  1 : (g0 g1 g2 ...)
			temp_row0 = temp[nrh - 1]; //  E    <= (e0 e1 e2 ...)
			out_row1 = ppOutput[nrh + 0];
			for (col = ncl; col < nch + 1; col++) {
				out_row1[col] = temp_row0[ col]  | scalar_or3(row2, col) | scalar_or3(row3, col); // E | F | G	
			}
			break;
		default:
			break;
	}
	free_ui8matrix(temp, nrl + s->nrl, nrh + s->nrh, ncl + s->ncl, nch + s->nch);
}
void ui8matrix_dilation_pipeline_LU3x3_O1xO3 (uint8** ppInput, long nrl, long nrh, long ncl, long nch, p_struct_elem_dim s, uint8 **ppOutput)
{
	const long order = 3;
	long row = nrl, col = ncl, x, y, r = 0;
	uint8 **temp = ui8matrix(nrl + s->nrl, nrh + s->nrh, ncl + s->ncl, nch + s->nch);
	uint8 *row0, *row1, *row2, *row3, x0, x1, x2, x3, x4, x5, y0, y1, y2;
	uint8 *temp_row0, *temp_row1, *temp_row2, *temp_row3, *temp_row4;
	uint8 *out_row0, *out_row1, *out_row2, *out_row3;
	
	r = (nch + 1) % order;
	// Prologue	
	row0 = ppInput[row - 1];
	row1 = ppInput[row + 0];
	row2 = ppInput[row + 1];
	
	// temp_row1 = temp[row - 1];
	// temp_row2 = temp[row + 0];
	// temp_row3 = temp[row + 1];
	for (row = nrl; row < nrh + 1; row++) {
		row0 = ppInput[row - 1];
		row1 = ppInput[row + 0];
		row2 = ppInput[row + 1];
		temp_row1 = temp[row + 0];

		temp_row1[ncl - 1] = row0[ncl - 1] | row1[ncl - 1] | row2[ncl - 1];
		temp_row1[ncl + 0] = row0[ncl + 0] | row1[ncl + 0] | row2[ncl + 0];
	}

	for (row = nrl; row < nrh + 1; row ++){
		temp_row1 = temp[row + 0];
		
		row0 = ppInput[row - 1];
		row1 = ppInput[row + 0];
		row2 = ppInput[row + 1];

		out_row1 = ppOutput[row]; 
		
		for (col = ncl; col < nch + 1 - r; col += order) {
			x0 = temp_row1[col - 1];		// LOAD => A
			x1 = temp_row1[col - 0];		// LOAD => B
			x2 = row0[col + 1] | row1[col + 1] | row2[col + 1];
			x3 = row0[col + 2] | row1[col + 2] | row2[col + 2];
			x4 = row0[col + 3] | row1[col + 3] | row2[col + 3];
			// x5 = row0[col + 3] | row1[col + 3] | row2[col + 3];
			y0 = x0 | x1 | x2;
			y1 = x1 | x2 | x3;
			y2 = x2 | x3 | x4;
			
			out_row1 [col + 0] = y0; // A | B | C
			out_row1 [col + 1] = y1; // B | C | D
			out_row1 [col + 2] = y2; // C | D | E

			// temp_row1[col + 2] = x2;		  // STORE => c
			temp_row1[col + 2] = x3;		  // STORE => D
			temp_row1[col + 3] = x4;		  // STORE => E 
										// NEXT  :  F 
										// NEXT  :  G 
										// NEXT  :  H
		}
	}
	switch (r) {
		case 2: 
			for (row = nrl; row < nrh + 1; row ++){
				temp_row1 = temp[row + 0];
				
				row0 = ppInput[row - 1];
				row1 = ppInput[row + 0];
				row2 = ppInput[row + 1];

				out_row0 = ppOutput[row]; 
				
				x0 = temp_row1[nch - 2];		// LOAD => E
				x1 = row0[nch - 1] | row1[nch - 1] | row2[nch - 1];
				x2 = row0[nch + 0] | row1[nch + 0] | row2[nch + 0];
				x3 = row0[nch + 1] | row1[nch + 1] | row2[nch + 1];
				
				y0 = x0 | x1 | x2;
				y1 = x1 | x2 | x3;
				// y2 = x2 | x3 | x4;
				out_row0[nch - 1] = y0; 	// A | B | C
				out_row0[nch + 0] = y1; 	// B | C | D
				// out_row0[nch + 1] = y2; 	// C | D | E				
			}
			// ui8matrix_dilation_pipeline3x3(ppInput, nrl, nrh, nch - 1, nch, s, ppOutput);
			break;
		case 1: 
			for (row = nrl; row < nrh + 1; row ++){
				temp_row1 = temp[row + 0];
				
				row0 = ppInput[row - 1];
				row1 = ppInput[row + 0];
				row2 = ppInput[row + 1];

				out_row0 = ppOutput[row]; 
				
				x0 = temp_row1[nch - 1];		// LOAD => E
				// x1 = row0[nch - 1] | row1[nch - 1] | row2[nch - 1];
				x2 = row0[nch + 0] | row1[nch + 0] | row2[nch + 0];
				x3 = row0[nch + 1] | row1[nch + 1] | row2[nch + 1];
				
				// y0 = x0 | x1 | x2;
				y1 = x0 | x2 | x3;
				// y2 = x2 | x3 | x4;
				// out_row0[nch - 1] = y0; 	// A | B | C
				out_row0[nch + 0] = y1; 	// B | C | D
				// out_row0[nch + 1] = y2; 	// C | D | E				
			}
		
			break;
		default:
			break;
	}
	// display_ui8matrix(ppOutput, nrl - 1, nrh + 1, ncl - 1, nch + 1, "%u", "TEst");
	free_ui8matrix(temp, nrl + s->nrl, nrh + s->nrh, ncl + s->ncl, nch + s->nch);
	
}


void ui8matrix_dilation_pipeline_LU3x3_O1xO3_PJ (uint8** ppInput, long nrl, long nrh, long ncl, long nch, p_struct_elem_dim s, uint8 **ppOutput)
{
	const long order = 3;
	long row = nrl, col = ncl, x, y, r = 0;
	uint8 **temp = ui8matrix(nrl + s->nrl, nrh + s->nrh, ncl + s->ncl, nch + s->nch);
	uint8 *row0, *row1, *row2, *row3, y0, y1, y2, y3, y4, y5, z0, z1, z2;
	uint8 x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14;
	uint8 *temp_row0, *temp_row1, *temp_row2, *temp_row3, *temp_row4;
	uint8 *out_row0, *out_row1, *out_row2, *out_row3;
	
	r = (nch + 1) % order;
	// Prologue	
	row0 = ppInput[row - 1];
	row1 = ppInput[row + 0];
	row2 = ppInput[row + 1];
	
	// temp_row1 = temp[row - 1];
	// temp_row2 = temp[row + 0];
	// temp_row3 = temp[row + 1];
	for (row = nrl; row < nrh + 1; row++) {
		row0 = ppInput[row - 1];
		row1 = ppInput[row + 0];
		row2 = ppInput[row + 1];
		temp_row1 = temp[row + 0];

		temp_row1[ncl - 1] = row0[ncl - 1] | row1[ncl - 1] | row2[ncl - 1];
		temp_row1[ncl + 0] = row0[ncl + 0] | row1[ncl + 0] | row2[ncl + 0];
	}

	row0 = ppInput[ncl - 1];
	row1 = ppInput[ncl + 0];
	for (row = nrl; row < nrh + 1; row ++){
		temp_row1 = temp[row + 0];
		
		row2 = ppInput[row + 1];

		out_row1 = ppOutput[row]; 
		
		y0 = temp_row1[ncl - 1];		// LOAD => A
		y1 = temp_row1[ncl - 0];		// LOAD => B
		
		for (col = ncl; col < nch + 1 - r; col += order) {
			x0 = row0[col + 1]; x1 = row1[col + 1]; x2 = row2[col + 1];
			x3 = row0[col + 2]; x4 = row1[col + 2]; x5 = row2[col + 2];
			x6 = row0[col + 3]; x7 = row1[col + 3]; x8 = row2[col + 3];
			
			z0 = (x0 | x1 | x2);
			z1 = (x3 | x4 | x5);
			z2 = (x6 | x7 | x8);
			
			out_row1 [col + 0] = y0	| y1 | z0; // A | B | C
			out_row1 [col + 1] = y1	| z0 | z1; // B | C | D
			out_row1 [col + 2] = z0 | z1 | z2; // C | D | E

			// temp_row1[col + 1] = (x0 | x1 | x2);		  // STORE => C 
			temp_row1[col + 2] = z1;		  // STORE => D
			temp_row1[col + 3] = z2;		  // STORE => E 
														  // NEXT  :  F 
														  // NEXT  :  G 
														  // NEXT  :  H
			y0 = z1; 						  // STORE => D
			y1 = z2; 						  // STORE => E 
											  // NEXT  :  F 
											  // NEXT  :  G 
											  // NEXT  :  H
		}
		row0 = row1;
		row1 = row2;
	}
	switch (r) {
		case 2: 
			row0 = ppInput[nrl - 1];
			row1 = ppInput[nrl + 0];
			for (row = nrl; row < nrh + 1; row++) {
				row2 = ppInput[row + 1];
				temp_row1 = temp[row + 0];
				out_row1 = ppOutput[row]; 


				x0 = row0[nch - 2]; x1 = row1[nch - 2]; x2 = row2[nch - 2];
				x3 = row0[nch - 1]; x4 = row1[nch - 1]; x5 = row2[nch - 1];
				x6 = row0[nch + 0]; x7 = row1[nch + 0]; x8 = row2[nch + 0];
				x9 = row0[nch + 1]; x10 = row1[nch + 1]; x11 = row2[nch + 1];
			
				out_row1 [nch - 1] = x0 | x1 | x2 | x3 | x4 | x5 | x6 | x7  | x8; // A | B | C
				out_row1 [nch + 0] = x3 | x4 | x5 | x6 | x7 | x8 | x9 | x10 | x11; // B | C | D
				row0 = row1;
				row1 = row2;
			}
			break;
		case 1: 
			row0 = ppInput[nrl - 1];
			row1 = ppInput[nrl + 0];
			for (row = nrl; row < nrh + 1; row++) {
				row2 = ppInput[row + 1];
				// temp_row1 = temp[row + 0];
				x3 = row0[nch - 1]; x4 = row1[nch - 1]; x5 = row2[nch - 1];
				x6 = row0[nch + 0]; x7 = row1[nch + 0]; x8 = row2[nch + 0];
				x9 = row0[nch + 1]; x10 = row1[nch + 1]; x11 = row2[nch + 1];
				ppOutput[row][nch] =  x3 | x4 | x5 | x6 | x7 | x8 | x9 | x10 | x11; // A | B | C
				// out_row1 [nch + 0] = x1 | x2 | row0[nch + 1] | row1[nch + 0] | row2[nch + 0]; // B | C | D
				row0 = row1;
				row1 = row2;
			}
			break;
		default:
			break;
	}
	// display_ui8matrix(ppOutput, nrl - 1, nrh + 1, ncl - 1, nch + 1, "%u", "TEst");
	free_ui8matrix(temp, nrl + s->nrl, nrh + s->nrh, ncl + s->ncl, nch + s->nch);
	
}


void ui8matrix_dilation_pipeline_LU3x3_O3xO1_RR (uint8** ppInput, long nrl, long nrh, long ncl, long nch, p_struct_elem_dim s, uint8 **ppOutput)
{
	const long order = 3;
	long row = nrl, col = ncl, r = 0;
	uint8 **temp = ui8matrix(nrl + s->nrl, nrh + s->nrh, ncl + s->ncl, nch + s->nch);
	uint8 *row0, *row1, *row2, *row3, x0, x1, x2, x3;
	uint8 y0, y1, y2, y3, y4, y5, y6, y7, y8, y9, y10, y11, y12, y13, y14;
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
	y0 = row0[ncl-1]; y1 = row0[ncl+0];
	y3 = row1[ncl-1]; y4 = row1[ncl+0];
	for (col = ncl; col < nch + 1; col++) {
		y2 = row0[col + 1];
		y5 = row1[col + 1];
		x1 = y0|y1|y2;//scalar_or3(row1, col); // OR3  => C
		x2 = y3|y4|y5;//scalar_or3(row2, col); // OR3  => D
		temp_row1[col] = y0|y1|y2;//scalar_or3(row0, col);
		temp_row2[col] = y3|y4|y5;//scalar_or3(row1, col);
		y0 = y1; y1 = y2;
		y3 = y4; y4 = y5;
			
	}

	// row0 = ppInput[row + 0];   //  0 : b0 b1 b2 ...
	temp_row0 = temp[row - 1]; // -1 : A    <= (a0 a1 a2 ...) Prologued
	temp_row1 = temp[row - 0]; //  0 : B    <= (b0 b1 b2 ...) Prologued
	for (row = row + 1; row < nrh + 1 - r; row += order){
		row1 = ppInput[row + 0];   //  1 : c0 c1 c2 ...
		row2 = ppInput[row + 1];   //  2 : d0 d1 d2 ...
		row3 = ppInput[row + 2];   //  3 : e0 e1 e2 ...

		temp_row2 = temp[row + 0]; //  1 : NULL <= (c0 c1 c2 ...)
		temp_row3 = temp[row + 1]; //  2 : NULL <= (d0 d1 d2 ...)
		temp_row4 = temp[row + 2]; //  3 : NULL <= (e0 e1 e2 ...)
		out_row0 = ppOutput[row - 1]; 
		out_row1 = ppOutput[row + 0];
		out_row2 = ppOutput[row + 1];

		y0 = row1[ncl-1]; y1 = row1[ncl+0];
		y3 = row2[ncl-1]; y4 = row2[ncl+0];
		y6 = row3[ncl-1]; y7 = row3[ncl+0];
		for (col = ncl; col < nch + 1; col++) {
			//  temp_row0[ col]			   LOAD => A										
			x0 = temp_row1[col];		// LOAD => B

			
			// y = row0[col - 1];  y  = row0[col + 0];  y2  = row0[col + 1];
			y2 = row1[col + 1];
			y5 = row2[col + 1];
			y8 = row3[col + 1];

			x1 = y0|y1|y2;//scalar_or3(row1, col); // OR3  => C
			x2 = y3|y4|y5;//scalar_or3(row2, col); // OR3  => D
			x3 = y6|y7|y8;//scalar_or3(row3, col); // OR3  => E
			out_row0[col] = temp_row0[ col] | x0 | x1; // A | B | C
			out_row1[col] = 			 x0 | x1 | x2; // B | C | D
			out_row2[col] = 			 x1 | x2 | x3; // C | D | E
			temp_row2[col] = x1;		// STORE => C 
			temp_row3[col] = x2;		// STORE => D
			temp_row4[col] = x3;		// STORE => E 
										// NEXT  :  F 
										// NEXT  :  G 
										// NEXT  :  H
			y0 = y1; y1 = y2;
			y3 = y4; y4 = y5;
			y6 = y7; y7 = y8;			
		}
		temp_row0 = temp_row3;
		temp_row1 = temp_row4;
	}
	
	
	switch (r) {
		case 2: 
			row1 = ppInput[nrh - 1];   // -1 : (f0 f1 f2 ...)
			row2 = ppInput[nrh + 0];   //  0 : (g0 g1 g2 ...)
			row3 = ppInput[nrh + 1];   //  1 : (h0 h1 h2 ...)
			temp_row0 = temp[nrh - 2]; //  E    <= (e0 e1 e2 ...)
			out_row0 = ppOutput[nrh - 1];
			out_row1 = ppOutput[nrh + 0];
			y0 = row1[ncl - 1]; y1 = row1[ncl + 0];
			y3 = row2[ncl - 1]; y4 = row2[ncl + 0];
			y6 = row3[ncl - 1]; y7 = row3[ncl + 0];
			for (col = ncl; col < nch + 1; col++) {
				//  temp_row0[ col]		   		 LOAD => E	
				y2 = row1[col + 1];
				y5 = row2[col + 1];
				y8 = row3[col + 1];

				x0 = y0|y1|y2;//scalar_or3(row1, col); // OR3  => C
				x1 = y3|y4|y5;//scalar_or3(row2, col); // OR3  => D

				out_row0[col] = temp_row0[ col] | x0 | x1; // A | B | C
				out_row1[col] = 			 x0 | x1 | y6|y7|y8; // B | C | D	
				y0 = y1; y1 = y2;
				y3 = y4; y4 = y5;
				y6 = y7; y7 = y8;
			}
			break;
		case 1: 
			row2 = ppInput[nrh + 0];   //  0 : (f0 f1 f2 ...)
			row3 = ppInput[nrh + 1];   //  1 : (g0 g1 g2 ...)
			temp_row0 = temp[nrh - 1]; //  E    <= (e0 e1 e2 ...)
			out_row1 = ppOutput[nrh + 0];
			y3 = row2[ncl - 1]; y4 = row2[ncl + 0];
			y6 = row3[ncl - 1]; y7 = row3[ncl + 0];
			for (col = ncl; col < nch + 1; col++) {
				y5 = row2[col + 1];
				y8 = row3[col + 1];
				out_row1[col] = temp_row0[ col] | y3 | y4 | y5 | y6 | y7 | y8; // E | F | G	
				y3 = y4; y4 = y5;
				y6 = y7; y7 = y8;
			}

			
			break;
		default:
			break;
	}
	free_ui8matrix(temp, nrl + s->nrl, nrh + s->nrh, ncl + s->ncl, nch + s->nch);
}


// DEPRECATED VERSIONS OF ui8matrix_dilation_LU3x3_O1xO3_RR
void ui8matrix_dilation_LU3x3_O1xO3_RR_ver2 (uint8** ppInput, long nrl, long nrh, long ncl, long nch, p_struct_elem_dim s, uint8 **ppOutput)
{
	const long order = 3;
	long row = nrl, col = ncl, x, y, r;
	uint8 y0, y1, y2, y3, y4, y5;
	uint8 x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16, x17, x18, x19;
	uint8 *row0, *row1, *row2, *row3, *row4;

	r = (nch + 1) % order;
	row0 = ppInput[row - 1];
	row1 = ppInput[row + 0];
	for (row = nrl; row < nrh + 1; row++) {
		row2 = ppInput[row + 1];

		x0  = row0[ncl - 1]; x1  = row0[ncl + 0]; 
		x5  = row1[ncl - 1]; x6  = row1[ncl + 0]; 
		x10 = row2[ncl - 1]; x11 = row2[ncl + 0]; 
		for (col = ncl; col < nch + 1 - r; col += order) {
			x2  = row0[col + 1]; x3  = row0[col + 2]; x4  = row0[col + 3];
			x7  = row1[col + 1]; x8  = row1[col + 2]; x9  = row1[col + 3];
			x12 = row2[col + 1]; x13 = row2[col + 2]; x14 = row2[col + 3];
			// y0 = x1  | x2; 
			// y1 = x2  | x3;
			// y2 = x6  | x7;
			// y3 = x7  | x8;
			// y4 = x11 | x12; 
			// y5 = x12 | x13;
			// ppOutput[row][col + 0] = (x0 | x1 | x2 ) | (x5 | x6 | x7) | (x10 | x11 | x12);
			// ppOutput[row][col + 1] = (x1 | x2 | x3 ) | (x6 | x7 | x8) | (x11 | x12 | x13);
			// ppOutput[row][col + 2] = (x2 | x3 | x4 ) | (x7 | x8 | x9) | (x12 | x13 | x14);
			
			ppOutput[row][col + 0] = x0 | x1 | x2 |           x5 | x6 | x7 |           x10 | x11 | x12; 			//8 calc -> 3 calc
			ppOutput[row][col + 1] =      x1 | x2 | x3 |           x6 | x7 | x8 |            x11 | x12 | x13; 		//8 calc -> 3 calc
			ppOutput[row][col + 2] =           x2 | x3 | x4 |           x7 | x8 | x9 |             x12 | x13 | x14; //8 calc -> 3 calc
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
			x0  = row0[col - 1]; x1  = row0[col +  0]; x2  = row0[col + 1]; x3  = row0[col + 2];
			x5  = row1[col - 1]; x6  = row1[col +  0]; x7  = row1[col + 1]; x8  = row1[col + 2];
			x10 = row2[col - 1]; x11 = row2[col +  0]; x12 = row2[col + 1]; x13 = row2[col + 2];
			ppOutput[row][col + 0] = (x0 | x1 | x2 ) | (x5 | x6 | x7) | (x10 | x11 | x12);
			ppOutput[row][col + 1] = (x1 | x2 | x3 ) | (x6 | x7 | x8) | (x11 | x12 | x13);
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


void ui8matrix_dilation_LU3x3_O1xO3_RR_ver3 (uint8** ppInput, long nrl, long nrh, long ncl, long nch, p_struct_elem_dim s, uint8 **ppOutput)
{
	const long order = 3;
	long row = nrl, col = ncl, x, y, r;
	uint8 y0, y1, y2, y3, y4, y5;
	uint8 x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16, x17, x18, x19;
	uint8 *row0, *row1, *row2, *row3, *row4;

	r = (nch + 1) % order;
	row0 = ppInput[row - 1];
	row1 = ppInput[row + 0];
	for (row = nrl; row < nrh + 1; row++) {
		row2 = ppInput[row + 1];

		x0  = row0[ncl - 1]; x1  = row0[ncl + 0]; 
		x5  = row1[ncl - 1]; x6  = row1[ncl + 0]; 
		x10 = row2[ncl - 1]; x11 = row2[ncl + 0]; 
		for (col = ncl; col < nch + 1 - r; col += order) {
			x2  = row0[col + 1]; x3  = row0[col + 2]; x4  = row0[col + 3];
			x7  = row1[col + 1]; x8  = row1[col + 2]; x9  = row1[col + 3];
			x12 = row2[col + 1]; x13 = row2[col + 2]; x14 = row2[col + 3];
			
			y1 = x1 | x2 | x6 | x7 | x11 | x12; // 5 calcs
			y2 = x3 | x8 | x13;					// 2 calcs
			ppOutput[row][col + 0] = x0 | x5 | y1 | x10; //3 calcs
			ppOutput[row][col + 1] =      y2 | y1; 		 //1 calc
			ppOutput[row][col + 2] =      y2 | x2 | x4          | x7 | x9             | x12 | x14; //6 calcs
			// 5 + 2 + 3 + 1 + 6 = 17						

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
			x0  = row0[col - 1]; x1  = row0[col +  0]; x2  = row0[col + 1]; x3  = row0[col + 2];
			x5  = row1[col - 1]; x6  = row1[col +  0]; x7  = row1[col + 1]; x8  = row1[col + 2];
			x10 = row2[col - 1]; x11 = row2[col +  0]; x12 = row2[col + 1]; x13 = row2[col + 2];
			y1 = x1 | x2 | x6 | x7 | x11 | x12; // 5 calcs
			y2 = x3 | x8 | x13;					// 2 calcs
			ppOutput[row][col + 0] = x0 | x5 | y1 | x10; //3 calcs
			ppOutput[row][col + 1] =      y2 | y1; 		 //1 calc
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
/*
void ui8matrix_dilation_LU3x3_O3xO3_RR(uint8** ppInput, long nrl, long nrh, long ncl, long nch, p_struct_elem_dim s, uint8 **ppOutput)
{
		const long order = 3;
	long row = nrl, col = ncl, x, y, rr, cr;
	uint8 y00, y10, y20, y30, y40, y01, y11, y21, y31, y41, y02, y12, y22, y32, y42;
	uint8 x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16, x17, x18, x19, x20, x21, x22, x23, x24, x25;
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
		out_row0 = ppOutput[row + 0];
		out_row1 = ppOutput[row + 1];
		out_row2 = ppOutput[row + 2];
		x0  = row0[ncl-1]; x1  = row0[ncl+ 0]; 
		x5  = row1[ncl-1]; x6  = row1[ncl+ 0]; 
		x10 = row2[ncl-1]; x11 = row2[ncl+ 0]; 
		x15 = row3[ncl-1]; x16 = row3[ncl+ 0]; 
		x20 = row4[ncl-1]; x21 = row4[ncl+ 0]; 
		for (col = ncl; col < nch + 1 - cr; col += order) {

			x2  = row0[col + 1]; x3  = row0[col + 2]; x4  = row0[col + 3];
			x7  = row1[col + 1]; x8  = row1[col + 2]; x9  = row1[col + 3];
			x12 = row2[col + 1]; x13 = row2[col + 2]; x14 = row2[col + 3];
			x17 = row3[col + 1]; x18 = row3[col + 2]; x19 = row3[col + 3];
			x22 = row4[col + 1]; x23 = row4[col + 2]; x24 = row4[col + 3];

			y00 = x0  | x1  | x2 ; y01 = x1  | x2  | x3 ; y02 = x2  | x3  | x4 ;
			y10 = x5  | x6  | x7 ; y11 = x6  | x7  | x8 ; y12 = x7  | x8  | x9 ;
			y20 = x10 | x11 | x12; y21 = x11 | x12 | x13; y22 = x12 | x13 | x14;
			y30 = x15 | x16 | x17; y31 = x16 | x17 | x18; y32 = x17 | x18 | x19 ;
			y40 = x20 | x21 | x22; y41 = x21 | x22 | x23; y42 = x22 | x23 | x24;

			out_row0[col + 0] = y00 | y10 | y20;
			out_row1[col + 0] = y10 | y20 | y30;
			out_row2[col + 0] = y20 | y30 | y40;
			// out_row2[col + 0] = y30 | y40 | y50;
			// y00 = y30;
			// y10 = y40;
			// y3 = 

			out_row0[col + 1] = y01 | y11 | y21;
			out_row1[col + 1] = y11 | y21 | y31;
			out_row2[col + 1] = y21 | y31 | y41;
			// y1 = y13;


			out_row0[col + 2] = y02 | y12 | y22;
			out_row1[col + 2] = y12 | y22 | y32;
			out_row2[col + 2] = y22 | y32 | y42;
			// y2 = y14;
			
			x0  = x3 ; x1  = x4 ; 
			x5  = x8 ; x6  = x9 ; 
			x10 = x13; x11 = x14; 
			x15 = x18; x16 = x19; 
			x20 = x23; x21 = x24; 

		}
	}
	switch (rr) {
		case 2 :
			switch(cr){
				case 2 :
					ui8matrix_dilation_LU3x3_O3xO1_RR(ppInput, nrl, nrh, nch - 1, nch, s, ppOutput);
					ui8matrix_dilation_LU3x3_O1xO3_RR(ppInput, nrh - 1, nrh, ncl, nch, s, ppOutput);
				break;
				case 1: 
					ui8matrix_dilation_LU3x3_O3xO1_RR(ppInput, nrl, nrh, nch, nch, s, ppOutput);		
					ui8matrix_dilation_LU3x3_O1xO3_RR(ppInput, nrh - 1, nrh, ncl, nch, s, ppOutput);	
				break;
				default :
					ui8matrix_dilation_LU3x3_O1xO3_RR(ppInput, nrh - 1, nrh, ncl, nch, s, ppOutput);
				break;
			}
			break;
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
*/

/*
const long order = 3;
	long row = nrl, col = ncl, x, y, r;
	uint8 y0, y1, y2, x0, x1, x2, x3, x4, x5, x6, x7, x8;
	uint8 *row0, *row1, *row2, *row3, *row4;
	uint8 *out_row0;
	r = (nch + 1)  % order;
	
	for (row = nrl; row < nrh + 1; row ++) {
		row0 = ppInput[row - 1];
		row1 = ppInput[row + 0];
		row2 = ppInput[row + 1];
		out_row0 = ppOutput[row + 0];
		for (col = ncl; col < nch + 1 - r; col += order) {
			// x0 = scalar_or3(row0, col); x3 = scalar_or3(row0, col + 1); x6 = scalar_or3(row0, col + 2);
			// x1 = scalar_or3(row1, col); x4 = scalar_or3(row1, col + 1); x7 = scalar_or3(row1, col + 2);
			// x2 = scalar_or3(row2, col); x5 = scalar_or3(row2, col + 1); x8 = scalar_or3(row2, col + 2);

			out_row0[col + 0] = scalar_or3(row0, col + 0)|
								scalar_or3(row1, col + 0)|
								scalar_or3(row2, col + 0);
			out_row0[col + 1] = scalar_or3(row0, col + 1)|
 								scalar_or3(row1, col + 1)|
 								scalar_or3(row2, col + 1);
			out_row0[col + 2] = scalar_or3(row0, col + 2)|
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
				out_row0 = ppOutput[row + 0];

				out_row0[nch - 1] = scalar_or3(row0, nch - 1)|
 									scalar_or3(row1, nch - 1)|
 									scalar_or3(row2, nch - 1);
				out_row0[nch + 0] = scalar_or3(row0, nch + 0)|
									scalar_or3(row1, nch + 0)|
									scalar_or3(row2, nch + 0);
			}
			break;
		case 1:
			for (row = nrl; row < nrh + 1; row ++) {
				ppOutput[row + 0][nch + 0] = scalar_or3(ppInput[row - 1], nch) | 
										  	 scalar_or3(ppInput[row + 0], nch) | 
										  	 scalar_or3(ppInput[row + 1], nch);
			}
			// ui8matrix_dilation_LU3x3_O1xO1(ppInput, nrl, nrh, nch, nch, s, ppOutput);
			break;
		default:
			break;
	}*/
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

