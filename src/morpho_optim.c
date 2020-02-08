#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include "nrdef.h"
#include "nrutil.h"
#include "vnrdef.h"
#include "vnrutil.h"
#include <img_SIMD.h>
#include "util.h"
#include <img.h>
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
#include "mynrutil.h"
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

#define HANDLE_EDGE_OF_ROW(x) 	x[ncl - 1] = x[ncl];\
							  	x[nch + 1] = x[nch];
#define HANDLE_EDGE(x, row)		x[row][ncl - 1] = x[row][ncl];\
							  	x[row][nch + 1] = x[row][nch];
#define ZERO_EDGE_OF_ROW(x) 	x[ncl - 1] = 0;\
							  	x[nch + 1] = 0;
#define ZERO_EDGE(x, row)		x[row][ncl - 1] = 0;\
							  	x[row][nch + 1] = 0;

void ui8matrix_sequence_drnc_fo(uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **Y, uint8 **Z)
{
	uint8  *in_row, *out_row;
	uint8 **in, **mid, **out;
	long row, col, nrow, r;
	
	in = X; mid = Z; out = Y;
	for (long row = nrl; row < nrh + 1; row++){HANDLE_EDGE_OF_ROW(in[row]);	ZERO_EDGE_OF_ROW(out[row]);}
	ui8matrix_erosion_divide_row_and_conquer(in, nrl, nrh, ncl, nch, mid, out);

	in = Y; mid = Z; out = X;
	ui8matrix_dilation5_divide_row_and_conquer(in, nrl, nrh, ncl, nch, mid, out);

	in = X; mid = Y; out = Z;
	for (long row = nrl; row < nrh + 1; row++) {HANDLE_EDGE_OF_ROW(in[row]);}
	ui8matrix_erosion_divide_row_and_conquer(X, nrl, nrh, ncl, nch, Y, Z);
}

void ui8matrix_sequence_drnc_fo_pipeline2(uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **Y, uint8 **Z)
{
	uint8 *temp_row, *temp_row0, *temp_row1, *temp_row2, *temp_row3, *temp_row4, *out_row, *in_row;
	uint8 **in, **mid, **out;
	long row, col, nrow, r;
	long d5_row, e3_row, nrl_prime, nrh_prime;
	const long PRE_D5_NROW = 5;
	const long PRE_E3_NROW = 3;
	
	in = X; mid = Z; out = Y;
	for (long row = nrl; row < nrh + 1; row++){
		in[row][ncl - 1] = in[row][ncl]; out[row][ncl - 1] = 0;
		in[row][nch + 1] = in[row][nch]; out[row][nch + 1] = 0;
	}



	nrl_prime = (nrl - 2) + PRE_D5_NROW ;
	nrh_prime = (nrl - 2) + PRE_D5_NROW + PRE_E3_NROW ;
	ui8matrix_erosion_divide_row_and_conquer(in, nrl, nrl_prime, ncl, nch, mid, out);
	
	d5_row = nrl;
	e3_row = nrl;
	for (row = nrl_prime + 1; row < nrh_prime + 1; row ++) {
		in = X; mid = Z; out = Y;
		ui8matrix_erosion_divide_row_and_conquer(in, row, row, ncl, nch, mid, out);
		in = Y; mid = Z; out = X;
		ui8matrix_dilation5_divide_row_and_conquer(in, d5_row, d5_row, ncl, nch, mid, out);
		out_row = out[d5_row];
		out_row[ncl - 1] = out_row[ncl];
		out_row[nch + 1] = out_row[nch];
		d5_row++;
	}
	
	for (row = nrh_prime + 1; row < nrh + 1; row ++) {
		in = X; mid = Z; out = Y;
		ui8matrix_erosion_divide_row_and_conquer(in, row, row, ncl, nch, mid, out);
		in = Y; mid = Z; out = X;
		ui8matrix_dilation5_divide_row_and_conquer(in, d5_row, d5_row, ncl, nch, mid, out);
		out_row = out[d5_row];
		out_row[ncl - 1] = out_row[ncl];
		out_row[nch + 1] = out_row[nch];

		in = X; mid = Y; out = Z;
		memset_ui8matrix(&Z[e3_row], 0, 0, 0, ncl - 2, nch + 2);
		ui8matrix_erosion_divide_row_and_conquer(in, e3_row, e3_row, ncl, nch, mid, out);
		d5_row++;
		e3_row++;
	}
	
	// d5_row--;
	e3_row--;
	for (row = d5_row; row < nrh + 1; row ++) {
		in = Y; mid = Z; out = X;
		ui8matrix_dilation5_divide_row_and_conquer(in, row, nrh, ncl, nch, mid, out);
		// display_ui8matrix(out, nrl-2, nrh+2, ncl-2, nch+2, "%4u", "E3-D5-E3(Pipeline)");
		out[row][ncl - 1] = out[row][ncl];
		out[row][nch + 1] = out[row][nch];
		in = X; mid = Y; out = Z;
		// memset_ui8matrix(&Z[e3_row], 0, 0, 0, ncl - 2, nch + 2);
		ui8matrix_erosion_divide_row_and_conquer(in, e3_row, e3_row, ncl, nch, mid, out);
		e3_row++;
	}
	in = X; mid = Y; out = Z;
	for (row = e3_row; row < nrh + 1; row ++) {
		// memset_ui8matrix(&Z[row], 0, 0, 0, ncl - 2, nch + 2);
		ui8matrix_erosion_divide_row_and_conquer(in, row, row, ncl, nch, mid, out);
	}
	// display_ui8matrix(out, nrl-2, nrh+2, ncl-2, nch+2, "%4u", "E3-D5-E3(Pipeline)");

	
}

void ui8matrix_sequence_drnc_fo_pipeline(uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **Y, uint8 **Z)
{
	
	long row = nrl, col = ncl, x, y;
	uint8 **in, **mid_e3, **mid_d5, **out;
	uint8 *mid_e3_row, *mid_e3_row0, *mid_e3_row1, *mid_e3_row2;
	uint8 *mid_d5_row, *mid_d5_row0, *mid_d5_row1, *mid_d5_row2, *mid_d5_row3, *mid_d5_row4;
	uint8 *in_row, *mid_row, *out_row;
	const long PRE_D5_NROW = 5;
	const long PRE_E3_NROW = 3;
	
	
	// memset_ui8matrix(Y, 0, nrl - 2, nrh + 2, ncl - 2, nch + 2);
	in     = X; 
	mid_e3 = Y;
	out    = Z;

	for (row = nrl; row < (nrh + 1); row++) {	
		mid_e3_row 		= mid_e3[row]; 
		mid_e3_row[ncl - 1] = 0; // Empty buffer, just in case
		mid_e3_row[nch + 1] = 0; // Empty buffer, just in case
		in_row     		=      X[row];
		in_row[ncl - 1] = in_row[ncl];
		in_row[nch + 1] = in_row[nch];
		for (col = ncl; col < nch + 1; col++) {
			mid_e3_row[col] = scalar_and3(in_row, col);
		}
	}

	// Pre E3
	mid_e3_row0 = mid_e3[nrl - 1];
	mid_e3_row1 = mid_e3[nrl + 0];
	mid_e3_row2 = mid_e3[nrl + 1];
	out_row		= 	 out[nrl + 0];
	for (col = ncl; col < nch + 1; col++) {
		out_row[col] = mid_e3_row0[col] &
					   mid_e3_row1[col] &
					   mid_e3_row2[col];
	}
	mid_d5 = X;
	// E3 & Pre D5 
	for (row = nrl + 1; row < nrh + 1; row++) {
		mid_e3_row0 = mid_e3[row - 1];
		mid_e3_row1 = mid_e3[row + 0];
		mid_e3_row2 = mid_e3[row + 1];
		out_row		= 	 out[row + 0];

		mid_d5_row = mid_d5[row - 1]; 
		mid_row    =    out[row - 1]; 
		for (col = ncl; col < nch + 1; col++) {
			out_row[col] = mid_e3_row0[col] &
					 	   mid_e3_row1[col] &
					 	   mid_e3_row2[col] ;
			mid_d5_row[col] = scalar_or5(mid_row, col);
		}
	}
	// E3 Done 
	// display_ui8matrix(out, nrl - 2, nrh + 2, ncl - 2, nch + 2, "%4u", "E3");
	// Pre D5 Epilogue & D5 Prologue
	mid_d5_row0 = mid_d5[nrl - 2];
	mid_d5_row1 = mid_d5[nrl - 1];
	mid_d5_row2 = mid_d5[nrl + 0];
	mid_d5_row3 = mid_d5[nrl + 1];
	mid_d5_row4 = mid_d5[nrl + 2];
	out_row		= 	 out[nrl + 0];
	mid_d5_row  = mid_d5[nrh    ]; 
	mid_row     =    out[nrh    ]; 
	for (col = ncl; col < nch + 1; col++) {
		out_row[col] = mid_d5_row0[col] |
					   mid_d5_row1[col] |			
					   mid_d5_row2[col] |
					   mid_d5_row3[col] |
					   mid_d5_row4[col];
		mid_d5_row[col] = scalar_or5(mid_row, col);
	}
	// D5 & Pre E3
	for (row = nrl + 1; row < nrh + 1; row++) {
		mid_d5_row0 = mid_d5[row - 2];
		mid_d5_row1 = mid_d5[row - 1];
		mid_d5_row2 = mid_d5[row + 0];
		mid_d5_row3 = mid_d5[row + 1];
		mid_d5_row4 = mid_d5[row + 2];
		out_row		= 	 out[row + 0];

		mid_e3_row  = mid_e3[row - 1]; 
		mid_row     =    out[row - 1]; 
		mid_row[ncl - 1] = mid_row[ncl];
		mid_row[nch + 1] = mid_row[nch];
		for (col = ncl; col < nch + 1; col++) {
			out_row[col] = mid_d5_row0[col] |
						   mid_d5_row1[col] |			
						   mid_d5_row2[col] |
						   mid_d5_row3[col] |
						   mid_d5_row4[col];
			mid_e3_row[col] = scalar_and3(mid_row, col);
		}
	}
	mid_e3_row  = mid_e3[nrh]; 
	mid_row     =    out[nrh]; 
	mid_row[ncl - 1] = mid_row[ncl];
	mid_row[nch + 1] = mid_row[nch];
	// display_ui8matrix(out, nrl - 2, nrh + 2, ncl - 2, nch + 2, "%4u", "E3-D5");
	// Pre E3 Epilogue
	mid_e3_row  = mid_e3[nrh    ]; // Epilogue
	mid_row     =    out[nrh    ]; // Epilogue
	for (col = ncl; col < nch + 1; col++) {
		mid_e3_row[col] = scalar_and3(mid_row, col);
	}
	
	// memset_ui8matrix(X, 0, nrl - 2, nrh + 2, ncl - 2, nch + 2);
	for (row = nrl; row < nrh + 1; row++) {
		mid_e3_row0 = mid_e3[row - 1];
		mid_e3_row1 = mid_e3[row + 0];
		mid_e3_row2 = mid_e3[row + 1];
		out_row		= 	 out[row + 0];
		for (col = ncl; col < nch + 1; col++) {
			out_row[col] = mid_e3_row0[col] &
					 	   mid_e3_row1[col] &
					 	   mid_e3_row2[col] ;
		}
	}
	
	// display_ui8matrix(Z, nrl - 2, nrh + 2, ncl - 2, nch + 2, "%4u", "E3-D5-E3");
	// getchar();
}
/*******************************************/
/******* Optimisation : Loop Unroll ********/
/*******************************************/

void ui8matrix_dilation_LU3x3_O1xO1(uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **temp_buffer, uint8 **Y)
{
	long row = nrl, col = ncl, x, y;
	// dilate
	for (row = nrl; row < nrh + 1; row++)
		for (col = ncl; col < nch + 1; col++)
            Y[row][col] = scalar_or3x3(&X[row], col);
}

void ui8matrix_dilation_LU3x3_InLU_O3 (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **temp_buffer, uint8 **Y)
{
	const long order = 3;
	long row = nrl, col = ncl, x, y, r;
	uint8 y0, y1, y2, x0, x1, x2, x3, x4, x5, x6, x7, x8;
	uint8 *row0, *row1, *row2, *row3, *row4;
	uint8 *out_row0;
	r = (nch + 1)  % order;
	
	for (row = nrl; row < nrh + 1; row ++) {
		row0 = X[row - 1];
		row1 = X[row + 0];
		row2 = X[row + 1];
		out_row0 = Y[row + 0];
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
				row0 = X[row - 1];
				row1 = X[row + 0];
				row2 = X[row + 1];
				out_row0 = Y[row + 0];

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
				Y[row + 0][nch + 0] = scalar_or3(X[row - 1], nch) | 
										  	 scalar_or3(X[row + 0], nch) | 
										  	 scalar_or3(X[row + 1], nch);
			}
			// ui8matrix_dilation_LU3x3_O1xO1(X, nrl, nrh, nch, nch, temp_buffer, Y);
			break;
		default:
			break;
	}
}

void ui8matrix_dilation_LU3x3_ExLU_O3 (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **temp_buffer, uint8 **Y)
{
	const long order = 3;
	long row = nrl, col = ncl, x, y;
	uint8 y0, y1, y2, x0, x1, x2, x3, r;
	uint8 *row0, *row1, *row2, *row3, *row4;
	uint8 *out_row0, *out_row1, *out_row2;

	r = (nrh + 1)  % order;
	
	for (row = nrl; row < nrh + 1 - r; row += order) {
		row0 = X[row - 1];
		row1 = X[row + 0];
		row2 = X[row + 1];
		row3 = X[row + 2];
		row4 = X[row + 3];

		out_row0 = Y[row + 0];
		out_row1 = Y[row + 1];
		out_row2 = Y[row + 2];
		
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
			row0 = X[row - 1];
			row1 = X[row + 0];
			row2 = X[row + 1];
			row3 = X[row + 2];
			out_row0 = Y[row + 0];
			out_row1 = Y[row + 1];
			for (col = ncl; col < nch + 1; col++) {
				x0 = scalar_or3(row1, col);
				x1 = scalar_or3(row2, col);
				out_row0[col] = scalar_or3(row0, col) | x0 | x1;
				out_row1[col] = 				     x0 | x1 | scalar_or3(row3, col);
			}
			break;
		case 1:
			row0 = X[row - 1];
			row1 = X[row + 0];
			row2 = X[row + 1];
			out_row0 = Y[row + 0];
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

void ui8matrix_dilation_LU3x3_ComLU_O3 (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **temp_buffer, uint8 **Y)
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
		row0 = X[row - 1];
		row1 = X[row + 0];
		row2 = X[row + 1];
		row3 = X[row + 2];
		row4 = X[row + 3];
		out_row0 = Y[row + 0];
		out_row1 = Y[row + 1];
		out_row2 = Y[row + 2];
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
					ui8matrix_dilation_LU3x3_ExLU_O3(X, nrl, nrh, nch - 1, nch, temp_buffer, Y);
					ui8matrix_dilation_LU3x3_InLU_O3(X, nrh - 1, nrh, ncl, nch, temp_buffer, Y);
				break;
				case 1: 
					ui8matrix_dilation_LU3x3_ExLU_O3(X, nrl, nrh, nch, nch, temp_buffer, Y);		
					ui8matrix_dilation_LU3x3_InLU_O3(X, nrh - 1, nrh, ncl, nch, temp_buffer, Y);	
				break;
				default :
					ui8matrix_dilation_LU3x3_InLU_O3(X, nrh - 1, nrh, ncl, nch, temp_buffer, Y);
				break;
			}
			break;
		case 1:
			switch(cr){
				case 2 :
					ui8matrix_dilation_LU3x3_ExLU_O3(X, nrl, nrh, nch - 1, nch, temp_buffer, Y);
					ui8matrix_dilation_LU3x3_InLU_O3(X, nrh, nrh, ncl, nch, temp_buffer, Y);
				break;
				case 1: 
					ui8matrix_dilation_LU3x3_ExLU_O3(X, nrl, nrh, nch, nch, temp_buffer, Y);		
					ui8matrix_dilation_LU3x3_InLU_O3(X, nrh, nrh, ncl, nch, temp_buffer, Y);	
				break;
				default :
					ui8matrix_dilation_LU3x3_InLU_O3(X, nrh, nrh, ncl, nch, temp_buffer, Y);
				break;
			}
		break;
		default:
			switch(cr){
				case 2 :
					ui8matrix_dilation_LU3x3_ExLU_O3(X, nrl, nrh, nch - 1, nch, temp_buffer, Y);
				break;
				case 1: 
					ui8matrix_dilation_LU3x3_ExLU_O3(X, nrl, nrh, nch, nch, temp_buffer, Y);	
				break;
				default :
				break;
			}
		break;
	}
}
/**************************************************/
/******* End of optimisation : Loop Unroll ********/
/**************************************************/

/**************************************************************************************/
/******* Optimisation : Loop Unroll +  Register Rotation of Values / Addresses ********/
/**************************************************************************************/
void ui8matrix_dilation_LU3x3_InLU_O3_ValAddrRR (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **temp_buffer, uint8 **Y)
{
	const long order = 3;
	long row = nrl, col = ncl, r;
	uint8 y0, y1, y2;
	uint8 x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14;
	uint8 *row0, *row1, *row2, *row3, *row4;
	uint8 *out_row0;

	r = (nch + 1) % order;
	row0 = X[nrl - 1];
	row1 = X[nrl + 0];
	for (row = nrl; row < nrh + 1; row++) {
		row2 = X[row + 1];

		out_row0 = Y[row + 0];
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
			row0 = X[row - 1];
			row1 = X[row + 0];
			row2 = X[row + 1];

			out_row0 = Y[row + 0];
		
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
			row0 = X[row - 1];
			row1 = X[row + 0];
			row2 = X[row + 1];

			out_row0 = Y[row + 0];
			out_row0[nch + 0] = (row0[nch - 1] |  row0[nch + 0] | row0[nch + 1])| 
			  					(row1[nch - 1] |  row1[nch + 0] | row1[nch + 1])| 
								(row2[nch - 1] |  row2[nch + 0] | row2[nch + 1]);
			
		}
		break;
	default:
		break;
	}
}
void ui8matrix_dilation_LU3x3_ExLU_O3_ValAddrRR(uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **temp_buffer, uint8 **Y)
{
	
	const long order = 3;
	long row = nrl, col = ncl, x, y;
	uint8 y0, y1, y2, y3, y5, r;
	uint8 x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14;
	uint8 *row0, *row1, *row2, *row3, *row4;
	uint8 *out_row0, *out_row1, *out_row2;

	r = (nrh + 1)  % order;
	
	row0 = X[row - 1];
	row1 = X[row + 0];
	for (row = nrl; row < nrh + 1 - r; row += order) {
		row2 = X[row + 1]; 
		row3 = X[row + 2]; 
		row4 = X[row + 3];
		x0 = row0[ncl - 1]; x5 = row0[ncl + 0];
		x1 = row1[ncl - 1]; x6 = row1[ncl + 0];
		x2 = row2[ncl - 1]; x7 = row2[ncl + 0];
		x3 = row3[ncl - 1]; x8 = row3[ncl + 0];
		x4 = row4[ncl - 1]; x9 = row4[ncl + 0];
		out_row0 = Y[row + 0];
		out_row1 = Y[row + 1];
		out_row2 = Y[row + 2];
		for (col = ncl; col < nch + 1; col++) {
			x10 = row0[col + 1]; //scalar_or3(row1, col);
			x11 = row1[col + 1]; //scalar_or3(row2, col);
			x12 = row2[col + 1]; //scalar_or3(row3, col);
			x13 = row3[col + 1]; //scalar_or3(row2, col);
			x14 = row4[col + 1]; //scalar_or3(row3, col);
			
			y1 = (x1 | x6 | x11); // 2 calc
			y2 = (x2 | x7 | x12); // 2 calc
			y3 = (x3 | x8 | x13); // 2 calc
			out_row0[col] = (x0 | x5 | x10) | y1 | y2; 				
			out_row1[col] = 				  y1 | y2 | y3;				
			out_row2[col] = 		          y2 | y3 | (x4 | x9 | x14); 
			// 6 + 12 = 18 calc
			x0 = x5; x5 = x10;
			x1 = x6; x6 = x11;
			x2 = x7; x7 = x12;
			x3 = x8; x8 = x13;
			x4 = x9; x9 = x14;
		}
		row0 = row3;
		row1 = row4;
	}
	row2 = X[row + 1];
	row3 = X[row + 2];
	row4 = X[row + 3];
	out_row0 = Y[nrh - 1];
	out_row1 = Y[nrh + 0];
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
void ui8matrix_dilation_LU3x3_ComLU_O3_ValAddrRR(uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **temp_buffer, uint8 **Y)
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
	
	row0 = X[row - 1];
	row1 = X[row + 0];

	for (row = nrl; row < nrh + 1 - rr; row += order) {
		row2 = X[row + 1];
		row3 = X[row + 2];
		row4 = X[row + 3];
		out_row0 = Y[row + 0];
		out_row1 = Y[row + 1];
		out_row2 = Y[row + 2];


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
					ui8matrix_dilation_LU3x3_ExLU_O3_ValAddrRR(X, nrl, nrh, nch - 1, nch, temp_buffer, Y);
					ui8matrix_dilation_LU3x3_InLU_O3_ValAddrRR(X, nrh - 1, nrh, ncl, nch, temp_buffer, Y);
				break;
				case 1: 
					ui8matrix_dilation_LU3x3_ExLU_O3_ValAddrRR(X, nrl, nrh, nch, nch, temp_buffer, Y);		
					ui8matrix_dilation_LU3x3_InLU_O3_ValAddrRR(X, nrh - 1, nrh, ncl, nch, temp_buffer, Y);	
				break;
				default :
					ui8matrix_dilation_LU3x3_InLU_O3_ValAddrRR(X, nrh - 1, nrh, ncl, nch, temp_buffer, Y);
				break;
			}
			break;
		case 1:
			switch(cr){
				case 2 :
					ui8matrix_dilation_LU3x3_ExLU_O3_ValAddrRR(X, nrl, nrh, nch - 1, nch, temp_buffer, Y);
					ui8matrix_dilation_LU3x3_InLU_O3_ValAddrRR(X, nrh, nrh, ncl, nch, temp_buffer, Y);
				break;
				case 1: 
					ui8matrix_dilation_LU3x3_ExLU_O3_ValAddrRR(X, nrl, nrh, nch, nch, temp_buffer, Y);		
					ui8matrix_dilation_LU3x3_InLU_O3_ValAddrRR(X, nrh, nrh, ncl, nch, temp_buffer, Y);	
					// display_ui8matrix(Y, nrl, nrh, ncl, nch, "%u", "rr=1; cr=1");
				break;
				default :
					ui8matrix_dilation_LU3x3_InLU_O3_ValAddrRR(X, nrh, nrh, ncl, nch, temp_buffer, Y);
				break;
			}
		break;
		default:
			switch(cr){
				case 2 :
					ui8matrix_dilation_LU3x3_ExLU_O3_ValAddrRR(X, nrl, nrh, nch - 1, nch, temp_buffer, Y);
				break;
				case 1: 
					ui8matrix_dilation_LU3x3_ExLU_O3_ValAddrRR(X, nrl, nrh, nch, nch, temp_buffer, Y);	
				break;
				default :
				break;
			}
		break;
	}
		
}
/*********************************************************************************************/
/******* End of optimisation : Loop Unroll +  Register Rotation of Values / Addresses ********/
/*********************************************************************************************/

/*****************************************************************************/
/******** Optimisation : Loop Unroll + Register Rotation of Addresses ********/
/*****************************************************************************/

void ui8matrix_dilation_LU3x3_InLU_O3_AddrRR(uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **temp_buffer, uint8 **Y)
{
	const long order = 3;
	long row = nrl, col = ncl, r;
	uint8 y0, y1, y2;
	uint8 x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14;
	uint8 *row0, *row1, *row2, *row3, *row4;
	uint8 *out_row0;

	r = (nch + 1) % order;
	row0 = X[nrl - 1];
	row1 = X[nrl + 0];
	for (row = nrl; row < nrh + 1; row++) {
		row2 = X[row + 1];

		out_row0 = Y[row + 0];
		for (col = ncl; col < nch + 1 - r; col += order) {
			x0  = row0[col - 1]; x1  = row0[col + 0]; x2  = row0[col + 1]; x3  = row0[col + 2]; x4  = row0[col + 3];
			x5  = row1[col - 1]; x6  = row1[col + 0]; x7  = row1[col + 1]; x8  = row1[col + 2]; x9  = row1[col + 3];
			x10 = row2[col - 1]; x11 = row2[col + 0]; x12 = row2[col + 1]; x13 = row2[col + 2]; x14 = row2[col + 3];
			

			out_row0[col + 0] = (x0  | x1  | x2 ) | 
								(x5  | x6  | x7 ) | 
								(x10 | x11 | x12);

			out_row0[col + 1] = (x1  | x2  | x3) | 
								(x6  | x7  | x8) | 
								(x11 | x12 | x13);
								
			out_row0[col + 2] = (x2  | x3  | x4) | 
								(x7  | x8  | x9) | 
								(x12 | x13 | x14); 			
		}
		row0 = row1;
		row1 = row2;
	}

	switch(r) {
		case 2:

		
		for (row = nrl; row < nrh + 1; row++) {
			row0 = X[row - 1];
			row1 = X[row + 0];
			row2 = X[row + 1];

			out_row0 = Y[row + 0];
		
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
			row0 = X[row - 1];
			row1 = X[row + 0];
			row2 = X[row + 1];

			out_row0 = Y[row + 0];
		
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
void ui8matrix_dilation_LU3x3_ExLU_O3_AddrRR(uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **temp_buffer, uint8 **Y)
{
	
	const long order = 3;
	long row = nrl, col = ncl;
	uint8 y0, y1, y2, y3, y5, r;
	uint8 x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14;
	uint8 *row0, *row1, *row2, *row3, *row4;
	uint8 *out_row0, *out_row1, *out_row2;

	r = (nrh + 1)  % order;
	
	row0 = X[row - 1];
	row1 = X[row + 0];
	for (row = nrl; row < nrh + 1 - r; row += order) {
		row2 = X[row + 1]; 
		row3 = X[row + 2]; 
		row4 = X[row + 3];
		
		
		out_row0 = Y[row + 0];
		out_row1 = Y[row + 1];
		out_row2 = Y[row + 2];
		for (col = ncl; col < nch + 1; col++) {
			x0 = row0[col - 1]; x5 = row0[col]; x10 = row0[col + 1];
			x1 = row1[col - 1]; x6 = row1[col]; x11 = row1[col + 1];
			x2 = row2[col - 1]; x7 = row2[col]; x12 = row2[col + 1];
			x3 = row3[col - 1]; x8 = row3[col]; x13 = row3[col + 1];
			x4 = row4[col - 1]; x9 = row4[col]; x14 = row4[col + 1];
			

			out_row0[col] = (x0 | x5 | x10) | (x1 | x6 | x11) | (x2 | x7 | x12); 				
			out_row1[col] = (x1 | x6 | x11) | (x2 | x7 | x12) | (x3 | x8 | x13);				
			out_row2[col] = (x2 | x7 | x12) | (x3 | x8 | x13) | (x4 | x9 | x14); 

		}
		row0 = row3;
		row1 = row4;
	}
	row2 = X[row + 1];
	row3 = X[row + 2];
	row4 = X[row + 3];
	out_row0 = Y[nrh - 1];
	out_row1 = Y[nrh + 0];
	switch(r) {
		case 2:
			for (col = ncl; col < nch + 1; col++) {
				x0 = row1[col - 1]; x3 = row1[col + 0]; x6 = row1[col + 1]; //scalar_or3(row1, col);
				x1 = row2[col - 1]; x4 = row2[col + 0]; x7 = row2[col + 1]; //scalar_or3(row2, col);
				x2 = row3[col - 1]; x5 = row3[col + 0]; x8 = row3[col + 1]; //scalar_or3(row3, col);
				out_row0[col] = scalar_or3(row0, col) | (x0 | x3 | x6) | (x1 | x4 | x7);
				out_row1[col] =        (x0 | x3 | x6) | (x1 | x4 | x7) | (x2 | x5 | x8);
			}
			break;
		case 1:
			for (col = ncl; col < nch + 1; col++) {
				x0 = row1[col - 1]; x3 = row1[col + 0]; x6 = row1[col + 1]; //scalar_or3(row1, col);
				x1 = row2[col - 1]; x4 = row2[col + 0]; x7 = row2[col + 1]; //scalar_or3(row2, col);
				out_row1[col] = scalar_or3(row0, col) | (x0 | x3 | x6) | (x1 | x4 | x7);
			}				
			break;
		default:
			break;
	}
}
void ui8matrix_dilation_LU3x3_ComLU_O3_AddrRR(uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **temp_buffer, uint8 **Y)
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
	
	row0 = X[row - 1];
	row1 = X[row + 0];

	for (row = nrl; row < nrh + 1 - rr; row += order) {
		row2 = X[row + 1];
		row3 = X[row + 2];
		row4 = X[row + 3];
		out_row0 = Y[row + 0];
		out_row1 = Y[row + 1];
		out_row2 = Y[row + 2];



		for (col = ncl; col < nch + 1 - cr; col += order) {
			x0  = row0[col - 1]; x1  = row0[col + 0]; x2  = row0[col + 1]; x3  = row0[col + 2]; x4  = row0[col + 3];
			x5  = row1[col - 1]; x6  = row1[col + 0]; x7  = row1[col + 1]; x8  = row1[col + 2]; x9  = row1[col + 3];
			x10 = row2[col - 1]; x11 = row2[col + 0]; x12 = row2[col + 1]; x13 = row2[col + 2]; x14 = row2[col + 3];
			x15 = row3[col - 1]; x16 = row3[col + 0]; x17 = row3[col + 1]; x18 = row3[col + 2]; x19 = row3[col + 3];
			x20 = row4[col - 1]; x21 = row4[col + 0]; x22 = row4[col + 1]; x23 = row4[col + 2]; x24 = row4[col + 3];


			y0  = x0  | x5  | x10; y1  = x1  | x6  | x11; y2  = x2  | x7  | x12; y3  = x3  | x8  | x13; y4  = x4  | x9  | x14; 
			y5  = x5  | x10 | x15; y6  = x6  | x11 | x16; y7  = x7  | x12 | x17; y8  = x8  | x13 | x18; y9  = x9  | x14 | x19; 
			y10 = x10 | x15 | x20; y11 = x11 | x16 | x21; y12 = x12 | x17 | x22; y13 = x13 | x18 | x23; y14 = x14 | x19 | x24; 
 
			out_row0[col + 0] = y0 | y1 | y2; 
			out_row0[col + 1] = y1 | y2 | y3; 
			out_row0[col + 2] = y2 | y3 | y4;

			out_row1[col + 0] = y5 | y6 | y7;
			out_row1[col + 1] = y6 | y7 | y8;
			out_row1[col + 2] = y7 | y8 | y9;
			
			out_row2[col + 0] = y10 | y11 | y12;
			out_row2[col + 1] = y11 | y12 | y13;
			out_row2[col + 2] = y12 | y13 | y14;
		}
		row0 = row3;
		row1 = row4;
	}
	switch (rr) {
		case 2 :
			switch(cr){
				case 2 :
					ui8matrix_dilation_LU3x3_ExLU_O3_AddrRR(X, nrl, nrh, nch - 1, nch, temp_buffer, Y);
					ui8matrix_dilation_LU3x3_InLU_O3_AddrRR(X, nrh - 1, nrh, ncl, nch, temp_buffer, Y);
				break;
				case 1: 
					ui8matrix_dilation_LU3x3_ExLU_O3_AddrRR(X, nrl, nrh, nch, nch, temp_buffer, Y);		
					ui8matrix_dilation_LU3x3_InLU_O3_AddrRR(X, nrh - 1, nrh, ncl, nch, temp_buffer, Y);	
				break;
				default :
					ui8matrix_dilation_LU3x3_InLU_O3_AddrRR(X, nrh - 1, nrh, ncl, nch, temp_buffer, Y);
				break;
			}
			break;
		case 1:
			switch(cr){
				case 2 :
					ui8matrix_dilation_LU3x3_ExLU_O3_AddrRR(X, nrl, nrh, nch - 1, nch, temp_buffer, Y);
					ui8matrix_dilation_LU3x3_InLU_O3_AddrRR(X, nrh, nrh, ncl, nch, temp_buffer, Y);
				break;
				case 1: 
					ui8matrix_dilation_LU3x3_ExLU_O3_AddrRR(X, nrl, nrh, nch, nch, temp_buffer, Y);		
					ui8matrix_dilation_LU3x3_InLU_O3_AddrRR(X, nrh, nrh, ncl, nch, temp_buffer, Y);	
					// display_ui8matrix(Y, nrl, nrh, ncl, nch, "%u", "rr=1; cr=1");
				break;
				default :
					ui8matrix_dilation_LU3x3_InLU_O3_AddrRR(X, nrh, nrh, ncl, nch, temp_buffer, Y);
				break;
			}
		break;
		default:
			switch(cr){
				case 2 :
					ui8matrix_dilation_LU3x3_ExLU_O3_AddrRR(X, nrl, nrh, nch - 1, nch, temp_buffer, Y);
				break;
				case 1: 
					ui8matrix_dilation_LU3x3_ExLU_O3_AddrRR(X, nrl, nrh, nch, nch, temp_buffer, Y);	
				break;
				default :
				break;
			}
		break;
	}
		
}
/************************************************************************************/
/******** End of optimisation : Loop Unroll + Register Rotation of Addresses ********/
/************************************************************************************/

/********************************************************/
/******* Optimisation :  Pipelining ********/
/********************************************************/
void ui8matrix_dilation_row_pipeline(uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **temp_buffer, uint8 **Y)
{
	long row = nrl, col = ncl, x, y;
	// uint8 **temp = ui8matrix(nrl + (-1), nrh + 1, ncl + (-1), nch + 1);
	uint8 *temp_row0, *temp_row1, *temp_row2, *out_row0;
	// Prologue	
	temp_row0 = temp_buffer[nrl - 1];
	temp_row1 = temp_buffer[nrl + 0];
	for (col = ncl; col < nch + 1; col++) {
    	temp_row0[col] = scalar_or3(X[row - 1], col);
		temp_row1[col] = scalar_or3(X[row + 0], col);
	}	

	for (row = nrl + 1; row < nrh + 1; row++){
		out_row0 = Y[row - 1];
		temp_row0 = temp_buffer[row - 2];
		temp_row1 = temp_buffer[row - 1];
		temp_row2 = temp_buffer[row - 0];
		for (col = ncl; col < nch + 1; col++) {
			temp_row2[col] = scalar_or3(X[row + 0], col);
			out_row0[col] = temp_row0[ col] |
							temp_row1[ col] |
							temp_row2[ col];
		}
	}
	temp_row0 = temp_buffer[nrh - 1];
	temp_row1 = temp_buffer[nrh + 0];
	temp_row2 = X[nrh + 1];
	out_row0 = Y[nrh + 0];
	for (col = ncl; col < nch + 1; col++) {

    	out_row0[col] = temp_row0[ col] |
				  	    temp_row1[ col] |
		     scalar_or3(temp_row2, col);
	}
	// free_ui8matrix(temp, nrl + (-1), nrh + 1, ncl + (-1), nch + 1);
}

void ui8matrix_dilation_col_pipeline(uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **temp_buffer, uint8 **Y)
{
	long row = nrl, col = ncl, x, y;
	// uint8 **temp = ui8matrix(nrl + (-1), nrh + 1, ncl + (-1), nch + 1);
	uint8 *in_row0, *in_row1, *in_row2, *temp_row, *out_row0;
	uint8 x0, x1, x2;
	// Prologue	
	for (row = nrl; row < nrh + 1; row++) {
		in_row0 = X[row - 1];
		in_row1 = X[row + 0];
		in_row2 = X[row + 1];
		temp_row = temp_buffer[row];
    	temp_row[ncl - 1] = in_row0[ncl - 1] | 
						    in_row1[ncl - 1] | 
							in_row2[ncl - 1];
    	temp_row[ncl + 0] = in_row0[ncl + 0] | 
							in_row1[ncl + 0] | 
							in_row2[ncl + 0];
	}	
	

	for (row = nrl; row < nrh + 1; row++){
		temp_row = temp_buffer[row];
		
		in_row0 = X[row - 1];
		in_row1 = X[row + 0];
		in_row2 = X[row + 1];

		out_row0 = Y[row]; 
		
		for (col = ncl; col < nch + 1; col ++) {
			x0 = temp_row[col - 1];		// LOAD => A
			x1 = temp_row[col - 0];		// LOAD => B
			x2 = in_row0[col + 1] | in_row1[col + 1] | in_row2[col + 1];
			
			out_row0[col + 0] = x0 | x1 | x2; // A | B | C
			temp_row[col + 1] = x2;		  
		}
	}
	// free_ui8matrix(temp, nrl + (-1), nrh + 1, ncl + (-1), nch + 1);
}

void ui8matrix_dilation_col_pipeline_RR(uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **temp_buffer, uint8 **Y)
{
	long row = nrl, col = ncl, x, y;
	// uint8 **temp = ui8matrix(nrl + (-1), nrh + 1, ncl + (-1), nch + 1);
	uint8 *in_row0, *in_row1, *in_row2, *temp_row, *out_row0;
	uint8 x0, x1, x2;
	// Prologue	
	for (row = nrl; row < nrh + 1; row++) {
		in_row0 = X[row - 1];
		in_row1 = X[row + 0];
		in_row2 = X[row + 1];
		temp_row = temp_buffer[row];
    	temp_row[ncl - 1] = in_row0[ncl - 1] | 
						    in_row1[ncl - 1] | 
							in_row2[ncl - 1];
    	temp_row[ncl + 0] = in_row0[ncl + 0] | 
							in_row1[ncl + 0] | 
							in_row2[ncl + 0];
	}	
	

	for (row = nrl; row < nrh + 1; row++){
		temp_row = temp_buffer[row];
		
		in_row0 = X[row - 1];
		in_row1 = X[row + 0];
		in_row2 = X[row + 1];

		out_row0 = Y[row]; 
		x0 = temp_row[ncl - 1];		// LOAD => A
		x1 = temp_row[ncl - 0];		// LOAD => B
		for (col = ncl; col < nch + 1; col ++) {
			x2 = in_row0[col + 1] | in_row1[col + 1] | in_row2[col + 1];
			out_row0[col + 0] = x0 | x1 | x2; // A | B | C
			temp_row[col + 1] = x2;	
			x0 = x1;
			x1 = x2;	  
		}
	}
	// for (row = nrl; row < nrh + 1; row++){
	// 	temp_row = temp_buffer[row];
	// 	Y[row][nch] = temp_row[nch - 1] | 
	// 						 temp_row[nch - 0] |
	// 				 X[row - 1][nch + 1] | 
	// 				 X[row + 0][nch + 1] |
	// 				 X[row + 1][nch + 1];
	// }
	// free_ui8matrix(temp, nrl + (-1), nrh + 1, ncl + (-1), nch + 1);
}
/*************************************************/
/******* End of optimisation : Pipelining ********/
/*************************************************/


/********************************************************/
/******* Optimisation : Loop Unroll + Pipelining ********/
/********************************************************/

void ui8matrix_dilation_pipeline_LU3x3_ExLU_O3 (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **temp_buffer, uint8 **Y)
{
	const long order = 3;
	long row = nrl, col = ncl, x, y, r = 0;
	
	uint8 *row0, *row1, *row2, *row3, x0, x1, x2, x3, x4;
	uint8 *temp_row0, *temp_row1, *temp_row2, *temp_row3, *temp_row4;
	uint8 *out_row0, *out_row1, *out_row2, *out_row3;
	
	r = (nrh + 1) % order;
	// Prologue	
	for (row = nrl - 1; row < nrh - r + 1; row += order){
		row0 = X[row + 0];
		row1 = X[row + 1];
		temp_row1 = temp_buffer[row + 0];
		temp_row2 = temp_buffer[row + 1];
		for (col = ncl; col < nch + 1; col++) {
			temp_row1[col] = scalar_or3(row0, col);
			temp_row2[col] = scalar_or3(row1, col);
		}
	}

	row0 = X[row - 1];
	row1 = X[row + 0];
	row2 = X[row + 1];
	
	temp_row2 = temp_buffer[row - 1];
	temp_row3 = temp_buffer[row + 0];
	temp_row4 = temp_buffer[row + 1];
	

	for (row = nrl; row < nrh + 1 - r; row += order){
		// row0 = X[row - 1];   //  0 : b0 b1 b2 ...
		row1 = X[row + 1];   //  1 : c0 c1 c2 ...
		// row2 = X[row + 1];   //  2 : d0 d1 d2 ...
		// row3 = X[row + 2];   //  3 : e0 e1 e2 ...

		temp_row0 = temp_buffer[row - 1]; // -1 : A    <= (a0 a1 a2 ...) Prologued
		temp_row1 = temp_buffer[row + 0]; //  0 : B    <= (b0 b1 b2 ...) Prologued
		temp_row2 = temp_buffer[row + 1]; //  1 : temp_buffer <= (c0 c1 c2 ...)
		temp_row3 = temp_buffer[row + 2]; //  2 : D	<= (d0 d1 d2 ...) Prologued
		temp_row4 = temp_buffer[row + 3]; //  3 : E	<= (e0 e1 e2 ...) Prologued

		out_row0 = Y[row + 0]; 
		out_row1 = Y[row + 1];
		out_row2 = Y[row + 2];

		for (col = ncl; col < nch + 1; col++) {
			x0 = temp_row0[col];		// LOAD => A										
			x1 = temp_row1[col];		// LOAD => B
			x2 = scalar_or3(row1, col); // OR3  => C
			x3 = temp_row3[col]; 		// LOAD => D
			x4 = temp_row4[col]; 		// LOAD => E
			
			out_row0[col] = x0 | x1 | x2; // A | B | C
			out_row1[col] = x1 | x2 | x3; // B | C | D
			out_row2[col] = x2 | x3 | x4; // C | D | E

			temp_row2[col] = x2;		// STORE => C 
										// NEXT  :  F 
										// NEXT  :  G 
										// NEXT  :  H
		}
	}
	// epilogue	
	for (row = nrh - r; row < nrh + 1; row++){
		row0 = X[row + 0];
		temp_row3 = temp_buffer[row + 0];
		for (col = ncl; col < nch + 1; col++) 
			temp_row3[col] = scalar_or3(row0, col);
	}

	switch(r) {
		case 1: 
		temp_row1 = temp_buffer[nrh - 1]; //  0 : B    <= (b0 b1 b2 ...) Prologued
		temp_row2 = temp_buffer[nrh + 0]; //  1 : temp_buffer <= (c0 c1 c2 ...)
		// temp_row3 = temp_buffer[nrh + 1]; //  2 : D	<= (d0 d1 d2 ...) Prologued
		row3 = X[nrh + 1];
		out_row1 = Y[nrh + 0];
		
		for (col = ncl; col < nch + 1; col++) {											
			out_row1[col] = temp_row1[col] | temp_row2[col] | scalar_or3(row3, col); 		
		}
		
		break;
		case 2: 
		temp_row0 = temp_buffer[nrh - 2]; // -1 : A    <= (a0 a1 a2 ...) Prologued
		temp_row1 = temp_buffer[nrh - 1]; //  0 : B    <= (b0 b1 b2 ...) Prologued
		temp_row2 = temp_buffer[nrh - 0]; //  1 : temp_buffer <= (c0 c1 c2 ...)
		row3 = X[nrh + 1];
		
		out_row0 = Y[nrh - 1]; 
		out_row1 = Y[nrh + 0];
		out_row2 = Y[nrh + 1];

		for (col = ncl; col < nch + 1; col++) {									
			out_row0[col] = temp_row0[col] | temp_row1[col] | temp_row2[col]; 		
			out_row1[col] = temp_row1[col] | temp_row2[col] | scalar_or3(row3, col); 	
		}
		break;
		default:
		break;

	}
	// display_ui8matrix(X,nrl + (-1), nrh + 1, ncl + (-1), nch + 1, "%u", "In");
	// display_ui8matrix(Y,nrl, nrh , ncl, nch , "%03u ", "Out");
	// display_ui8matrix(temp,nrl + (-1), nrh + 1, ncl, nch, "%u ", "Out");
	// free_ui8matrix(temp,nrl + (-1), nrh + 1, ncl, nch);
}
void ui8matrix_dilation_pipeline_LU3x3_InLU_O3 (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **temp_buffer, uint8 **Y)
{
	const long order = 3;
	long row = nrl, col = ncl, x, y, r = 0;
	// uint8 **temp = ui8matrix(nrl, nrh, ncl + (-1), nch + 1);
	uint8 *row0, *row1, *row2, *row3, x0, x1, x2, x3, x4, x5;
	uint8 *temp_row0, *temp_row1, *temp_row2, *temp_row3, *temp_row4;
	uint8 *out_row0, *out_row1, *out_row2, *out_row3;
	
	r = (nch + 1) % order;
	// Prologue	
		
	for (row = nrl; row < nrh + 1; row++) {
		row0 = X[row - 1];
		row1 = X[row + 0];
		row2 = X[row + 1];

		temp_row1 = temp_buffer[row + 0];
		for (col = ncl - 1; col < nch - r + 1 ; col += order) {
			temp_row1[col + 0] = row0[col + 0] | row1[col + 0] | row2[col + 0];
			temp_row1[col + 1] = row0[col + 1] | row1[col + 1] | row2[col + 1];
			
		}
	}

	for (row = nrl; row < nrh + 1; row ++){
		temp_row1 = temp_buffer[row + 0];
		row0 = X[row - 1];
		row1 = X[row + 0];
		row2 = X[row + 1];

		out_row1 = Y[row]; 

		for (col = ncl; col < nch + 1 - r; col += order) {
			x0 = temp_row1[col - 1];
			x1 = temp_row1[col + 0];
			x2 = row0[col + 1] | row1[col + 1] | row2[col + 1];
			x3 = temp_row1[col + 2];
			x4 = temp_row1[col + 3];
			
			out_row1 [col + 0] = x0 | x1 | x2; // A | B | C
			out_row1 [col + 1] = x1 | x2 | x3; // B | C | D
			out_row1 [col + 2] = x2 | x3 | x4; // C | D | E
			temp_row1[col] = x2;

		}
	}
	// epilogue	
	for (row = nrl; row < nrh + 1; row++){
		row0 = X[row - 1];
		row1 = X[row + 0];
		row2 = X[row + 1];
		temp_row3 = temp_buffer[row + 0];
		for (col = nch - r; col < nch + 1; col++) 
			temp_row3[col] = row0[col] | row1[col] | row2[col];
	}

	switch(r){
		case 2:
		for (row = nrl; row < nrh + 1; row++) {											
			temp_row0 = temp_buffer[row];
			row0 = X[row - 1];
			row1 = X[row + 0];
			row2 = X[row + 1];
			out_row0 = Y[row];
			out_row0[nch - 1] = temp_row0[nch - 2] | 
								temp_row0[nch - 1] | 
								temp_row0[nch - 0] ;

			out_row0[nch + 0] = temp_row0[nch - 1] | 
								temp_row0[nch - 0] |
								     row0[nch + 1] | 	
									 row1[nch + 1] | 	
									 row2[nch + 1]; 	
		}
		break;
		case 1:
		
		for (row = nrl; row < nrh + 1; row++) {											
			temp_row0 = temp_buffer[row];
			row0 = X[row - 1];
			row1 = X[row + 0];
			row2 = X[row + 1];
			out_row0 = Y[row];
			out_row0[nch + 0] = temp_row0[nch - 1] | 
								temp_row0[nch - 0] |
								     row0[nch + 1] | 	
									 row1[nch + 1] | 	
									 row2[nch + 1]; 	
		}
		break;
		default:
		break;
	}
	// display_ui8matrix(X,nrl + (-1), nrh + 1, ncl + (-1), nch + 1, "%u", "In");
	// display_ui8matrix(Y,nrl, nrh, ncl, nch , "%03u ", "Out");
	// display_ui8matrix(temp,nrl, nrh, ncl + (-1), nch + 1, "%u ", "Temp");

	// ui8matrix_dilation_LU3x3_O1xO1(X, nrl, nrh, ncl, nch, temp_buffer, Y);
	// free_ui8matrix(temp, nrl, nrh, ncl + (-1), nch + 1);
	
}

/***************************************************************/
/******* End of optimisation : Loop Unroll + Pipelining ********/
/***************************************************************/

/****************************************************************************/
/******* Optimisation : Loop Unroll + Register Rotation + Pipelining ********/
/****************************************************************************/

void ui8matrix_dilation_pipeline_LU3x3_ExLU_O3_RR (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **temp_buffer, uint8 **Y)
{
	const long order = 3;
	long row = nrl, col = ncl, x, y, r = 0;
	
	uint8 *row0, *row1, *row2, *row3;
	uint8 y0, y1, y2, y3, y4, x0, x1, x2, x3, x4, x5;
	uint8 *temp_row0, *temp_row1, *temp_row2, *temp_row3, *temp_row4;
	uint8 *out_row0, *out_row1, *out_row2, *out_row3;
	
	r = (nrh + 1) % order;
	// Prologue	
	for (row = nrl - 1; row < nrh - r + 1; row += order){
		row0 = X[row + 0];
		row1 = X[row + 1];
		temp_row1 = temp_buffer[row + 0];
		temp_row2 = temp_buffer[row + 1];
		x0  = row0[ncl - 1]; x1 = row0[ncl + 0];
		x3  = row1[ncl - 1]; x4 = row1[ncl + 0];
		for (col = ncl; col < nch + 1; col++) {
			x2  = row0[col + 1];
			x5  = row1[col + 1];
			temp_row1[col] = scalar_or3(row0, col);
			temp_row2[col] = scalar_or3(row1, col);
			x0 = x1; x1 = x2;
			x3 = x4; x4 = x5;
		}
	}

	row0 = X[row - 1];
	row1 = X[row + 0];
	row2 = X[row + 1];
	
	temp_row2 = temp_buffer[row - 1];
	temp_row3 = temp_buffer[row + 0];
	temp_row4 = temp_buffer[row + 1];
	

	for (row = nrl; row < nrh + 1 - r; row += order){
		// row0 = X[row - 1];   //  0 : b0 b1 b2 ...
		row1 = X[row + 1];   //  1 : c0 c1 c2 ...
		// row2 = X[row + 1];   //  2 : d0 d1 d2 ...
		// row3 = X[row + 2];   //  3 : e0 e1 e2 ...

		temp_row0 = temp_buffer[row - 1]; // -1 : A    <= (a0 a1 a2 ...) Prologued
		temp_row1 = temp_buffer[row + 0]; //  0 : B    <= (b0 b1 b2 ...) Prologued
		temp_row2 = temp_buffer[row + 1]; //  1 : temp_buffer <= (c0 c1 c2 ...)
		temp_row3 = temp_buffer[row + 2]; //  2 : D	<= (d0 d1 d2 ...) Prologued
		temp_row4 = temp_buffer[row + 3]; //  3 : E	<= (e0 e1 e2 ...) Prologued

		out_row0 = Y[row + 0]; 
		out_row1 = Y[row + 1];
		out_row2 = Y[row + 2];

		x3  = row1[ncl - 1]; x4  = row1[ncl + 0]; 
		for (col = ncl; col < nch + 1; col++) {
			x5  = row1[col + 1];
			y0 = temp_row0[col];		// LOAD => A										
			y1 = temp_row1[col];		// LOAD => B
			y2 = x3 | x4 | x5;
			y3 = temp_row3[col]; 		// LOAD => D
			y4 = temp_row4[col]; 		// LOAD => E
			

			out_row0[col] = y0 | y1 | y2; // A | B | C
			out_row1[col] = y1 | y2 | y3; // B | C | D
			out_row2[col] = y2 | y3 | y4; // C | D | E

			temp_row2[col] = y2;		// STORE => C 
										// NEXT  :  F 
										// NEXT  :  G 
										// NEXT  :  H
			x3 = x4; x4 = x5;
		}
	}
	// epilogue	
	for (row = nrh - r; row < nrh + 1; row++){
		row0 = X[row + 0];
		temp_row3 = temp_buffer[row + 0];
		for (col = ncl; col < nch + 1; col++) 
			temp_row3[col] = scalar_or3(row0, col);
	}

	switch(r) {
		case 1: 
		temp_row1 = temp_buffer[nrh - 1]; //  0 : B    <= (b0 b1 b2 ...) Prologued
		temp_row2 = temp_buffer[nrh + 0]; //  1 : temp_buffer <= (c0 c1 c2 ...)
		// temp_row3 = temp_buffer[nrh + 1]; //  2 : D	<= (d0 d1 d2 ...) Prologued
		row3 = X[nrh + 1];
		out_row1 = Y[nrh + 0];
		
		for (col = ncl; col < nch + 1; col++) {											
			out_row1[col] = temp_row1[col] | temp_row2[col] | scalar_or3(row3, col); 		
		}
		
		break;
		case 2: 
		temp_row0 = temp_buffer[nrh - 2]; // -1 : A    <= (a0 a1 a2 ...) Prologued
		temp_row1 = temp_buffer[nrh - 1]; //  0 : B    <= (b0 b1 b2 ...) Prologued
		temp_row2 = temp_buffer[nrh - 0]; //  1 : temp_buffer <= (c0 c1 c2 ...)
		row3 = X[nrh + 1];
		
		out_row0 = Y[nrh - 1]; 
		out_row1 = Y[nrh + 0];
		out_row2 = Y[nrh + 1];

		for (col = ncl; col < nch + 1; col++) {									
			out_row0[col] = temp_row0[col] | temp_row1[col] | temp_row2[col]; 		
			out_row1[col] = temp_row1[col] | temp_row2[col] | scalar_or3(row3, col); 	
		}
		break;
		default:
		break;

	}
	// free_ui8matrix(temp,nrl + (-1), nrh + 1, ncl, nch);
}
void ui8matrix_dilation_pipeline_LU3x3_InLU_O3_RR (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **temp_buffer, uint8 **Y)
{
	const long order = 3;
	long row = nrl, col = ncl, x, y, r = 0;
	// uint8 **temp = ui8matrix(nrl, nrh, ncl + (-1), nch + 1);
	uint8 *row0, *row1, *row2, *row3, x0, x1, x2, x3, x4, x5;
	uint8 *temp_row0, *temp_row1, *temp_row2, *temp_row3, *temp_row4;
	uint8 *out_row0, *out_row1, *out_row2, *out_row3;
	
	r = (nch + 1) % order;
	// Prologue	
		
	row0 = X[nrl - 1];
	row1 = X[nrl + 0];
	for (row = nrl; row < nrh + 1; row++) {
		row2 = X[row + 1];

		temp_row1 = temp_buffer[row + 0];
		for (col = ncl - 1; col < nch - r + 1 ; col += order) {
			temp_row1[col + 0] = row0[col + 0] | row1[col + 0] | row2[col + 0];
			temp_row1[col + 1] = row0[col + 1] | row1[col + 1] | row2[col + 1];
			
		}
		row0 = row1;
		row1 = row2;
	}

	row0 = X[nrl - 1];
	row1 = X[nrl + 0];
	for (row = nrl; row < nrh + 1; row ++){
		temp_row1 = temp_buffer[row + 0];
		row2 = X[row + 1];

		out_row1 = Y[row]; 

		x0 = temp_row1[ncl - 1];
		x1 = temp_row1[ncl + 0];
		for (col = ncl; col < nch + 1 - r; col += order) {
			x2 = row0[col + 1] | row1[col + 1] | row2[col + 1];
			x3 = temp_row1[col + 2];
			x4 = temp_row1[col + 3];
			
			out_row1 [col + 0] = x0 | x1 | x2; // A | B | C
			out_row1 [col + 1] = x1 | x2 | x3; // B | C | D
			out_row1 [col + 2] = x2 | x3 | x4; // C | D | E
			temp_row1[col] = x2;

			x0 = x3;
			x1 = x4;

		}
		row0 = row1;
		row1 = row2;
	}
	// epilogue	
	for (row = nrl; row < nrh + 1; row++){
		row0 = X[row - 1];
		row1 = X[row + 0];
		row2 = X[row + 1];
		temp_row3 = temp_buffer[row + 0];
		for (col = nch - r; col < nch + 1; col++) 
			temp_row3[col] = row0[col] | row1[col] | row2[col];
	}

	switch(r){
		case 2:
		row0 = X[nrl - 1];
		row1 = X[nrl + 0];
		for (row = nrl; row < nrh + 1; row++) {											
			temp_row0 = temp_buffer[row];
			row2 = X[row + 1];
			out_row0 = Y[row];
			out_row0[nch - 1] = temp_row0[nch - 2] | 
								temp_row0[nch - 1] | 
								temp_row0[nch - 0] ;

			out_row0[nch + 0] = temp_row0[nch - 1] | 
								temp_row0[nch - 0] |
								     row0[nch + 1] | 	
									 row1[nch + 1] | 	
									 row2[nch + 1]; 
			row0 = row1;
			row1 = row2;	
		}
		break;
		case 1:
		
		row0 = X[nrl - 1];
		row1 = X[nrl + 0];
		for (row = nrl; row < nrh + 1; row++) {											
			temp_row0 = temp_buffer[row];
			row0 = X[row - 1];
			row1 = X[row + 0];
			row2 = X[row + 1];
			out_row0 = Y[row];
			out_row0[nch + 0] = temp_row0[nch - 1] | 
								temp_row0[nch - 0] |
								     row0[nch + 1] | 	
									 row1[nch + 1] | 	
									 row2[nch + 1]; 	
			row0 = row1;
			row1 = row2;
		}
		break;
		default:
		break;
	}
	// display_ui8matrix(X,nrl + (-1), nrh + 1, ncl + (-1), nch + 1, "%u", "In");
	// display_ui8matrix(Y,nrl, nrh, ncl, nch , "%03u ", "Out");
	// display_ui8matrix(temp,nrl, nrh, ncl + (-1), nch + 1, "%u ", "Temp");

	// ui8matrix_dilation_LU3x3_O1xO1(X, nrl, nrh, ncl, nch, temp_buffer, Y);
	// free_ui8matrix(temp, nrl, nrh, ncl + (-1), nch + 1);
	
}

void ui8matrix_dilation5_divide_row_and_conquer (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **Y, uint8 **Z)
{

	uint8 *temp_row, *temp_row0, *temp_row1, *temp_row2, *temp_row3, *temp_row4, *out_row, *in_row;
	long row, col, nrow, r;

	for (row = nrl - 2; row < (nrh + 1) + 2; row ++) {
		temp_row = Y[row];
		in_row   = X[row];
		for(col = ncl - 1; col < (nch + 1) + 1; col++) {
			temp_row[col] = scalar_or5(in_row, col);
		}
	}

	for (row = nrl; row < nrh + 1; row ++)	{
		temp_row0 = Y[row - 2];
		temp_row1 = Y[row - 1];
		temp_row2 = Y[row + 0];
		temp_row3 = Y[row + 1];
		temp_row4 = Y[row + 2];

		out_row   = Z[row + 0];
		for(col = ncl; col < nch + 1; col++) {
			out_row[col] = temp_row0[col]|
						temp_row1[col]|
						temp_row2[col]|
						temp_row3[col]|
						temp_row4[col];
		}
	}
}
void ui8matrix_dilation_divide_row_and_conquer (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **temp_buffer, uint8 **Y)
{
	const long order = 3;
	
	uint8 *temp_row, *temp_row0, *temp_row1, *temp_row2, *out_row, *in_row;
	long row, col, nrow, r;
	int nb_threads;

	nrow = nrh - nrl + 1;
	row = 0;
	
	for (row = nrl - 1; row < (nrh + 1) + 1; row ++)
	{
		temp_row =    temp_buffer[row];
		in_row   = X[row];
		for(col = ncl; col < nch + 1; col++) {
			temp_row[col] = scalar_or3(in_row, col);
		}
	}
	
	for (row = nrl; row < nrh + 1; row ++)
	{
		temp_row0 =     temp_buffer[row - 1];
		temp_row1 =     temp_buffer[row + 0];
		temp_row2 =     temp_buffer[row + 1];
		out_row   = Y[row];
		for(col = ncl; col < nch + 1; col++) {
			out_row[col] = temp_row0[col] |
						   temp_row1[col] |
						   temp_row2[col];
		}
	}
}

/** Deprcated **/
void ui8matrix_dilation_unrolled_divide_row_and_conquer (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **temp_buffer, uint8 **Y)
{
	const long order = 3;
	uint8 *temp_row, *temp_row0, *temp_row1, *temp_row2, *out_row, *in_row;
	long row, col, nrow, r;
	int nb_threads;

	nrow = nrh - nrl + 1;
	row = 0;
	
	for (row = nrl - 1; row < nrh + 1; row += order)
	{
		temp_row =    temp_buffer[row];
		in_row   = X[row];
		for(col = ncl; col < nch + 1; col++) {
			temp_row[col] = scalar_or3(in_row, col);
		}
	}
	for (row = nrl + 0; row < nrh + 1; row += order)
	{

		temp_row =    temp_buffer[row];
		in_row   = X[row];
		for(col = ncl; col < nch + 1; col++) {
			temp_row[col] = scalar_or3(in_row, col);
		}
	}
	for (row = nrl + 1; row < (nrh + 1) + 1; row += order)
	{
		temp_row =    temp_buffer[row];
		in_row   = X[row];
		for(col = ncl; col < nch + 1; col++) {
			temp_row[col] = scalar_or3(in_row, col);
		}
	}
	// Epilogue
	temp_row = temp_buffer[nrh + 1];
	in_row = X[nrh + 1];
	for(col = ncl; col < nch + 1; col++) {
		temp_row[col] = scalar_or3(in_row, col);
	}
	
	for (row = nrl; row < nrh + 1; row ++)
	{
		temp_row0 =     temp_buffer[row - 1];
		temp_row1 =     temp_buffer[row + 0];
		temp_row2 =     temp_buffer[row + 1];
		out_row   = Y[row];
		for(col = ncl; col < nch + 1; col++) {
			out_row[col] = temp_row0[col]|
						   temp_row1[col]|
						   temp_row2[col];
		}
	}
}
void ui8matrix_dilation_divide_row_and_conquer_ExLU_O3 (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **temp_buffer, uint8 **Y)
{
	const long order = 3;
	
	uint8 *temp_row, *temp_row0, *temp_row1, *temp_row2, *temp_row3, *temp_row4, *in_row, *in_row0, *in_row1, *in_row2, *out_row;
	long row, col, nrow, r;
	int nb_threads;

	nrow = nrh - nrl + 1;
	row = 0;

	r = (nrh + 1) % order;
	temp_row0 =    temp_buffer[nrl - 1];
	in_row0   = X[nrl - 1];
	for(col = ncl; col < nch + 1; col++) {
		temp_row0[col] = scalar_or3(in_row0, col);
	}
	for (row = nrl; row < nrh + 1 - r; row += order)
	{
		temp_row0 =    temp_buffer[row + 0];
		temp_row1 =    temp_buffer[row + 1];
		temp_row2 =    temp_buffer[row + 2];
		in_row0    = X[row + 0];
		in_row1    = X[row + 1];
		in_row2    = X[row + 2];
		for(col = ncl; col < nch + 1; col++) {
			temp_row0[col] = scalar_or3(in_row0, col);
			temp_row1[col] = scalar_or3(in_row1, col);
			temp_row2[col] = scalar_or3(in_row2, col);
		}
	}
	// Epilogue
	switch(r) {
		case 2 :
		temp_row0  = temp_buffer[nrh - 1];
		in_row0 = X[nrh - 1];
		temp_row1  = temp_buffer[nrh + 0];
		in_row1 = X[nrh + 0];
		temp_row2  = temp_buffer[nrh + 1];
		in_row2 = X[nrh + 1];
		for(col = ncl; col < nch + 1; col++) {
			temp_row0[col] = scalar_or3(in_row0, col);
			temp_row1[col] = scalar_or3(in_row1, col);
			temp_row2[col] = scalar_or3(in_row2, col);
		}
		break;
		
		case 1 :
		temp_row1  = temp_buffer[nrh + 0];
		in_row1 = X[nrh + 0];
		temp_row2  = temp_buffer[nrh + 1];
		in_row2 = X[nrh + 1];
		for(col = ncl; col < nch + 1; col++) {
			temp_row1[col] = scalar_or3(in_row1, col);
			temp_row2[col] = scalar_or3(in_row2, col);
		}
		break;
		case 0 :
		temp_row2  = temp_buffer[nrh + 1];
		in_row2 = X[nrh + 1];
		for(col = ncl; col < nch + 1; col++) {
			temp_row2[col] = scalar_or3(in_row2, col);
		}
		break;
	}
	// display_ui8matrix(temp, nrl + (-1), nrh + 1, ncl, nch, "%u", "Temp");
	for (row = nrl; row < nrh + 1; row ++)
	{
		temp_row0 =     temp_buffer[row - 1];
		temp_row1 =     temp_buffer[row + 0];
		temp_row2 =     temp_buffer[row + 1];
		out_row   = Y[row];
		for(col = ncl; col < nch + 1; col++) {
			out_row[col] = temp_row0[col]|
						   temp_row1[col]|
						   temp_row2[col];
		}
	}

	
	// free_ui8matrix(temp, nrl + (-1), nrh + 1, ncl, nch);
}
void ui8matrix_dilation_divide_row_and_conquer_InLU_O3 (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **temp_buffer, uint8 **Y)
{
	const long order = 3;
	
	uint8 *temp_row, *temp_row0, *temp_row1, *temp_row2, *temp_row3, *temp_row4, *out_row0, *out_row1, *out_row2, *in_row, *out_row;
	long row, col, nrow, r;
	int nb_threads;

	nrow = nrh - nrl + 1;
	row = 0;
	r = (nch + 1) % order;
	for (row = nrl - 1; row < (nrh + 1) + 1; row ++)
	{
		temp_row =    temp_buffer[row];
		in_row   = X[row];
		for(col = ncl; col < nch + 1 - r; col += order) {
			temp_row[col + 0] = scalar_or3(in_row, col + 0);
			temp_row[col + 1] = scalar_or3(in_row, col + 1);
			temp_row[col + 2] = scalar_or3(in_row, col + 2);
		}
		
	}
	switch(r) {
		case 2:
		for (row = nrl - 1; row < (nrh + 1) + 1; row ++)
		{
			temp_row =    temp_buffer[row];
			in_row   = X[row];
			temp_row[nch - 1] = scalar_or3(in_row, nch - 1);
			temp_row[nch + 0] = scalar_or3(in_row, nch + 0);
		}
		break;
		case 1:
		for (row = nrl - 1; row < (nrh + 1) + 1; row ++)
		{
			temp_row =    temp_buffer[row];
			in_row   = X[row];
			temp_row[nch + 0] = scalar_or3(in_row, nch + 0);
		}
		break;
		default: 
		break;
	}
	// Epilogue
	// temp_row = temp_buffer[nrh + 1];
	// in_row = X[nrh + 1];
	// for(col = ncl; col < nch + 1; col++) {
	// 	temp_row[col] = scalar_or3(in_row, col);
	// }
	for (row = nrl; row < nrh + 1; row ++)
	{
		temp_row0 =     temp_buffer[row - 1];
		temp_row1 =     temp_buffer[row + 0];
		temp_row2 =     temp_buffer[row + 1];
		out_row   = Y[row];
		for(col = ncl; col < nch + 1; col++) {
			out_row[col] = temp_row0[col]|
						   temp_row1[col]|
						   temp_row2[col];
		}
	}
	
		
	// free_ui8matrix(temp, nrl + (-1), nrh + 1, ncl, nch);
}
void ui8matrix_dilation_divide_col_and_conquer (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **temp_buffer, uint8 **Y)
{
	const long order = 3;
	// uint8 **temp = ui8matrix(nrl, nrh, ncl + (-1), nch + 1);
	uint8 *temp_row, *temp_row0, *temp_row1, *temp_row2, *out_row, *in_row0, *in_row1, *in_row2;
	long row, col, nrow, r;
	int nb_threads;

	nrow = nrh - nrl + 1;
	row = 0;
	
	for (row = nrl; row < nrh + 1; row ++)
	{
		temp_row =    temp_buffer[row];
		in_row0   = X[row - 1];
		in_row1   = X[row + 0];
		in_row2   = X[row + 1];
		for(col = ncl + (-1); col < nch + 1 + 1; col++) 
			temp_row[col] = in_row0[col] | in_row1[col] | in_row2[col];
	}
	
	for (row = nrl; row < nrh + 1; row ++)
	{
		temp_row = temp_buffer[row];
		out_row  = Y[row];
		for(col = ncl; col < nch + 1; col++) 
			out_row[col] = scalar_or3(temp_row, col);
	}
	// free_ui8matrix(temp, nrl, nrh, ncl + (-1), nch + 1);
}
void ui8matrix_dilation_divide_col_and_conquer_ExLU_O3 (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **temp_buffer, uint8 **Y)
{
	const long order = 3;
	// uint8 **temp = ui8matrix(nrl, nrh, ncl + (-1), nch + 1);
	uint8 *temp_row, *temp_row0, *temp_row1, *temp_row2, *temp_row3, *out_row0,*out_row1,*out_row2, *in_row0, *in_row1, *in_row2;
	long row, col, nrow, r;
	int nb_threads;

	nrow = nrh - nrl + 1;
	row = 0;
	for (row = nrl; row < nrh + 1; row ++)
	{
		temp_row =    temp_buffer[row];
		in_row0   = X[row - 1];
		in_row1   = X[row + 0];
		in_row2   = X[row + 1];
		for(col = ncl + (-1); col < nch + 1 + 1; col++) 
			temp_row[col] = in_row0[col] | in_row1[col] | in_row2[col];
	}
	
	r = (nrh + 1) % order;
	for (row = nrl; row < nrh + 1 - r; row += order)
	{
		temp_row0 =     temp_buffer[row + 0];
		temp_row1 =     temp_buffer[row + 1];
		temp_row2 =     temp_buffer[row + 2];
		out_row0   = Y[row + 0];
		out_row1   = Y[row + 1];
		out_row2   = Y[row + 2];

		for(col = ncl; col < nch + 1; col++) {
			out_row0[col] = scalar_or3(temp_row0, col);
			out_row1[col] = scalar_or3(temp_row1, col);
			out_row2[col] = scalar_or3(temp_row2, col);
		}
	}

	switch(r){
		case 2:
		temp_row1 =     temp_buffer[nrh - 1];
		temp_row2 =     temp_buffer[nrh + 0];

		out_row0  = Y[nrh - 1];
		out_row1  = Y[nrh - 0];

		for(col = ncl; col < nch + 1; col++) {
			out_row0[col] = scalar_or3(temp_row1, col);
			out_row1[col] = scalar_or3(temp_row2, col);
		}
		break;
		case 1:
		temp_row2 =     temp_buffer[nrh + 0];
		out_row1  = Y[nrh - 0];

		for(col = ncl; col < nch + 1; col++) {
			out_row1[col] = scalar_or3(temp_row2, col);
		}
		break;
		case 0:
		break;
	}
	// free_ui8matrix(temp, nrl, nrh, ncl + (-1), nch + 1);
}
void ui8matrix_dilation_divide_col_and_conquer_InLU_O3 (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **temp_buffer, uint8 **Y)
{
	const long order = 3;
	// uint8 **temp = ui8matrix(nrl, nrh, ncl + (-1), nch + 1);
	uint8 *temp_row, *temp_row0, *temp_row1, *temp_row2, *out_row, *in_row0, *in_row1, *in_row2;
	long row, col, nrow, r;
	int nb_threads;

	// for (row = nrl; row < nrh + 1; row ++)
	nrow = nrh - nrl + 1;
	row = 0;
	for (row = nrl; row < nrh + 1; row ++)
	{
		temp_row =    temp_buffer[row];
		in_row0   = X[row - 1];
		in_row1   = X[row + 0];
		in_row2   = X[row + 1];
		for(col = ncl + (-1); col < nch + 1 + 1; col++) 
			temp_row[col] = in_row0[col] | in_row1[col] | in_row2[col];
	}
	r = (nch + 1) % order;
	for (row = nrl; row < nrh + 1; row ++)
	{
		temp_row = temp_buffer[row];
		out_row  = Y[row];
		for(col = ncl; col < nch + 1 - r; col+= order) {
			out_row[col + 0] = scalar_or3(temp_row, col + 0);
			out_row[col + 1] = scalar_or3(temp_row, col + 1);
			out_row[col + 2] = scalar_or3(temp_row, col + 2);
		}
	}

	switch(r){
		case 2:
		for (row = nrl; row < nrh + 1; row ++) {
			temp_row = temp_buffer[row];
			out_row  = Y[row];
			out_row[nch - 1] = scalar_or3(temp_row, nch - 1);
			out_row[nch + 0] = scalar_or3(temp_row, nch + 0);
		}
		break;
		case 1:
		for (row = nrl; row < nrh + 1; row ++) {
			temp_row = temp_buffer[row];
			Y[row][nch + 0] = scalar_or3(temp_row, nch + 0);
		}
		break;
		case 0:
		break;
	}
	// display_ui8matrix(temp, nrl, nrh, ncl + (-1), nch + 1, "%u", "temp");
	// free_ui8matrix(temp, nrl, nrh, ncl + (-1), nch + 1);
}
void ui8matrix_dilation_pipeline2_LU3x3_ExLU_O3 (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **temp_buffer, uint8 **Y)
{
	const long order = 3;
	long row = nrl, col = ncl, x, y, r = 0;
	// uint8 **temp = ui8matrix(nrl + (-1), nrh + 1, ncl + (-1), nch + 1);
	uint8 *row0, *row1, *row2, *row3, x0, x1, x2, x3;
	uint8 *temp_row0, *temp_row1, *temp_row2, *temp_row3, *temp_row4;
	uint8 *out_row0, *out_row1, *out_row2, *out_row3;
	
	r = (nrh + 1) % order;
	// Prologue	
	row0 = X[row - 1];
	row1 = X[row + 0];
	row2 = X[row + 1];
	
	temp_row1 = temp_buffer[row - 1];
	temp_row2 = temp_buffer[row + 0];
	temp_row3 = temp_buffer[row + 1];
	for (col = ncl; col < nch + 1; col++) {
		temp_row1[col] = scalar_or3(row0, col);
		temp_row2[col] = scalar_or3(row1, col);
	}

	for (row = row + 1; row < nrh + 1 - r; row += order){
		row0 = X[row - 1];   //  0 : b0 b1 b2 ...
		row1 = X[row + 0];   //  1 : c0 c1 c2 ...
		row2 = X[row + 1];   //  2 : d0 d1 d2 ...
		row3 = X[row + 2];   //  3 : e0 e1 e2 ...

		temp_row0 = temp_buffer[row - 2]; // -1 : A    <= (a0 a1 a2 ...) Prologued
		temp_row1 = temp_buffer[row - 1]; //  0 : B    <= (b0 b1 b2 ...) Prologued
		temp_row2 = temp_buffer[row + 0]; //  1 : temp_buffer <= (c0 c1 c2 ...)
		temp_row3 = temp_buffer[row + 1]; //  2 : temp_buffer <= (d0 d1 d2 ...)
		temp_row4 = temp_buffer[row + 2]; //  3 : temp_buffer <= (e0 e1 e2 ...)
		out_row0 = Y[row - 1]; 
		out_row1 = Y[row + 0];
		out_row2 = Y[row + 1];

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
			row1 = X[nrh - 1];   // -1 : (f0 f1 f2 ...)
			row2 = X[nrh + 0];   //  0 : (g0 g1 g2 ...)
			row3 = X[nrh + 1];   //  1 : (h0 h1 h2 ...)
			temp_row0 = temp_buffer[nrh - 2]; //  E    <= (e0 e1 e2 ...)
			out_row0 = Y[nrh - 1];
			out_row1 = Y[nrh + 0];
			for (col = ncl; col < nch + 1; col++) {
				//  temp_row0[ col]		   		 LOAD => E		
				x0 = scalar_or3(row1, col);	  // OR3  => F	
				x1 = scalar_or3(row2, col);	  // OR3  => G
				out_row0[col] = temp_row0[ col] | x0 | x1; // A | B | C
				out_row1[col] = 			 x0 | x1 | scalar_or3(row3, col); // B | C | D	
			}
			break;
		case 1: 
			row2 = X[nrh + 0];   //  0 : (f0 f1 f2 ...)
			row3 = X[nrh + 1];   //  1 : (g0 g1 g2 ...)
			temp_row0 = temp_buffer[nrh - 1]; //  E    <= (e0 e1 e2 ...)
			out_row1 = Y[nrh + 0];
			for (col = ncl; col < nch + 1; col++) {
				out_row1[col] = temp_row0[ col]  | scalar_or3(row2, col) | scalar_or3(row3, col); // E | F | G	
			}
			break;
		default:
			break;
	}
	// free_ui8matrix(temp, nrl + (-1), nrh + 1, ncl + (-1), nch + 1);
}
void ui8matrix_dilation_pipeline2_LU3x3_InLU_O3 (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **temp_buffer, uint8 **Y)
{
	const long order = 3;
	long row = nrl, col = ncl, x, y, r = 0;
	// uint8 **temp = ui8matrix(nrl + (-1), nrh + 1, ncl + (-1), nch + 1);
	uint8 *row0, *row1, *row2, *row3, x0, x1, x2, x3, x4, x5, y0, y1, y2;
	uint8 *temp_row0, *temp_row1, *temp_row2, *temp_row3, *temp_row4;
	uint8 *out_row0, *out_row1, *out_row2, *out_row3;
	
	r = (nch + 1) % order;
	// Prologue	
	row0 = X[row - 1];
	row1 = X[row + 0];
	row2 = X[row + 1];
	
	// temp_row1 = temp_buffer[row - 1];
	// temp_row2 = temp_buffer[row + 0];
	// temp_row3 = temp_buffer[row + 1];
	for (row = nrl; row < nrh + 1; row++) {
		row0 = X[row - 1];
		row1 = X[row + 0];
		row2 = X[row + 1];
		temp_row1 = temp_buffer[row + 0];

		temp_row1[ncl - 1] = row0[ncl - 1] | row1[ncl - 1] | row2[ncl - 1];
		temp_row1[ncl + 0] = row0[ncl + 0] | row1[ncl + 0] | row2[ncl + 0];
	}

	for (row = nrl; row < nrh + 1; row ++){
		temp_row1 = temp_buffer[row + 0];
		
		row0 = X[row - 1];
		row1 = X[row + 0];
		row2 = X[row + 1];

		out_row1 = Y[row]; 
		
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
				temp_row1 = temp_buffer[row + 0];
				
				row0 = X[row - 1];
				row1 = X[row + 0];
				row2 = X[row + 1];

				out_row0 = Y[row]; 
				
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
			// ui8matrix_dilation_pipeline3x3(X, nrl, nrh, nch - 1, nch, temp_buffer, Y);
			break;
		case 1: 
			for (row = nrl; row < nrh + 1; row ++){
				temp_row1 = temp_buffer[row + 0];
				
				row0 = X[row - 1];
				row1 = X[row + 0];
				row2 = X[row + 1];

				out_row0 = Y[row]; 
				
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
	// display_ui8matrix(Y, nrl - 1, nrh + 1, ncl - 1, nch + 1, "%u", "TEst");
	// free_ui8matrix(temp, nrl + (-1), nrh + 1, ncl + (-1), nch + 1);
	
}


void ui8matrix_dilation_pipeline2_LU3x3_ExLU_O3_RR (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **temp_buffer, uint8 **Y)
{
	const long order = 3;
	long row = nrl, col = ncl, r = 0;
	// uint8 **temp = ui8matrix(nrl + (-1), nrh + 1, ncl + (-1), nch + 1);
	uint8 *row0, *row1, *row2, *row3, x0, x1, x2, x3;
	uint8 y0, y1, y2, y3, y4, y5, y6, y7, y8, y9, y10, y11, y12, y13, y14;
	uint8 *temp_row0, *temp_row1, *temp_row2, *temp_row3, *temp_row4;
	uint8 *out_row0, *out_row1, *out_row2, *out_row3;
	
	r = (nrh + 1) % order;
	// Prologue	
	row0 = X[row - 1];
	row1 = X[row + 0];
	row2 = X[row + 1];
	
	temp_row1 = temp_buffer[row - 1];
	temp_row2 = temp_buffer[row + 0];
	temp_row3 = temp_buffer[row + 1];
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

	// row0 = X[row + 0];   //  0 : b0 b1 b2 ...
	temp_row0 = temp_buffer[row - 1]; // -1 : A    <= (a0 a1 a2 ...) Prologued
	temp_row1 = temp_buffer[row - 0]; //  0 : B    <= (b0 b1 b2 ...) Prologued
	for (row = row + 1; row < nrh + 1 - r; row += order){
		row1 = X[row + 0];   //  1 : c0 c1 c2 ...
		row2 = X[row + 1];   //  2 : d0 d1 d2 ...
		row3 = X[row + 2];   //  3 : e0 e1 e2 ...

		temp_row2 = temp_buffer[row + 0]; //  1 : temp_buffer <= (c0 c1 c2 ...)
		temp_row3 = temp_buffer[row + 1]; //  2 : temp_buffer <= (d0 d1 d2 ...)
		temp_row4 = temp_buffer[row + 2]; //  3 : temp_buffer <= (e0 e1 e2 ...)
		out_row0 = Y[row - 1]; 
		out_row1 = Y[row + 0];
		out_row2 = Y[row + 1];

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
			row1 = X[nrh - 1];   // -1 : (f0 f1 f2 ...)
			row2 = X[nrh + 0];   //  0 : (g0 g1 g2 ...)
			row3 = X[nrh + 1];   //  1 : (h0 h1 h2 ...)
			temp_row0 = temp_buffer[nrh - 2]; //  E    <= (e0 e1 e2 ...)
			out_row0 = Y[nrh - 1];
			out_row1 = Y[nrh + 0];
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
			row2 = X[nrh + 0];   //  0 : (f0 f1 f2 ...)
			row3 = X[nrh + 1];   //  1 : (g0 g1 g2 ...)
			temp_row0 = temp_buffer[nrh - 1]; //  E    <= (e0 e1 e2 ...)
			out_row1 = Y[nrh + 0];
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
	// free_ui8matrix(temp, nrl + (-1), nrh + 1, ncl + (-1), nch + 1);
}
void ui8matrix_dilation_pipeline2_LU3x3_InLU_O3_RR (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **temp_buffer, uint8 **Y)
{
	const long order = 3;
	long row = nrl, col = ncl, x, y, r = 0;
	// uint8 **temp = ui8matrix(nrl + (-1), nrh + 1, ncl + (-1), nch + 1);
	uint8 *row0, *row1, *row2, *row3, x0, x1, x2, x3, x4, x5, y0, y1, y2;
	uint8 *temp_row0, *temp_row1, *temp_row2, *temp_row3, *temp_row4;
	uint8 *out_row0, *out_row1, *out_row2, *out_row3;
	
	r = (nch + 1) % order;
	// Prologue	
	row0 = X[row - 1];
	row1 = X[row + 0];
	row2 = X[row + 1];
	
	// temp_row1 = temp_buffer[row - 1];
	// temp_row2 = temp_buffer[row + 0];
	// temp_row3 = temp_buffer[row + 1];
	for (row = nrl; row < nrh + 1; row++) {
		row0 = X[row - 1];
		row1 = X[row + 0];
		row2 = X[row + 1];
		temp_row1 = temp_buffer[row + 0];

		temp_row1[ncl - 1] = row0[ncl - 1] | row1[ncl - 1] | row2[ncl - 1];
		temp_row1[ncl + 0] = row0[ncl + 0] | row1[ncl + 0] | row2[ncl + 0];
	}

	for (row = nrl; row < nrh + 1; row ++){
		temp_row1 = temp_buffer[row + 0];
		
		row0 = X[row - 1];
		row1 = X[row + 0];
		row2 = X[row + 1];

		out_row1 = Y[row]; 
		
		x0 = temp_row1[ncl - 1];		// LOAD => A
		x1 = temp_row1[ncl - 0];		// LOAD => B
		for (col = ncl; col < nch + 1 - r; col += order) {
			x2 = row0[col + 1] | row1[col + 1] | row2[col + 1];
			x3 = row0[col + 2] | row1[col + 2] | row2[col + 2];
			x4 = row0[col + 3] | row1[col + 3] | row2[col + 3];
			
			out_row1 [col + 0] = x0 | x1 | x2; // A | B | C
			out_row1 [col + 1] = x1 | x2 | x3; // B | C | D
			out_row1 [col + 2] = x2 | x3 | x4; // C | D | E

			x0 = x3; x1 = x4;
		}
	}
	switch (r) {
		case 2: 
			row0 = X[nrl - 1];
			row1 = X[nrl + 0];
			for (row = nrl; row < nrh + 1; row++) {
				row2 = X[row + 1];
				out_row1 = Y[row]; 

				x1 = row0[nch - 1] | row1[nch - 1] | row2[nch - 1]|
					 row0[nch + 0] | row1[nch + 0] | row2[nch + 0];

				out_row1 [nch - 1] = row0[nch - 2] | row1[nch - 2] | row2[nch - 2] | x1 ; // A | B | C
				out_row1 [nch + 0] = x1 | row0[nch + 1] | row1[nch + 1] | row2[nch + 1]; // B | C | D
				row0 = row1;
				row1 = row2;
			}
			break;
		case 1: 
			row0 = X[nrl - 1];
			row1 = X[nrl + 0];
			for (row = nrl; row < nrh + 1; row++) {
				row2 = X[row + 1];
				Y[row][nch] = scalar_or3(row0, nch) |
									 scalar_or3(row1, nch) |
									 scalar_or3(row2, nch);
				row0 = row1;
				row1 = row2;
			}
			break;
		default:
			break;
	}
	// display_ui8matrix(Y, nrl - 1, nrh + 1, ncl - 1, nch + 1, "%u", "TEst");
	// free_ui8matrix(temp, nrl + (-1), nrh + 1, ncl + (-1), nch + 1);
	
}
/***********************************************************************************/
/******* End of optimisation : Loop Unroll + Register Rotation + Pipelining ********/
/***********************************************************************************/

void ui8matrix_dilation_LU3x3_ExLU_O3_NS (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **temp_buffer, uint8 **Y)
{
	const long order = 3;
	long row = nrl, col = ncl, x, y;
	uint8 y0, y1, y2, x0, x1, x2, x3, r;
	uint8 *row0, *row1, *row2, *row3, *row4;

	r = (nrh + 1)  % order;
	
	for (row = nrl; row < nrh + 1 - r; row += order) {
		for (col = ncl; col < nch + 1; col++) {
			x0 = scalar_or3(X[row + 0], col);
			x1 = scalar_or3(X[row + 1], col);
			x2 = scalar_or3(X[row + 2], col);

			Y[row + 0][col] = scalar_or3(X[row - 1], col) | x0 | x1;
			Y[row + 1][col] = 				     			x0 | x1 | x2;
			Y[row + 2][col] = 				     			x1 | x2 | scalar_or3(X[row + 3], col);
		}
	}
	
	switch(r) {
		case 2:
			for (col = ncl; col < nch + 1; col++) {
				x0 = scalar_or3(X[row + 0], col);
				x1 = scalar_or3(X[row + 1], col);
				Y[row + 0][col] = scalar_or3(X[row - 1], col) | x0 | x1;
				Y[row + 1][col] = 				     		    x0 | x1 | scalar_or3(X[row + 2], col);
			}
			break;
		case 1:
			for (col = ncl; col < nch + 1; col++) {
				Y[row + 0][col] = scalar_or3(X[row - 1], col) | 
										 scalar_or3(X[row + 0], col) |
										 scalar_or3(X[row + 1], col);
				
			}
			break;
		default:
			break;
	}
}
void ui8matrix_dilation_LU3x3_InLU_O3_NS (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **temp_buffer, uint8 **Y)
{
	const long order = 3;
	long row = nrl, col = ncl, x, y, r;
	uint8 y0, y1, y2, x0, x1, x2, x3, x4, x5, x6, x7, x8;
	uint8 *row0, *row1, *row2, *row3, *row4;

	r = (nch + 1)  % order;
	
	for (row = nrl; row < nrh + 1; row ++) {
		row0 = X[row - 1];
		row1 = X[row + 0];
		row2 = X[row + 1];
		for (col = ncl; col < nch + 1 - r; col += order) {
			// x0 = scalar_or3(row0, col); x3 = scalar_or3(row0, col + 1); x6 = scalar_or3(row0, col + 2);
			// x1 = scalar_or3(row1, col); x4 = scalar_or3(row1, col + 1); x7 = scalar_or3(row1, col + 2);
			// x2 = scalar_or3(row2, col); x5 = scalar_or3(row2, col + 1); x8 = scalar_or3(row2, col + 2);

			Y[row + 0][col + 0] = scalar_or3(row0, col + 0)|
										 scalar_or3(row1, col + 0)|
										 scalar_or3(row2, col + 0);
			Y[row + 0][col + 1] = scalar_or3(row0, col + 1)|
 										 scalar_or3(row1, col + 1)|
 										 scalar_or3(row2, col + 1);
			Y[row + 0][col + 2] = scalar_or3(row0, col + 2)|
										 scalar_or3(row1, col + 2)|
										 scalar_or3(row2, col + 2);
		}
	}

	switch(r) {
		case 2:
			for (row = nrl; row < nrh + 1; row ++) {
				row0 = X[row - 1];
				row1 = X[row + 0];
				row2 = X[row + 1];

				Y[row + 0][col + 0] = scalar_or3(row0, col + 0)|
 											 scalar_or3(row1, col + 0)|
 											 scalar_or3(row2, col + 0);
				Y[row + 0][col + 1] = scalar_or3(row0, col + 1)|
											 scalar_or3(row1, col + 1)|
											 scalar_or3(row2, col + 1);
			}
			break;
		case 1:
			for (row = nrl; row < nrh + 1; row ++) {
				Y[row + 0][col + 0] = scalar_or3(X[row - 1], col) | 
											 scalar_or3(X[row + 0], col) | 
											 scalar_or3(X[row + 1], col);
			}
			break;
		default:
			break;
	}
}
void ui8matrix_dilation_LU3x3_ComLU_O3_NS (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **temp_buffer, uint8 **Y)
{
	const long order = 3;
	long row = nrl, col = ncl, x, y;
	uint8 y0, y1, y2, rr, cr;
	uint8 x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14;
	uint8 *row0, *row1, *row2, *row3, *row4;

	rr = (nrh + 1) % order;
	cr = (nch + 1) % order;
	
	for (row = nrl; row < nrh + 1 - rr; row += order) {
		for (col = ncl; col < nch + 1 - cr; col += order) {
			x0 = scalar_or3(X[row + 0], col); x1 = scalar_or3(X[row + 0], col + 1); x2 = scalar_or3(X[row + 0], col + 2);
			x3 = scalar_or3(X[row + 1], col); x4 = scalar_or3(X[row + 1], col + 1); x5 = scalar_or3(X[row + 1], col + 2);
			x6 = scalar_or3(X[row + 2], col); x7 = scalar_or3(X[row + 2], col + 1); x8 = scalar_or3(X[row + 2], col + 2);

			Y[row + 0][col + 0] = scalar_or3(X[row - 1], col + 0) | x0 | x3;
			Y[row + 1][col + 0] = 				         			x0 | x3 | x6;
			Y[row + 2][col + 0] = 				         			x3 | x6 | scalar_or3(X[row + 3], col + 0);


			Y[row + 0][col + 1] = scalar_or3(X[row - 1], col + 1) | x1 | x4;
			Y[row + 1][col + 1] =			  			            x1 | x4 | x7;
			Y[row + 2][col + 1] =			  			            x4 | x7 | scalar_or3(X[row + 3], col + 1);


			Y[row + 0][col + 2] = scalar_or3(X[row - 1], col + 2) | x2 | x5;
			Y[row + 1][col + 2] = 				                    x2 | x5 | x8;
			Y[row + 2][col + 2] = 				                    x5 | x8 | scalar_or3(X[row + 3], col + 2);

		}
	}
	switch (rr) {
		case 2 :
			switch(cr){
				case 2 :
					ui8matrix_dilation_LU3x3_ExLU_O3_NS(X, nrl, nrh, nch - 1, nch, temp_buffer, Y);
					ui8matrix_dilation_LU3x3_InLU_O3_NS(X, nrh - 1, nrh, ncl, nch, temp_buffer, Y);
				break;
				case 1: 
					ui8matrix_dilation_LU3x3_ExLU_O3_NS(X, nrl, nrh, nch, nch, temp_buffer, Y);		
					ui8matrix_dilation_LU3x3_InLU_O3_NS(X, nrh - 1, nrh, ncl, nch, temp_buffer, Y);	
				break;
				default :
					ui8matrix_dilation_LU3x3_InLU_O3_NS(X, nrh - 1, nrh, ncl, nch, temp_buffer, Y);
				break;
			}
		break;
		case 1:
			switch(cr){
				case 2 :
					ui8matrix_dilation_LU3x3_ExLU_O3_NS(X, nrl, nrh, nch - 1, nch, temp_buffer, Y);
					ui8matrix_dilation_LU3x3_InLU_O3_NS(X, nrh, nrh, ncl, nch, temp_buffer, Y);
				break;
				case 1: 
					ui8matrix_dilation_LU3x3_ExLU_O3_NS(X, nrl, nrh, nch, nch, temp_buffer, Y);		
					ui8matrix_dilation_LU3x3_InLU_O3_NS(X, nrh, nrh, ncl, nch, temp_buffer, Y);	
				break;
				default :
					ui8matrix_dilation_LU3x3_InLU_O3_NS(X, nrh, nrh, ncl, nch, temp_buffer, Y);
				break;
			}
		break;
		default:
			switch(cr){
				case 2 :
					ui8matrix_dilation_LU3x3_ExLU_O3_NS(X, nrl, nrh, nch - 1, nch, temp_buffer, Y);
				break;
				case 1: 
					ui8matrix_dilation_LU3x3_ExLU_O3_NS(X, nrl, nrh, nch, nch, temp_buffer, Y);	
				break;
				default :
				break;
			}
		break;
	}
	// ui8matrix_dilation_LU3x3(X, row, nrh, col, nrh, Y);
}
void ui8matrix_dilation_LU3x3_ExLU_O3_RR_NS(uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **temp_buffer, uint8 **Y)
{
	const long order = 3;
	long row = nrl, col = ncl, x, y;
	uint8 y0, y1, y2, r;
	uint8 x0, x1, x2, x3, x4, x5, x6, x7, x8;
	uint8 *row0, *row1, *row2, *row3, *row4;

	r = (nrh + 1)  % order;
	
	row0 = X[row - 1];
	row1 = X[row + 0];
	for (row = nrl; row < nrh + 1 - r; row += order) {
		row2 = X[row + 1]; 
		row3 = X[row + 2]; 
		row4 = X[row + 3];
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

			Y[row + 0][col] = scalar_or3(row0, col) | y0 | y1; 						// 3 + 2    =  5
			Y[row + 1][col] =         			 y0 | y1 | y2;						//     2    =  2
			Y[row + 2][col] = 		 			 y1 | y2 | scalar_or3(row4, col);  //     2 + 3=  5
			// 6 + 12 = 18 calc
			x0 = x3; x3 = x6;
			x1 = x4; x4 = x7;
			x2 = x5; x5 = x8;
		}
		row0 = row3;
		row1 = row4;
	}
	row2 = X[row + 1];
	row3 = X[row + 2];
	row4 = X[row + 3];
	switch(r) {
		case 2:
			x0 = row1[-1]; x3 = row1[ 0];
			x1 = row2[-1]; x4 = row2[ 0];
			x2 = row3[-1]; x5 = row3[ 0];
			for (col = ncl; col < nch + 1; col++) {
				x6 = row1[col + 1]; //scalar_or3(row1, col);
				x7 = row2[col + 1]; //scalar_or3(row2, col);
				x8 = row3[col + 1]; //scalar_or3(row3, col);
				Y[row + 0][col] = scalar_or3(row0, col) | (x0 | x3 | x6) | (x1 | x4 | x7);;
				Y[row + 1][col] =         (x0 | x3 | x6) | (x1 | x4 | x7) | (x2 | x5 | x8);;
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
				Y[row + 0][col] = scalar_or3(row0, col) | (x0 | x3 | x6) | (x1 | x4 | x7);
				x0 = x3; x3 = x6;
				x1 = x4; x4 = x7;
			}				
			break;
		default:
			break;
	}
}
void ui8matrix_dilation_LU3x3_InLU_O3_RR_NS (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **temp_buffer, uint8 **Y)
{
	const long order = 3;
	long row = nrl, col = ncl, x, y, r;
	uint8 y0, y1, y2, y3, y4, y5;
	uint8 x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16, x17, x18, x19;
	uint8 *row0, *row1, *row2, *row3, *row4;

	r = (nch + 1) % order;
	row0 = X[row - 1];
	row1 = X[row + 0];
	for (row = nrl; row < nrh + 1; row++) {
		row2 = X[row + 1];

		x0  = row0[ncl - 1]; x1  = row0[ncl + 0]; 
		x5  = row1[ncl - 1]; x6  = row1[ncl + 0]; 
		x10 = row2[ncl - 1]; x11 = row2[ncl + 0]; 
		for (col = ncl; col < nch + 1 - r; col += order) {
			x2  = row0[col + 1]; x3  = row0[col + 2]; x4  = row0[col + 3];
			x7  = row1[col + 1]; x8  = row1[col + 2]; x9  = row1[col + 3];
			x12 = row2[col + 1]; x13 = row2[col + 2]; x14 = row2[col + 3];
			
			y1 = x1 | x2 | x6 | x7 | x11 | x12; // 5 calcs
			y2 = x3 | x8 | x13;					// 2 calcs
			Y[row][col + 0] = x0 | x5 | y1 | x10; //3 calcs
			Y[row][col + 1] =      y2 | y1; 		 //1 calc
			Y[row][col + 2] =      y2 | x2 | x4          | x7 | x9             | x12 | x14; //6 calcs
			// 5 + 2 + 3 + 1 + 6 = 17						

			x0  = x3;  x1  = x4;
			x5  = x8;  x6  = x9;
			x10 = x13; x11 = x14;
			
		}
		row0 = row1;
		row1 = row2;
	}

		
	row = nrl;
	row0 = X[nrl - 1];
	row1 = X[nrl + 0];
	
	switch(r) {
	case 2:
		for (row; row < nrh + 1; row++) {
			row2 = X[row + 1];
			x0  = row0[col - 1]; x1  = row0[col +  0]; x2  = row0[col + 1]; x3  = row0[col + 2];
			x5  = row1[col - 1]; x6  = row1[col +  0]; x7  = row1[col + 1]; x8  = row1[col + 2];
			x10 = row2[col - 1]; x11 = row2[col +  0]; x12 = row2[col + 1]; x13 = row2[col + 2];
			y1 = x1 | x2 | x6 | x7 | x11 | x12; // 5 calcs
			y2 = x3 | x8 | x13;					// 2 calcs
			Y[row][col + 0] = x0 | x5 | y1 | x10; //3 calcs
			Y[row][col + 1] =      y2 | y1; 		 //1 calc
			row0 = row1;
			row1 = row2;
		}
		break;
	case 1:
		for (row; row < nrh + 1; row++) {			
			row2 = X[row + 1];
			x0  = row0[col - 1]; x1  = row0[col +  0]; x2  = row0[col + 1];
			x5  = row1[col - 1]; x6  = row1[col +  0]; x7  = row1[col + 1];
			x10 = row2[col - 1]; x11 = row2[col +  0]; x12 = row2[col + 1];
			Y[row][col + 0] = (x0 | x1 | x2 ) | (x5 | x6 | x7) | (x10 | x11 | x12);
			row0 = row1;
			row1 = row2;
		}
		break;
	default:
		break;
	}
}


/*******************************************/
/******* Optimisation : Loop Unroll ********/
/*******************************************/

void ui8matrix_erosion_LU3x3_O1xO1(uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **temp_buffer, uint8 **Y)
{
	long row = nrl, col = ncl, x, y;
	// dilate
	for (row = nrl; row < nrh + 1; row++)
		for (col = ncl; col < nch + 1; col++)
            Y[row][col] = scalar_and3x3(&X[row], col);
}

void ui8matrix_erosion_LU3x3_InLU_O3 (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **temp_buffer, uint8 **Y)
{
	const long order = 3;
	long row = nrl, col = ncl, x, y, r;
	uint8 y0, y1, y2, x0, x1, x2, x3, x4, x5, x6, x7, x8;
	uint8 *row0, *row1, *row2, *row3, *row4;
	uint8 *out_row0;
	r = (nch + 1)  % order;
	
	for (row = nrl; row < nrh + 1; row ++) {
		row0 = X[row - 1];
		row1 = X[row + 0];
		row2 = X[row + 1];
		out_row0 = Y[row + 0];
		for (col = ncl; col < nch + 1 - r; col += order) {
			// x0 = scalar_and3(row0, col); x3 = scalar_and3(row0, col + 1); x6 = scalar_and3(row0, col + 2);
			// x1 = scalar_and3(row1, col); x4 = scalar_and3(row1, col + 1); x7 = scalar_and3(row1, col + 2);
			// x2 = scalar_and3(row2, col); x5 = scalar_and3(row2, col + 1); x8 = scalar_and3(row2, col + 2);

			out_row0[col + 0] = scalar_and3(row0, col + 0)&
								scalar_and3(row1, col + 0)&
								scalar_and3(row2, col + 0);
			out_row0[col + 1] = scalar_and3(row0, col + 1)&
 								scalar_and3(row1, col + 1)&
 								scalar_and3(row2, col + 1);
			out_row0[col + 2] = scalar_and3(row0, col + 2)&
								scalar_and3(row1, col + 2)&
								scalar_and3(row2, col + 2);
		}
	}

	switch(r) {
		case 2:
			for (row = nrl; row < nrh + 1; row ++) {
				row0 = X[row - 1];
				row1 = X[row + 0];
				row2 = X[row + 1];
				out_row0 = Y[row + 0];

				out_row0[nch - 1] = scalar_and3(row0, nch - 1)&
 									scalar_and3(row1, nch - 1)&
 									scalar_and3(row2, nch - 1);
				out_row0[nch + 0] = scalar_and3(row0, nch + 0)&
									scalar_and3(row1, nch + 0)&
									scalar_and3(row2, nch + 0);
			}
			break;
		case 1:
			for (row = nrl; row < nrh + 1; row ++) {
				Y[row + 0][nch + 0] = scalar_and3(X[row - 1], nch) & 
										  	 scalar_and3(X[row + 0], nch) & 
										  	 scalar_and3(X[row + 1], nch);
			}
			// ui8matrix_erosion_LU3x3_O1xO1(X, nrl, nrh, nch, nch, temp_buffer, Y);
			break;
		default:
			break;
	}
}

void ui8matrix_erosion_LU3x3_ExLU_O3 (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **temp_buffer, uint8 **Y)
{
	const long order = 3;
	long row = nrl, col = ncl, x, y;
	uint8 y0, y1, y2, x0, x1, x2, x3, r;
	uint8 *row0, *row1, *row2, *row3, *row4;
	uint8 *out_row0, *out_row1, *out_row2;

	r = (nrh + 1)  % order;
	
	for (row = nrl; row < nrh + 1 - r; row += order) {
		row0 = X[row - 1];
		row1 = X[row + 0];
		row2 = X[row + 1];
		row3 = X[row + 2];
		row4 = X[row + 3];

		out_row0 = Y[row + 0];
		out_row1 = Y[row + 1];
		out_row2 = Y[row + 2];
		
		for (col = ncl; col < nch + 1; col++) {
			x0 = scalar_and3(row1, col);
			x1 = scalar_and3(row2, col);
			x2 = scalar_and3(row3, col);

			out_row0[col] = scalar_and3(row0, col) & x0 & x1;
			out_row1[col] = 			       x0 & x1 & x2;
			out_row2[col] = 			       x1 & x2 & scalar_and3(row4, col);
		}
	}
	
	switch(r) {
		case 2:
			row0 = X[row - 1];
			row1 = X[row + 0];
			row2 = X[row + 1];
			row3 = X[row + 2];
			out_row0 = Y[row + 0];
			out_row1 = Y[row + 1];
			for (col = ncl; col < nch + 1; col++) {
				x0 = scalar_and3(row1, col);
				x1 = scalar_and3(row2, col);
				out_row0[col] = scalar_and3(row0, col) & x0 & x1;
				out_row1[col] = 				     x0 & x1 & scalar_and3(row3, col);
			}
			break;
		case 1:
			row0 = X[row - 1];
			row1 = X[row + 0];
			row2 = X[row + 1];
			out_row0 = Y[row + 0];
			for (col = ncl; col < nch + 1; col++) {
				out_row0[col] = scalar_and3(row0, col) & 
										 scalar_and3(row1, col) &
										 scalar_and3(row2, col);
				
			}
			break;
		default:
			break;
	}
}

void ui8matrix_erosion_LU3x3_ComLU_O3 (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **temp_buffer, uint8 **Y)
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
		row0 = X[row - 1];
		row1 = X[row + 0];
		row2 = X[row + 1];
		row3 = X[row + 2];
		row4 = X[row + 3];
		out_row0 = Y[row + 0];
		out_row1 = Y[row + 1];
		out_row2 = Y[row + 2];
		for (col = ncl; col < nch + 1 - cr; col += order) {
			x0 = scalar_and3(row1, col); x1 = scalar_and3(row1, col + 1); x2 = scalar_and3(row1, col + 2);
			x3 = scalar_and3(row2, col); x4 = scalar_and3(row2, col + 1); x5 = scalar_and3(row2, col + 2);
			x6 = scalar_and3(row3, col); x7 = scalar_and3(row3, col + 1); x8 = scalar_and3(row3, col + 2);

			out_row0[col + 0] = scalar_and3(row0, col + 0) & x0 & x3;
			out_row1[col + 0] = 				         x0 & x3 & x6;
			out_row2[col + 0] = 				         x3 & x6 & scalar_and3(row4, col + 0);


			out_row0[col + 1] = scalar_and3(row0, col + 1) & x1 & x4;
			out_row1[col + 1] = 				         x1 & x4 & x7;
			out_row2[col + 1] = 				         x4 & x7 & scalar_and3(row4, col + 1);


			out_row0[col + 2] = scalar_and3(row0, col + 2) & x2 & x5;
			out_row1[col + 2] = 				         x2 & x5 & x8;
			out_row2[col + 2] = 				         x5 & x8 & scalar_and3(row4, col + 2);

		}
	}
	switch (rr) {
		case 2 :
			switch(cr){
				case 2 :
					ui8matrix_erosion_LU3x3_ExLU_O3(X, nrl, nrh, nch - 1, nch, temp_buffer, Y);
					ui8matrix_erosion_LU3x3_InLU_O3(X, nrh - 1, nrh, ncl, nch, temp_buffer, Y);
				break;
				case 1: 
					ui8matrix_erosion_LU3x3_ExLU_O3(X, nrl, nrh, nch, nch, temp_buffer, Y);		
					ui8matrix_erosion_LU3x3_InLU_O3(X, nrh - 1, nrh, ncl, nch, temp_buffer, Y);	
				break;
				default :
					ui8matrix_erosion_LU3x3_InLU_O3(X, nrh - 1, nrh, ncl, nch, temp_buffer, Y);
				break;
			}
			break;
		case 1:
			switch(cr){
				case 2 :
					ui8matrix_erosion_LU3x3_ExLU_O3(X, nrl, nrh, nch - 1, nch, temp_buffer, Y);
					ui8matrix_erosion_LU3x3_InLU_O3(X, nrh, nrh, ncl, nch, temp_buffer, Y);
				break;
				case 1: 
					ui8matrix_erosion_LU3x3_ExLU_O3(X, nrl, nrh, nch, nch, temp_buffer, Y);		
					ui8matrix_erosion_LU3x3_InLU_O3(X, nrh, nrh, ncl, nch, temp_buffer, Y);	
				break;
				default :
					ui8matrix_erosion_LU3x3_InLU_O3(X, nrh, nrh, ncl, nch, temp_buffer, Y);
				break;
			}
		break;
		default:
			switch(cr){
				case 2 :
					ui8matrix_erosion_LU3x3_ExLU_O3(X, nrl, nrh, nch - 1, nch, temp_buffer, Y);
				break;
				case 1: 
					ui8matrix_erosion_LU3x3_ExLU_O3(X, nrl, nrh, nch, nch, temp_buffer, Y);	
				break;
				default :
				break;
			}
		break;
	}
}
/**************************************************/
/******* End of optimisation : Loop Unroll ********/
/**************************************************/

/**************************************************************************************/
/******* Optimisation : Loop Unroll +  Register Rotation of Values / Addresses ********/
/**************************************************************************************/
void ui8matrix_erosion_LU3x3_InLU_O3_ValAddrRR (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **temp_buffer, uint8 **Y)
{
	const long order = 3;
	long row = nrl, col = ncl, r;
	uint8 y0, y1, y2;
	uint8 x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14;
	uint8 *row0, *row1, *row2, *row3, *row4;
	uint8 *out_row0;

	r = (nch + 1) % order;
	row0 = X[nrl - 1];
	row1 = X[nrl + 0];
	for (row = nrl; row < nrh + 1; row++) {
		row2 = X[row + 1];

		out_row0 = Y[row + 0];
		x0  = row0[ncl - 1]; x1  = row0[ncl + 0];
		x5  = row1[ncl - 1]; x6  = row1[ncl + 0];
		x10 = row2[ncl - 1]; x11 = row2[ncl + 0];
		for (col = ncl; col < nch + 1 - r; col += order) {
			x2  = row0[col + 1]; x3  = row0[col + 2]; x4  = row0[col + 3];
			x7  = row1[col + 1]; x8  = row1[col + 2]; x9  = row1[col + 3];
			x12 = row2[col + 1]; x13 = row2[col + 2]; x14 = row2[col + 3];
			
			out_row0[col + 0] = (x0 & x1 & x2) & (x5 & x6 & x7) & (x10 & x11 & x12);
			out_row0[col + 1] = (x1 & x2 & x3) & (x6 & x7 & x8) & (x11 & x12 & x13);
			out_row0[col + 2] = (x2 & x3 & x4) & (x7 & x8 & x9) & (x12 & x13 & x14); 

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
			row0 = X[row - 1];
			row1 = X[row + 0];
			row2 = X[row + 1];

			out_row0 = Y[row + 0];
		
			// y1 = x1 & x2 & x6 & x7 & x11 & x12; // 5
			x0  = row0[nch - 2]; x1  = row0[nch - 1]; x2  = row0[nch + 0]; x3  = row0[nch + 1]; 
			x5  = row1[nch - 2]; x6  = row1[nch - 1]; x7  = row1[nch + 0]; x8  = row1[nch + 1]; 
			x10 = row2[nch - 2]; x11 = row2[nch - 1]; x12 = row2[nch + 0]; x13 = row2[nch + 1]; 

			out_row0[nch - 1] = (x0  & x1  & x2) & 
								(x5  & x6  & x7) & 
								(x10 & x11 & x12);

			out_row0[nch + 0] = (x1  & x2  & x3) & 
			  					(x6  & x7  & x8) & 
								(x11 & x12 & x13);
			
		}
		break;
	case 1:
		for (row = nrl; row < nrh + 1; row++) {
			row0 = X[row - 1];
			row1 = X[row + 0];
			row2 = X[row + 1];

			out_row0 = Y[row + 0];
			out_row0[nch + 0] = (row0[nch - 1] &  row0[nch + 0] & row0[nch + 1])& 
			  					(row1[nch - 1] &  row1[nch + 0] & row1[nch + 1])& 
								(row2[nch - 1] &  row2[nch + 0] & row2[nch + 1]);
			
		}
		break;
	default:
		break;
	}
}
void ui8matrix_erosion_LU3x3_ExLU_O3_ValAddrRR(uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **temp_buffer, uint8 **Y)
{
	
	const long order = 3;
	long row = nrl, col = ncl, x, y;
	uint8 y0, y1, y2, y3, y5, r;
	uint8 x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14;
	uint8 *row0, *row1, *row2, *row3, *row4;
	uint8 *out_row0, *out_row1, *out_row2;

	r = (nrh + 1)  % order;
	
	row0 = X[row - 1];
	row1 = X[row + 0];
	for (row = nrl; row < nrh + 1 - r; row += order) {
		row2 = X[row + 1]; 
		row3 = X[row + 2]; 
		row4 = X[row + 3];
		x0 = row0[ncl - 1]; x5 = row0[ncl + 0];
		x1 = row1[ncl - 1]; x6 = row1[ncl + 0];
		x2 = row2[ncl - 1]; x7 = row2[ncl + 0];
		x3 = row3[ncl - 1]; x8 = row3[ncl + 0];
		x4 = row4[ncl - 1]; x9 = row4[ncl + 0];
		out_row0 = Y[row + 0];
		out_row1 = Y[row + 1];
		out_row2 = Y[row + 2];
		for (col = ncl; col < nch + 1; col++) {
			x10 = row0[col + 1]; //scalar_and3(row1, col);
			x11 = row1[col + 1]; //scalar_and3(row2, col);
			x12 = row2[col + 1]; //scalar_and3(row3, col);
			x13 = row3[col + 1]; //scalar_and3(row2, col);
			x14 = row4[col + 1]; //scalar_and3(row3, col);
			
			y1 = (x1 & x6 & x11); // 2 calc
			y2 = (x2 & x7 & x12); // 2 calc
			y3 = (x3 & x8 & x13); // 2 calc
			out_row0[col] = (x0 & x5 & x10) & y1 & y2; 				
			out_row1[col] = 				  y1 & y2 & y3;				
			out_row2[col] = 		          y2 & y3 & (x4 & x9 & x14); 
			// 6 + 12 = 18 calc
			x0 = x5; x5 = x10;
			x1 = x6; x6 = x11;
			x2 = x7; x7 = x12;
			x3 = x8; x8 = x13;
			x4 = x9; x9 = x14;
		}
		row0 = row3;
		row1 = row4;
	}
	row2 = X[row + 1];
	row3 = X[row + 2];
	row4 = X[row + 3];
	out_row0 = Y[nrh - 1];
	out_row1 = Y[nrh + 0];
	switch(r) {
		case 2:
			x0 = row1[ncl - 1]; x3 = row1[ncl + 0];
			x1 = row2[ncl - 1]; x4 = row2[ncl + 0];
			x2 = row3[ncl - 1]; x5 = row3[ncl + 0];
			for (col = ncl; col < nch + 1; col++) {
				x6 = row1[col + 1]; //scalar_and3(row1, col);
				x7 = row2[col + 1]; //scalar_and3(row2, col);
				x8 = row3[col + 1]; //scalar_and3(row3, col);
				out_row0[col] = scalar_and3(row0, col) & (x0 & x3 & x6) & (x1 & x4 & x7);;
				out_row1[col] =         (x0 & x3 & x6) & (x1 & x4 & x7) & (x2 & x5 & x8);;
				x0 = x3; x3 = x6;
				x1 = x4; x4 = x7;
				x2 = x5; x5 = x8;
			}
			break;
		case 1:
			x0 = row1[ncl - 1]; x3 = row1[ncl + 0];
			x1 = row2[ncl - 1]; x4 = row2[ncl + 0];
			for (col = ncl; col < nch + 1; col++) {
				x6 = row1[col + 1]; //scalar_and3(row1, col);
				x7 = row2[col + 1]; //scalar_and3(row2, col);
				out_row1[col] = scalar_and3(row0, col) & (x0 & x3 & x6) & (x1 & x4 & x7);
				x0 = x3; x3 = x6;
				x1 = x4; x4 = x7;
			}				
			break;
		default:
			break;
	}
}
void ui8matrix_erosion_LU3x3_ComLU_O3_ValAddrRR(uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **temp_buffer, uint8 **Y)
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
	
	row0 = X[row - 1];
	row1 = X[row + 0];

	for (row = nrl; row < nrh + 1 - rr; row += order) {
		row2 = X[row + 1];
		row3 = X[row + 2];
		row4 = X[row + 3];
		out_row0 = Y[row + 0];
		out_row1 = Y[row + 1];
		out_row2 = Y[row + 2];


		x0  = row0[ncl - 1]; x1  = row0[ncl + 0]; x2  = row0[ncl + 1]; x3  = row0[ncl + 2]; x4  = row0[ncl + 3];
		x5  = row1[ncl - 1]; x6  = row1[ncl + 0]; x7  = row1[ncl + 1]; x8  = row1[ncl + 2]; x9  = row1[ncl + 3];
		x10 = row2[ncl - 1]; x11 = row2[ncl + 0]; x12 = row2[ncl + 1]; x13 = row2[ncl + 2]; x14 = row2[ncl + 3];
		x15 = row3[ncl - 1]; x16 = row3[ncl + 0]; x17 = row3[ncl + 1]; x18 = row3[ncl + 2]; x19 = row3[ncl + 3];
		x20 = row4[ncl - 1]; x21 = row4[ncl + 0]; x22 = row4[ncl + 1]; x23 = row4[ncl + 2]; x24 = row4[ncl + 3];
		y0  = x0  & x5  & x10; y1  = x1  & x6  & x11;
		y5  = x5  & x10 & x15; y6  = x6  & x11 & x16;
		y10 = x10 & x15 & x20; y11 = x11 & x16 & x21;

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


			y2  = x2  & x7  & x12; y3  = x3  & x8  & x13; y4  = x4  & x9  & x14; 
			y7  = x7  & x12 & x17; y8  = x8  & x13 & x18; y9  = x9  & x14 & x19; 
			y12 = x12 & x17 & x22; y13 = x13 & x18 & x23; y14 = x14 & x19 & x24; 
 
			out_row0[col + 0] = y0 & y1 & y2; 
			out_row0[col + 1] = y1 & y2 & y3; 
			out_row0[col + 2] = y2 & y3 & y4;
			y0 = y3; y1 = y4;

			out_row1[col + 0] = y5 & y6 & y7;
			out_row1[col + 1] = y6 & y7 & y8;
			out_row1[col + 2] = y7 & y8 & y9;
			y5 = y8; y6 = y9;
			
			out_row2[col + 0] = y10 & y11 & y12;
			out_row2[col + 1] = y11 & y12 & y13;
			out_row2[col + 2] = y12 & y13 & y14;
			y10 = y13; y11 = y14;
		}
		row0 = row3;
		row1 = row4;
	}
	switch (rr) {
		case 2 :
			switch(cr){
				case 2 :
					ui8matrix_erosion_LU3x3_ExLU_O3_ValAddrRR(X, nrl, nrh, nch - 1, nch, temp_buffer, Y);
					ui8matrix_erosion_LU3x3_InLU_O3_ValAddrRR(X, nrh - 1, nrh, ncl, nch, temp_buffer, Y);
				break;
				case 1: 
					ui8matrix_erosion_LU3x3_ExLU_O3_ValAddrRR(X, nrl, nrh, nch, nch, temp_buffer, Y);		
					ui8matrix_erosion_LU3x3_InLU_O3_ValAddrRR(X, nrh - 1, nrh, ncl, nch, temp_buffer, Y);	
				break;
				default :
					ui8matrix_erosion_LU3x3_InLU_O3_ValAddrRR(X, nrh - 1, nrh, ncl, nch, temp_buffer, Y);
				break;
			}
			break;
		case 1:
			switch(cr){
				case 2 :
					ui8matrix_erosion_LU3x3_ExLU_O3_ValAddrRR(X, nrl, nrh, nch - 1, nch, temp_buffer, Y);
					ui8matrix_erosion_LU3x3_InLU_O3_ValAddrRR(X, nrh, nrh, ncl, nch, temp_buffer, Y);
				break;
				case 1: 
					ui8matrix_erosion_LU3x3_ExLU_O3_ValAddrRR(X, nrl, nrh, nch, nch, temp_buffer, Y);		
					ui8matrix_erosion_LU3x3_InLU_O3_ValAddrRR(X, nrh, nrh, ncl, nch, temp_buffer, Y);	
					// display_ui8matrix(Y, nrl, nrh, ncl, nch, "%u", "rr=1; cr=1");
				break;
				default :
					ui8matrix_erosion_LU3x3_InLU_O3_ValAddrRR(X, nrh, nrh, ncl, nch, temp_buffer, Y);
				break;
			}
		break;
		default:
			switch(cr){
				case 2 :
					ui8matrix_erosion_LU3x3_ExLU_O3_ValAddrRR(X, nrl, nrh, nch - 1, nch, temp_buffer, Y);
				break;
				case 1: 
					ui8matrix_erosion_LU3x3_ExLU_O3_ValAddrRR(X, nrl, nrh, nch, nch, temp_buffer, Y);	
				break;
				default :
				break;
			}
		break;
	}
		
}
/*********************************************************************************************/
/******* End of optimisation : Loop Unroll +  Register Rotation of Values / Addresses ********/
/*********************************************************************************************/

/*****************************************************************************/
/******** Optimisation : Loop Unroll + Register Rotation of Addresses ********/
/*****************************************************************************/

void ui8matrix_erosion_LU3x3_InLU_O3_AddrRR(uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **temp_buffer, uint8 **Y)
{
	const long order = 3;
	long row = nrl, col = ncl, r;
	uint8 y0, y1, y2;
	uint8 x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14;
	uint8 *row0, *row1, *row2, *row3, *row4;
	uint8 *out_row0;

	r = (nch + 1) % order;
	row0 = X[nrl - 1];
	row1 = X[nrl + 0];
	for (row = nrl; row < nrh + 1; row++) {
		row2 = X[row + 1];

		out_row0 = Y[row + 0];
		for (col = ncl; col < nch + 1 - r; col += order) {
			x0  = row0[col - 1]; x1  = row0[col + 0]; x2  = row0[col + 1]; x3  = row0[col + 2]; x4  = row0[col + 3];
			x5  = row1[col - 1]; x6  = row1[col + 0]; x7  = row1[col + 1]; x8  = row1[col + 2]; x9  = row1[col + 3];
			x10 = row2[col - 1]; x11 = row2[col + 0]; x12 = row2[col + 1]; x13 = row2[col + 2]; x14 = row2[col + 3];
			

			out_row0[col + 0] = (x0  & x1  & x2 ) & 
								(x5  & x6  & x7 ) & 
								(x10 & x11 & x12);

			out_row0[col + 1] = (x1  & x2  & x3) & 
								(x6  & x7  & x8) & 
								(x11 & x12 & x13);
								
			out_row0[col + 2] = (x2  & x3  & x4) & 
								(x7  & x8  & x9) & 
								(x12 & x13 & x14); 			
		}
		row0 = row1;
		row1 = row2;
	}

	switch(r) {
		case 2:

		
		for (row = nrl; row < nrh + 1; row++) {
			row0 = X[row - 1];
			row1 = X[row + 0];
			row2 = X[row + 1];

			out_row0 = Y[row + 0];
		
			// y1 = x1 & x2 & x6 & x7 & x11 & x12; // 5
			x0  = row0[nch - 2]; x1  = row0[nch - 1]; x2  = row0[nch + 0]; x3  = row0[nch + 1]; 
			x5  = row1[nch - 2]; x6  = row1[nch - 1]; x7  = row1[nch + 0]; x8  = row1[nch + 1]; 
			x10 = row2[nch - 2]; x11 = row2[nch - 1]; x12 = row2[nch + 0]; x13 = row2[nch + 1]; 

			out_row0[nch - 1] = (x0  & x1  & x2) & 
								(x5  & x6  & x7) & 
								(x10 & x11 & x12);

			out_row0[nch + 0] = (x1  & x2  & x3) & 
			  					(x6  & x7  & x8) & 
								(x11 & x12 & x13);
			
		}
		break;
	case 1:
		for (row = nrl; row < nrh + 1; row++) {
			row0 = X[row - 1];
			row1 = X[row + 0];
			row2 = X[row + 1];

			out_row0 = Y[row + 0];
		
			// x1  = row0[nch - 1] & x2  = row0[nch + 0] & x3  = row0[nch + 1]; 
			// x6  = row1[nch - 1] & x7  = row1[nch + 0] & x8  = row1[nch + 1]; 
			// x11 = row2[nch - 1] & x12 = row2[nch + 0] & x13 = row2[nch + 1]; 

			out_row0[nch + 0] = (row0[nch - 1] &  row0[nch + 0] & row0[nch + 1])& 
			  					(row1[nch - 1] &  row1[nch + 0] & row1[nch + 1])& 
								(row2[nch - 1] &  row2[nch + 0] & row2[nch + 1]);
			
		}
		break;
	default:
		break;
	}
}
void ui8matrix_erosion_LU3x3_ExLU_O3_AddrRR(uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **temp_buffer, uint8 **Y)
{
	
	const long order = 3;
	long row = nrl, col = ncl;
	uint8 y0, y1, y2, y3, y5, r;
	uint8 x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14;
	uint8 *row0, *row1, *row2, *row3, *row4;
	uint8 *out_row0, *out_row1, *out_row2;

	r = (nrh + 1)  % order;
	
	row0 = X[row - 1];
	row1 = X[row + 0];
	for (row = nrl; row < nrh + 1 - r; row += order) {
		row2 = X[row + 1]; 
		row3 = X[row + 2]; 
		row4 = X[row + 3];
		
		
		out_row0 = Y[row + 0];
		out_row1 = Y[row + 1];
		out_row2 = Y[row + 2];
		for (col = ncl; col < nch + 1; col++) {
			x0 = row0[col - 1]; x5 = row0[col]; x10 = row0[col + 1];
			x1 = row1[col - 1]; x6 = row1[col]; x11 = row1[col + 1];
			x2 = row2[col - 1]; x7 = row2[col]; x12 = row2[col + 1];
			x3 = row3[col - 1]; x8 = row3[col]; x13 = row3[col + 1];
			x4 = row4[col - 1]; x9 = row4[col]; x14 = row4[col + 1];
			

			out_row0[col] = (x0 & x5 & x10) & (x1 & x6 & x11) & (x2 & x7 & x12); 				
			out_row1[col] = (x1 & x6 & x11) & (x2 & x7 & x12) & (x3 & x8 & x13);				
			out_row2[col] = (x2 & x7 & x12) & (x3 & x8 & x13) & (x4 & x9 & x14); 

		}
		row0 = row3;
		row1 = row4;
	}
	row2 = X[row + 1];
	row3 = X[row + 2];
	row4 = X[row + 3];
	out_row0 = Y[nrh - 1];
	out_row1 = Y[nrh + 0];
	switch(r) {
		case 2:
			for (col = ncl; col < nch + 1; col++) {
				x0 = row1[col - 1]; x3 = row1[col + 0]; x6 = row1[col + 1]; //scalar_and3(row1, col);
				x1 = row2[col - 1]; x4 = row2[col + 0]; x7 = row2[col + 1]; //scalar_and3(row2, col);
				x2 = row3[col - 1]; x5 = row3[col + 0]; x8 = row3[col + 1]; //scalar_and3(row3, col);
				out_row0[col] = scalar_and3(row0, col) & (x0 & x3 & x6) & (x1 & x4 & x7);
				out_row1[col] =        (x0 & x3 & x6) & (x1 & x4 & x7) & (x2 & x5 & x8);
			}
			break;
		case 1:
			for (col = ncl; col < nch + 1; col++) {
				x0 = row1[col - 1]; x3 = row1[col + 0]; x6 = row1[col + 1]; //scalar_and3(row1, col);
				x1 = row2[col - 1]; x4 = row2[col + 0]; x7 = row2[col + 1]; //scalar_and3(row2, col);
				out_row1[col] = scalar_and3(row0, col) & (x0 & x3 & x6) & (x1 & x4 & x7);
			}				
			break;
		default:
			break;
	}
}
void ui8matrix_erosion_LU3x3_ComLU_O3_AddrRR(uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **temp_buffer, uint8 **Y)
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
	
	row0 = X[row - 1];
	row1 = X[row + 0];

	for (row = nrl; row < nrh + 1 - rr; row += order) {
		row2 = X[row + 1];
		row3 = X[row + 2];
		row4 = X[row + 3];
		out_row0 = Y[row + 0];
		out_row1 = Y[row + 1];
		out_row2 = Y[row + 2];



		for (col = ncl; col < nch + 1 - cr; col += order) {
			x0  = row0[col - 1]; x1  = row0[col + 0]; x2  = row0[col + 1]; x3  = row0[col + 2]; x4  = row0[col + 3];
			x5  = row1[col - 1]; x6  = row1[col + 0]; x7  = row1[col + 1]; x8  = row1[col + 2]; x9  = row1[col + 3];
			x10 = row2[col - 1]; x11 = row2[col + 0]; x12 = row2[col + 1]; x13 = row2[col + 2]; x14 = row2[col + 3];
			x15 = row3[col - 1]; x16 = row3[col + 0]; x17 = row3[col + 1]; x18 = row3[col + 2]; x19 = row3[col + 3];
			x20 = row4[col - 1]; x21 = row4[col + 0]; x22 = row4[col + 1]; x23 = row4[col + 2]; x24 = row4[col + 3];


			y0  = x0  & x5  & x10; y1  = x1  & x6  & x11; y2  = x2  & x7  & x12; y3  = x3  & x8  & x13; y4  = x4  & x9  & x14; 
			y5  = x5  & x10 & x15; y6  = x6  & x11 & x16; y7  = x7  & x12 & x17; y8  = x8  & x13 & x18; y9  = x9  & x14 & x19; 
			y10 = x10 & x15 & x20; y11 = x11 & x16 & x21; y12 = x12 & x17 & x22; y13 = x13 & x18 & x23; y14 = x14 & x19 & x24; 
 
			out_row0[col + 0] = y0 & y1 & y2; 
			out_row0[col + 1] = y1 & y2 & y3; 
			out_row0[col + 2] = y2 & y3 & y4;

			out_row1[col + 0] = y5 & y6 & y7;
			out_row1[col + 1] = y6 & y7 & y8;
			out_row1[col + 2] = y7 & y8 & y9;
			
			out_row2[col + 0] = y10 & y11 & y12;
			out_row2[col + 1] = y11 & y12 & y13;
			out_row2[col + 2] = y12 & y13 & y14;
		}
		row0 = row3;
		row1 = row4;
	}
	switch (rr) {
		case 2 :
			switch(cr){
				case 2 :
					ui8matrix_erosion_LU3x3_ExLU_O3_AddrRR(X, nrl, nrh, nch - 1, nch, temp_buffer, Y);
					ui8matrix_erosion_LU3x3_InLU_O3_AddrRR(X, nrh - 1, nrh, ncl, nch, temp_buffer, Y);
				break;
				case 1: 
					ui8matrix_erosion_LU3x3_ExLU_O3_AddrRR(X, nrl, nrh, nch, nch, temp_buffer, Y);		
					ui8matrix_erosion_LU3x3_InLU_O3_AddrRR(X, nrh - 1, nrh, ncl, nch, temp_buffer, Y);	
				break;
				default :
					ui8matrix_erosion_LU3x3_InLU_O3_AddrRR(X, nrh - 1, nrh, ncl, nch, temp_buffer, Y);
				break;
			}
			break;
		case 1:
			switch(cr){
				case 2 :
					ui8matrix_erosion_LU3x3_ExLU_O3_AddrRR(X, nrl, nrh, nch - 1, nch, temp_buffer, Y);
					ui8matrix_erosion_LU3x3_InLU_O3_AddrRR(X, nrh, nrh, ncl, nch, temp_buffer, Y);
				break;
				case 1: 
					ui8matrix_erosion_LU3x3_ExLU_O3_AddrRR(X, nrl, nrh, nch, nch, temp_buffer, Y);		
					ui8matrix_erosion_LU3x3_InLU_O3_AddrRR(X, nrh, nrh, ncl, nch, temp_buffer, Y);	
					// display_ui8matrix(Y, nrl, nrh, ncl, nch, "%u", "rr=1; cr=1");
				break;
				default :
					ui8matrix_erosion_LU3x3_InLU_O3_AddrRR(X, nrh, nrh, ncl, nch, temp_buffer, Y);
				break;
			}
		break;
		default:
			switch(cr){
				case 2 :
					ui8matrix_erosion_LU3x3_ExLU_O3_AddrRR(X, nrl, nrh, nch - 1, nch, temp_buffer, Y);
				break;
				case 1: 
					ui8matrix_erosion_LU3x3_ExLU_O3_AddrRR(X, nrl, nrh, nch, nch, temp_buffer, Y);	
				break;
				default :
				break;
			}
		break;
	}
		
}
/************************************************************************************/
/******** End of optimisation : Loop Unroll + Register Rotation of Addresses ********/
/************************************************************************************/

/********************************************************/
/******* Optimisation :  Pipelining ********/
/********************************************************/
void ui8matrix_erosion_row_pipeline(uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **temp_buffer, uint8 **Y)
{
	long row = nrl, col = ncl, x, y;
	// uint8 **temp = ui8matrix(nrl + (-1), nrh + 1, ncl + (-1), nch + 1);
	uint8 *temp_row0, *temp_row1, *temp_row2, *out_row0;
	// Prologue	
	temp_row0 = temp_buffer[nrl - 1];
	temp_row1 = temp_buffer[nrl + 0];
	for (col = ncl; col < nch + 1; col++) {
    	temp_row0[col] = scalar_and3(X[row - 1], col);
		temp_row1[col] = scalar_and3(X[row + 0], col);
	}	

	for (row = nrl + 1; row < nrh + 1; row++){
		out_row0 = Y[row - 1];
		temp_row0 = temp_buffer[row - 2];
		temp_row1 = temp_buffer[row - 1];
		temp_row2 = temp_buffer[row - 0];
		for (col = ncl; col < nch + 1; col++) {
			temp_row2[col] = scalar_and3(X[row + 0], col);
			out_row0[col] = temp_row0[ col] &
							temp_row1[ col] &
							temp_row2[ col];
		}
	}
	temp_row0 = temp_buffer[nrh - 1];
	temp_row1 = temp_buffer[nrh + 0];
	temp_row2 = X[nrh + 1];
	out_row0 = Y[nrh + 0];
	for (col = ncl; col < nch + 1; col++) {

    	out_row0[col] = temp_row0[ col] &
				  	    temp_row1[ col] &
		     scalar_and3(temp_row2, col);
	}
	// free_ui8matrix(temp, nrl + (-1), nrh + 1, ncl + (-1), nch + 1);
}

void ui8matrix_erosion_col_pipeline(uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **temp_buffer, uint8 **Y)
{
	long row = nrl, col = ncl, x, y;
	// uint8 **temp = ui8matrix(nrl + (-1), nrh + 1, ncl + (-1), nch + 1);
	uint8 *in_row0, *in_row1, *in_row2, *temp_row, *out_row0;
	uint8 x0, x1, x2;
	// Prologue	
	for (row = nrl; row < nrh + 1; row++) {
		in_row0 = X[row - 1];
		in_row1 = X[row + 0];
		in_row2 = X[row + 1];
		temp_row = temp_buffer[row];
    	temp_row[ncl - 1] = in_row0[ncl - 1] & 
						    in_row1[ncl - 1] & 
							in_row2[ncl - 1];
    	temp_row[ncl + 0] = in_row0[ncl + 0] & 
							in_row1[ncl + 0] & 
							in_row2[ncl + 0];
	}	
	

	for (row = nrl; row < nrh + 1; row++){
		temp_row = temp_buffer[row];
		
		in_row0 = X[row - 1];
		in_row1 = X[row + 0];
		in_row2 = X[row + 1];

		out_row0 = Y[row]; 
		
		for (col = ncl; col < nch + 1; col ++) {
			x0 = temp_row[col - 1];		// LOAD => A
			x1 = temp_row[col - 0];		// LOAD => B
			x2 = in_row0[col + 1] & in_row1[col + 1] & in_row2[col + 1];
			
			out_row0[col + 0] = x0 & x1 & x2; // A & B & C
			temp_row[col + 1] = x2;		  
		}
	}
	// free_ui8matrix(temp, nrl + (-1), nrh + 1, ncl + (-1), nch + 1);
}

void ui8matrix_erosion_col_pipeline_RR(uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **temp_buffer, uint8 **Y)
{
	long row = nrl, col = ncl, x, y;
	// uint8 **temp = ui8matrix(nrl + (-1), nrh + 1, ncl + (-1), nch + 1);
	uint8 *in_row0, *in_row1, *in_row2, *temp_row, *out_row0;
	uint8 x0, x1, x2;
	// Prologue	
	for (row = nrl; row < nrh + 1; row++) {
		in_row0 = X[row - 1];
		in_row1 = X[row + 0];
		in_row2 = X[row + 1];
		temp_row = temp_buffer[row];
    	temp_row[ncl - 1] = in_row0[ncl - 1] & 
						    in_row1[ncl - 1] & 
							in_row2[ncl - 1];
    	temp_row[ncl + 0] = in_row0[ncl + 0] & 
							in_row1[ncl + 0] & 
							in_row2[ncl + 0];
	}	
	

	for (row = nrl; row < nrh + 1; row++){
		temp_row = temp_buffer[row];
		
		in_row0 = X[row - 1];
		in_row1 = X[row + 0];
		in_row2 = X[row + 1];

		out_row0 = Y[row]; 
		x0 = temp_row[ncl - 1];		// LOAD => A
		x1 = temp_row[ncl - 0];		// LOAD => B
		for (col = ncl; col < nch + 1; col ++) {
			x2 = in_row0[col + 1] & in_row1[col + 1] & in_row2[col + 1];
			out_row0[col + 0] = x0 & x1 & x2; // A & B & C
			temp_row[col + 1] = x2;	
			x0 = x1;
			x1 = x2;	  
		}
	}
	// for (row = nrl; row < nrh + 1; row++){
	// 	temp_row = temp_buffer[row];
	// 	Y[row][nch] = temp_row[nch - 1] & 
	// 						 temp_row[nch - 0] &
	// 				 X[row - 1][nch + 1] & 
	// 				 X[row + 0][nch + 1] &
	// 				 X[row + 1][nch + 1];
	// }
	// free_ui8matrix(temp, nrl + (-1), nrh + 1, ncl + (-1), nch + 1);
}
/*************************************************/
/******* End of optimisation : Pipelining ********/
/*************************************************/


/********************************************************/
/******* Optimisation : Loop Unroll + Pipelining ********/
/********************************************************/

void ui8matrix_erosion_pipeline_LU3x3_ExLU_O3 (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **temp_buffer, uint8 **Y)
{
	const long order = 3;
	long row = nrl, col = ncl, x, y, r = 0;
	
	uint8 *row0, *row1, *row2, *row3, x0, x1, x2, x3, x4;
	uint8 *temp_row0, *temp_row1, *temp_row2, *temp_row3, *temp_row4;
	uint8 *out_row0, *out_row1, *out_row2, *out_row3;
	
	r = (nrh + 1) % order;
	// Prologue	
	for (row = nrl - 1; row < nrh - r + 1; row += order){
		row0 = X[row + 0];
		row1 = X[row + 1];
		temp_row1 = temp_buffer[row + 0];
		temp_row2 = temp_buffer[row + 1];
		for (col = ncl; col < nch + 1; col++) {
			temp_row1[col] = scalar_and3(row0, col);
			temp_row2[col] = scalar_and3(row1, col);
		}
	}

	row0 = X[row - 1];
	row1 = X[row + 0];
	row2 = X[row + 1];
	
	temp_row2 = temp_buffer[row - 1];
	temp_row3 = temp_buffer[row + 0];
	temp_row4 = temp_buffer[row + 1];
	

	for (row = nrl; row < nrh + 1 - r; row += order){
		// row0 = X[row - 1];   //  0 : b0 b1 b2 ...
		row1 = X[row + 1];   //  1 : c0 c1 c2 ...
		// row2 = X[row + 1];   //  2 : d0 d1 d2 ...
		// row3 = X[row + 2];   //  3 : e0 e1 e2 ...

		temp_row0 = temp_buffer[row - 1]; // -1 : A    <= (a0 a1 a2 ...) Prologued
		temp_row1 = temp_buffer[row + 0]; //  0 : B    <= (b0 b1 b2 ...) Prologued
		temp_row2 = temp_buffer[row + 1]; //  1 : temp_buffer <= (c0 c1 c2 ...)
		temp_row3 = temp_buffer[row + 2]; //  2 : D	<= (d0 d1 d2 ...) Prologued
		temp_row4 = temp_buffer[row + 3]; //  3 : E	<= (e0 e1 e2 ...) Prologued

		out_row0 = Y[row + 0]; 
		out_row1 = Y[row + 1];
		out_row2 = Y[row + 2];

		for (col = ncl; col < nch + 1; col++) {
			x0 = temp_row0[col];		// LOAD => A										
			x1 = temp_row1[col];		// LOAD => B
			x2 = scalar_and3(row1, col); // OR3  => C
			x3 = temp_row3[col]; 		// LOAD => D
			x4 = temp_row4[col]; 		// LOAD => E
			
			out_row0[col] = x0 & x1 & x2; // A & B & C
			out_row1[col] = x1 & x2 & x3; // B & C & D
			out_row2[col] = x2 & x3 & x4; // C & D & E

			temp_row2[col] = x2;		// STORE => C 
										// NEXT  :  F 
										// NEXT  :  G 
										// NEXT  :  H
		}
	}
	// epilogue	
	for (row = nrh - r; row < nrh + 1; row++){
		row0 = X[row + 0];
		temp_row3 = temp_buffer[row + 0];
		for (col = ncl; col < nch + 1; col++) 
			temp_row3[col] = scalar_and3(row0, col);
	}

	switch(r) {
		case 1: 
		temp_row1 = temp_buffer[nrh - 1]; //  0 : B    <= (b0 b1 b2 ...) Prologued
		temp_row2 = temp_buffer[nrh + 0]; //  1 : temp_buffer <= (c0 c1 c2 ...)
		// temp_row3 = temp_buffer[nrh + 1]; //  2 : D	<= (d0 d1 d2 ...) Prologued
		row3 = X[nrh + 1];
		out_row1 = Y[nrh + 0];
		
		for (col = ncl; col < nch + 1; col++) {											
			out_row1[col] = temp_row1[col] & temp_row2[col] & scalar_and3(row3, col); 		
		}
		
		break;
		case 2: 
		temp_row0 = temp_buffer[nrh - 2]; // -1 : A    <= (a0 a1 a2 ...) Prologued
		temp_row1 = temp_buffer[nrh - 1]; //  0 : B    <= (b0 b1 b2 ...) Prologued
		temp_row2 = temp_buffer[nrh - 0]; //  1 : temp_buffer <= (c0 c1 c2 ...)
		row3 = X[nrh + 1];
		
		out_row0 = Y[nrh - 1]; 
		out_row1 = Y[nrh + 0];
		out_row2 = Y[nrh + 1];

		for (col = ncl; col < nch + 1; col++) {									
			out_row0[col] = temp_row0[col] & temp_row1[col] & temp_row2[col]; 		
			out_row1[col] = temp_row1[col] & temp_row2[col] & scalar_and3(row3, col); 	
		}
		break;
		default:
		break;

	}
	// display_ui8matrix(X,nrl + (-1), nrh + 1, ncl + (-1), nch + 1, "%u", "In");
	// display_ui8matrix(Y,nrl, nrh , ncl, nch , "%03u ", "Out");
	// display_ui8matrix(temp,nrl + (-1), nrh + 1, ncl, nch, "%u ", "Out");
	// free_ui8matrix(temp,nrl + (-1), nrh + 1, ncl, nch);
}
void ui8matrix_erosion_pipeline_LU3x3_InLU_O3 (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **temp_buffer, uint8 **Y)
{
	const long order = 3;
	long row = nrl, col = ncl, x, y, r = 0;
	// uint8 **temp = ui8matrix(nrl, nrh, ncl + (-1), nch + 1);
	uint8 *row0, *row1, *row2, *row3, x0, x1, x2, x3, x4, x5;
	uint8 *temp_row0, *temp_row1, *temp_row2, *temp_row3, *temp_row4;
	uint8 *out_row0, *out_row1, *out_row2, *out_row3;
	
	r = (nch + 1) % order;
	// Prologue	
		
	for (row = nrl; row < nrh + 1; row++) {
		row0 = X[row - 1];
		row1 = X[row + 0];
		row2 = X[row + 1];

		temp_row1 = temp_buffer[row + 0];
		for (col = ncl - 1; col < nch - r + 1 ; col += order) {
			temp_row1[col + 0] = row0[col + 0] & row1[col + 0] & row2[col + 0];
			temp_row1[col + 1] = row0[col + 1] & row1[col + 1] & row2[col + 1];
			
		}
	}

	for (row = nrl; row < nrh + 1; row ++){
		temp_row1 = temp_buffer[row + 0];
		row0 = X[row - 1];
		row1 = X[row + 0];
		row2 = X[row + 1];

		out_row1 = Y[row]; 

		for (col = ncl; col < nch + 1 - r; col += order) {
			x0 = temp_row1[col - 1];
			x1 = temp_row1[col + 0];
			x2 = row0[col + 1] & row1[col + 1] & row2[col + 1];
			x3 = temp_row1[col + 2];
			x4 = temp_row1[col + 3];
			
			out_row1 [col + 0] = x0 & x1 & x2; // A & B & C
			out_row1 [col + 1] = x1 & x2 & x3; // B & C & D
			out_row1 [col + 2] = x2 & x3 & x4; // C & D & E
			temp_row1[col] = x2;

		}
	}
	// epilogue	
	for (row = nrl; row < nrh + 1; row++){
		row0 = X[row - 1];
		row1 = X[row + 0];
		row2 = X[row + 1];
		temp_row3 = temp_buffer[row + 0];
		for (col = nch - r; col < nch + 1; col++) 
			temp_row3[col] = row0[col] & row1[col] & row2[col];
	}

	switch(r){
		case 2:
		for (row = nrl; row < nrh + 1; row++) {											
			temp_row0 = temp_buffer[row];
			row0 = X[row - 1];
			row1 = X[row + 0];
			row2 = X[row + 1];
			out_row0 = Y[row];
			out_row0[nch - 1] = temp_row0[nch - 2] & 
								temp_row0[nch - 1] & 
								temp_row0[nch - 0] ;

			out_row0[nch + 0] = temp_row0[nch - 1] & 
								temp_row0[nch - 0] &
								     row0[nch + 1] & 	
									 row1[nch + 1] & 	
									 row2[nch + 1]; 	
		}
		break;
		case 1:
		
		for (row = nrl; row < nrh + 1; row++) {											
			temp_row0 = temp_buffer[row];
			row0 = X[row - 1];
			row1 = X[row + 0];
			row2 = X[row + 1];
			out_row0 = Y[row];
			out_row0[nch + 0] = temp_row0[nch - 1] & 
								temp_row0[nch - 0] &
								     row0[nch + 1] & 	
									 row1[nch + 1] & 	
									 row2[nch + 1]; 	
		}
		break;
		default:
		break;
	}
	// display_ui8matrix(X,nrl + (-1), nrh + 1, ncl + (-1), nch + 1, "%u", "In");
	// display_ui8matrix(Y,nrl, nrh, ncl, nch , "%03u ", "Out");
	// display_ui8matrix(temp,nrl, nrh, ncl + (-1), nch + 1, "%u ", "Temp");

	// ui8matrix_erosion_LU3x3_O1xO1(X, nrl, nrh, ncl, nch, temp_buffer, Y);
	// free_ui8matrix(temp, nrl, nrh, ncl + (-1), nch + 1);
	
}

/***************************************************************/
/******* End of optimisation : Loop Unroll + Pipelining ********/
/***************************************************************/

/****************************************************************************/
/******* Optimisation : Loop Unroll + Register Rotation + Pipelining ********/
/****************************************************************************/

void ui8matrix_erosion_pipeline_LU3x3_ExLU_O3_RR (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **temp_buffer, uint8 **Y)
{
	const long order = 3;
	long row = nrl, col = ncl, x, y, r = 0;
	
	uint8 *row0, *row1, *row2, *row3;
	uint8 y0, y1, y2, y3, y4, x0, x1, x2, x3, x4, x5;
	uint8 *temp_row0, *temp_row1, *temp_row2, *temp_row3, *temp_row4;
	uint8 *out_row0, *out_row1, *out_row2, *out_row3;
	
	r = (nrh + 1) % order;
	// Prologue	
	for (row = nrl - 1; row < nrh - r + 1; row += order){
		row0 = X[row + 0];
		row1 = X[row + 1];
		temp_row1 = temp_buffer[row + 0];
		temp_row2 = temp_buffer[row + 1];
		x0  = row0[ncl - 1]; x1 = row0[ncl + 0];
		x3  = row1[ncl - 1]; x4 = row1[ncl + 0];
		for (col = ncl; col < nch + 1; col++) {
			x2  = row0[col + 1];
			x5  = row1[col + 1];
			temp_row1[col] = scalar_and3(row0, col);
			temp_row2[col] = scalar_and3(row1, col);
			x0 = x1; x1 = x2;
			x3 = x4; x4 = x5;
		}
	}

	row0 = X[row - 1];
	row1 = X[row + 0];
	row2 = X[row + 1];
	
	temp_row2 = temp_buffer[row - 1];
	temp_row3 = temp_buffer[row + 0];
	temp_row4 = temp_buffer[row + 1];
	

	for (row = nrl; row < nrh + 1 - r; row += order){
		// row0 = X[row - 1];   //  0 : b0 b1 b2 ...
		row1 = X[row + 1];   //  1 : c0 c1 c2 ...
		// row2 = X[row + 1];   //  2 : d0 d1 d2 ...
		// row3 = X[row + 2];   //  3 : e0 e1 e2 ...

		temp_row0 = temp_buffer[row - 1]; // -1 : A    <= (a0 a1 a2 ...) Prologued
		temp_row1 = temp_buffer[row + 0]; //  0 : B    <= (b0 b1 b2 ...) Prologued
		temp_row2 = temp_buffer[row + 1]; //  1 : temp_buffer <= (c0 c1 c2 ...)
		temp_row3 = temp_buffer[row + 2]; //  2 : D	<= (d0 d1 d2 ...) Prologued
		temp_row4 = temp_buffer[row + 3]; //  3 : E	<= (e0 e1 e2 ...) Prologued

		out_row0 = Y[row + 0]; 
		out_row1 = Y[row + 1];
		out_row2 = Y[row + 2];

		x3  = row1[ncl - 1]; x4  = row1[ncl + 0]; 
		for (col = ncl; col < nch + 1; col++) {
			x5  = row1[col + 1];
			y0 = temp_row0[col];		// LOAD => A										
			y1 = temp_row1[col];		// LOAD => B
			y2 = x3 & x4 & x5;
			y3 = temp_row3[col]; 		// LOAD => D
			y4 = temp_row4[col]; 		// LOAD => E
			

			out_row0[col] = y0 & y1 & y2; // A & B & C
			out_row1[col] = y1 & y2 & y3; // B & C & D
			out_row2[col] = y2 & y3 & y4; // C & D & E

			temp_row2[col] = y2;		// STORE => C 
										// NEXT  :  F 
										// NEXT  :  G 
										// NEXT  :  H
			x3 = x4; x4 = x5;
		}
	}
	// epilogue	
	for (row = nrh - r; row < nrh + 1; row++){
		row0 = X[row + 0];
		temp_row3 = temp_buffer[row + 0];
		for (col = ncl; col < nch + 1; col++) 
			temp_row3[col] = scalar_and3(row0, col);
	}

	switch(r) {
		case 1: 
		temp_row1 = temp_buffer[nrh - 1]; //  0 : B    <= (b0 b1 b2 ...) Prologued
		temp_row2 = temp_buffer[nrh + 0]; //  1 : temp_buffer <= (c0 c1 c2 ...)
		// temp_row3 = temp_buffer[nrh + 1]; //  2 : D	<= (d0 d1 d2 ...) Prologued
		row3 = X[nrh + 1];
		out_row1 = Y[nrh + 0];
		
		for (col = ncl; col < nch + 1; col++) {											
			out_row1[col] = temp_row1[col] & temp_row2[col] & scalar_and3(row3, col); 		
		}
		
		break;
		case 2: 
		temp_row0 = temp_buffer[nrh - 2]; // -1 : A    <= (a0 a1 a2 ...) Prologued
		temp_row1 = temp_buffer[nrh - 1]; //  0 : B    <= (b0 b1 b2 ...) Prologued
		temp_row2 = temp_buffer[nrh - 0]; //  1 : temp_buffer <= (c0 c1 c2 ...)
		row3 = X[nrh + 1];
		
		out_row0 = Y[nrh - 1]; 
		out_row1 = Y[nrh + 0];
		out_row2 = Y[nrh + 1];

		for (col = ncl; col < nch + 1; col++) {									
			out_row0[col] = temp_row0[col] & temp_row1[col] & temp_row2[col]; 		
			out_row1[col] = temp_row1[col] & temp_row2[col] & scalar_and3(row3, col); 	
		}
		break;
		default:
		break;

	}
	// free_ui8matrix(temp,nrl + (-1), nrh + 1, ncl, nch);
}
void ui8matrix_erosion_pipeline_LU3x3_InLU_O3_RR (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **temp_buffer, uint8 **Y)
{
	const long order = 3;
	long row = nrl, col = ncl, x, y, r = 0;
	// uint8 **temp = ui8matrix(nrl, nrh, ncl + (-1), nch + 1);
	uint8 *row0, *row1, *row2, *row3, x0, x1, x2, x3, x4, x5;
	uint8 *temp_row0, *temp_row1, *temp_row2, *temp_row3, *temp_row4;
	uint8 *out_row0, *out_row1, *out_row2, *out_row3;
	
	r = (nch + 1) % order;
	// Prologue	
		
	row0 = X[nrl - 1];
	row1 = X[nrl + 0];
	for (row = nrl; row < nrh + 1; row++) {
		row2 = X[row + 1];

		temp_row1 = temp_buffer[row + 0];
		for (col = ncl - 1; col < nch - r + 1 ; col += order) {
			temp_row1[col + 0] = row0[col + 0] & row1[col + 0] & row2[col + 0];
			temp_row1[col + 1] = row0[col + 1] & row1[col + 1] & row2[col + 1];
			
		}
		row0 = row1;
		row1 = row2;
	}

	row0 = X[nrl - 1];
	row1 = X[nrl + 0];
	for (row = nrl; row < nrh + 1; row ++){
		temp_row1 = temp_buffer[row + 0];
		row2 = X[row + 1];

		out_row1 = Y[row]; 

		x0 = temp_row1[ncl - 1];
		x1 = temp_row1[ncl + 0];
		for (col = ncl; col < nch + 1 - r; col += order) {
			x2 = row0[col + 1] & row1[col + 1] & row2[col + 1];
			x3 = temp_row1[col + 2];
			x4 = temp_row1[col + 3];
			
			out_row1 [col + 0] = x0 & x1 & x2; // A & B & C
			out_row1 [col + 1] = x1 & x2 & x3; // B & C & D
			out_row1 [col + 2] = x2 & x3 & x4; // C & D & E
			temp_row1[col] = x2;

			x0 = x3;
			x1 = x4;

		}
		row0 = row1;
		row1 = row2;
	}
	// epilogue	
	for (row = nrl; row < nrh + 1; row++){
		row0 = X[row - 1];
		row1 = X[row + 0];
		row2 = X[row + 1];
		temp_row3 = temp_buffer[row + 0];
		for (col = nch - r; col < nch + 1; col++) 
			temp_row3[col] = row0[col] & row1[col] & row2[col];
	}

	switch(r){
		case 2:
		row0 = X[nrl - 1];
		row1 = X[nrl + 0];
		for (row = nrl; row < nrh + 1; row++) {											
			temp_row0 = temp_buffer[row];
			row2 = X[row + 1];
			out_row0 = Y[row];
			out_row0[nch - 1] = temp_row0[nch - 2] & 
								temp_row0[nch - 1] & 
								temp_row0[nch - 0] ;

			out_row0[nch + 0] = temp_row0[nch - 1] & 
								temp_row0[nch - 0] &
								     row0[nch + 1] & 	
									 row1[nch + 1] & 	
									 row2[nch + 1]; 
			row0 = row1;
			row1 = row2;	
		}
		break;
		case 1:
		
		row0 = X[nrl - 1];
		row1 = X[nrl + 0];
		for (row = nrl; row < nrh + 1; row++) {											
			temp_row0 = temp_buffer[row];
			row0 = X[row - 1];
			row1 = X[row + 0];
			row2 = X[row + 1];
			out_row0 = Y[row];
			out_row0[nch + 0] = temp_row0[nch - 1] & 
								temp_row0[nch - 0] &
								     row0[nch + 1] & 	
									 row1[nch + 1] & 	
									 row2[nch + 1]; 	
			row0 = row1;
			row1 = row2;
		}
		break;
		default:
		break;
	}
	// display_ui8matrix(X,nrl + (-1), nrh + 1, ncl + (-1), nch + 1, "%u", "In");
	// display_ui8matrix(Y,nrl, nrh, ncl, nch , "%03u ", "Out");
	// display_ui8matrix(temp,nrl, nrh, ncl + (-1), nch + 1, "%u ", "Temp");

	// ui8matrix_erosion_LU3x3_O1xO1(X, nrl, nrh, ncl, nch, temp_buffer, Y);
	// free_ui8matrix(temp, nrl, nrh, ncl + (-1), nch + 1);
	
}
void ui8matrix_erosion_divide_row_and_conquer (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **temp_buffer, uint8 **Y)
{
	const long order = 3;
	
	uint8 *temp_row, *temp_row0, *temp_row1, *temp_row2, *out_row, *in_row;
	long row, col, nrow, r;

	nrow = nrh - nrl + 1;
	row = 0;
	
	for (row = nrl - 1; row < (nrh + 1) + 1; row ++)
	{
		temp_row =    temp_buffer[row];
		in_row   = X[row];
		for(col = ncl; col < nch + 1; col++) {
			temp_row[col] = scalar_and3(in_row, col);
		}
	}
	
	for (row = nrl; row < nrh + 1; row ++)
	{
		temp_row0 =     temp_buffer[row - 1];
		temp_row1 =     temp_buffer[row + 0];
		temp_row2 =     temp_buffer[row + 1];
		out_row   = Y[row];
		for(col = ncl; col < nch + 1; col++) {
			out_row[col] = temp_row0[col]&
						   temp_row1[col]&
						   temp_row2[col];
		}
	}
	// free_ui8matrix(temp, nrl + (-1), nrh + 1, ncl, nch);
}
/** Deprcated **/
void ui8matrix_erosion_unrolled_divide_row_and_conquer (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **temp_buffer, uint8 **Y)
{
	const long order = 3;
	
	uint8 *temp_row, *temp_row0, *temp_row1, *temp_row2, *out_row, *in_row;
	long row, col, nrow, r;
	int nb_threads;

	nrow = nrh - nrl + 1;
	row = 0;
	
	for (row = nrl - 1; row < nrh + 1; row += order)
	{
		temp_row =    temp_buffer[row];
		in_row   = X[row];
		for(col = ncl; col < nch + 1; col++) {
			temp_row[col] = scalar_and3(in_row, col);
		}
	}
	for (row = nrl + 0; row < nrh + 1; row += order)
	{

		temp_row =    temp_buffer[row];
		in_row   = X[row];
		for(col = ncl; col < nch + 1; col++) {
			temp_row[col] = scalar_and3(in_row, col);
		}
	}
	for (row = nrl + 1; row < (nrh + 1) + 1; row += order)
	{
		temp_row =    temp_buffer[row];
		in_row   = X[row];
		for(col = ncl; col < nch + 1; col++) {
			temp_row[col] = scalar_and3(in_row, col);
		}
	}
	// Epilogue
	temp_row = temp_buffer[nrh + 1];
	in_row = X[nrh + 1];
	for(col = ncl; col < nch + 1; col++) {
		temp_row[col] = scalar_and3(in_row, col);
	}
	
	for (row = nrl; row < nrh + 1; row ++)
	{
		temp_row0 =     temp_buffer[row - 1];
		temp_row1 =     temp_buffer[row + 0];
		temp_row2 =     temp_buffer[row + 1];
		out_row   = Y[row];
		for(col = ncl; col < nch + 1; col++) {
			out_row[col] = temp_row0[col]&
						   temp_row1[col]&
						   temp_row2[col];
		}
	}
	// free_ui8matrix(temp, nrl + (-1), nrh + 1, ncl, nch);
}
void ui8matrix_erosion_divide_row_and_conquer_ExLU_O3 (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **temp_buffer, uint8 **Y)
{
	const long order = 3;
	
	uint8 *temp_row, *temp_row0, *temp_row1, *temp_row2, *temp_row3, *temp_row4, *in_row, *in_row0, *in_row1, *in_row2, *out_row;
	long row, col, nrow, r;
	int nb_threads;

	nrow = nrh - nrl + 1;
	row = 0;

	r = (nrh + 1) % order;
	temp_row0 =    temp_buffer[nrl - 1];
	in_row0   = X[nrl - 1];
	for(col = ncl; col < nch + 1; col++) {
		temp_row0[col] = scalar_and3(in_row0, col);
	}
	for (row = nrl; row < nrh + 1 - r; row += order)
	{
		temp_row0 =    temp_buffer[row + 0];
		temp_row1 =    temp_buffer[row + 1];
		temp_row2 =    temp_buffer[row + 2];
		in_row0    = X[row + 0];
		in_row1    = X[row + 1];
		in_row2    = X[row + 2];
		for(col = ncl; col < nch + 1; col++) {
			temp_row0[col] = scalar_and3(in_row0, col);
			temp_row1[col] = scalar_and3(in_row1, col);
			temp_row2[col] = scalar_and3(in_row2, col);
		}
	}
	// Epilogue
	switch(r) {
		case 2 :
		temp_row0  = temp_buffer[nrh - 1];
		in_row0 = X[nrh - 1];
		temp_row1  = temp_buffer[nrh + 0];
		in_row1 = X[nrh + 0];
		temp_row2  = temp_buffer[nrh + 1];
		in_row2 = X[nrh + 1];
		for(col = ncl; col < nch + 1; col++) {
			temp_row0[col] = scalar_and3(in_row0, col);
			temp_row1[col] = scalar_and3(in_row1, col);
			temp_row2[col] = scalar_and3(in_row2, col);
		}
		break;
		
		case 1 :
		temp_row1  = temp_buffer[nrh + 0];
		in_row1 = X[nrh + 0];
		temp_row2  = temp_buffer[nrh + 1];
		in_row2 = X[nrh + 1];
		for(col = ncl; col < nch + 1; col++) {
			temp_row1[col] = scalar_and3(in_row1, col);
			temp_row2[col] = scalar_and3(in_row2, col);
		}
		break;
		case 0 :
		temp_row2  = temp_buffer[nrh + 1];
		in_row2 = X[nrh + 1];
		for(col = ncl; col < nch + 1; col++) {
			temp_row2[col] = scalar_and3(in_row2, col);
		}
		break;
	}
	// display_ui8matrix(temp, nrl + (-1), nrh + 1, ncl, nch, "%u", "Temp");
	for (row = nrl; row < nrh + 1; row ++)
	{
		temp_row0 =     temp_buffer[row - 1];
		temp_row1 =     temp_buffer[row + 0];
		temp_row2 =     temp_buffer[row + 1];
		out_row   = Y[row];
		for(col = ncl; col < nch + 1; col++) {
			out_row[col] = temp_row0[col]&
						   temp_row1[col]&
						   temp_row2[col];
		}
	}

	
	// free_ui8matrix(temp, nrl + (-1), nrh + 1, ncl, nch);
}
void ui8matrix_erosion_divide_row_and_conquer_InLU_O3 (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **temp_buffer, uint8 **Y)
{
	const long order = 3;
	
	uint8 *temp_row, *temp_row0, *temp_row1, *temp_row2, *temp_row3, *temp_row4, *out_row0, *out_row1, *out_row2, *in_row, *out_row;
	long row, col, nrow, r;
	int nb_threads;

	nrow = nrh - nrl + 1;
	row = 0;
	r = (nch + 1) % order;
	for (row = nrl - 1; row < (nrh + 1) + 1; row ++)
	{
		temp_row =    temp_buffer[row];
		in_row   = X[row];
		for(col = ncl; col < nch + 1 - r; col += order) {
			temp_row[col + 0] = scalar_and3(in_row, col + 0);
			temp_row[col + 1] = scalar_and3(in_row, col + 1);
			temp_row[col + 2] = scalar_and3(in_row, col + 2);
		}
		
	}
	switch(r) {
		case 2:
		for (row = nrl - 1; row < (nrh + 1) + 1; row ++)
		{
			temp_row =    temp_buffer[row];
			in_row   = X[row];
			temp_row[nch - 1] = scalar_and3(in_row, nch - 1);
			temp_row[nch + 0] = scalar_and3(in_row, nch + 0);
		}
		break;
		case 1:
		for (row = nrl - 1; row < (nrh + 1) + 1; row ++)
		{
			temp_row =    temp_buffer[row];
			in_row   = X[row];
			temp_row[nch + 0] = scalar_and3(in_row, nch + 0);
		}
		break;
		default: 
		break;
	}
	// Epilogue
	// temp_row = temp_buffer[nrh + 1];
	// in_row = X[nrh + 1];
	// for(col = ncl; col < nch + 1; col++) {
	// 	temp_row[col] = scalar_and3(in_row, col);
	// }
	for (row = nrl; row < nrh + 1; row ++)
	{
		temp_row0 =     temp_buffer[row - 1];
		temp_row1 =     temp_buffer[row + 0];
		temp_row2 =     temp_buffer[row + 1];
		out_row   = Y[row];
		for(col = ncl; col < nch + 1; col++) {
			out_row[col] = temp_row0[col]&
						   temp_row1[col]&
						   temp_row2[col];
		}
	}
	
		
	// free_ui8matrix(temp, nrl + (-1), nrh + 1, ncl, nch);
}
void ui8matrix_erosion_divide_col_and_conquer (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **temp_buffer, uint8 **Y)
{
	const long order = 3;
	// uint8 **temp = ui8matrix(nrl, nrh, ncl + (-1), nch + 1);
	uint8 *temp_row, *temp_row0, *temp_row1, *temp_row2, *out_row, *in_row0, *in_row1, *in_row2;
	long row, col, nrow, r;
	int nb_threads;

	nrow = nrh - nrl + 1;
	row = 0;
	
	for (row = nrl; row < nrh + 1; row ++)
	{
		temp_row =    temp_buffer[row];
		in_row0   = X[row - 1];
		in_row1   = X[row + 0];
		in_row2   = X[row + 1];
		for(col = ncl + (-1); col < nch + 1 + 1; col++) 
			temp_row[col] = in_row0[col] & in_row1[col] & in_row2[col];
	}
	
	for (row = nrl; row < nrh + 1; row ++)
	{
		temp_row = temp_buffer[row];
		out_row  = Y[row];
		for(col = ncl; col < nch + 1; col++) 
			out_row[col] = scalar_and3(temp_row, col);
	}
	// free_ui8matrix(temp, nrl, nrh, ncl + (-1), nch + 1);
}
void ui8matrix_erosion_divide_col_and_conquer_ExLU_O3 (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **temp_buffer, uint8 **Y)
{
	const long order = 3;
	// uint8 **temp = ui8matrix(nrl, nrh, ncl + (-1), nch + 1);
	uint8 *temp_row, *temp_row0, *temp_row1, *temp_row2, *temp_row3, *out_row0,*out_row1,*out_row2, *in_row0, *in_row1, *in_row2;
	long row, col, nrow, r;
	int nb_threads;

	nrow = nrh - nrl + 1;
	row = 0;
	for (row = nrl; row < nrh + 1; row ++)
	{
		temp_row =    temp_buffer[row];
		in_row0   = X[row - 1];
		in_row1   = X[row + 0];
		in_row2   = X[row + 1];
		for(col = ncl + (-1); col < nch + 1 + 1; col++) 
			temp_row[col] = in_row0[col] & in_row1[col] & in_row2[col];
	}
	
	r = (nrh + 1) % order;
	for (row = nrl; row < nrh + 1 - r; row += order)
	{
		temp_row0 =     temp_buffer[row + 0];
		temp_row1 =     temp_buffer[row + 1];
		temp_row2 =     temp_buffer[row + 2];
		out_row0   = Y[row + 0];
		out_row1   = Y[row + 1];
		out_row2   = Y[row + 2];

		for(col = ncl; col < nch + 1; col++) {
			out_row0[col] = scalar_and3(temp_row0, col);
			out_row1[col] = scalar_and3(temp_row1, col);
			out_row2[col] = scalar_and3(temp_row2, col);
		}
	}

	switch(r){
		case 2:
		temp_row1 =     temp_buffer[nrh - 1];
		temp_row2 =     temp_buffer[nrh + 0];

		out_row0  = Y[nrh - 1];
		out_row1  = Y[nrh - 0];

		for(col = ncl; col < nch + 1; col++) {
			out_row0[col] = scalar_and3(temp_row1, col);
			out_row1[col] = scalar_and3(temp_row2, col);
		}
		break;
		case 1:
		temp_row2 =     temp_buffer[nrh + 0];
		out_row1  = Y[nrh - 0];

		for(col = ncl; col < nch + 1; col++) {
			out_row1[col] = scalar_and3(temp_row2, col);
		}
		break;
		case 0:
		break;
	}
	// free_ui8matrix(temp, nrl, nrh, ncl + (-1), nch + 1);
}
void ui8matrix_erosion_divide_col_and_conquer_InLU_O3 (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **temp_buffer, uint8 **Y)
{
	const long order = 3;
	// uint8 **temp = ui8matrix(nrl, nrh, ncl + (-1), nch + 1);
	uint8 *temp_row, *temp_row0, *temp_row1, *temp_row2, *out_row, *in_row0, *in_row1, *in_row2;
	long row, col, nrow, r;
	int nb_threads;

	// for (row = nrl; row < nrh + 1; row ++)
	nrow = nrh - nrl + 1;
	row = 0;
	for (row = nrl; row < nrh + 1; row ++)
	{
		temp_row =    temp_buffer[row];
		in_row0   = X[row - 1];
		in_row1   = X[row + 0];
		in_row2   = X[row + 1];
		for(col = ncl + (-1); col < nch + 1 + 1; col++) 
			temp_row[col] = in_row0[col] & in_row1[col] & in_row2[col];
	}
	r = (nch + 1) % order;
	for (row = nrl; row < nrh + 1; row ++)
	{
		temp_row = temp_buffer[row];
		out_row  = Y[row];
		for(col = ncl; col < nch + 1 - r; col+= order) {
			out_row[col + 0] = scalar_and3(temp_row, col + 0);
			out_row[col + 1] = scalar_and3(temp_row, col + 1);
			out_row[col + 2] = scalar_and3(temp_row, col + 2);
		}
	}

	switch(r){
		case 2:
		for (row = nrl; row < nrh + 1; row ++) {
			temp_row = temp_buffer[row];
			out_row  = Y[row];
			out_row[nch - 1] = scalar_and3(temp_row, nch - 1);
			out_row[nch + 0] = scalar_and3(temp_row, nch + 0);
		}
		break;
		case 1:
		for (row = nrl; row < nrh + 1; row ++) {
			temp_row = temp_buffer[row];
			Y[row][nch + 0] = scalar_and3(temp_row, nch + 0);
		}
		break;
		case 0:
		break;
	}
	// display_ui8matrix(temp, nrl, nrh, ncl + (-1), nch + 1, "%u", "temp");
	// free_ui8matrix(temp, nrl, nrh, ncl + (-1), nch + 1);
}
void ui8matrix_erosion_pipeline2_LU3x3_ExLU_O3 (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **temp_buffer, uint8 **Y)
{
	const long order = 3;
	long row = nrl, col = ncl, x, y, r = 0;
	// uint8 **temp = ui8matrix(nrl + (-1), nrh + 1, ncl + (-1), nch + 1);
	uint8 *row0, *row1, *row2, *row3, x0, x1, x2, x3;
	uint8 *temp_row0, *temp_row1, *temp_row2, *temp_row3, *temp_row4;
	uint8 *out_row0, *out_row1, *out_row2, *out_row3;
	
	r = (nrh + 1) % order;
	// Prologue	
	row0 = X[row - 1];
	row1 = X[row + 0];
	row2 = X[row + 1];
	
	temp_row1 = temp_buffer[row - 1];
	temp_row2 = temp_buffer[row + 0];
	temp_row3 = temp_buffer[row + 1];
	for (col = ncl; col < nch + 1; col++) {
		temp_row1[col] = scalar_and3(row0, col);
		temp_row2[col] = scalar_and3(row1, col);
	}

	for (row = row + 1; row < nrh + 1 - r; row += order){
		row0 = X[row - 1];   //  0 : b0 b1 b2 ...
		row1 = X[row + 0];   //  1 : c0 c1 c2 ...
		row2 = X[row + 1];   //  2 : d0 d1 d2 ...
		row3 = X[row + 2];   //  3 : e0 e1 e2 ...

		temp_row0 = temp_buffer[row - 2]; // -1 : A    <= (a0 a1 a2 ...) Prologued
		temp_row1 = temp_buffer[row - 1]; //  0 : B    <= (b0 b1 b2 ...) Prologued
		temp_row2 = temp_buffer[row + 0]; //  1 : temp_buffer <= (c0 c1 c2 ...)
		temp_row3 = temp_buffer[row + 1]; //  2 : temp_buffer <= (d0 d1 d2 ...)
		temp_row4 = temp_buffer[row + 2]; //  3 : temp_buffer <= (e0 e1 e2 ...)
		out_row0 = Y[row - 1]; 
		out_row1 = Y[row + 0];
		out_row2 = Y[row + 1];

		for (col = ncl; col < nch + 1; col++) {
			//  temp_row0[ col]			   LOAD => A										
			x0 = temp_row1[col];		// LOAD => B
			x1 = scalar_and3(row1, col); // OR3  => C
			x2 = scalar_and3(row2, col); // OR3  => D
			x3 = scalar_and3(row3, col); // OR3  => E
			out_row0[col] = temp_row0[ col] & x0 & x1; // A & B & C
			out_row1[col] = 			 x0 & x1 & x2; // B & C & D
			out_row2[col] = 			 x1 & x2 & x3; // C & D & E
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
			row1 = X[nrh - 1];   // -1 : (f0 f1 f2 ...)
			row2 = X[nrh + 0];   //  0 : (g0 g1 g2 ...)
			row3 = X[nrh + 1];   //  1 : (h0 h1 h2 ...)
			temp_row0 = temp_buffer[nrh - 2]; //  E    <= (e0 e1 e2 ...)
			out_row0 = Y[nrh - 1];
			out_row1 = Y[nrh + 0];
			for (col = ncl; col < nch + 1; col++) {
				//  temp_row0[ col]		   		 LOAD => E		
				x0 = scalar_and3(row1, col);	  // OR3  => F	
				x1 = scalar_and3(row2, col);	  // OR3  => G
				out_row0[col] = temp_row0[ col] & x0 & x1; // A & B & C
				out_row1[col] = 			 x0 & x1 & scalar_and3(row3, col); // B & C & D	
			}
			break;
		case 1: 
			row2 = X[nrh + 0];   //  0 : (f0 f1 f2 ...)
			row3 = X[nrh + 1];   //  1 : (g0 g1 g2 ...)
			temp_row0 = temp_buffer[nrh - 1]; //  E    <= (e0 e1 e2 ...)
			out_row1 = Y[nrh + 0];
			for (col = ncl; col < nch + 1; col++) {
				out_row1[col] = temp_row0[ col]  & scalar_and3(row2, col) & scalar_and3(row3, col); // E & F & G	
			}
			break;
		default:
			break;
	}
	// free_ui8matrix(temp, nrl + (-1), nrh + 1, ncl + (-1), nch + 1);
}
void ui8matrix_erosion_pipeline2_LU3x3_InLU_O3 (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **temp_buffer, uint8 **Y)
{
	const long order = 3;
	long row = nrl, col = ncl, x, y, r = 0;
	// uint8 **temp = ui8matrix(nrl + (-1), nrh + 1, ncl + (-1), nch + 1);
	uint8 *row0, *row1, *row2, *row3, x0, x1, x2, x3, x4, x5, y0, y1, y2;
	uint8 *temp_row0, *temp_row1, *temp_row2, *temp_row3, *temp_row4;
	uint8 *out_row0, *out_row1, *out_row2, *out_row3;
	
	r = (nch + 1) % order;
	// Prologue	
	row0 = X[row - 1];
	row1 = X[row + 0];
	row2 = X[row + 1];
	
	// temp_row1 = temp_buffer[row - 1];
	// temp_row2 = temp_buffer[row + 0];
	// temp_row3 = temp_buffer[row + 1];
	for (row = nrl; row < nrh + 1; row++) {
		row0 = X[row - 1];
		row1 = X[row + 0];
		row2 = X[row + 1];
		temp_row1 = temp_buffer[row + 0];

		temp_row1[ncl - 1] = row0[ncl - 1] & row1[ncl - 1] & row2[ncl - 1];
		temp_row1[ncl + 0] = row0[ncl + 0] & row1[ncl + 0] & row2[ncl + 0];
	}

	for (row = nrl; row < nrh + 1; row ++){
		temp_row1 = temp_buffer[row + 0];
		
		row0 = X[row - 1];
		row1 = X[row + 0];
		row2 = X[row + 1];

		out_row1 = Y[row]; 
		
		for (col = ncl; col < nch + 1 - r; col += order) {
			x0 = temp_row1[col - 1];		// LOAD => A
			x1 = temp_row1[col - 0];		// LOAD => B
			x2 = row0[col + 1] & row1[col + 1] & row2[col + 1];
			x3 = row0[col + 2] & row1[col + 2] & row2[col + 2];
			x4 = row0[col + 3] & row1[col + 3] & row2[col + 3];
			// x5 = row0[col + 3] & row1[col + 3] & row2[col + 3];
			y0 = x0 & x1 & x2;
			y1 = x1 & x2 & x3;
			y2 = x2 & x3 & x4;
			
			out_row1 [col + 0] = y0; // A & B & C
			out_row1 [col + 1] = y1; // B & C & D
			out_row1 [col + 2] = y2; // C & D & E

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
				temp_row1 = temp_buffer[row + 0];
				
				row0 = X[row - 1];
				row1 = X[row + 0];
				row2 = X[row + 1];

				out_row0 = Y[row]; 
				
				x0 = temp_row1[nch - 2];		// LOAD => E
				x1 = row0[nch - 1] & row1[nch - 1] & row2[nch - 1];
				x2 = row0[nch + 0] & row1[nch + 0] & row2[nch + 0];
				x3 = row0[nch + 1] & row1[nch + 1] & row2[nch + 1];
				
				y0 = x0 & x1 & x2;
				y1 = x1 & x2 & x3;
				// y2 = x2 & x3 & x4;
				out_row0[nch - 1] = y0; 	// A & B & C
				out_row0[nch + 0] = y1; 	// B & C & D
				// out_row0[nch + 1] = y2; 	// C & D & E				
			}
			// ui8matrix_erosion_pipeline3x3(X, nrl, nrh, nch - 1, nch, temp_buffer, Y);
			break;
		case 1: 
			for (row = nrl; row < nrh + 1; row ++){
				temp_row1 = temp_buffer[row + 0];
				
				row0 = X[row - 1];
				row1 = X[row + 0];
				row2 = X[row + 1];

				out_row0 = Y[row]; 
				
				x0 = temp_row1[nch - 1];		// LOAD => E
				// x1 = row0[nch - 1] & row1[nch - 1] & row2[nch - 1];
				x2 = row0[nch + 0] & row1[nch + 0] & row2[nch + 0];
				x3 = row0[nch + 1] & row1[nch + 1] & row2[nch + 1];
				
				// y0 = x0 & x1 & x2;
				y1 = x0 & x2 & x3;
				// y2 = x2 & x3 & x4;
				// out_row0[nch - 1] = y0; 	// A & B & C
				out_row0[nch + 0] = y1; 	// B & C & D
				// out_row0[nch + 1] = y2; 	// C & D & E				
			}
		
			break;
		default:
			break;
	}
	// display_ui8matrix(Y, nrl - 1, nrh + 1, ncl - 1, nch + 1, "%u", "TEst");
	// free_ui8matrix(temp, nrl + (-1), nrh + 1, ncl + (-1), nch + 1);
	
}


void ui8matrix_erosion_pipeline2_LU3x3_ExLU_O3_RR (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **temp_buffer, uint8 **Y)
{
	const long order = 3;
	long row = nrl, col = ncl, r = 0;
	// uint8 **temp = ui8matrix(nrl + (-1), nrh + 1, ncl + (-1), nch + 1);
	uint8 *row0, *row1, *row2, *row3, x0, x1, x2, x3;
	uint8 y0, y1, y2, y3, y4, y5, y6, y7, y8, y9, y10, y11, y12, y13, y14;
	uint8 *temp_row0, *temp_row1, *temp_row2, *temp_row3, *temp_row4;
	uint8 *out_row0, *out_row1, *out_row2, *out_row3;
	
	r = (nrh + 1) % order;
	// Prologue	
	row0 = X[row - 1];
	row1 = X[row + 0];
	row2 = X[row + 1];
	
	temp_row1 = temp_buffer[row - 1];
	temp_row2 = temp_buffer[row + 0];
	temp_row3 = temp_buffer[row + 1];
	y0 = row0[ncl-1]; y1 = row0[ncl+0];
	y3 = row1[ncl-1]; y4 = row1[ncl+0];
	for (col = ncl; col < nch + 1; col++) {
		y2 = row0[col + 1];
		y5 = row1[col + 1];
		x1 = y0&y1&y2;//scalar_and3(row1, col); // OR3  => C
		x2 = y3&y4&y5;//scalar_and3(row2, col); // OR3  => D
		temp_row1[col] = y0&y1&y2;//scalar_and3(row0, col);
		temp_row2[col] = y3&y4&y5;//scalar_and3(row1, col);
		y0 = y1; y1 = y2;
		y3 = y4; y4 = y5;
			
	}

	// row0 = X[row + 0];   //  0 : b0 b1 b2 ...
	temp_row0 = temp_buffer[row - 1]; // -1 : A    <= (a0 a1 a2 ...) Prologued
	temp_row1 = temp_buffer[row - 0]; //  0 : B    <= (b0 b1 b2 ...) Prologued
	for (row = row + 1; row < nrh + 1 - r; row += order){
		row1 = X[row + 0];   //  1 : c0 c1 c2 ...
		row2 = X[row + 1];   //  2 : d0 d1 d2 ...
		row3 = X[row + 2];   //  3 : e0 e1 e2 ...

		temp_row2 = temp_buffer[row + 0]; //  1 : temp_buffer <= (c0 c1 c2 ...)
		temp_row3 = temp_buffer[row + 1]; //  2 : temp_buffer <= (d0 d1 d2 ...)
		temp_row4 = temp_buffer[row + 2]; //  3 : temp_buffer <= (e0 e1 e2 ...)
		out_row0 = Y[row - 1]; 
		out_row1 = Y[row + 0];
		out_row2 = Y[row + 1];

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

			x1 = y0&y1&y2;//scalar_and3(row1, col); // OR3  => C
			x2 = y3&y4&y5;//scalar_and3(row2, col); // OR3  => D
			x3 = y6&y7&y8;//scalar_and3(row3, col); // OR3  => E
			out_row0[col] = temp_row0[ col] & x0 & x1; // A & B & C
			out_row1[col] = 			 x0 & x1 & x2; // B & C & D
			out_row2[col] = 			 x1 & x2 & x3; // C & D & E
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
			row1 = X[nrh - 1];   // -1 : (f0 f1 f2 ...)
			row2 = X[nrh + 0];   //  0 : (g0 g1 g2 ...)
			row3 = X[nrh + 1];   //  1 : (h0 h1 h2 ...)
			temp_row0 = temp_buffer[nrh - 2]; //  E    <= (e0 e1 e2 ...)
			out_row0 = Y[nrh - 1];
			out_row1 = Y[nrh + 0];
			y0 = row1[ncl - 1]; y1 = row1[ncl + 0];
			y3 = row2[ncl - 1]; y4 = row2[ncl + 0];
			y6 = row3[ncl - 1]; y7 = row3[ncl + 0];
			for (col = ncl; col < nch + 1; col++) {
				//  temp_row0[ col]		   		 LOAD => E	
				y2 = row1[col + 1];
				y5 = row2[col + 1];
				y8 = row3[col + 1];

				x0 = y0&y1&y2;//scalar_and3(row1, col); // OR3  => C
				x1 = y3&y4&y5;//scalar_and3(row2, col); // OR3  => D

				out_row0[col] = temp_row0[ col] & x0 & x1; // A & B & C
				out_row1[col] = 			 x0 & x1 & y6&y7&y8; // B & C & D	
				y0 = y1; y1 = y2;
				y3 = y4; y4 = y5;
				y6 = y7; y7 = y8;
			}
			break;
		case 1: 
			row2 = X[nrh + 0];   //  0 : (f0 f1 f2 ...)
			row3 = X[nrh + 1];   //  1 : (g0 g1 g2 ...)
			temp_row0 = temp_buffer[nrh - 1]; //  E    <= (e0 e1 e2 ...)
			out_row1 = Y[nrh + 0];
			y3 = row2[ncl - 1]; y4 = row2[ncl + 0];
			y6 = row3[ncl - 1]; y7 = row3[ncl + 0];
			for (col = ncl; col < nch + 1; col++) {
				y5 = row2[col + 1];
				y8 = row3[col + 1];
				out_row1[col] = temp_row0[ col] & y3 & y4 & y5 & y6 & y7 & y8; // E & F & G	
				y3 = y4; y4 = y5;
				y6 = y7; y7 = y8;
			}

			
			break;
		default:
			break;
	}
	// free_ui8matrix(temp, nrl + (-1), nrh + 1, ncl + (-1), nch + 1);
}
void ui8matrix_erosion_pipeline2_LU3x3_InLU_O3_RR (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **temp_buffer, uint8 **Y)
{
	const long order = 3;
	long row = nrl, col = ncl, x, y, r = 0;
	// uint8 **temp = ui8matrix(nrl + (-1), nrh + 1, ncl + (-1), nch + 1);
	uint8 *row0, *row1, *row2, *row3, x0, x1, x2, x3, x4, x5, y0, y1, y2;
	uint8 *temp_row0, *temp_row1, *temp_row2, *temp_row3, *temp_row4;
	uint8 *out_row0, *out_row1, *out_row2, *out_row3;
	
	r = (nch + 1) % order;
	// Prologue	
	row0 = X[row - 1];
	row1 = X[row + 0];
	row2 = X[row + 1];
	
	// temp_row1 = temp_buffer[row - 1];
	// temp_row2 = temp_buffer[row + 0];
	// temp_row3 = temp_buffer[row + 1];
	for (row = nrl; row < nrh + 1; row++) {
		row0 = X[row - 1];
		row1 = X[row + 0];
		row2 = X[row + 1];
		temp_row1 = temp_buffer[row + 0];

		temp_row1[ncl - 1] = row0[ncl - 1] & row1[ncl - 1] & row2[ncl - 1];
		temp_row1[ncl + 0] = row0[ncl + 0] & row1[ncl + 0] & row2[ncl + 0];
	}

	for (row = nrl; row < nrh + 1; row ++){
		temp_row1 = temp_buffer[row + 0];
		
		row0 = X[row - 1];
		row1 = X[row + 0];
		row2 = X[row + 1];

		out_row1 = Y[row]; 
		
		x0 = temp_row1[ncl - 1];		// LOAD => A
		x1 = temp_row1[ncl - 0];		// LOAD => B
		for (col = ncl; col < nch + 1 - r; col += order) {
			x2 = row0[col + 1] & row1[col + 1] & row2[col + 1];
			x3 = row0[col + 2] & row1[col + 2] & row2[col + 2];
			x4 = row0[col + 3] & row1[col + 3] & row2[col + 3];
			
			out_row1 [col + 0] = x0 & x1 & x2; // A & B & C
			out_row1 [col + 1] = x1 & x2 & x3; // B & C & D
			out_row1 [col + 2] = x2 & x3 & x4; // C & D & E

			x0 = x3; x1 = x4;
		}
	}
	switch (r) {
		case 2: 
			row0 = X[nrl - 1];
			row1 = X[nrl + 0];
			for (row = nrl; row < nrh + 1; row++) {
				row2 = X[row + 1];
				out_row1 = Y[row]; 

				x1 = row0[nch - 1] & row1[nch - 1] & row2[nch - 1]&
					 row0[nch + 0] & row1[nch + 0] & row2[nch + 0];

				out_row1 [nch - 1] = row0[nch - 2] & row1[nch - 2] & row2[nch - 2] & x1 ; // A & B & C
				out_row1 [nch + 0] = x1 & row0[nch + 1] & row1[nch + 1] & row2[nch + 1]; // B & C & D
				row0 = row1;
				row1 = row2;
			}
			break;
		case 1: 
			row0 = X[nrl - 1];
			row1 = X[nrl + 0];
			for (row = nrl; row < nrh + 1; row++) {
				row2 = X[row + 1];
				Y[row][nch] = scalar_and3(row0, nch) &
									 scalar_and3(row1, nch) &
									 scalar_and3(row2, nch);
				row0 = row1;
				row1 = row2;
			}
			break;
		default:
			break;
	}
	// display_ui8matrix(Y, nrl - 1, nrh + 1, ncl - 1, nch + 1, "%u", "TEst");
	// free_ui8matrix(temp, nrl + (-1), nrh + 1, ncl + (-1), nch + 1);
	
}
/***********************************************************************************/
/******* End of optimisation : Loop Unroll + Register Rotation + Pipelining ********/
/***********************************************************************************/

void ui8matrix_erosion_LU3x3_ExLU_O3_NS (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **temp_buffer, uint8 **Y)
{
	const long order = 3;
	long row = nrl, col = ncl, x, y;
	uint8 y0, y1, y2, x0, x1, x2, x3, r;
	uint8 *row0, *row1, *row2, *row3, *row4;

	r = (nrh + 1)  % order;
	
	for (row = nrl; row < nrh + 1 - r; row += order) {
		for (col = ncl; col < nch + 1; col++) {
			x0 = scalar_and3(X[row + 0], col);
			x1 = scalar_and3(X[row + 1], col);
			x2 = scalar_and3(X[row + 2], col);

			Y[row + 0][col] = scalar_and3(X[row - 1], col) & x0 & x1;
			Y[row + 1][col] = 				     			x0 & x1 & x2;
			Y[row + 2][col] = 				     			x1 & x2 & scalar_and3(X[row + 3], col);
		}
	}
	
	switch(r) {
		case 2:
			for (col = ncl; col < nch + 1; col++) {
				x0 = scalar_and3(X[row + 0], col);
				x1 = scalar_and3(X[row + 1], col);
				Y[row + 0][col] = scalar_and3(X[row - 1], col) & x0 & x1;
				Y[row + 1][col] = 				     		    x0 & x1 & scalar_and3(X[row + 2], col);
			}
			break;
		case 1:
			for (col = ncl; col < nch + 1; col++) {
				Y[row + 0][col] = scalar_and3(X[row - 1], col) & 
										 scalar_and3(X[row + 0], col) &
										 scalar_and3(X[row + 1], col);
				
			}
			break;
		default:
			break;
	}
}
void ui8matrix_erosion_LU3x3_InLU_O3_NS (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **temp_buffer, uint8 **Y)
{
	const long order = 3;
	long row = nrl, col = ncl, x, y, r;
	uint8 y0, y1, y2, x0, x1, x2, x3, x4, x5, x6, x7, x8;
	uint8 *row0, *row1, *row2, *row3, *row4;

	r = (nch + 1)  % order;
	
	for (row = nrl; row < nrh + 1; row ++) {
		row0 = X[row - 1];
		row1 = X[row + 0];
		row2 = X[row + 1];
		for (col = ncl; col < nch + 1 - r; col += order) {
			// x0 = scalar_and3(row0, col); x3 = scalar_and3(row0, col + 1); x6 = scalar_and3(row0, col + 2);
			// x1 = scalar_and3(row1, col); x4 = scalar_and3(row1, col + 1); x7 = scalar_and3(row1, col + 2);
			// x2 = scalar_and3(row2, col); x5 = scalar_and3(row2, col + 1); x8 = scalar_and3(row2, col + 2);

			Y[row + 0][col + 0] = scalar_and3(row0, col + 0)&
										 scalar_and3(row1, col + 0)&
										 scalar_and3(row2, col + 0);
			Y[row + 0][col + 1] = scalar_and3(row0, col + 1)&
 										 scalar_and3(row1, col + 1)&
 										 scalar_and3(row2, col + 1);
			Y[row + 0][col + 2] = scalar_and3(row0, col + 2)&
										 scalar_and3(row1, col + 2)&
										 scalar_and3(row2, col + 2);
		}
	}

	switch(r) {
		case 2:
			for (row = nrl; row < nrh + 1; row ++) {
				row0 = X[row - 1];
				row1 = X[row + 0];
				row2 = X[row + 1];

				Y[row + 0][col + 0] = scalar_and3(row0, col + 0)&
 											 scalar_and3(row1, col + 0)&
 											 scalar_and3(row2, col + 0);
				Y[row + 0][col + 1] = scalar_and3(row0, col + 1)&
											 scalar_and3(row1, col + 1)&
											 scalar_and3(row2, col + 1);
			}
			break;
		case 1:
			for (row = nrl; row < nrh + 1; row ++) {
				Y[row + 0][col + 0] = scalar_and3(X[row - 1], col) & 
											 scalar_and3(X[row + 0], col) & 
											 scalar_and3(X[row + 1], col);
			}
			break;
		default:
			break;
	}
}
void ui8matrix_erosion_LU3x3_ComLU_O3_NS (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **temp_buffer, uint8 **Y)
{
	const long order = 3;
	long row = nrl, col = ncl, x, y;
	uint8 y0, y1, y2, rr, cr;
	uint8 x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14;
	uint8 *row0, *row1, *row2, *row3, *row4;

	rr = (nrh + 1) % order;
	cr = (nch + 1) % order;
	
	for (row = nrl; row < nrh + 1 - rr; row += order) {
		for (col = ncl; col < nch + 1 - cr; col += order) {
			x0 = scalar_and3(X[row + 0], col); x1 = scalar_and3(X[row + 0], col + 1); x2 = scalar_and3(X[row + 0], col + 2);
			x3 = scalar_and3(X[row + 1], col); x4 = scalar_and3(X[row + 1], col + 1); x5 = scalar_and3(X[row + 1], col + 2);
			x6 = scalar_and3(X[row + 2], col); x7 = scalar_and3(X[row + 2], col + 1); x8 = scalar_and3(X[row + 2], col + 2);

			Y[row + 0][col + 0] = scalar_and3(X[row - 1], col + 0) & x0 & x3;
			Y[row + 1][col + 0] = 				         			x0 & x3 & x6;
			Y[row + 2][col + 0] = 				         			x3 & x6 & scalar_and3(X[row + 3], col + 0);


			Y[row + 0][col + 1] = scalar_and3(X[row - 1], col + 1) & x1 & x4;
			Y[row + 1][col + 1] =			  			            x1 & x4 & x7;
			Y[row + 2][col + 1] =			  			            x4 & x7 & scalar_and3(X[row + 3], col + 1);


			Y[row + 0][col + 2] = scalar_and3(X[row - 1], col + 2) & x2 & x5;
			Y[row + 1][col + 2] = 				                    x2 & x5 & x8;
			Y[row + 2][col + 2] = 				                    x5 & x8 & scalar_and3(X[row + 3], col + 2);

		}
	}
	switch (rr) {
		case 2 :
			switch(cr){
				case 2 :
					ui8matrix_erosion_LU3x3_ExLU_O3_NS(X, nrl, nrh, nch - 1, nch, temp_buffer, Y);
					ui8matrix_erosion_LU3x3_InLU_O3_NS(X, nrh - 1, nrh, ncl, nch, temp_buffer, Y);
				break;
				case 1: 
					ui8matrix_erosion_LU3x3_ExLU_O3_NS(X, nrl, nrh, nch, nch, temp_buffer, Y);		
					ui8matrix_erosion_LU3x3_InLU_O3_NS(X, nrh - 1, nrh, ncl, nch, temp_buffer, Y);	
				break;
				default :
					ui8matrix_erosion_LU3x3_InLU_O3_NS(X, nrh - 1, nrh, ncl, nch, temp_buffer, Y);
				break;
			}
		break;
		case 1:
			switch(cr){
				case 2 :
					ui8matrix_erosion_LU3x3_ExLU_O3_NS(X, nrl, nrh, nch - 1, nch, temp_buffer, Y);
					ui8matrix_erosion_LU3x3_InLU_O3_NS(X, nrh, nrh, ncl, nch, temp_buffer, Y);
				break;
				case 1: 
					ui8matrix_erosion_LU3x3_ExLU_O3_NS(X, nrl, nrh, nch, nch, temp_buffer, Y);		
					ui8matrix_erosion_LU3x3_InLU_O3_NS(X, nrh, nrh, ncl, nch, temp_buffer, Y);	
				break;
				default :
					ui8matrix_erosion_LU3x3_InLU_O3_NS(X, nrh, nrh, ncl, nch, temp_buffer, Y);
				break;
			}
		break;
		default:
			switch(cr){
				case 2 :
					ui8matrix_erosion_LU3x3_ExLU_O3_NS(X, nrl, nrh, nch - 1, nch, temp_buffer, Y);
				break;
				case 1: 
					ui8matrix_erosion_LU3x3_ExLU_O3_NS(X, nrl, nrh, nch, nch, temp_buffer, Y);	
				break;
				default :
				break;
			}
		break;
	}
	// ui8matrix_erosion_LU3x3(X, row, nrh, col, nrh, Y);
}
void ui8matrix_erosion_LU3x3_ExLU_O3_RR_NS(uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **temp_buffer, uint8 **Y)
{
	const long order = 3;
	long row = nrl, col = ncl, x, y;
	uint8 y0, y1, y2, r;
	uint8 x0, x1, x2, x3, x4, x5, x6, x7, x8;
	uint8 *row0, *row1, *row2, *row3, *row4;

	r = (nrh + 1)  % order;
	
	row0 = X[row - 1];
	row1 = X[row + 0];
	for (row = nrl; row < nrh + 1 - r; row += order) {
		row2 = X[row + 1]; 
		row3 = X[row + 2]; 
		row4 = X[row + 3];
		x0 = row1[-1]; x3 = row1[ 0];
		x1 = row2[-1]; x4 = row2[ 0];
		x2 = row3[-1]; x5 = row3[ 0];
		for (col = ncl; col < nch + 1; col++) {
			x6 = row1[col + 1]; //scalar_and3(row1, col);
			x7 = row2[col + 1]; //scalar_and3(row2, col);
			x8 = row3[col + 1]; //scalar_and3(row3, col);
			y0 = (x0 & x3 & x6); // 2 calc
			y1 = (x1 & x4 & x7); // 2 calc
			y2 = (x2 & x5 & x8); // 2 calc

			Y[row + 0][col] = scalar_and3(row0, col) & y0 & y1; 						// 3 + 2    =  5
			Y[row + 1][col] =         			 y0 & y1 & y2;						//     2    =  2
			Y[row + 2][col] = 		 			 y1 & y2 & scalar_and3(row4, col);  //     2 + 3=  5
			// 6 + 12 = 18 calc
			x0 = x3; x3 = x6;
			x1 = x4; x4 = x7;
			x2 = x5; x5 = x8;
		}
		row0 = row3;
		row1 = row4;
	}
	row2 = X[row + 1];
	row3 = X[row + 2];
	row4 = X[row + 3];
	switch(r) {
		case 2:
			x0 = row1[-1]; x3 = row1[ 0];
			x1 = row2[-1]; x4 = row2[ 0];
			x2 = row3[-1]; x5 = row3[ 0];
			for (col = ncl; col < nch + 1; col++) {
				x6 = row1[col + 1]; //scalar_and3(row1, col);
				x7 = row2[col + 1]; //scalar_and3(row2, col);
				x8 = row3[col + 1]; //scalar_and3(row3, col);
				Y[row + 0][col] = scalar_and3(row0, col) & (x0 & x3 & x6) & (x1 & x4 & x7);;
				Y[row + 1][col] =         (x0 & x3 & x6) & (x1 & x4 & x7) & (x2 & x5 & x8);;
				x0 = x3; x3 = x6;
				x1 = x4; x4 = x7;
				x2 = x5; x5 = x8;
			}
			break;
		case 1:
			x0 = row1[-1]; x3 = row1[ 0];
			x1 = row2[-1]; x4 = row2[ 0];
			for (col = ncl; col < nch + 1; col++) {
				x6 = row1[col + 1]; //scalar_and3(row1, col);
				x7 = row2[col + 1]; //scalar_and3(row2, col);
				Y[row + 0][col] = scalar_and3(row0, col) & (x0 & x3 & x6) & (x1 & x4 & x7);
				x0 = x3; x3 = x6;
				x1 = x4; x4 = x7;
			}				
			break;
		default:
			break;
	}
}
void ui8matrix_erosion_LU3x3_InLU_O3_RR_NS (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **temp_buffer, uint8 **Y)
{
	const long order = 3;
	long row = nrl, col = ncl, x, y, r;
	uint8 y0, y1, y2, y3, y4, y5;
	uint8 x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16, x17, x18, x19;
	uint8 *row0, *row1, *row2, *row3, *row4;

	r = (nch + 1) % order;
	row0 = X[row - 1];
	row1 = X[row + 0];
	for (row = nrl; row < nrh + 1; row++) {
		row2 = X[row + 1];

		x0  = row0[ncl - 1]; x1  = row0[ncl + 0]; 
		x5  = row1[ncl - 1]; x6  = row1[ncl + 0]; 
		x10 = row2[ncl - 1]; x11 = row2[ncl + 0]; 
		for (col = ncl; col < nch + 1 - r; col += order) {
			x2  = row0[col + 1]; x3  = row0[col + 2]; x4  = row0[col + 3];
			x7  = row1[col + 1]; x8  = row1[col + 2]; x9  = row1[col + 3];
			x12 = row2[col + 1]; x13 = row2[col + 2]; x14 = row2[col + 3];
			
			y1 = x1 & x2 & x6 & x7 & x11 & x12; // 5 calcs
			y2 = x3 & x8 & x13;					// 2 calcs
			Y[row][col + 0] = x0 & x5 & y1 & x10; //3 calcs
			Y[row][col + 1] =      y2 & y1; 		 //1 calc
			Y[row][col + 2] =      y2 & x2 & x4          & x7 & x9             & x12 & x14; //6 calcs
			// 5 + 2 + 3 + 1 + 6 = 17						

			x0  = x3;  x1  = x4;
			x5  = x8;  x6  = x9;
			x10 = x13; x11 = x14;
			
		}
		row0 = row1;
		row1 = row2;
	}

		
	row = nrl;
	row0 = X[nrl - 1];
	row1 = X[nrl + 0];
	
	switch(r) {
	case 2:
		for (row; row < nrh + 1; row++) {
			row2 = X[row + 1];
			x0  = row0[col - 1]; x1  = row0[col +  0]; x2  = row0[col + 1]; x3  = row0[col + 2];
			x5  = row1[col - 1]; x6  = row1[col +  0]; x7  = row1[col + 1]; x8  = row1[col + 2];
			x10 = row2[col - 1]; x11 = row2[col +  0]; x12 = row2[col + 1]; x13 = row2[col + 2];
			y1 = x1 & x2 & x6 & x7 & x11 & x12; // 5 calcs
			y2 = x3 & x8 & x13;					// 2 calcs
			Y[row][col + 0] = x0 & x5 & y1 & x10; //3 calcs
			Y[row][col + 1] =      y2 & y1; 		 //1 calc
			row0 = row1;
			row1 = row2;
		}
		break;
	case 1:
		for (row; row < nrh + 1; row++) {			
			row2 = X[row + 1];
			x0  = row0[col - 1]; x1  = row0[col +  0]; x2  = row0[col + 1];
			x5  = row1[col - 1]; x6  = row1[col +  0]; x7  = row1[col + 1];
			x10 = row2[col - 1]; x11 = row2[col +  0]; x12 = row2[col + 1];
			Y[row][col + 0] = (x0 & x1 & x2 ) & (x5 & x6 & x7) & (x10 & x11 & x12);
			row0 = row1;
			row1 = row2;
		}
		break;
	default:
		break;
	}
}
