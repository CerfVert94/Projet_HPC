
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include "nrdef.h"
#include "nrutil.h"
#include <util.h>
#include <img.h>
#include <omp.h>
//#include <malloc.h>
#include <morpho.h>
// #include <test_morpho.h>

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


void ui8matrix_sequence_divide_row_and_conquer_OMP(uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **temp_buffer, uint8 **Y)
{
	// uint8 **ppPreOutput0, **ppPreOutput1;
	// ppPreOutput0 = ui8matrix(nrl + SE_NRL, nrh + SE_NRH, ncl + SE_NCL, nch + SE_NCH);
	// ppPreOutput1 = ui8matrix(nrl + SE_NRL, nrh + SE_NRH, ncl + SE_NCL, nch + SE_NCH);
	// memset_ui8matrix(ppPreOutput0, 0, nrl + SE_NRL, nrh + SE_NRH, ncl + SE_NCL, nch + SE_NCH);
	// ui8matrix_erosion_divide_row_and_conquer_OMP (X     , nrl, nrh, ncl, nch, ppPreOutput0);
	// memset_ui8matrix(ppPreOutput1, 0, nrl + SE_NRL, nrh + SE_NRH, ncl + SE_NCL, nch + SE_NCH);
	// ui8matrix_dilation_divide_row_and_conquer_OMP(ppPreOutput0, nrl, nrh, ncl, nch, ppPreOutput1);
	// memset_ui8matrix(ppPreOutput0, 0, nrl + SE_NRL, nrh + SE_NRH, ncl + SE_NCL, nch + SE_NCH);
	// ui8matrix_dilation_divide_row_and_conquer_OMP(ppPreOutput1, nrl, nrh, ncl, nch, ppPreOutput0);
	// ui8matrix_erosion_divide_row_and_conquer_OMP (ppPreOutput0, nrl, nrh, ncl, nch, temp_buffer, Y);
	// free_ui8matrix(ppPreOutput0, nrl + SE_NRL, nrh + SE_NRH, ncl + SE_NCL, nch + SE_NCH);
	// free_ui8matrix(ppPreOutput1, nrl + SE_NRL, nrh + SE_NRH, ncl + SE_NCL, nch + SE_NCH);

}
/*******************************************/
/******* Optimisation : Loop Unroll ********/
/*******************************************/

void ui8matrix_dilation_LU3x3_O1xO1_OMP(uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **temp_buffer, uint8 **Y)
{
	long row = nrl, col = ncl, x, y;
	// dilate
	omp_set_num_threads(omp_get_max_threads());
	#pragma omp parallel for default(none) shared(nrl, ncl, nrh, nch, Y, X) private(row, col) 
	for (row = nrl; row < nrh + 1; row++)
		for (col = ncl; col < nch + 1; col++)
            Y[row][col] = scalar_or3x3(&X[row], col);
}

void ui8matrix_dilation_LU3x3_InLU_O3_OMP(uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **temp_buffer, uint8 **Y)
{
	const long order = 3;
	long row = nrl, col = ncl, x, y, r;
	uint8 y0, y1, y2, x0, x1, x2, x3, x4, x5, x6, x7, x8;
	uint8 *row0, *row1, *row2, *row3, *row4;
	uint8 *out_row0;
	r = (nch + 1)  % order;
	
	omp_set_num_threads(omp_get_max_threads());

	#pragma omp parallel for default(none) shared(nrl, ncl, nrh, nch, Y, X, r) private(row, col, row0, row1, row2, out_row0) 		
	for (row = nrl; row < nrh + 1; row ++) {
		row0 = X[row - 1];
		row1 = X[row + 0];
		row2 = X[row + 1];
		out_row0 = Y[row + 0];
		for (col = ncl; col < nch + 1 - r; col += order) {
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
			#pragma omp parallel for default(none) shared(nrl, ncl, nrh, nch, Y, X, r) private(row, col, row0, row1, row2, out_row0) 		
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
			#pragma omp parallel for default(none) shared(nrl, ncl, nrh, nch, Y, X, r) private(row, col, row0, row1, row2, out_row0) 		
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

void ui8matrix_dilation_LU3x3_ExLU_O3_OMP(uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **temp_buffer, uint8 **Y)
{
	const long order = 3;
	long row = nrl, col = ncl, x, y;
	uint8 y0, y1, y2, x0, x1, x2, x3, r;
	uint8 *row0, *row1, *row2, *row3, *row4;
	uint8 *out_row0, *out_row1, *out_row2;

	r = (nrh + 1)  % order;
	
	omp_set_num_threads(omp_get_max_threads());
	#pragma omp parallel for default(none) shared(nrl, ncl, nrh, nch, Y, X, r) private(y0, y1, y2, x0, x1, x2, x3,row, col, row0, row1, row2,row3,row4, out_row0,out_row1, out_row2)
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
	// printf("remainder :%ld\n", r);
	switch(r) {
		case 2:
			row0 = X[nrh - 2];
			row1 = X[nrh - 1];
			row2 = X[nrh + 0];
			row3 = X[nrh + 1];
			out_row0 = Y[nrh - 1];
			out_row1 = Y[nrh + 0];
			#pragma omp parallel for default(none) shared(nrl, ncl, nrh, nch, Y, X, r, row0, row1, row2,row3,row4, out_row0,out_row1, out_row2) private(y0, y1, y2, x0, x1, x2, x3,row, col)		
			for (col = ncl; col < nch + 1; col++) {
				x0 = scalar_or3(row1, col);
				x1 = scalar_or3(row2, col);
				out_row0[col] = scalar_or3(row0, col) | x0 | x1;
				out_row1[col] = 				        x0 | x1 | scalar_or3(row3, col);
			}
			break;
		case 1:
			row0 = X[nrh - 1];
			row1 = X[nrh + 0];
			row2 = X[nrh + 1];
			out_row0 = Y[nrh + 0];
			#pragma omp parallel for default(none) shared(nrl, ncl, nrh, nch, Y, X, r, row0, row1, row2,row3,row4, out_row0,out_row1, out_row2) private(y0, y1, y2, x0, x1, x2, x3,row, col)		
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

void ui8matrix_dilation_LU3x3_ComLU_O3_OMP(uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **temp_buffer, uint8 **Y)
{
	const long order = 3;
	long row = nrl, col = ncl, x, y;
	uint8 y0, y1, y2, rr, cr;
	uint8 x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14;
	uint8 *row0, *row1, *row2, *row3, *row4;
	uint8 *out_row0, *out_row1, *out_row2;

	rr = (nrh + 1) % order;
	cr = (nch + 1) % order;
	
	omp_set_num_threads(omp_get_max_threads());
	#pragma omp parallel for default(none) shared(nrl, ncl, nrh, nch, Y, X, rr, cr) private(y0, y1, y2, x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14,row, col, row0, row1, row2,row3,row4, out_row0,out_row1, out_row2)	schedule(dynamic)
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
					ui8matrix_dilation_LU3x3_ExLU_O3_OMP(X, nrl, nrh, nch - 1, nch, temp_buffer, Y);
					ui8matrix_dilation_LU3x3_InLU_O3_OMP(X, nrh - 1, nrh, ncl, nch, temp_buffer, Y);
				break;
				case 1: 
					ui8matrix_dilation_LU3x3_ExLU_O3_OMP(X, nrl, nrh, nch, nch, temp_buffer, Y);		
					ui8matrix_dilation_LU3x3_InLU_O3_OMP(X, nrh - 1, nrh, ncl, nch, temp_buffer, Y);	
				break;
				default :
					ui8matrix_dilation_LU3x3_InLU_O3_OMP(X, nrh - 1, nrh, ncl, nch, temp_buffer, Y);
				break;
			}
			break;
		case 1:
			switch(cr){
				case 2 :
					ui8matrix_dilation_LU3x3_ExLU_O3_OMP(X, nrl, nrh, nch - 1, nch, temp_buffer, Y);
					ui8matrix_dilation_LU3x3_InLU_O3_OMP(X, nrh, nrh, ncl, nch, temp_buffer, Y);
				break;
				case 1: 
					ui8matrix_dilation_LU3x3_ExLU_O3_OMP(X, nrl, nrh, nch, nch, temp_buffer, Y);		
					ui8matrix_dilation_LU3x3_InLU_O3_OMP(X, nrh, nrh, ncl, nch, temp_buffer, Y);	
				break;
				default :
					ui8matrix_dilation_LU3x3_InLU_O3_OMP(X, nrh, nrh, ncl, nch, temp_buffer, Y);
				break;
			}
		break;
		default:
			switch(cr){
				case 2 :
					ui8matrix_dilation_LU3x3_ExLU_O3_OMP(X, nrl, nrh, nch - 1, nch, temp_buffer, Y);
				break;
				case 1: 
					ui8matrix_dilation_LU3x3_ExLU_O3_OMP(X, nrl, nrh, nch, nch, temp_buffer, Y);	
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


/*****************************************************************************/
/******** Optimisation : Loop Unroll + Register Rotation of Addresses ********/
/*****************************************************************************/

void ui8matrix_dilation_LU3x3_InLU_O3_AddrRR_OMP(uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **temp_buffer, uint8 **Y)
{
	const long order = 3;
	long row = nrl, col = ncl, r, nrow;
	uint8 y0, y1, y2, x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14;
	uint8 *row0, *row1, *row2, *row3, *row4;
	uint8 *out_row0;
	int nb_threads;


 
	r = (nch + 1) % order;
	nrow = nrh - nrl + 1;
	
	omp_set_num_threads(omp_get_max_threads());
	#pragma omp parallel //shared(X, ncl, nch, Y, nb_threads, r, nrow)
	{
		row0 = X[nrl - 1];
		row1 = X[nrl + 0];
		#pragma omp parallel for firstprivate(row0, row1) private(out_row0, row2, y0, y1, y2, x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, row, col) 
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
	}

	switch(r) {
		case 2:

		#pragma omp parallel for default(none) firstprivate(nrl, nrh) private(out_row0, row0, row1, row2, row3, x0, x1, x2, x3, x5, x6, x7, x8, x10, x11, x12, x13, row, col) shared(X, ncl, nch, Y)
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
		#pragma omp parallel for default(none) firstprivate(nrl, nrh) private(out_row0, row0, row1, row2, row3, row) shared(X, ncl, nch, Y)
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
void ui8matrix_dilation_LU3x3_ExLU_O3_AddrRR_OMP(uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **temp_buffer, uint8 **Y)
{
	
	const long order = 3;
	long row = nrl, col = ncl, x, y, col0, col1, col2, nrow;
	uint8 y0, y1, y2, y3, y5, r;
	uint8 x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14;
	uint8 *row0, *row1, *row2, *row3, *row4;
	uint8 *out_row0, *out_row1, *out_row2;
	int nb_threads;
	r = (nrh + 1)  % order;
	nrow = nrh - nrl + 1;
	omp_set_num_threads(omp_get_max_threads());
	#pragma omp parallel 
	{
		row0 = X[nrl - 1];
		row1 = X[nrl + 0];
		#pragma omp parallel for firstprivate(row0, row1) private(out_row0,out_row1, out_row2,row2, row3, row4,y0, y1, y2, y3, y5, x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, row, col) 
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
	}
	
	out_row0 = Y[nrh - 1];
	out_row1 = Y[nrh + 0];
	switch(r) {
		case 2:
			row0 = X[nrh - 2];
			row1 = X[nrh - 1];
			row2 = X[nrh + 0];
			row3 = X[nrh + 1];
			#pragma omp parallel for default(none) private(row,x0,x1,x2,x3,x4,x5,x6,x7,x8,col) shared(X, ncl, nch, Y, out_row0,out_row1, row0, row1, row2, row3)
			for (col = ncl; col < nch + 1; col++) {
				x0 = row1[col - 1]; x3 = row1[col + 0]; x6 = row1[col + 1]; //scalar_or3(row1, col);
				x1 = row2[col - 1]; x4 = row2[col + 0]; x7 = row2[col + 1]; //scalar_or3(row2, col);
				x2 = row3[col - 1]; x5 = row3[col + 0]; x8 = row3[col + 1]; //scalar_or3(row3, col);
				out_row0[col] = scalar_or3(row0, col) | (x0 | x3 | x6) | (x1 | x4 | x7);
				out_row1[col] =        (x0 | x3 | x6) | (x1 | x4 | x7) | (x2 | x5 | x8);
			}
			break;
		case 1:
			row0 = X[nrh - 1];
			row1 = X[nrh + 0];
			row2 = X[nrh + 1];
			#pragma omp parallel for default(none) private(row,x0,x1,x2,x3,x4,x5,x6,x7,x8,col) shared(X, ncl, nch, Y, out_row0,out_row1, row0, row1, row2, row3)
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
void ui8matrix_dilation_LU3x3_ComLU_O3_AddrRR_OMP(uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **temp_buffer, uint8 **Y)
{
	const long order = 3;
	long row = nrl, col = ncl, rr, cr, nrow;
	uint8 r0c0, r1c0, r2c0, r0c1,r1c1, r2c1, r0c2, r1c2, r2c2;
	uint8 *row0, *row1, *row2, *row3, *row4;
	uint8 *out_row0, *out_row1, *out_row2;
	int nb_threads;

	rr = (nrh + 1) % order;
	cr = (nch + 1) % order;
	nrow = nrh - nrl + 1;
	omp_set_num_threads(omp_get_max_threads());
	#pragma omp parallel 
	{
		uint8 y0, y1, y2, y3, y4, y5, y6, y7, y8, y9, y10, y11, y12, y13, y14;
		uint8 x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16, x17, x18, x19, x20, x21, x22, x23, x24, x25;
		
		row0 = X[nrl - 1];
		row1 = X[nrl + 0];
		#pragma omp parallel for default(none) firstprivate(row0, row1) private(out_row0, out_row1, out_row2, row2, row3, row4, row, col, y0, y1, y2, y3, y4, y5, y6, y7, y8, y9, y10, y11, y12, y13, y14, x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16, x17, x18, x19, x20, x21, x22, x23, x24, x25) shared(X, nrl, nrh, ncl, nch, Y, nb_threads, cr, rr, nrow)
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
	}
	switch (rr) {
		case 2 :
			switch(cr){
				case 2 :
					ui8matrix_dilation_LU3x3_ExLU_O3_AddrRR_OMP(X, nrl, nrh, nch - 1, nch, temp_buffer, Y);
					ui8matrix_dilation_LU3x3_InLU_O3_AddrRR_OMP(X, nrh - 1, nrh, ncl, nch, temp_buffer, Y);
				break;
				case 1: 
					ui8matrix_dilation_LU3x3_ExLU_O3_AddrRR_OMP(X, nrl, nrh, nch, nch, temp_buffer, Y);		
					ui8matrix_dilation_LU3x3_InLU_O3_AddrRR_OMP(X, nrh - 1, nrh, ncl, nch, temp_buffer, Y);	
				break;
				default :
					ui8matrix_dilation_LU3x3_InLU_O3_AddrRR_OMP(X, nrh - 1, nrh, ncl, nch, temp_buffer, Y);
				break;
			}
			break;
		case 1:
			switch(cr){
				case 2 :
					ui8matrix_dilation_LU3x3_ExLU_O3_AddrRR_OMP(X, nrl, nrh, nch - 1, nch, temp_buffer, Y);
					ui8matrix_dilation_LU3x3_InLU_O3_AddrRR_OMP(X, nrh, nrh, ncl, nch, temp_buffer, Y);
				break;
				case 1: 
					ui8matrix_dilation_LU3x3_ExLU_O3_AddrRR_OMP(X, nrl, nrh, nch, nch, temp_buffer, Y);		
					ui8matrix_dilation_LU3x3_InLU_O3_AddrRR_OMP(X, nrh, nrh, ncl, nch, temp_buffer, Y);	
					// display_ui8matrix(Y, nrl, nrh, ncl, nch, "%u", "rr=1; cr=1");
				break;
				default :
					ui8matrix_dilation_LU3x3_InLU_O3_AddrRR_OMP(X, nrh, nrh, ncl, nch, temp_buffer, Y);
				break;
			}
		break;
		default:
			switch(cr){
				case 2 :
					ui8matrix_dilation_LU3x3_ExLU_O3_AddrRR_OMP(X, nrl, nrh, nch - 1, nch, temp_buffer, Y);
				break;
				case 1: 
					ui8matrix_dilation_LU3x3_ExLU_O3_AddrRR_OMP(X, nrl, nrh, nch, nch, temp_buffer, Y);	
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
/*
void ui8matrix_dilation_row_pipeline_OMP(uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **temp_buffer, uint8 **Y)
{
	const long order = 3;
	// uint8 **temp = ui8matrix(nrl + SE_NRL, nrh + SE_NRH, ncl, nch);
	uint8 *temp_row, *temp_row0, *temp_row1, *temp_row2, *out_row, *in_row;
	long row, col, nrow, r;
	int nb_threads;

	nrow = nrh - nrl + 1;
	#pragma omp parallel default(none) private(row, col, temp_row0, temp_row1) shared(nrow, ncl, nch, nrl, nrh, X, Y, temp) 
	{
		uint8 *temp_row, *in_row, *out_row;
		row = 0;
		// Preprocess the epilogue (last border)
		#pragma omp master 
		{
			temp_row = temp_buffer[nrh + 1];
			in_row = X[nrh + 1];
			for(col = ncl; col < nch + 1; col++) {
				temp_row[col] = scalar_or3(in_row, col);
			}
		}

		#pragma omp for nowait private(row, col, temp_row, in_row)
		for (row = nrl - 1; row < nrh + 1; row += order)
		{
			temp_row =    temp_buffer[row];
			in_row   = X[row];
			for(col = ncl; col < nch + 1; col++) {
				temp_row[col] = scalar_or3(in_row, col);
			}
		}
		row = 0;
		#pragma omp for nowait private(row, col,temp_row, in_row)
		for (row = nrl + 0; row < nrh + 1; row += order)
		{

			temp_row =    temp_buffer[row];
			in_row   = X[row];
			for(col = ncl; col < nch + 1; col++) {
				temp_row[col] = scalar_or3(in_row, col);
			}
		}
		row = 0;
		#pragma omp for nowait private(row, col, temp_row, in_row)
		for (row = nrl + 1; row < (nrh + 1) + 1; row += order)
		{
			temp_row =    temp_buffer[row];
			in_row   = X[row];
			for(col = ncl; col < nch + 1; col++) {
				temp_row[col] = scalar_or3(in_row, col);
			}
		}
	}
	
	for (row = nrl; row < nrh + 1; row ++)
	{
		temp_row0 =     temp_buffer[row - 1];
		temp_row1 =     temp_buffer[row + 0];
		temp_row2 =     temp_buffer[row + 1];
		out_row   = Y[row];
		for(col = nrl; col < nch + 1; col++) {
			out_row[col] = temp_row0[col]|
						   temp_row1[col]|
						   temp_row2[col];
		}
	}
	// free_ui8matrix(temp, nrl + SE_NRL, nrh + SE_NRH, ncl, nch);
}
*/
/*
void ui8matrix_dilation_divide_row_and_conquer_OMP (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **temp_buffer, uint8 **Y)
{
	const long order = 3;
	// uint8 **temp = ui8matrix(nrl + SE_NRL, nrh + SE_NRH, ncl, nch);
	uint8 *temp_row, *temp_row0, *temp_row1, *temp_row2, *out_row, *in_row;
	long row, col, nrow, r;
	int nb_threads;

	nrow = nrh - nrl + 1;
	row = 0;
	omp_set_num_threads(omp_get_max_threads());
	#pragma omp parallel shared(temp_buffer, X, Y, ncl, nrl, nch, nrh) 
	{
		#pragma omp master 
		{
			uint8 *temp_row, *in_row;
			long col;

			temp_row = temp_buffer[nrh + 1];
			in_row = X[nrh + 1];
			for(col = ncl; col < nch + 1; col++) {
				temp_row[col] = scalar_or3(in_row, col);
			}
		}

		#pragma omp for nowait private(temp_row, in_row, row, col)
		for (row = nrl - 1; row < nrh + 1; row += order)
		{
			temp_row =    temp_buffer[row];
			in_row   = X[row];
			for(col = ncl; col < nch + 1; col++) {
				temp_row[col] = scalar_or3(in_row, col);
			}
		}
		#pragma omp for nowait private(temp_row, in_row, row, col)
		for (row = nrl + 0; row < nrh + 1; row += order)
		{

			temp_row =    temp_buffer[row];
			in_row   = X[row];
			for(col = ncl; col < nch + 1; col++) {
				temp_row[col] = scalar_or3(in_row, col);
			}
		}
		#pragma omp for nowait private(temp_row, in_row, row, col)
		for (row = nrl + 1; row < (nrh + 1) + 1; row += order)
		{
			temp_row =    temp_buffer[row];
			in_row   = X[row];
			for(col = ncl; col < nch + 1; col++) {
				temp_row[col] = scalar_or3(in_row, col);
			}
		}
		#pragma omp barrier
		#pragma omp for nowait private(temp_row0, temp_row1, temp_row2, out_row, row, col) schedule(auto)
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
	// free_ui8matrix(temp, nrl + SE_NRL, nrh + SE_NRH, ncl, nch);
}
*/
void ui8matrix_dilation_divide_row_and_conquer_OMP (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **temp_buffer, uint8 **Y)
{
	const long order = 3;
	// uint8 **temp = ui8matrix(nrl + SE_NRL, nrh + SE_NRH, ncl, nch);
	uint8 *temp_row, *temp_row0, *temp_row1, *temp_row2, *out_row, *in_row;
	long row, col, nrow, r;
	int nb_threads;

	nrow = nrh - nrl + 1;
	row = 0;
	omp_set_num_threads(omp_get_max_threads());
	#pragma omp parallel shared(temp_buffer, X, Y, ncl, nrl, nch, nrh) 
	{
		#pragma omp for nowait private(temp_row, in_row, row, col) schedule(auto)
		for (row = nrl - 1; row < nrh + 1 + SE_NRH; row ++)
		{
			temp_row =    temp_buffer[row];
			in_row   = X[row];
			for(col = ncl; col < nch + 1; col++) {
				temp_row[col] = scalar_or3(in_row, col);
			}
		}
		#pragma omp barrier
		#pragma omp for nowait private(temp_row0, temp_row1, temp_row2, out_row, row, col) schedule(auto)
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
	// free_ui8matrix(temp, nrl + SE_NRL, nrh + SE_NRH, ncl, nch);
}

void ui8matrix_dilation_divide_row_and_conquer_ExLU_O3_OMP (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **temp_buffer, uint8 **Y)
{
	const long order = 3;
	// uint8 **temp = ui8matrix(nrl + SE_NRL, nrh + SE_NRH, ncl, nch);
	uint8 *temp_row, *temp_row0, *temp_row1, *temp_row2, *temp_row3, *temp_row4, *in_row, *in_row0, *in_row1, *in_row2, *out_row;
	long row, col, nrow, r;
	int nb_threads;

	nrow = nrh - nrl + 1;
	row = 0;

	r = (nrh + 1) % order;
	#pragma omp parallel shared(temp_buffer, X, Y, ncl, nrl, nch, nrh, r) 
	{
	
		temp_row0 =    temp_buffer[nrl - 1];
		in_row0   = X[nrl - 1];
		#pragma omp for nowait firstprivate(temp_row0, in_row0) private(col)
		for(col = ncl; col < nch + 1; col++) {
			temp_row0[col] = scalar_or3(in_row0, col);
		}
		#pragma omp for nowait private(temp_row0, temp_row1, temp_row2, in_row0, in_row1, in_row2, row, col)
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
			#pragma omp for nowait firstprivate(temp_row0, temp_row1, temp_row2, in_row0, in_row1, in_row2) private(row, col)
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
			#pragma omp for nowait firstprivate(temp_row1, temp_row2, in_row1, in_row2) private(row, col)
			for(col = ncl; col < nch + 1; col++) {
				temp_row1[col] = scalar_or3(in_row1, col);
				temp_row2[col] = scalar_or3(in_row2, col);
			}
			break;
			case 0 :
			temp_row2  = temp_buffer[nrh + 1];
			in_row2 = X[nrh + 1];
			#pragma omp for nowait firstprivate(temp_row2, in_row2) private(row, col)
			for(col = ncl; col < nch + 1; col++) {
				temp_row2[col] = scalar_or3(in_row2, col);
			}
			break;
		}
		#pragma omp barrier
		#pragma omp for nowait private(temp_row0, temp_row1, temp_row2, out_row, row, col) schedule(auto)
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
	
	// free_ui8matrix(temp, nrl + SE_NRL, nrh + SE_NRH, ncl, nch);
}


void ui8matrix_dilation_divide_row_and_conquer_InLU_O3_OMP (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **temp_buffer, uint8 **Y)
{
	const long order = 3;
	// uint8 **temp = ui8matrix(nrl + SE_NRL, nrh + SE_NRH, ncl, nch);
	uint8 *temp_row, *temp_row0, *temp_row1, *temp_row2, *temp_row3, *temp_row4, *out_row0, *out_row1, *out_row2, *in_row, *out_row;
	long row, col, nrow, r;
	int nb_threads;

	nrow = nrh - nrl + 1;
	row = 0;
	r = (nch + 1) % order;
	#pragma omp parallel shared(temp_buffer, X, Y, ncl, nrl, nch, nrh, r) 
	{
		
		#pragma omp for nowait private(temp_row, in_row, row, col)
		for (row = nrl - 1; row < (nrh + 1) + SE_NRH; row ++)
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
			#pragma omp for nowait private(temp_row, in_row, row)
			for (row = nrl - 1; row < (nrh + 1) + SE_NRH; row ++)
			{
				temp_row =    temp_buffer[row];
				in_row   = X[row];
				temp_row[nch - 1] = scalar_or3(in_row, nch - 1);
				temp_row[nch + 0] = scalar_or3(in_row, nch + 0);
			}
			break;
			case 1:
			#pragma omp for nowait private(temp_row, in_row, row)
			for (row = nrl - 1; row < (nrh + 1) + SE_NRH; row ++)
			{
				temp_row =    temp_buffer[row];
				in_row   = X[row];
				temp_row[nch + 0] = scalar_or3(in_row, nch + 0);
			}
			break;
			default: 
			break;
		}

		#pragma omp barrier
		#pragma omp for nowait private(temp_row0, temp_row1, temp_row2, out_row, row, col) schedule(auto)
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
	// free_ui8matrix(temp, nrl + SE_NRL, nrh + SE_NRH, ncl, nch);
}

void ui8matrix_dilation_pipeline2_LU3x3_InLU_O3_RR_OMP (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **temp_buffer, uint8 **Y)
{
	const long order = 3;
	long row = nrl, col = ncl, x, y, r = 0;
	// uint8 **temp = ui8matrix(nrl + SE_NRL, nrh + SE_NRH, ncl + SE_NCL, nch + SE_NCH);
	uint8 *row0, *row1, *row2, *row3, x0, x1, x2, x3, x4, x5, y0, y1, y2;
	uint8 *temp_row0, *temp_row1, *temp_row2, *temp_row3, *temp_row4;
	uint8 *out_row0, *out_row1, *out_row2, *out_row3;
	
	r = (nch + 1) % order;
	// Prologue	
	#pragma omp parallel for private(row0, row1, row2, temp_row1, row) 
	for (row = nrl; row < nrh + 1; row++) {
		row0 = X[row - 1];
		row1 = X[row + 0];
		row2 = X[row + 1];
		temp_row1 = temp_buffer[row + 0];

		temp_row1[ncl - 1] = row0[ncl - 1] | row1[ncl - 1] | row2[ncl - 1];
		temp_row1[ncl + 0] = row0[ncl + 0] | row1[ncl + 0] | row2[ncl + 0];
	}
	#pragma omp parallel for private(row0, row1, row2, temp_row1, row, x0, x1, x2, x3, x4) 
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
			// x5 = row0[col + 3] | row1[col + 3] | row2[col + 3];
			
			out_row1 [col + 0] = x0 | x1 | x2; // A | B | C
			out_row1 [col + 1] = x1 | x2 | x3; // B | C | D
			out_row1 [col + 2] = x2 | x3 | x4; // C | D | E

			x0 = x3; x1 = x4;
		}
	}
	switch (r) {
		case 2: 
			#pragma omp parallel 	
			{
				row0 = X[nrl - 1];
				row1 = X[nrl + 0];
				#pragma omp parallel for default(none) firstprivate(row0, row1) private(out_row1, row, row2, x1) shared(X, Y, nrl, nrh, ncl, nch)
				for (row = nrl; row < nrh + 1; row++) {
					row2 = X[row + 1];
					out_row1 = Y[row]; 

					x1 = row0[nch - 1] | row1[nch - 1] | row2[nch - 1]|
						row0[nch + 0] | row1[nch + 0] | row2[nch + 0];

					out_row1 [nch - 1] = x1 | row0[nch - 2] | row1[nch - 2] | row2[nch - 2]; // A | B | C
					out_row1 [nch + 0] = x1 | row0[nch + 1] | row1[nch + 1] | row2[nch + 1]; // B | C | D
					row0 = row1;
					row1 = row2;
				}
			}
			break;
		case 1: 
			#pragma omp parallel 	
			{
				row0 = X[nrl - 1];
				row1 = X[nrl + 0];
				#pragma omp parallel for default(none) firstprivate(row0, row1) private(out_row1, row, row2, x1) shared(X, Y, nrl, nrh, ncl, nch)
				for (row = nrl; row < nrh + 1; row++) {
					row2 = X[row + 1];
					Y[row][nch] = scalar_or3(row0, nch) |
										scalar_or3(row1, nch) |
										scalar_or3(row2, nch);
					row0 = row1;
					row1 = row2;
				}
			}
			break;
		default:
			break;
	}
	// display_ui8matrix(Y, nrl - 1, nrh + 1, ncl - 1, nch + 1, "%u", "TEst");
	// free_ui8matrix(temp, nrl + SE_NRL, nrh + SE_NRH, ncl + SE_NCL, nch + SE_NCH);
	
}
void ui8matrix_dilation_pipeline_LU3x3_ExLU_O3_RR_OMP (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **temp_buffer, uint8 **Y)
{
	const long order = 3;
	long row = nrl, col = ncl, x, y, r = 0;
	// uint8 **temp = ui8matrix(nrl + SE_NRL, nrh + SE_NRH, ncl, nch);
	uint8 *row0, *row1, *row2, *row3;
	uint8 y0, y1, y2, y3, y4, x0, x1, x2, x3, x4, x5;
	uint8 *temp_row0, *temp_row1, *temp_row2, *temp_row3, *temp_row4;
	uint8 *out_row0, *out_row1, *out_row2, *out_row3;
	
	r = (nrh + 1) % order;
	// Prologue	
	#pragma omp parallel for private(row0, row1, row2, temp_row1, temp_row2, row, x0, x1,x2,x3,x4,x5)
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
	
	
	#pragma omp parallel for private(row0, row1, row2, temp_row0, temp_row1, temp_row2, temp_row3, temp_row4, row, x0, x1, x2, x3, x4, x5, y0, y1, y2, y3, y4)
	for (row = nrl; row < nrh + 1 - r; row += order){
		// row0 = X[row - 1];   //  0 : b0 b1 b2 ...
		row1 = X[row + 1];   //  1 : c0 c1 c2 ...
		// row2 = X[row + 1];   //  2 : d0 d1 d2 ...
		// row3 = X[row + 2];   //  3 : e0 e1 e2 ...

		temp_row0 = temp_buffer[row - 1]; // -1 : A    <= (a0 a1 a2 ...) Prologued
		temp_row1 = temp_buffer[row + 0]; //  0 : B    <= (b0 b1 b2 ...) Prologued
		temp_row2 = temp_buffer[row + 1]; //  1 : NULL <= (c0 c1 c2 ...)
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
		temp_row2 = temp_buffer[nrh + 0]; //  1 : NULL <= (c0 c1 c2 ...)
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
		temp_row2 = temp_buffer[nrh - 0]; //  1 : NULL <= (c0 c1 c2 ...)
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
	// free_ui8matrix(temp,nrl + SE_NRL, nrh + SE_NRH, ncl, nch);
}



/*******************************************/
/******* Optimisation : Loop Unroll ********/
/*******************************************/

/*******************************************/
/******* Optimisation : Loop Unroll ********/
/*******************************************/

void ui8matrix_erosion_LU3x3_O1xO1_OMP(uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **temp_buffer, uint8 **Y)
{
	long row = nrl, col = ncl, x, y;
	// dilate
	omp_set_num_threads(omp_get_max_threads());
	#pragma omp parallel for default(none) shared(nrl, ncl, nrh, nch, Y, X) private(row, col) 
	for (row = nrl; row < nrh + 1; row++)
		for (col = ncl; col < nch + 1; col++)
            Y[row][col] = scalar_and3x3(&X[row], col);
}

void ui8matrix_erosion_LU3x3_InLU_O3_OMP(uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **temp_buffer, uint8 **Y)
{
	const long order = 3;
	long row = nrl, col = ncl, x, y, r;
	uint8 y0, y1, y2, x0, x1, x2, x3, x4, x5, x6, x7, x8;
	uint8 *row0, *row1, *row2, *row3, *row4;
	uint8 *out_row0;
	r = (nch + 1)  % order;
	
	omp_set_num_threads(omp_get_max_threads());

	#pragma omp parallel for default(none) shared(nrl, ncl, nrh, nch, Y, X, r) private(row, col, row0, row1, row2, out_row0) 		
	for (row = nrl; row < nrh + 1; row ++) {
		row0 = X[row - 1];
		row1 = X[row + 0];
		row2 = X[row + 1];
		out_row0 = Y[row + 0];
		for (col = ncl; col < nch + 1 - r; col += order) {
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
			#pragma omp parallel for default(none) shared(nrl, ncl, nrh, nch, Y, X, r) private(row, col, row0, row1, row2, out_row0) 		
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
			#pragma omp parallel for default(none) shared(nrl, ncl, nrh, nch, Y, X, r) private(row, col, row0, row1, row2, out_row0) 		
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

void ui8matrix_erosion_LU3x3_ExLU_O3_OMP(uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **temp_buffer, uint8 **Y)
{
	const long order = 3;
	long row = nrl, col = ncl, x, y;
	uint8 y0, y1, y2, x0, x1, x2, x3, r;
	uint8 *row0, *row1, *row2, *row3, *row4;
	uint8 *out_row0, *out_row1, *out_row2;

	r = (nrh + 1)  % order;
	
	omp_set_num_threads(omp_get_max_threads());
	#pragma omp parallel for default(none) shared(nrl, ncl, nrh, nch, Y, X, r) private(y0, y1, y2, x0, x1, x2, x3,row, col, row0, row1, row2,row3,row4, out_row0,out_row1, out_row2)
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
	// printf("remainder :%ld\n", r);
	switch(r) {
		case 2:
			row0 = X[nrh - 2];
			row1 = X[nrh - 1];
			row2 = X[nrh + 0];
			row3 = X[nrh + 1];
			out_row0 = Y[nrh - 1];
			out_row1 = Y[nrh + 0];
			#pragma omp parallel for default(none) shared(nrl, ncl, nrh, nch, Y, X, r, row0, row1, row2,row3,row4, out_row0,out_row1, out_row2) private(y0, y1, y2, x0, x1, x2, x3,row, col)		
			for (col = ncl; col < nch + 1; col++) {
				x0 = scalar_and3(row1, col);
				x1 = scalar_and3(row2, col);
				out_row0[col] = scalar_and3(row0, col) & x0 & x1;
				out_row1[col] = 				        x0 & x1 & scalar_and3(row3, col);
			}
			break;
		case 1:
			row0 = X[nrh - 1];
			row1 = X[nrh + 0];
			row2 = X[nrh + 1];
			out_row0 = Y[nrh + 0];
			#pragma omp parallel for default(none) shared(nrl, ncl, nrh, nch, Y, X, r, row0, row1, row2,row3,row4, out_row0,out_row1, out_row2) private(y0, y1, y2, x0, x1, x2, x3,row, col)		
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

void ui8matrix_erosion_LU3x3_ComLU_O3_OMP(uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **temp_buffer, uint8 **Y)
{
	const long order = 3;
	long row = nrl, col = ncl, x, y;
	uint8 y0, y1, y2, rr, cr;
	uint8 x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14;
	uint8 *row0, *row1, *row2, *row3, *row4;
	uint8 *out_row0, *out_row1, *out_row2;

	rr = (nrh + 1) % order;
	cr = (nch + 1) % order;
	
	omp_set_num_threads(omp_get_max_threads());
	#pragma omp parallel for default(none) shared(nrl, ncl, nrh, nch, Y, X, rr, cr) private(y0, y1, y2, x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14,row, col, row0, row1, row2,row3,row4, out_row0,out_row1, out_row2)	schedule(dynamic)
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
					ui8matrix_erosion_LU3x3_ExLU_O3_OMP(X, nrl, nrh, nch - 1, nch, temp_buffer, Y);
					ui8matrix_erosion_LU3x3_InLU_O3_OMP(X, nrh - 1, nrh, ncl, nch, temp_buffer, Y);
				break;
				case 1: 
					ui8matrix_erosion_LU3x3_ExLU_O3_OMP(X, nrl, nrh, nch, nch, temp_buffer, Y);		
					ui8matrix_erosion_LU3x3_InLU_O3_OMP(X, nrh - 1, nrh, ncl, nch, temp_buffer, Y);	
				break;
				default :
					ui8matrix_erosion_LU3x3_InLU_O3_OMP(X, nrh - 1, nrh, ncl, nch, temp_buffer, Y);
				break;
			}
			break;
		case 1:
			switch(cr){
				case 2 :
					ui8matrix_erosion_LU3x3_ExLU_O3_OMP(X, nrl, nrh, nch - 1, nch, temp_buffer, Y);
					ui8matrix_erosion_LU3x3_InLU_O3_OMP(X, nrh, nrh, ncl, nch, temp_buffer, Y);
				break;
				case 1: 
					ui8matrix_erosion_LU3x3_ExLU_O3_OMP(X, nrl, nrh, nch, nch, temp_buffer, Y);		
					ui8matrix_erosion_LU3x3_InLU_O3_OMP(X, nrh, nrh, ncl, nch, temp_buffer, Y);	
				break;
				default :
					ui8matrix_erosion_LU3x3_InLU_O3_OMP(X, nrh, nrh, ncl, nch, temp_buffer, Y);
				break;
			}
		break;
		default:
			switch(cr){
				case 2 :
					ui8matrix_erosion_LU3x3_ExLU_O3_OMP(X, nrl, nrh, nch - 1, nch, temp_buffer, Y);
				break;
				case 1: 
					ui8matrix_erosion_LU3x3_ExLU_O3_OMP(X, nrl, nrh, nch, nch, temp_buffer, Y);	
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


/*****************************************************************************/
/******** Optimisation : Loop Unroll + Register Rotation of Addresses ********/
/*****************************************************************************/

void ui8matrix_erosion_LU3x3_InLU_O3_AddrRR_OMP(uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **temp_buffer, uint8 **Y)
{
	const long order = 3;
	long row = nrl, col = ncl, r, nrow;
	uint8 y0, y1, y2, x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14;
	uint8 *row0, *row1, *row2, *row3, *row4;
	uint8 *out_row0;
	int nb_threads;


 
	r = (nch + 1) % order;
	nrow = nrh - nrl + 1;
	
	omp_set_num_threads(omp_get_max_threads());
	#pragma omp parallel //shared(X, ncl, nch, Y, nb_threads, r, nrow)
	{
		row0 = X[nrl - 1];
		row1 = X[nrl + 0];
		#pragma omp parallel for firstprivate(row0, row1) private(out_row0, row2, y0, y1, y2, x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, row, col) 
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
	}

	switch(r) {
		case 2:

		#pragma omp parallel for default(none) firstprivate(nrl, nrh) private(out_row0, row0, row1, row2, row3, x0, x1, x2, x3, x5, x6, x7, x8, x10, x11, x12, x13, row, col) shared(X, ncl, nch, Y)
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
		#pragma omp parallel for default(none) firstprivate(nrl, nrh) private(out_row0, row0, row1, row2, row3, row) shared(X, ncl, nch, Y)
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
void ui8matrix_erosion_LU3x3_ExLU_O3_AddrRR_OMP(uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **temp_buffer, uint8 **Y)
{
	
	const long order = 3;
	long row = nrl, col = ncl, x, y, col0, col1, col2, nrow;
	uint8 y0, y1, y2, y3, y5, r;
	uint8 x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14;
	uint8 *row0, *row1, *row2, *row3, *row4;
	uint8 *out_row0, *out_row1, *out_row2;
	int nb_threads;
	r = (nrh + 1)  % order;
	nrow = nrh - nrl + 1;
	omp_set_num_threads(omp_get_max_threads());
	#pragma omp parallel 
	{
		row0 = X[nrl - 1];
		row1 = X[nrl + 0];
		#pragma omp parallel for firstprivate(row0, row1) private(out_row0,out_row1, out_row2,row2, row3, row4,y0, y1, y2, y3, y5, x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, row, col) 
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
	}
	
	out_row0 = Y[nrh - 1];
	out_row1 = Y[nrh + 0];
	switch(r) {
		case 2:
			row0 = X[nrh - 2];
			row1 = X[nrh - 1];
			row2 = X[nrh + 0];
			row3 = X[nrh + 1];
			#pragma omp parallel for default(none) private(row,x0,x1,x2,x3,x4,x5,x6,x7,x8,col) shared(X, ncl, nch, Y, out_row0,out_row1, row0, row1, row2, row3)
			for (col = ncl; col < nch + 1; col++) {
				x0 = row1[col - 1]; x3 = row1[col + 0]; x6 = row1[col + 1]; //scalar_and3(row1, col);
				x1 = row2[col - 1]; x4 = row2[col + 0]; x7 = row2[col + 1]; //scalar_and3(row2, col);
				x2 = row3[col - 1]; x5 = row3[col + 0]; x8 = row3[col + 1]; //scalar_and3(row3, col);
				out_row0[col] = scalar_and3(row0, col) & (x0 & x3 & x6) & (x1 & x4 & x7);
				out_row1[col] =        (x0 & x3 & x6) & (x1 & x4 & x7) & (x2 & x5 & x8);
			}
			break;
		case 1:
			row0 = X[nrh - 1];
			row1 = X[nrh + 0];
			row2 = X[nrh + 1];
			#pragma omp parallel for default(none) private(row,x0,x1,x2,x3,x4,x5,x6,x7,x8,col) shared(X, ncl, nch, Y, out_row0,out_row1, row0, row1, row2, row3)
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
void ui8matrix_erosion_LU3x3_ComLU_O3_AddrRR_OMP(uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **temp_buffer, uint8 **Y)
{
	const long order = 3;
	long row = nrl, col = ncl, rr, cr, nrow;
	uint8 r0c0, r1c0, r2c0, r0c1,r1c1, r2c1, r0c2, r1c2, r2c2;
	uint8 *row0, *row1, *row2, *row3, *row4;
	uint8 *out_row0, *out_row1, *out_row2;
	int nb_threads;

	rr = (nrh + 1) % order;
	cr = (nch + 1) % order;
	nrow = nrh - nrl + 1;
	omp_set_num_threads(omp_get_max_threads());
	#pragma omp parallel 
	{
		uint8 y0, y1, y2, y3, y4, y5, y6, y7, y8, y9, y10, y11, y12, y13, y14;
		uint8 x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16, x17, x18, x19, x20, x21, x22, x23, x24, x25;
		
		row0 = X[nrl - 1];
		row1 = X[nrl + 0];
		#pragma omp parallel for default(none) firstprivate(row0, row1) private(out_row0, out_row1, out_row2, row2, row3, row4, row, col, y0, y1, y2, y3, y4, y5, y6, y7, y8, y9, y10, y11, y12, y13, y14, x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16, x17, x18, x19, x20, x21, x22, x23, x24, x25) shared(X, nrl, nrh, ncl, nch, Y, nb_threads, cr, rr, nrow)
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
	}
	switch (rr) {
		case 2 :
			switch(cr){
				case 2 :
					ui8matrix_erosion_LU3x3_ExLU_O3_AddrRR_OMP(X, nrl, nrh, nch - 1, nch, temp_buffer, Y);
					ui8matrix_erosion_LU3x3_InLU_O3_AddrRR_OMP(X, nrh - 1, nrh, ncl, nch, temp_buffer, Y);
				break;
				case 1: 
					ui8matrix_erosion_LU3x3_ExLU_O3_AddrRR_OMP(X, nrl, nrh, nch, nch, temp_buffer, Y);		
					ui8matrix_erosion_LU3x3_InLU_O3_AddrRR_OMP(X, nrh - 1, nrh, ncl, nch, temp_buffer, Y);	
				break;
				default :
					ui8matrix_erosion_LU3x3_InLU_O3_AddrRR_OMP(X, nrh - 1, nrh, ncl, nch, temp_buffer, Y);
				break;
			}
			break;
		case 1:
			switch(cr){
				case 2 :
					ui8matrix_erosion_LU3x3_ExLU_O3_AddrRR_OMP(X, nrl, nrh, nch - 1, nch, temp_buffer, Y);
					ui8matrix_erosion_LU3x3_InLU_O3_AddrRR_OMP(X, nrh, nrh, ncl, nch, temp_buffer, Y);
				break;
				case 1: 
					ui8matrix_erosion_LU3x3_ExLU_O3_AddrRR_OMP(X, nrl, nrh, nch, nch, temp_buffer, Y);		
					ui8matrix_erosion_LU3x3_InLU_O3_AddrRR_OMP(X, nrh, nrh, ncl, nch, temp_buffer, Y);	
					// display_ui8matrix(Y, nrl, nrh, ncl, nch, "%u", "rr=1; cr=1");
				break;
				default :
					ui8matrix_erosion_LU3x3_InLU_O3_AddrRR_OMP(X, nrh, nrh, ncl, nch, temp_buffer, Y);
				break;
			}
		break;
		default:
			switch(cr){
				case 2 :
					ui8matrix_erosion_LU3x3_ExLU_O3_AddrRR_OMP(X, nrl, nrh, nch - 1, nch, temp_buffer, Y);
				break;
				case 1: 
					ui8matrix_erosion_LU3x3_ExLU_O3_AddrRR_OMP(X, nrl, nrh, nch, nch, temp_buffer, Y);	
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
/*
void ui8matrix_erosion_row_pipeline_OMP(uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **temp_buffer, uint8 **Y)
{
	const long order = 3;
	// uint8 **temp = ui8matrix(nrl + SE_NRL, nrh + SE_NRH, ncl, nch);
	uint8 *temp_row, *temp_row0, *temp_row1, *temp_row2, *out_row, *in_row;
	long row, col, nrow, r;
	int nb_threads;

	nrow = nrh - nrl + 1;
	#pragma omp parallel default(none) private(row, col, temp_row0, temp_row1) shared(nrow, ncl, nch, nrl, nrh, X, Y, temp) 
	{
		uint8 *temp_row, *in_row, *out_row;
		row = 0;
		// Preprocess the epilogue (last border)
		#pragma omp master 
		{
			temp_row = temp_buffer[nrh + 1];
			in_row = X[nrh + 1];
			for(col = ncl; col < nch + 1; col++) {
				temp_row[col] = scalar_and3(in_row, col);
			}
		}

		#pragma omp for nowait private(row, col, temp_row, in_row)
		for (row = nrl - 1; row < nrh + 1; row += order)
		{
			temp_row =    temp_buffer[row];
			in_row   = X[row];
			for(col = ncl; col < nch + 1; col++) {
				temp_row[col] = scalar_and3(in_row, col);
			}
		}
		row = 0;
		#pragma omp for nowait private(row, col,temp_row, in_row)
		for (row = nrl + 0; row < nrh + 1; row += order)
		{

			temp_row =    temp_buffer[row];
			in_row   = X[row];
			for(col = ncl; col < nch + 1; col++) {
				temp_row[col] = scalar_and3(in_row, col);
			}
		}
		row = 0;
		#pragma omp for nowait private(row, col, temp_row, in_row)
		for (row = nrl + 1; row < (nrh + 1) + 1; row += order)
		{
			temp_row =    temp_buffer[row];
			in_row   = X[row];
			for(col = ncl; col < nch + 1; col++) {
				temp_row[col] = scalar_and3(in_row, col);
			}
		}
	}
	
	for (row = nrl; row < nrh + 1; row ++)
	{
		temp_row0 =     temp_buffer[row - 1];
		temp_row1 =     temp_buffer[row + 0];
		temp_row2 =     temp_buffer[row + 1];
		out_row   = Y[row];
		for(col = nrl; col < nch + 1; col++) {
			out_row[col] = temp_row0[col]&
						   temp_row1[col]&
						   temp_row2[col];
		}
	}
	// free_ui8matrix(temp, nrl + SE_NRL, nrh + SE_NRH, ncl, nch);
}
*/
/*
void ui8matrix_erosion_divide_row_and_conquer_OMP (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **temp_buffer, uint8 **Y)
{
	const long order = 3;
	// uint8 **temp = ui8matrix(nrl + SE_NRL, nrh + SE_NRH, ncl, nch);
	uint8 *temp_row, *temp_row0, *temp_row1, *temp_row2, *out_row, *in_row;
	long row, col, nrow, r;
	int nb_threads;

	nrow = nrh - nrl + 1;
	row = 0;
	omp_set_num_threads(omp_get_max_threads());
	#pragma omp parallel shared(temp_buffer, X, Y, ncl, nrl, nch, nrh) 
	{
		#pragma omp master 
		{
			uint8 *temp_row, *in_row;
			long col;

			temp_row = temp_buffer[nrh + 1];
			in_row = X[nrh + 1];
			for(col = ncl; col < nch + 1; col++) {
				temp_row[col] = scalar_and3(in_row, col);
			}
		}

		#pragma omp for nowait private(temp_row, in_row, row, col)
		for (row = nrl - 1; row < nrh + 1; row += order)
		{
			temp_row =    temp_buffer[row];
			in_row   = X[row];
			for(col = ncl; col < nch + 1; col++) {
				temp_row[col] = scalar_and3(in_row, col);
			}
		}
		#pragma omp for nowait private(temp_row, in_row, row, col)
		for (row = nrl + 0; row < nrh + 1; row += order)
		{

			temp_row =    temp_buffer[row];
			in_row   = X[row];
			for(col = ncl; col < nch + 1; col++) {
				temp_row[col] = scalar_and3(in_row, col);
			}
		}
		#pragma omp for nowait private(temp_row, in_row, row, col)
		for (row = nrl + 1; row < (nrh + 1) + 1; row += order)
		{
			temp_row =    temp_buffer[row];
			in_row   = X[row];
			for(col = ncl; col < nch + 1; col++) {
				temp_row[col] = scalar_and3(in_row, col);
			}
		}
		#pragma omp barrier
		#pragma omp for nowait private(temp_row0, temp_row1, temp_row2, out_row, row, col) schedule(auto)
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
	}
	// free_ui8matrix(temp, nrl + SE_NRL, nrh + SE_NRH, ncl, nch);
}
*/
void ui8matrix_erosion_divide_row_and_conquer_OMP (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **temp_buffer, uint8 **Y)
{
	const long order = 3;
	// uint8 **temp = ui8matrix(nrl + SE_NRL, nrh + SE_NRH, ncl, nch);
	uint8 *temp_row, *temp_row0, *temp_row1, *temp_row2, *out_row, *in_row;
	long row, col, nrow, r;
	int nb_threads;

	nrow = nrh - nrl + 1;
	row = 0;
	omp_set_num_threads(omp_get_max_threads());
	#pragma omp parallel shared(temp_buffer, X, Y, ncl, nrl, nch, nrh) 
	{
		#pragma omp for nowait private(temp_row, in_row, row, col) schedule(auto)
		for (row = nrl - 1; row < nrh + 1 + SE_NRH; row ++)
		{
			temp_row =    temp_buffer[row];
			in_row   = X[row];
			for(col = ncl; col < nch + 1; col++) {
				temp_row[col] = scalar_and3(in_row, col);
			}
		}
		#pragma omp barrier
		#pragma omp for nowait private(temp_row0, temp_row1, temp_row2, out_row, row, col) schedule(auto)
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
	}
	// free_ui8matrix(temp, nrl + SE_NRL, nrh + SE_NRH, ncl, nch);
}

void ui8matrix_erosion_divide_row_and_conquer_ExLU_O3_OMP (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **temp_buffer, uint8 **Y)
{
	const long order = 3;
	// uint8 **temp = ui8matrix(nrl + SE_NRL, nrh + SE_NRH, ncl, nch);
	uint8 *temp_row, *temp_row0, *temp_row1, *temp_row2, *temp_row3, *temp_row4, *in_row, *in_row0, *in_row1, *in_row2, *out_row;
	long row, col, nrow, r;
	int nb_threads;

	nrow = nrh - nrl + 1;
	row = 0;

	r = (nrh + 1) % order;
	#pragma omp parallel shared(temp_buffer, X, Y, ncl, nrl, nch, nrh, r) 
	{
	
		temp_row0 =    temp_buffer[nrl - 1];
		in_row0   = X[nrl - 1];
		#pragma omp for nowait firstprivate(temp_row0, in_row0) private(col)
		for(col = ncl; col < nch + 1; col++) {
			temp_row0[col] = scalar_and3(in_row0, col);
		}
		#pragma omp for nowait private(temp_row0, temp_row1, temp_row2, in_row0, in_row1, in_row2, row, col)
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
			#pragma omp for nowait firstprivate(temp_row0, temp_row1, temp_row2, in_row0, in_row1, in_row2) private(row, col)
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
			#pragma omp for nowait firstprivate(temp_row1, temp_row2, in_row1, in_row2) private(row, col)
			for(col = ncl; col < nch + 1; col++) {
				temp_row1[col] = scalar_and3(in_row1, col);
				temp_row2[col] = scalar_and3(in_row2, col);
			}
			break;
			case 0 :
			temp_row2  = temp_buffer[nrh + 1];
			in_row2 = X[nrh + 1];
			#pragma omp for nowait firstprivate(temp_row2, in_row2) private(row, col)
			for(col = ncl; col < nch + 1; col++) {
				temp_row2[col] = scalar_and3(in_row2, col);
			}
			break;
		}
		#pragma omp barrier
		#pragma omp for nowait private(temp_row0, temp_row1, temp_row2, out_row, row, col) schedule(auto)
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
	}
	
	// free_ui8matrix(temp, nrl + SE_NRL, nrh + SE_NRH, ncl, nch);
}


void ui8matrix_erosion_divide_row_and_conquer_InLU_O3_OMP (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **temp_buffer, uint8 **Y)
{
	const long order = 3;
	// uint8 **temp = ui8matrix(nrl + SE_NRL, nrh + SE_NRH, ncl, nch);
	uint8 *temp_row, *temp_row0, *temp_row1, *temp_row2, *temp_row3, *temp_row4, *out_row0, *out_row1, *out_row2, *in_row, *out_row;
	long row, col, nrow, r;
	int nb_threads;

	nrow = nrh - nrl + 1;
	row = 0;
	r = (nch + 1) % order;
	#pragma omp parallel shared(temp_buffer, X, Y, ncl, nrl, nch, nrh, r) 
	{
		
		#pragma omp for nowait private(temp_row, in_row, row, col)
		for (row = nrl - 1; row < (nrh + 1) + SE_NRH; row ++)
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
			#pragma omp for nowait private(temp_row, in_row, row)
			for (row = nrl - 1; row < (nrh + 1) + SE_NRH; row ++)
			{
				temp_row =    temp_buffer[row];
				in_row   = X[row];
				temp_row[nch - 1] = scalar_and3(in_row, nch - 1);
				temp_row[nch + 0] = scalar_and3(in_row, nch + 0);
			}
			break;
			case 1:
			#pragma omp for nowait private(temp_row, in_row, row)
			for (row = nrl - 1; row < (nrh + 1) + SE_NRH; row ++)
			{
				temp_row =    temp_buffer[row];
				in_row   = X[row];
				temp_row[nch + 0] = scalar_and3(in_row, nch + 0);
			}
			break;
			default: 
			break;
		}

		#pragma omp barrier
		#pragma omp for nowait private(temp_row0, temp_row1, temp_row2, out_row, row, col) schedule(auto)
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
	
	}
	// free_ui8matrix(temp, nrl + SE_NRL, nrh + SE_NRH, ncl, nch);
}

void ui8matrix_erosion_pipeline2_LU3x3_InLU_O3_RR_OMP (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **temp_buffer, uint8 **Y)
{
	const long order = 3;
	long row = nrl, col = ncl, x, y, r = 0;
	// uint8 **temp = ui8matrix(nrl + SE_NRL, nrh + SE_NRH, ncl + SE_NCL, nch + SE_NCH);
	uint8 *row0, *row1, *row2, *row3, x0, x1, x2, x3, x4, x5, y0, y1, y2;
	uint8 *temp_row0, *temp_row1, *temp_row2, *temp_row3, *temp_row4;
	uint8 *out_row0, *out_row1, *out_row2, *out_row3;
	
	r = (nch + 1) % order;
	// Prologue	
	#pragma omp parallel for private(row0, row1, row2, temp_row1, row) 
	for (row = nrl; row < nrh + 1; row++) {
		row0 = X[row - 1];
		row1 = X[row + 0];
		row2 = X[row + 1];
		temp_row1 = temp_buffer[row + 0];

		temp_row1[ncl - 1] = row0[ncl - 1] & row1[ncl - 1] & row2[ncl - 1];
		temp_row1[ncl + 0] = row0[ncl + 0] & row1[ncl + 0] & row2[ncl + 0];
	}
	#pragma omp parallel for private(row0, row1, row2, temp_row1, row, x0, x1, x2, x3, x4) 
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
			// x5 = row0[col + 3] & row1[col + 3] & row2[col + 3];
			
			out_row1 [col + 0] = x0 & x1 & x2; // A & B & C
			out_row1 [col + 1] = x1 & x2 & x3; // B & C & D
			out_row1 [col + 2] = x2 & x3 & x4; // C & D & E

			x0 = x3; x1 = x4;
		}
	}
	switch (r) {
		case 2: 
			#pragma omp parallel 	
			{
				row0 = X[nrl - 1];
				row1 = X[nrl + 0];
				#pragma omp parallel for default(none) firstprivate(row0, row1) private(out_row1, row, row2, x1) shared(X, Y, nrl, nrh, ncl, nch)
				for (row = nrl; row < nrh + 1; row++) {
					row2 = X[row + 1];
					out_row1 = Y[row]; 

					x1 = row0[nch - 1] & row1[nch - 1] & row2[nch - 1]&
						row0[nch + 0] & row1[nch + 0] & row2[nch + 0];

					out_row1 [nch - 1] = x1 & row0[nch - 2] & row1[nch - 2] & row2[nch - 2]; // A & B & C
					out_row1 [nch + 0] = x1 & row0[nch + 1] & row1[nch + 1] & row2[nch + 1]; // B & C & D
					row0 = row1;
					row1 = row2;
				}
			}
			break;
		case 1: 
			#pragma omp parallel 	
			{
				row0 = X[nrl - 1];
				row1 = X[nrl + 0];
				#pragma omp parallel for default(none) firstprivate(row0, row1) private(out_row1, row, row2, x1) shared(X, Y, nrl, nrh, ncl, nch)
				for (row = nrl; row < nrh + 1; row++) {
					row2 = X[row + 1];
					Y[row][nch] = scalar_and3(row0, nch) &
										scalar_and3(row1, nch) &
										scalar_and3(row2, nch);
					row0 = row1;
					row1 = row2;
				}
			}
			break;
		default:
			break;
	}
	// display_ui8matrix(Y, nrl - 1, nrh + 1, ncl - 1, nch + 1, "%u", "TEst");
	// free_ui8matrix(temp, nrl + SE_NRL, nrh + SE_NRH, ncl + SE_NCL, nch + SE_NCH);
	
}
void ui8matrix_erosion_pipeline_LU3x3_ExLU_O3_RR_OMP (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **temp_buffer, uint8 **Y)
{
	const long order = 3;
	long row = nrl, col = ncl, x, y, r = 0;
	// uint8 **temp = ui8matrix(nrl + SE_NRL, nrh + SE_NRH, ncl, nch);
	uint8 *row0, *row1, *row2, *row3;
	uint8 y0, y1, y2, y3, y4, x0, x1, x2, x3, x4, x5;
	uint8 *temp_row0, *temp_row1, *temp_row2, *temp_row3, *temp_row4;
	uint8 *out_row0, *out_row1, *out_row2, *out_row3;
	
	r = (nrh + 1) % order;
	// Prologue	
	#pragma omp parallel for private(row0, row1, row2, temp_row1, temp_row2, row, x0, x1,x2,x3,x4,x5)
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
	
	
	#pragma omp parallel for private(row0, row1, row2, temp_row0, temp_row1, temp_row2, temp_row3, temp_row4, row, x0, x1, x2, x3, x4, x5, y0, y1, y2, y3, y4)
	for (row = nrl; row < nrh + 1 - r; row += order){
		// row0 = X[row - 1];   //  0 : b0 b1 b2 ...
		row1 = X[row + 1];   //  1 : c0 c1 c2 ...
		// row2 = X[row + 1];   //  2 : d0 d1 d2 ...
		// row3 = X[row + 2];   //  3 : e0 e1 e2 ...

		temp_row0 = temp_buffer[row - 1]; // -1 : A    <= (a0 a1 a2 ...) Prologued
		temp_row1 = temp_buffer[row + 0]; //  0 : B    <= (b0 b1 b2 ...) Prologued
		temp_row2 = temp_buffer[row + 1]; //  1 : NULL <= (c0 c1 c2 ...)
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
		temp_row2 = temp_buffer[nrh + 0]; //  1 : NULL <= (c0 c1 c2 ...)
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
		temp_row2 = temp_buffer[nrh - 0]; //  1 : NULL <= (c0 c1 c2 ...)
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
	// free_ui8matrix(temp,nrl + SE_NRL, nrh + SE_NRH, ncl, nch);
}



/*******************************************/
/******* Optimisation : Loop Unroll ********/
/*******************************************/
