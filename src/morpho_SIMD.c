/* ------------------------ */
/* ----- morpho_SIMD.c ---- */
/* ------------------------ */

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <emmintrin.h>
#include <tmmintrin.h>
#include <simd_macro.h>

#include "nrdef.h"
#include "nrutil.h"

#include "vnrdef.h"
#include "vnrutil.h"
#include "myvnrutil.h"

#include "mutil.h"

#include "mymacro.h"
#include "simd_macro.h"

#include "util.h"
#include "img.h"
#include "img_SIMD.h"
#include "morpho.h"
#include "morpho_SIMD.h"
#include <omp.h>
#define SE_NRL -1
#define SE_NRH  1
#define SE_NCL -1
#define SE_NCH  1

#define get_vec_edge(vX, i, v1, mask) _mm_and_si128(vec_right1(vX[i][v1], vX[i][v1]), mask);

void ui8matrix_sequence_SIMD_Pipeline_FO_InLU_O3_ValAddrRR_OMP(vuint8** X, int nrl, int nrh, long ncl, long nch, int v0, int v1, vuint8 **Y , vuint8 **Z) {
vuint8  *in_row, *out_row;
	vuint8 **in, **mid, **out;
	int d5_row, e3_row, row, col;
	int nrl_prime, nrh_prime;
	int thread_num, nb_threads, nrow;
	const int PRE_D5_NROW = 4;
	const int PRE_E3_NROW = 1;
	int card = card_vuint8(); 
	int last_v    = (nch) / card,  last_vcol = (nch) % card; 
	// Mask for edge
	vuint8 vMask = _mm_setzero_si128(); 
	uint8 *mask= (uint8*)&vMask;
	for (col = 0; col < last_vcol + 1; col++) mask[col] = 0xFF;

	nrow = (nrh - nrl) + 1;
	if ((v1 - v0) == 2) {
		printf("Start %d ~ %d\n", v0, v1);
		// ui8matrix_sequence_SIMD_Pipeline_FO_InLU_O3_ValAddrRR(X, nrl, nrh, ncl, nch, v0, v0, Y, Z);
		print_vui8matrix(X, nrl - 2, nrh + 2, v0, v1, "%4u", "Input 1\n");
		ui8matrix_sequence_SIMD_Pipeline_FO_InLU_O3_ValAddrRR(X, nrl, nrh, ncl, nch, 0, 0, Y, Z);
		print_vui8matrix(Z, nrl - 2, nrh + 2, v0, v1, "%4u", "Resultat 1\n");
		print_vui8matrix(X, nrl - 2, nrh + 2, v0, v1, "%4u", "Input 2\n");
		ui8matrix_sequence_SIMD_Pipeline_FO_InLU_O3_ValAddrRR(X, nrl, nrh, ncl, nch, 1, 1, Y, Z);
		print_vui8matrix(Z, nrl - 2, nrh + 2, v0, v1, "%4u", "Resultat 2\n");
		
	}
	else {
		ui8matrix_sequence_SIMD_Pipeline_FO_InLU_O3_ValAddrRR(X, nrl, nrh, ncl, nch, v0, v1, Y, Z);
	}
}

void ui8matrix_sequence_SIMD_FO_InLU_O3_ValAddrRR(vuint8** X, int nrl, int nrh, long ncl, long nch, int v0, int v1, vuint8 **Y , vuint8 **Z) {
	vuint8  *in_row, *out_row;
	vuint8 **in, **mid, **out;
	long row, col, nrow, r;
	int card = card_vuint8(); 
	int last_v    = (nch) / card; 
	int last_vcol = (nch) % card; 
	vuint8 vMask = _mm_setzero_si128(); 
	uint8 *mask= (uint8*)&vMask;
	for (col = 0; col < last_vcol + 1; col++) mask[col] = 0xFF;
	// printf("%d / %d\n", last_v, last_vcol);
	// print_vui8vector(&vMask, 0, 0, 0, 15, " %u", "Mask : "); printf("\n");
	

	in = X; mid = NULL; out = Y;
	ui8matrix_erosion_SIMD_InLU_O3_ValAddrRR(in, nrl, nrh, ncl, nch, v0, v1, mid, out);/// zero_vui8matrix(X, nrl, nrh, v0, v1);

	
	in = Y; mid = NULL; out = X;
	ui8matrix_dilation5_SIMD_InLU_O3_ValAddrRR(in, nrl, nrh, ncl, nch, v0, v1, mid, out); //zero_vui8matrix(Z, nrl, nrh, v0, v1);
	
	for (row = nrl - 2; row < nrh + 1 + 2; row++) {
		out[row][last_v] = _mm_and_si128(out[row][last_v], vMask);
		for (col = last_v + 1; col < v1 + 1; col++) 
			out[row][col] = _mm_setzero_si128();
	}

	in = X; mid = NULL; out = Z;
	ui8matrix_erosion_SIMD_InLU_O3_ValAddrRR(X, nrl, nrh, ncl, nch, v0, v1, Y, Z);
}

void ui8matrix_sequence_SIMD_FO_InLU_O3_ValAddrRR_OMP(vuint8** X, int nrl, int nrh, long ncl, long nch, int v0, int v1, vuint8 **Y , vuint8 **Z) {
	vuint8  *in_row, *out_row;
	vuint8 **in, **mid, **out;
	long row, col, nrow, r;
	int card = card_vuint8(); 
	int last_v    = (nch) / card; 
	int last_vcol = (nch) % card; 
	vuint8 vMask = _mm_setzero_si128(); 
	uint8 *mask= (uint8*)&vMask;
	for (col = 0; col < last_vcol + 1; col++) mask[col] = 0xFF;
	int thread_num, nb_threads; 
	omp_set_num_threads(omp_get_max_threads());

	nrow = (nrh - nrl) + 1;
	// printf("%d ~ %d (%d)\n", nrl, nrh, nrow / nb_threads);
	
	#pragma omp parallel firstprivate(nrl, nrh) private(in, mid, out) firstprivate(X, Y, Z) shared(nrow)
	{
		nb_threads = omp_get_num_threads();
		thread_num = omp_get_thread_num();
		nrl = nrl + (thread_num + 0) * (nrow / nb_threads); 
		if ( thread_num < nb_threads -1 )
			nrh = (thread_num + 1) * (nrow / nb_threads) - 1;   


		in = X; mid = NULL; out = Y;
		ui8matrix_erosion_SIMD_InLU_O3_ValAddrRR(in, nrl, nrh, ncl, nch, v0, v1, mid, out); 
		#pragma omp barrier
		// zero_vui8matrix(X, nrl, nrh, v0, v1);
		
		in = Y; mid = NULL; out = X;
		ui8matrix_dilation5_SIMD_InLU_O3_ValAddrRR(in, nrl, nrh, ncl, nch, v0, v1, mid, out); 

		#pragma omp barrier
		// zero_vui8matrix(Z, nrl, nrh, v0, v1);

		for (row = nrl - 2; row < nrh + 1 + 2; row++) {
			out[row][last_v] = _mm_and_si128(out[row][last_v], vMask);
			for (col = last_v + 1; col < v1 + 1; col++) 
				out[row][col] = _mm_setzero_si128();
		}

		in = X; mid = NULL; out = Z;
		ui8matrix_erosion_SIMD_InLU_O3_ValAddrRR(X, nrl, nrh, ncl, nch, v0, v1, Y, Z);
	}
		
	
	
	
	// print_vui8matrix(Z, nrl, nrh, v0, v1, "%4u", "Test:\n");
}

void ui8matrix_sequence_SIMD_Pipeline_FO_InLU_O3_ValAddrRR(vuint8** X, int nrl, int nrh, long ncl, long nch, int v0, int v1, vuint8 **Y , vuint8 **Z) {
	
	vuint8  *in_row, *out_row;
	vuint8 **in, **mid, **out;
	int d5_row, e3_row, row, col;
	int nrl_prime, nrh_prime;
	const int PRE_D5_NROW = 4;
	const int PRE_E3_NROW = 1;
	int card = card_vuint8(); 
	int last_v    = (nch) / card,  last_vcol = (nch) % card; 
	// Mask for edge
	vuint8 vMask = _mm_setzero_si128(); 
	uint8 *mask= (uint8*)&vMask;
	for (col = 0; col < last_vcol + 1; col++) mask[col] = 0xFF;

	// Prologue
	in = X; out = Y;
	nrl_prime = nrl + PRE_D5_NROW;
	ui8matrix_erosion_SIMD_InLU_O3_ValAddrRR(in, nrl, nrl_prime, ncl, nch, v0, v1, mid, out);
	in = Y; out = X;
	nrl_prime = nrl + PRE_E3_NROW;
	ui8matrix_dilation5_SIMD_InLU_O3_ValAddrRR(in, nrl, nrl_prime, ncl, nch, v0, v1, mid, out);
	for (d5_row = nrl; d5_row < nrl_prime + 1; d5_row++){
		// Handle edge
		out_row = out[d5_row];
		out_row[last_v] = _mm_and_si128(out_row[last_v], vMask);
		for (col = last_v + 1; col < v1 + 1; col++) out_row[col] = _mm_setzero_si128();
	}

	// Story
	nrl_prime = nrl + PRE_D5_NROW + 1;
	for (row = nrl_prime; row < nrh + 1; row ++) {
		e3_row = row - (PRE_D5_NROW + 1);
		d5_row = row - (PRE_D5_NROW - 1);


		in = X; out = Z;
		ui8matrix_erosion_SIMD_InLU_O3_ValAddrRR  (in, e3_row, e3_row, ncl, nch, v0, v1, mid, out);
		in = X; out = Y;  
		ui8matrix_erosion_SIMD_InLU_O3_ValAddrRR  (in,    row,    row, ncl, nch, v0, v1, mid, out);
		
		in = Y; out = X;
		ui8matrix_dilation5_SIMD_InLU_O3_ValAddrRR(in, d5_row, d5_row, ncl, nch, v0, v1, mid, out);		
		// Handle edge
		out_row = out[d5_row];
		out_row[last_v] = _mm_and_si128(out_row[last_v], vMask);
		for (col = last_v + 1; col < v1 + 1; col++) out_row[col] = _mm_setzero_si128();
	}
	// Epilogue
	in = Y; out = X;
	nrl_prime = nrh - 2;
	ui8matrix_dilation5_SIMD_InLU_O3_ValAddrRR(in, nrl_prime, nrh, ncl, nch, v0, v1, mid, out);		
	// Handle edge
	for (d5_row = nrl_prime - 1; d5_row < nrh + 1; d5_row++) {
		out_row = out[d5_row];
		out_row[last_v] = _mm_and_si128(out_row[last_v], vMask);
		for (col = last_v + 1; col < v1 + 1; col++) out_row[col] = _mm_setzero_si128();
	}
	in = X; out = Z;
	nrl_prime = nrh - 4;
	ui8matrix_erosion_SIMD_InLU_O3_ValAddrRR(in, nrl_prime, nrh, ncl, nch, v0, v1, mid, out);	
}

/*------------------------------------------------------------------------------------------*/
void ui8matrix_erosion_SIMD_naive(vuint8** X, int nrl, int nrh, long ncl, long nch, int v0, int v1, vuint8 **vTempBuffer , vuint8 **Y) {
/*------------------------------------------------------------------------------------------*/

	long row, col, x, y, i;
	long snrl = SE_NRL, snrh = SE_NRH, sncl = SE_NCL, snch = SE_NCH;


	vuint8 v0_0, v0_1, v0_2; //X[row-1][col-1], X[row-1][col  ], X[row-1][col+1]
	vuint8 v1_0, v1_1, v1_2; //X[row  ][col-1], X[row  ][col  ], X[row  ][col+1]
	vuint8 v2_0, v2_1, v2_2; //X[row+1][col-1], X[row+1][col  ], X[row+1][col+1]

	vuint8 vY0, vY1, vY2;
	vuint8 vTMP;

	// Erode
	for (row = nrl; row <= nrh; row++) {
		for (col = v0; col <= v1; col++) {

			v0_0 = _mm_load_si128((vuint8*) &X[row-1][col-1]);
			v1_0 = _mm_load_si128((vuint8*) &X[row  ][col-1]);
			v2_0 = _mm_load_si128((vuint8*) &X[row+1][col-1]);

			v0_1 = _mm_load_si128((vuint8*) &X[row-1][col  ]);
			v1_1 = _mm_load_si128((vuint8*) &X[row  ][col  ]);
			v2_1 = _mm_load_si128((vuint8*) &X[row+1][col  ]);

			v0_2 = _mm_load_si128((vuint8*) &X[row-1][col+1]);
			v1_2 = _mm_load_si128((vuint8*) &X[row  ][col+1]);
			v2_2 = _mm_load_si128((vuint8*) &X[row+1][col+1]);


			//Column -1 "and" operator
			vY0 = vector_and3(v0_0, v1_0, v2_0);
			//Column  0 "and" operator 
			vY1  = vector_and3(v0_1, v1_1, v2_1);
			//Column  1 "and" operator
			vY2  = vector_and3(v0_2, v1_2, v2_2);

			//Row operator
			vTMP = vector_and3_row1shift(vY0, vY1, vY2);

			_mm_store_si128((vuint8*) &Y[row][col], vTMP);	
		}
	}

} 
/*------------------------------------------------------------------------------------------*/
void ui8matrix_erosion_SIMD_divide_row_and_conquer(vuint8** X, int nrl, int nrh, long ncl, long nch, int v0, int v1, vuint8 **vTempBuffer , vuint8 **Y) {
/*------------------------------------------------------------------------------------------*/

	long row, col;
	vuint8 vec0, vec1, vec2;
	vuint8 *in_row, *out_row, *mid_row;
	vuint8 *mid_row0, *mid_row1, *mid_row2;

	// Erode
	for (row = nrl - 1; row <= nrh + 1; row++) {
		in_row  = X[row];
		mid_row = vTempBuffer[row];
		for (col = v0; col <= v1; col++) {
			vec0 = _mm_load_si128((vuint8*) &in_row[col - 1]);
			vec1 = _mm_load_si128((vuint8*) &in_row[col + 0]);
			vec2 = _mm_load_si128((vuint8*) &in_row[col + 1]);
			_mm_store_si128((vuint8*)&mid_row[col], vector_and3_row1shift(vec0, vec1, vec2));
			vec0 = vec1;
			vec1 = vec2;
		} 
	}
	
	for (row = nrl; row <= nrh; row ++) {
		mid_row0 = vTempBuffer[row - 1];
		mid_row1 = vTempBuffer[row + 0];
		mid_row2 = vTempBuffer[row + 1];
		out_row  = Y[row];
		for (col = v0; col <= v1; col++) {
			vec0 = mid_row0[col];
			vec1 = mid_row1[col];
			vec2 = mid_row2[col];
			_mm_store_si128((vuint8*)&out_row[col], vector_and3(vec0, vec1, vec2));
		} 
	}

} 

void ui8matrix_erosion_SIMD_pipeline2_LU3x3_InLU_O3_RR (vuint8** X, int nrl, int nrh, long ncl, long nch, int v0, int v1, vuint8 **vTempBuffer , vuint8 **Y) 
{
	const long order = 3;
	long row = nrl, col = v0, x, y, r = 0;
	vuint8 *row0, *row1, *row2, *row3, x0, x1, x2, x3, x4, x5;
	vuint8 *temp_row0, *temp_row1, *temp_row2, *temp_row3, *temp_row4;
	vuint8 *out_row0, *out_row1, *out_row2;
	
	r = (v1 + 1) % order;

	// Prologue
	for (row = nrl; row < nrh + 1; row++) {
		row0 = X[row - 1];
		row1 = X[row + 0];
		row2 = X[row + 1];
		temp_row1 = vTempBuffer[row + 0];

		temp_row1[v0 - 1] = vector_and3(row0[v0 - 1],  row1[v0 - 1], row2[v0 - 1]);
		temp_row1[v0 + 0] = vector_and3(row0[v0 - 0],  row1[v0 - 0], row2[v0 - 0]);
	}

	for (row = nrl; row < nrh + 1; row ++){
		temp_row1 = vTempBuffer[row + 0];
		row0 = X[row - 1];
		row1 = X[row + 0];
		row2 = X[row + 1];

		out_row1 = Y[row]; 
		
		x0 = temp_row1[v0 - 1];		// LOAD => A
		x1 = temp_row1[v0 - 0];		// LOAD => B
		for (col = v0; col < v1 + 1 - r; col += order) {
			x2 = vector_and3(row0[col + 1], row1[col + 1], row2[col + 1]);
			x3 = vector_and3(row0[col + 2], row1[col + 2], row2[col + 2]);
			x4 = vector_and3(row0[col + 3], row1[col + 3], row2[col + 3]);
			
			out_row1[col + 0] = vector_and3_row1shift(x0, x1, x2); // A | B | C
			out_row1[col + 1] = vector_and3_row1shift(x1, x2, x3); // B | C | D
			out_row1[col + 2] = vector_and3_row1shift(x2, x3, x4); // C | D | E

			x0 = x3;
			x1 = x4;
		}
	}
	
	switch (r) {
		case 2: 
			row0 = X[nrl - 1];
			row1 = X[nrl + 0];
			for (row = nrl; row < nrh + 1; row++) {
				row2 = X[row + 1];
				out_row1 = Y[row]; 

				x0 = vector_and3(row0[v1 - 2], row1[v1 - 2], row2[v1 - 2]);
				x1 = vector_and3(row0[v1 - 1], row1[v1 - 1], row2[v1 - 1]);
				x2 = vector_and3(row0[v1 + 0], row1[v1 + 0], row2[v1 + 0]);
				x3 = vector_and3(row0[v1 + 1], row1[v1 + 1], row2[v1 + 2]);
				out_row1[v1 - 1] = vector_and3_row1shift(x0, x1, x2);
				out_row1[v1 + 0] = vector_and3_row1shift(x1, x2, x3);
				row0 = row1;
				row1 = row2;
			}
			break;
		case 1: 
			row0 = X[nrl - 1];
			row1 = X[nrl + 0];
			for (row = nrl; row < nrh + 1; row++) {
				row2 = X[row + 1];

				x0 = vector_and3(row0[v1 - 1], row1[v1 - 1], row2[v1 - 1]);
				x1 = vector_and3(row0[v1 - 0], row1[v1 + 0], row2[v1 - 0]);
				x2 = vector_and3(row0[v1 + 1], row1[v1 + 1], row2[v1 + 1]);
				Y[row][v1] = vector_and3_row1shift(x0, x1, x2);
				row0 = row1;
				row1 = row2;
			}
			break;
		default:
			break;
	}

}

/*------------------------------------------------------------------------------------------*/
void ui8matrix_erosion_SIMD_divide_col_and_conquer(vuint8** X, int nrl, int nrh, long ncl, long nch, int v0, int v1, vuint8 **vTempBuffer , vuint8 **Y) {
/*------------------------------------------------------------------------------------------*/

	long row, col;
	vuint8 vec0, vec1, vec2;
	vuint8 *in_row, *out_row, *mid_row;
	vuint8 *in_row0, *in_row1, *in_row2;

	// Erode
	for (row = nrl; row < nrh + 1; row ++){
		
		in_row0 = X[row - 1];
		in_row1 = X[row + 0];
		in_row2 = X[row + 1];
		mid_row = vTempBuffer[row];
		for (col = v0 - 1; col <= v1 + 1; col++) {
			vec0 = _mm_load_si128((vuint8*) &in_row0[col]);
			vec1 = _mm_load_si128((vuint8*) &in_row1[col]);
			vec2 = _mm_load_si128((vuint8*) &in_row2[col]);
			_mm_store_si128((vuint8*)&mid_row[col], vector_and3(vec0, vec1, vec2));
		} 
	}
	
	for (row = nrl; row < nrh + 1; row ++){
		mid_row  = vTempBuffer[row];
		out_row  = Y[row];
		for (col = v0; col <= v1; col++) {
			vec0 = mid_row[col - 1];
			vec1 = mid_row[col + 0];
			vec2 = mid_row[col + 1];
			_mm_store_si128((vuint8*)&out_row[col], vector_and3_row1shift(vec0, vec1, vec2));
		} 
	}

} 

/*------------------------------------------------------------------------------------*/
void ui8matrix_dilation_SIMD_naive(vuint8** X, int nrl, int nrh, long ncl, long nch, int v0, int v1, vuint8 **vTempBuffer , vuint8 **Y) {
/*------------------------------------------------------------------------------------*/

	long row, col, x, y, i;
	long snrl = SE_NRL, snrh = SE_NRH, sncl = SE_NCL, snch = SE_NCH;


	vuint8 v0_0, v0_1, v0_2; 
	vuint8 v1_0, v1_1, v1_2; 
	vuint8 v2_0, v2_1, v2_2; 

	vuint8 vY0, vY1, vY2;
	vuint8 vTMP;
	vuint8 vV0[3], *vX0 = &vV0[1];
	vuint8 vV1[3], *vX1 = &vV1[1];
	vuint8 vV2[3], *vX2 = &vV2[1];

	// Dilation
	for (row = nrl; row <= nrh; row++) {
		for (col = v0; col <= v1; col++) {

			
			v0_0 = _mm_load_si128((vuint8*) &X[row-1][col-1]);
			v1_0 = _mm_load_si128((vuint8*) &X[row  ][col-1]);
			v2_0 = _mm_load_si128((vuint8*) &X[row+1][col-1]);

			v0_1 = _mm_load_si128((vuint8*) &X[row-1][col  ]);
			v1_1 = _mm_load_si128((vuint8*) &X[row  ][col  ]);
			v2_1 = _mm_load_si128((vuint8*) &X[row+1][col  ]);

			v0_2 = _mm_load_si128((vuint8*) &X[row-1][col+1]);
			v1_2 = _mm_load_si128((vuint8*) &X[row  ][col+1]);
			v2_2 = _mm_load_si128((vuint8*) &X[row+1][col+1]);


			//Column -1 "or" operator
			vY0 = vector_or3(v0_0, v1_0, v2_0);
			//Column  0 "or" operator 
			vY1 = vector_or3(v0_1, v1_1, v2_1);
			//Column  1 "or" operator
			vY2 = vector_or3(v0_2, v1_2, v2_2);

			vTMP = vector_or3_row1shift(vY0, vY1, vY2);

			_mm_store_si128((vuint8*) &Y[row][col], vTMP);	
		}
	}

}

/*---------------------------------------------------------------------------------------------------------*/
void ui8matrix_erosion_SIMD_RR_row (vuint8** X, int nrl, int nrh, long ncl, long nch, int v0, int v1, vuint8 **vTempBuffer , vuint8 **Y) {
/*---------------------------------------------------------------------------------------------------------*/
	
	long row, col, x, y, i;
	long snrl = SE_NRL, snrh = SE_NRH, sncl = SE_NCL, snch = SE_NCH;
	
	vuint8 v0_0, v0_1, v0_2; //X[row-1][col-1], X[row-1][col  ], X[row-1][col+1]
	vuint8 v1_0, v1_1, v1_2; //X[row  ][col-1], X[row  ][col  ], X[row  ][col+1]
	vuint8 v2_0, v2_1, v2_2; //X[row+1][col-1], X[row+1][col  ], X[row+1][col+1]

	vuint8 vY0, vY1, vY2;
	vuint8 vTMP;

	// Erosion
	for (row = nrl; row <= nrh; row++) {

		col = 0;
		v0_0 = _mm_load_si128((vuint8*) &X[row-1][col-1]);
		v1_0 = _mm_load_si128((vuint8*) &X[row  ][col-1]);
		v2_0 = _mm_load_si128((vuint8*) &X[row+1][col-1]);

		v0_1 = _mm_load_si128((vuint8*) &X[row-1][col  ]);
		v1_1 = _mm_load_si128((vuint8*) &X[row  ][col  ]);
		v2_1 = _mm_load_si128((vuint8*) &X[row+1][col  ]);

		//Column -1 "and" operator
		vY0 = vector_and3(v0_0, v1_0, v2_0);
		//Column  0 "and" operator 
		vY1  = vector_and3(v0_1, v1_1, v2_1);

		for (col = v0; col <= v1; col++) {	

			v0_2 = _mm_load_si128((vuint8*) &X[row-1][col+1]);
			v1_2 = _mm_load_si128((vuint8*) &X[row  ][col+1]);
			v2_2 = _mm_load_si128((vuint8*) &X[row+1][col+1]);

			//Column  1 "and" operator
			vY2  = vector_and3(v0_2, v1_2, v2_2);

			//Row operator
			vTMP = vector_and3_row1shift(vY0, vY1, vY2);

			_mm_store_si128((vuint8*) &Y[row][col], vTMP);

			vY0 = vY1;
			vY1 = vY2;

		}
	}
}
void ui8matrix_erosion_SIMD_InLU_O3_AddrRR_OMP (vuint8** X, int nrl, int nrh, long ncl, long nch, int v0, int v1, vuint8 **vTempBuffer , vuint8 **Y) 
{
	const long order = 3;
	long row = nrl, col = v0, r;
	vuint8 x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, y0, y1, y2, y3, y4;
	vuint8 *row0, *row1, *row2, *row3, *row4;
	vuint8 *out_row0;


	r = (v1 + 1) % order;
	row0 = X[nrl - 1];
	row1 = X[nrl + 0];
	#pragma omp parallel default(none) firstprivate(row0, row1, row2, row3, row4, out_row0) private(row, col, x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, y0, y1, y2, y3, y4) shared(X, Y, nrl, nrh, ncl, nch, v0, v1, r)
	for (row = nrl; row < nrh + 1; row++) {
		row2 = X[row + 1];

		out_row0 = Y[row + 0];
		for (col = v0; col < v1 + 1 - r; col += order) {
			x0 = _mm_load_si128(&row0[col - 1]); x5  = _mm_load_si128(&row1[col - 1]); x10 = _mm_load_si128(&row2[col - 1]); 
			x1 = _mm_load_si128(&row0[col + 0]); x6  = _mm_load_si128(&row1[col + 0]); x11 = _mm_load_si128(&row2[col + 0]); 
			x2 = _mm_load_si128(&row0[col + 1]); x7  = _mm_load_si128(&row1[col + 1]); x12 = _mm_load_si128(&row2[col + 1]); 
			x3 = _mm_load_si128(&row0[col + 2]); x8  = _mm_load_si128(&row1[col + 2]); x13 = _mm_load_si128(&row2[col + 2]);
			x4 = _mm_load_si128(&row0[col + 3]); x9  = _mm_load_si128(&row1[col + 3]); x14 = _mm_load_si128(&row2[col + 3]);
			y0 = vector_and3(x0, x5, x10);
			y1 = vector_and3(x1, x6, x11);
			y2 = vector_and3(x2, x7, x12);
			y3 = vector_and3(x3, x8, x13);
			y4 = vector_and3(x4, x9, x14);

			out_row0[col + 0] = vector_and3_row1shift(y0, y1, y2);
			out_row0[col + 1] = vector_and3_row1shift(y1, y2, y3);
			out_row0[col + 2] = vector_and3_row1shift(y2, y3, y4);
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
		
			x0 = _mm_load_si128(&row0[v1 - 2]); x5  = _mm_load_si128(&row1[v1 - 2]); x10 = _mm_load_si128(&row2[v1 - 2]); 
			x1 = _mm_load_si128(&row0[v1 - 1]); x6  = _mm_load_si128(&row1[v1 - 1]); x11 = _mm_load_si128(&row2[v1 - 1]); 
			x2 = _mm_load_si128(&row0[v1 + 0]); x7  = _mm_load_si128(&row1[v1 + 0]); x12 = _mm_load_si128(&row2[v1 + 0]); 
			x3 = _mm_load_si128(&row0[v1 + 1]); x8  = _mm_load_si128(&row1[v1 + 1]); x13 = _mm_load_si128(&row2[v1 + 1]);
			y0 = vector_and3(x0, x5, x10);
			y1 = vector_and3(x1, x6, x11);
			y2 = vector_and3(x2, x7, x12);
			y3 = vector_and3(x3, x8, x13);
			out_row0[v1 - 1] = vector_and3_row1shift(y0, y1, y2);
			out_row0[v1 + 0] = vector_and3_row1shift(y1, y2, y3);
			
		}
		break;
	case 1:
		for (row = nrl; row < nrh + 1; row++) {
			row0 = X[row - 1];
			row1 = X[row + 0];
			row2 = X[row + 1];

			out_row0 = Y[row + 0];
			x0 = _mm_load_si128(&row0[v1 - 1]); x5  = _mm_load_si128(&row1[v1 - 1]); x10 = _mm_load_si128(&row2[v1 - 1]); 
			x1 = _mm_load_si128(&row0[v1 + 0]); x6  = _mm_load_si128(&row1[v1 + 0]); x11 = _mm_load_si128(&row2[v1 + 0]); 
			x2 = _mm_load_si128(&row0[v1 + 1]); x7  = _mm_load_si128(&row1[v1 + 1]); x12 = _mm_load_si128(&row2[v1 + 1]); 
			y0 = vector_and3(x0, x5, x10);
			y1 = vector_and3(x1, x6, x11);
			y2 = vector_and3(x2, x7, x12);
			out_row0[v1 + 0] = vector_and3_row1shift(y0, y1, y2);;
			
		}
		break;
	default:
		break;
	}
}
void ui8matrix_dilation5_SIMD_InLU_O3_ValAddrRR_OMP(vuint8** X, int nrl, int nrh, long ncl, long nch, int v0, int v1, vuint8 **vTempBuffer , vuint8 **Y) 
{
	const long order = 3;
	long row = nrl, col = v0, r;
	vuint8 x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16, x17, x18, x19, x20, x21, x22, x23, x24;
	vuint8 y0, y1, y2, y3, y4;
	vuint8 *row0, *row1, *row2, *row3, *row4;
	vuint8 *out_row0;

	r = (v1 + 1) % order;
	row0 = X[nrl - 2];
	row1 = X[nrl - 1];
	for (row = nrl; row < nrh + 1; row++) {
		col = v0;
		row2 = X[row + 0];
		row3 = X[row + 1];
		row4 = X[row + 2];
		out_row0 = Y[row + 0];
		
		x0 = _mm_load_si128(&row0[col - 1]); x5  = _mm_load_si128(&row1[col - 1]); x10 = _mm_load_si128(&row2[col - 1]); x15 = _mm_load_si128(&row3[col - 1]); x20 = _mm_load_si128(&row4[col - 1]);
		x1 = _mm_load_si128(&row0[col + 0]); x6  = _mm_load_si128(&row1[col + 0]); x11 = _mm_load_si128(&row2[col + 0]); x16 = _mm_load_si128(&row3[col + 0]); x21 = _mm_load_si128(&row4[col + 0]);
		y0 = vector_or5(x0, x5, x10, x15, x20);
		y1 = vector_or5(x1, x6, x11, x16, x21);
		for (col = v0; col < v1 + 1 - r; col += order) { 
			x2 = _mm_load_si128(&row0[col + 1]); x7  = _mm_load_si128(&row1[col + 1]); x12 = _mm_load_si128(&row2[col + 1]); x17 = _mm_load_si128(&row3[col + 1]); x22 = _mm_load_si128(&row4[col + 1]); 
			x3 = _mm_load_si128(&row0[col + 2]); x8  = _mm_load_si128(&row1[col + 2]); x13 = _mm_load_si128(&row2[col + 2]); x18 = _mm_load_si128(&row3[col + 2]); x23 = _mm_load_si128(&row4[col + 2]);
			x4 = _mm_load_si128(&row0[col + 3]); x9  = _mm_load_si128(&row1[col + 3]); x14 = _mm_load_si128(&row2[col + 3]); x19 = _mm_load_si128(&row3[col + 3]); x24 = _mm_load_si128(&row4[col + 3]);
			
			y2 = vector_or5(x2, x7, x12, x17, x22);
			y3 = vector_or5(x3, x8, x13, x18, x23);
			y4 = vector_or5(x4, x9, x14, x19, x24);

			out_row0[col + 0] = _mm_or_si128(vector_or3_row1shift(y0, y1, y2), vector_or3_row2shift(y0, y1, y2));
			out_row0[col + 1] = _mm_or_si128(vector_or3_row1shift(y1, y2, y3), vector_or3_row2shift(y1, y2, y3));
			out_row0[col + 2] = _mm_or_si128(vector_or3_row1shift(y2, y3, y4), vector_or3_row2shift(y2, y3, y4));
			y0 = y3; y1 = y4;

		}
		row0 = row1;
		row1 = row2;
	}
	switch(r) {
		case 2:
		for (row = nrl; row < nrh + 1; row++) {
			row0 = X[row - 2];
			row1 = X[row - 1];
			row2 = X[row + 0];
			row3 = X[row + 1];
			row4 = X[row + 2];

			out_row0 = Y[row + 0];
		
			x0 = _mm_load_si128(&row0[v1 - 2]); x5  = _mm_load_si128(&row1[v1 - 2]); x10 = _mm_load_si128(&row2[v1 - 2]); x15 = _mm_load_si128(&row3[v1 - 2]); x20 = _mm_load_si128(&row4[v1 - 2]);
			x1 = _mm_load_si128(&row0[v1 - 1]); x6  = _mm_load_si128(&row1[v1 - 1]); x11 = _mm_load_si128(&row2[v1 - 1]); x16 = _mm_load_si128(&row3[v1 - 1]); x21 = _mm_load_si128(&row4[v1 - 1]);
			x2 = _mm_load_si128(&row0[v1 + 0]); x7  = _mm_load_si128(&row1[v1 + 0]); x12 = _mm_load_si128(&row2[v1 + 0]); x17 = _mm_load_si128(&row3[v1 + 0]); x22 = _mm_load_si128(&row4[v1 + 0]);
			x3 = _mm_load_si128(&row0[v1 + 1]); x8  = _mm_load_si128(&row1[v1 + 1]); x13 = _mm_load_si128(&row2[v1 + 1]); x18 = _mm_load_si128(&row3[v1 + 1]); x23 = _mm_load_si128(&row4[v1 + 1]);
			y0 = vector_or5(x0, x5, x10, x15, x20);
			y1 = vector_or5(x1, x6, x11, x16, x21);
			y2 = vector_or5(x2, x7, x12, x17, x22);
			y3 = vector_or5(x3, x8, x13, x18, x23);
			out_row0[v1 - 1] = _mm_or_si128(vector_or3_row1shift(y0, y1, y2), vector_or3_row2shift(y0, y1, y2));
			out_row0[v1 + 0] = _mm_or_si128(vector_or3_row1shift(y1, y2, y3), vector_or3_row2shift(y1, y2, y3));
			
		}
		break;
	case 1:
		for (row = nrl; row < nrh + 1; row++) {
			row0 = X[row - 2];
			row1 = X[row - 1];
			row2 = X[row + 0];
			row3 = X[row + 1];
			row4 = X[row + 2];

			out_row0 = Y[row + 0];
			x0 = _mm_load_si128(&row0[v1 - 1]); x5  = _mm_load_si128(&row1[v1 - 1]); x10 = _mm_load_si128(&row2[v1 - 1]); x15 = _mm_load_si128(&row3[v1 - 1]); x20 = _mm_load_si128(&row4[v1 - 1]);
			x1 = _mm_load_si128(&row0[v1 + 0]); x6  = _mm_load_si128(&row1[v1 + 0]); x11 = _mm_load_si128(&row2[v1 + 0]); x16 = _mm_load_si128(&row3[v1 + 0]); x21 = _mm_load_si128(&row4[v1 + 0]);
			x2 = _mm_load_si128(&row0[v1 + 1]); x7  = _mm_load_si128(&row1[v1 + 1]); x12 = _mm_load_si128(&row2[v1 + 1]); x17 = _mm_load_si128(&row3[v1 + 1]); x22 = _mm_load_si128(&row4[v1 + 1]);
			y0 = vector_or5(x0, x5, x10, x15, x20);
			y1 = vector_or5(x1, x6, x11, x16, x21);
			y2 = vector_or5(x2, x7, x12, x17, x22);
			out_row0[v1 + 0] = _mm_or_si128(vector_or3_row1shift(y0, y1, y2), vector_or3_row2shift(y0, y1, y2));
			
		}
		break;
	default:
		break;
	}
}
void ui8matrix_erosion_SIMD_InLU_O3_AddrRR (vuint8** X, int nrl, int nrh, long ncl, long nch, int v0, int v1, vuint8 **vTempBuffer , vuint8 **Y) 
{
	const long order = 3;
	long row = nrl, col = v0, r;
	vuint8 x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, y0, y1, y2, y3, y4;
	vuint8 *row0, *row1, *row2, *row3, *row4;
	vuint8 *out_row0;

	r = (v1 + 1) % order;
	row0 = X[nrl - 1];
	row1 = X[nrl + 0];
	for (row = nrl; row < nrh + 1; row++) {
		row2 = X[row + 1];

		out_row0 = Y[row + 0];
		for (col = v0; col < v1 + 1 - r; col += order) {
			x0 = _mm_load_si128(&row0[col - 1]); x5  = _mm_load_si128(&row1[col - 1]); x10 = _mm_load_si128(&row2[col - 1]); 
			x1 = _mm_load_si128(&row0[col + 0]); x6  = _mm_load_si128(&row1[col + 0]); x11 = _mm_load_si128(&row2[col + 0]); 
			x2 = _mm_load_si128(&row0[col + 1]); x7  = _mm_load_si128(&row1[col + 1]); x12 = _mm_load_si128(&row2[col + 1]); 
			x3 = _mm_load_si128(&row0[col + 2]); x8  = _mm_load_si128(&row1[col + 2]); x13 = _mm_load_si128(&row2[col + 2]);
			x4 = _mm_load_si128(&row0[col + 3]); x9  = _mm_load_si128(&row1[col + 3]); x14 = _mm_load_si128(&row2[col + 3]);
			y0 = vector_and3(x0, x5, x10);
			y1 = vector_and3(x1, x6, x11);
			y2 = vector_and3(x2, x7, x12);
			y3 = vector_and3(x3, x8, x13);
			y4 = vector_and3(x4, x9, x14);

			out_row0[col + 0] = vector_and3_row1shift(y0, y1, y2);
			out_row0[col + 1] = vector_and3_row1shift(y1, y2, y3);
			out_row0[col + 2] = vector_and3_row1shift(y2, y3, y4);
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
		
			x0 = _mm_load_si128(&row0[v1 - 2]); x5  = _mm_load_si128(&row1[v1 - 2]); x10 = _mm_load_si128(&row2[v1 - 2]); 
			x1 = _mm_load_si128(&row0[v1 - 1]); x6  = _mm_load_si128(&row1[v1 - 1]); x11 = _mm_load_si128(&row2[v1 - 1]); 
			x2 = _mm_load_si128(&row0[v1 + 0]); x7  = _mm_load_si128(&row1[v1 + 0]); x12 = _mm_load_si128(&row2[v1 + 0]); 
			x3 = _mm_load_si128(&row0[v1 + 1]); x8  = _mm_load_si128(&row1[v1 + 1]); x13 = _mm_load_si128(&row2[v1 + 1]);
			y0 = vector_and3(x0, x5, x10);
			y1 = vector_and3(x1, x6, x11);
			y2 = vector_and3(x2, x7, x12);
			y3 = vector_and3(x3, x8, x13);
			out_row0[v1 - 1] = vector_and3_row1shift(y0, y1, y2);
			out_row0[v1 + 0] = vector_and3_row1shift(y1, y2, y3);
			
		}
		break;
	case 1:
		for (row = nrl; row < nrh + 1; row++) {
			row0 = X[row - 1];
			row1 = X[row + 0];
			row2 = X[row + 1];

			out_row0 = Y[row + 0];
			x0 = _mm_load_si128(&row0[v1 - 1]); x5  = _mm_load_si128(&row1[v1 - 1]); x10 = _mm_load_si128(&row2[v1 - 1]); 
			x1 = _mm_load_si128(&row0[v1 + 0]); x6  = _mm_load_si128(&row1[v1 + 0]); x11 = _mm_load_si128(&row2[v1 + 0]); 
			x2 = _mm_load_si128(&row0[v1 + 1]); x7  = _mm_load_si128(&row1[v1 + 1]); x12 = _mm_load_si128(&row2[v1 + 1]); 
			y0 = vector_and3(x0, x5, x10);
			y1 = vector_and3(x1, x6, x11);
			y2 = vector_and3(x2, x7, x12);
			out_row0[v1 + 0] = vector_and3_row1shift(y0, y1, y2);;
			
		}
		break;
	default:
		break;
	}
}

void ui8matrix_dilation_SIMD_InLU_O3_AddrRR (vuint8** X, int nrl, int nrh, long ncl, long nch, int v0, int v1, vuint8 **vTempBuffer , vuint8 **Y) 
{
	const long order = 3;
	long row = nrl, col = v0, r;
	vuint8 x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, y0, y1, y2, y3, y4;
	vuint8 *row0, *row1, *row2, *row3, *row4;
	vuint8 *out_row0;

	r = (v1 + 1) % order;
	row0 = X[nrl - 1];
	row1 = X[nrl + 0];
	for (row = nrl; row < nrh + 1; row++) {
		row2 = X[row + 1];

		out_row0 = Y[row + 0];
		for (col = v0; col < v1 + 1 - r; col += order) {
			x0 = _mm_load_si128(&row0[col - 1]); x5  = _mm_load_si128(&row1[col - 1]); x10 = _mm_load_si128(&row2[col - 1]); 
			x1 = _mm_load_si128(&row0[col + 0]); x6  = _mm_load_si128(&row1[col + 0]); x11 = _mm_load_si128(&row2[col + 0]); 
			x2 = _mm_load_si128(&row0[col + 1]); x7  = _mm_load_si128(&row1[col + 1]); x12 = _mm_load_si128(&row2[col + 1]); 
			x3 = _mm_load_si128(&row0[col + 2]); x8  = _mm_load_si128(&row1[col + 2]); x13 = _mm_load_si128(&row2[col + 2]);
			x4 = _mm_load_si128(&row0[col + 3]); x9  = _mm_load_si128(&row1[col + 3]); x14 = _mm_load_si128(&row2[col + 3]);
			y0 = vector_or3(x0, x5, x10);
			y1 = vector_or3(x1, x6, x11);
			y2 = vector_or3(x2, x7, x12);
			y3 = vector_or3(x3, x8, x13);
			y4 = vector_or3(x4, x9, x14);

			out_row0[col + 0] = vector_or3_row1shift(y0, y1, y2);
			out_row0[col + 1] = vector_or3_row1shift(y1, y2, y3);
			out_row0[col + 2] = vector_or3_row1shift(y2, y3, y4);
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
		
			x0 = _mm_load_si128(&row0[v1 - 2]); x5  = _mm_load_si128(&row1[v1 - 2]); x10 = _mm_load_si128(&row2[v1 - 2]); 
			x1 = _mm_load_si128(&row0[v1 - 1]); x6  = _mm_load_si128(&row1[v1 - 1]); x11 = _mm_load_si128(&row2[v1 - 1]); 
			x2 = _mm_load_si128(&row0[v1 + 0]); x7  = _mm_load_si128(&row1[v1 + 0]); x12 = _mm_load_si128(&row2[v1 + 0]); 
			x3 = _mm_load_si128(&row0[v1 + 1]); x8  = _mm_load_si128(&row1[v1 + 1]); x13 = _mm_load_si128(&row2[v1 + 1]);
			y0 = vector_or3(x0, x5, x10);
			y1 = vector_or3(x1, x6, x11);
			y2 = vector_or3(x2, x7, x12);
			y3 = vector_or3(x3, x8, x13);
			out_row0[v1 - 1] = vector_or3_row1shift(y0, y1, y2);
			out_row0[v1 + 0] = vector_or3_row1shift(y1, y2, y3);
			
		}
		break;
	case 1:
		for (row = nrl; row < nrh + 1; row++) {
			row0 = X[row - 1];
			row1 = X[row + 0];
			row2 = X[row + 1];

			out_row0 = Y[row + 0];
			x0 = _mm_load_si128(&row0[v1 - 1]); x5  = _mm_load_si128(&row1[v1 - 1]); x10 = _mm_load_si128(&row2[v1 - 1]); 
			x1 = _mm_load_si128(&row0[v1 + 0]); x6  = _mm_load_si128(&row1[v1 + 0]); x11 = _mm_load_si128(&row2[v1 + 0]); 
			x2 = _mm_load_si128(&row0[v1 + 1]); x7  = _mm_load_si128(&row1[v1 + 1]); x12 = _mm_load_si128(&row2[v1 + 1]); 
			y0 = vector_or3(x0, x5, x10);
			y1 = vector_or3(x1, x6, x11);
			y2 = vector_or3(x2, x7, x12);
			out_row0[v1 + 0] = vector_or3_row1shift(y0, y1, y2);;
			
		}
		break;
	default:
		break;
	}
}
void ui8matrix_dilation5_SIMD_InLU_O3_ValAddrRR (vuint8** X, int nrl, int nrh, long ncl, long nch, int v0, int v1, vuint8 **vTempBuffer , vuint8 **Y) 
{
	const long order = 3;
	long row = nrl, col = v0, r;
	vuint8 x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16, x17, x18, x19, x20, x21, x22, x23, x24;
	vuint8 y0, y1, y2, y3, y4;
	vuint8 *row0, *row1, *row2, *row3, *row4;
	vuint8 *out_row0;

	r = (v1 + 1) % order;
	row0 = X[nrl - 2];
	row1 = X[nrl - 1];
	for (row = nrl; row < nrh + 1; row++) {
		col = v0;
		row2 = X[row + 0];
		row3 = X[row + 1];
		row4 = X[row + 2];
		out_row0 = Y[row + 0];
		
		x0 = _mm_load_si128(&row0[col - 1]); x5  = _mm_load_si128(&row1[col - 1]); x10 = _mm_load_si128(&row2[col - 1]); x15 = _mm_load_si128(&row3[col - 1]); x20 = _mm_load_si128(&row4[col - 1]);
		x1 = _mm_load_si128(&row0[col + 0]); x6  = _mm_load_si128(&row1[col + 0]); x11 = _mm_load_si128(&row2[col + 0]); x16 = _mm_load_si128(&row3[col + 0]); x21 = _mm_load_si128(&row4[col + 0]);
		y0 = vector_or5(x0, x5, x10, x15, x20);
		y1 = vector_or5(x1, x6, x11, x16, x21);
		for (col = v0; col < v1 + 1 - r; col += order) { 
			x2 = _mm_load_si128(&row0[col + 1]); x7  = _mm_load_si128(&row1[col + 1]); x12 = _mm_load_si128(&row2[col + 1]); x17 = _mm_load_si128(&row3[col + 1]); x22 = _mm_load_si128(&row4[col + 1]); 
			x3 = _mm_load_si128(&row0[col + 2]); x8  = _mm_load_si128(&row1[col + 2]); x13 = _mm_load_si128(&row2[col + 2]); x18 = _mm_load_si128(&row3[col + 2]); x23 = _mm_load_si128(&row4[col + 2]);
			x4 = _mm_load_si128(&row0[col + 3]); x9  = _mm_load_si128(&row1[col + 3]); x14 = _mm_load_si128(&row2[col + 3]); x19 = _mm_load_si128(&row3[col + 3]); x24 = _mm_load_si128(&row4[col + 3]);
			
			y2 = vector_or5(x2, x7, x12, x17, x22);
			y3 = vector_or5(x3, x8, x13, x18, x23);
			y4 = vector_or5(x4, x9, x14, x19, x24);

			out_row0[col + 0] = _mm_or_si128(vector_or3_row1shift(y0, y1, y2), vector_or3_row2shift(y0, y1, y2));
			out_row0[col + 1] = _mm_or_si128(vector_or3_row1shift(y1, y2, y3), vector_or3_row2shift(y1, y2, y3));
			out_row0[col + 2] = _mm_or_si128(vector_or3_row1shift(y2, y3, y4), vector_or3_row2shift(y2, y3, y4));
			y0 = y3; y1 = y4;

		}
		row0 = row1;
		row1 = row2;
	}
	switch(r) {
		case 2:
		for (row = nrl; row < nrh + 1; row++) {
			row0 = X[row - 2];
			row1 = X[row - 1];
			row2 = X[row + 0];
			row3 = X[row + 1];
			row4 = X[row + 2];

			out_row0 = Y[row + 0];
		
			x0 = _mm_load_si128(&row0[v1 - 2]); x5  = _mm_load_si128(&row1[v1 - 2]); x10 = _mm_load_si128(&row2[v1 - 2]); x15 = _mm_load_si128(&row3[v1 - 2]); x20 = _mm_load_si128(&row4[v1 - 2]);
			x1 = _mm_load_si128(&row0[v1 - 1]); x6  = _mm_load_si128(&row1[v1 - 1]); x11 = _mm_load_si128(&row2[v1 - 1]); x16 = _mm_load_si128(&row3[v1 - 1]); x21 = _mm_load_si128(&row4[v1 - 1]);
			x2 = _mm_load_si128(&row0[v1 + 0]); x7  = _mm_load_si128(&row1[v1 + 0]); x12 = _mm_load_si128(&row2[v1 + 0]); x17 = _mm_load_si128(&row3[v1 + 0]); x22 = _mm_load_si128(&row4[v1 + 0]);
			x3 = _mm_load_si128(&row0[v1 + 1]); x8  = _mm_load_si128(&row1[v1 + 1]); x13 = _mm_load_si128(&row2[v1 + 1]); x18 = _mm_load_si128(&row3[v1 + 1]); x23 = _mm_load_si128(&row4[v1 + 1]);
			y0 = vector_or5(x0, x5, x10, x15, x20);
			y1 = vector_or5(x1, x6, x11, x16, x21);
			y2 = vector_or5(x2, x7, x12, x17, x22);
			y3 = vector_or5(x3, x8, x13, x18, x23);
			out_row0[v1 - 1] = _mm_or_si128(vector_or3_row1shift(y0, y1, y2), vector_or3_row2shift(y0, y1, y2));
			out_row0[v1 + 0] = _mm_or_si128(vector_or3_row1shift(y1, y2, y3), vector_or3_row2shift(y1, y2, y3));
			
		}
		break;
	case 1:
		for (row = nrl; row < nrh + 1; row++) {
			row0 = X[row - 2];
			row1 = X[row - 1];
			row2 = X[row + 0];
			row3 = X[row + 1];
			row4 = X[row + 2];

			out_row0 = Y[row + 0];
			x0 = _mm_load_si128(&row0[v1 - 1]); x5  = _mm_load_si128(&row1[v1 - 1]); x10 = _mm_load_si128(&row2[v1 - 1]); x15 = _mm_load_si128(&row3[v1 - 1]); x20 = _mm_load_si128(&row4[v1 - 1]);
			x1 = _mm_load_si128(&row0[v1 + 0]); x6  = _mm_load_si128(&row1[v1 + 0]); x11 = _mm_load_si128(&row2[v1 + 0]); x16 = _mm_load_si128(&row3[v1 + 0]); x21 = _mm_load_si128(&row4[v1 + 0]);
			x2 = _mm_load_si128(&row0[v1 + 1]); x7  = _mm_load_si128(&row1[v1 + 1]); x12 = _mm_load_si128(&row2[v1 + 1]); x17 = _mm_load_si128(&row3[v1 + 1]); x22 = _mm_load_si128(&row4[v1 + 1]);
			y0 = vector_or5(x0, x5, x10, x15, x20);
			y1 = vector_or5(x1, x6, x11, x16, x21);
			y2 = vector_or5(x2, x7, x12, x17, x22);
			out_row0[v1 + 0] = _mm_or_si128(vector_or3_row1shift(y0, y1, y2), vector_or3_row2shift(y0, y1, y2));
			
		}
		break;
	default:
		break;
	}
}
void ui8matrix_dilation_SIMD_InLU_O3_ValAddrRR (vuint8** X, int nrl, int nrh, long ncl, long nch, int v0, int v1, vuint8 **vTempBuffer , vuint8 **Y) 
{
	const long order = 3;
	long row = nrl, col = v0, r;
	vuint8 x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, y0, y1, y2, y3, y4;
	vuint8 *row0, *row1, *row2, *row3, *row4;
	vuint8 *out_row0;

	r = (v1 + 1) % order;
	row0 = X[nrl - 1];
	row1 = X[nrl + 0];
	for (row = nrl; row < nrh + 1; row++) {
		col = v0;
		row2 = X[row + 1];
		out_row0 = Y[row + 0];

		x0 = _mm_load_si128(&row0[col - 1]); x5  = _mm_load_si128(&row1[col - 1]); x10 = _mm_load_si128(&row2[col - 1]);
		x1 = _mm_load_si128(&row0[col + 0]); x6  = _mm_load_si128(&row1[col + 0]); x11 = _mm_load_si128(&row2[col + 0]); 
		y0 = vector_or3(x0, x5, x10);
		y1 = vector_or3(x1, x6, x11);
		for (col = v0; col < v1 + 1 - r; col += order) { 
			x2 = _mm_load_si128(&row0[col + 1]); x7  = _mm_load_si128(&row1[col + 1]); x12 = _mm_load_si128(&row2[col + 1]); 
			x3 = _mm_load_si128(&row0[col + 2]); x8  = _mm_load_si128(&row1[col + 2]); x13 = _mm_load_si128(&row2[col + 2]);
			x4 = _mm_load_si128(&row0[col + 3]); x9  = _mm_load_si128(&row1[col + 3]); x14 = _mm_load_si128(&row2[col + 3]);
			
			y2 = vector_or3(x2, x7, x12);
			y3 = vector_or3(x3, x8, x13);
			y4 = vector_or3(x4, x9, x14);

			out_row0[col + 0] = vector_or3_row1shift(y0, y1, y2);
			out_row0[col + 1] = vector_or3_row1shift(y1, y2, y3);
			out_row0[col + 2] = vector_or3_row1shift(y2, y3, y4);
			y0 = y3; y1 = y4;

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
		
			x0 = _mm_load_si128(&row0[v1 - 2]); x5  = _mm_load_si128(&row1[v1 - 2]); x10 = _mm_load_si128(&row2[v1 - 2]); 
			x1 = _mm_load_si128(&row0[v1 - 1]); x6  = _mm_load_si128(&row1[v1 - 1]); x11 = _mm_load_si128(&row2[v1 - 1]); 
			x2 = _mm_load_si128(&row0[v1 + 0]); x7  = _mm_load_si128(&row1[v1 + 0]); x12 = _mm_load_si128(&row2[v1 + 0]); 
			x3 = _mm_load_si128(&row0[v1 + 1]); x8  = _mm_load_si128(&row1[v1 + 1]); x13 = _mm_load_si128(&row2[v1 + 1]);
			y0 = vector_or3(x0, x5, x10);
			y1 = vector_or3(x1, x6, x11);
			y2 = vector_or3(x2, x7, x12);
			y3 = vector_or3(x3, x8, x13);
			out_row0[v1 - 1] = vector_or3_row1shift(y0, y1, y2);
			out_row0[v1 + 0] = vector_or3_row1shift(y1, y2, y3);
			
		}
		break;
	case 1:
		for (row = nrl; row < nrh + 1; row++) {
			row0 = X[row - 1];
			row1 = X[row + 0];
			row2 = X[row + 1];

			out_row0 = Y[row + 0];
			x0 = _mm_load_si128(&row0[v1 - 1]); x5  = _mm_load_si128(&row1[v1 - 1]); x10 = _mm_load_si128(&row2[v1 - 1]); 
			x1 = _mm_load_si128(&row0[v1 + 0]); x6  = _mm_load_si128(&row1[v1 + 0]); x11 = _mm_load_si128(&row2[v1 + 0]); 
			x2 = _mm_load_si128(&row0[v1 + 1]); x7  = _mm_load_si128(&row1[v1 + 1]); x12 = _mm_load_si128(&row2[v1 + 1]); 
			y0 = vector_or3(x0, x5, x10);
			y1 = vector_or3(x1, x6, x11);
			y2 = vector_or3(x2, x7, x12);
			out_row0[v1 + 0] = vector_or3_row1shift(y0, y1, y2);;
			
		}
		break;
	default:
		break;
	}
}
void ui8matrix_erosion_SIMD_InLU_O3_ValAddrRR (vuint8** X, int nrl, int nrh, long ncl, long nch, int v0, int v1, vuint8 **vTempBuffer , vuint8 **Y) 
{
	const long order = 3;
	long row = nrl, col = v0, r;
	vuint8 x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, y0, y1, y2, y3, y4;
	vuint8 *row0, *row1, *row2, *row3, *row4;
	vuint8 *out_row0;

	r = (v1 + 1) % order;
	row0 = X[nrl - 1];
	row1 = X[nrl + 0];
	for (row = nrl; row < nrh + 1; row++) {
		col = v0;
		row2 = X[row + 1];
		out_row0 = Y[row + 0];

		x0 = _mm_load_si128(&row0[col - 1]); x5  = _mm_load_si128(&row1[col - 1]); x10 = _mm_load_si128(&row2[col - 1]);
		x1 = _mm_load_si128(&row0[col + 0]); x6  = _mm_load_si128(&row1[col + 0]); x11 = _mm_load_si128(&row2[col + 0]); 
		y0 = vector_and3(x0, x5, x10);
		y1 = vector_and3(x1, x6, x11);
		for (col = v0; col < v1 + 1 - r; col += order) { 
			x2 = _mm_load_si128(&row0[col + 1]); x7  = _mm_load_si128(&row1[col + 1]); x12 = _mm_load_si128(&row2[col + 1]); 
			x3 = _mm_load_si128(&row0[col + 2]); x8  = _mm_load_si128(&row1[col + 2]); x13 = _mm_load_si128(&row2[col + 2]);
			x4 = _mm_load_si128(&row0[col + 3]); x9  = _mm_load_si128(&row1[col + 3]); x14 = _mm_load_si128(&row2[col + 3]);
			
			y2 = vector_and3(x2, x7, x12);
			y3 = vector_and3(x3, x8, x13);
			y4 = vector_and3(x4, x9, x14);

			out_row0[col + 0] = vector_and3_row1shift(y0, y1, y2);
			out_row0[col + 1] = vector_and3_row1shift(y1, y2, y3);
			out_row0[col + 2] = vector_and3_row1shift(y2, y3, y4);
			y0 = y3; y1 = y4;

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
		
			x0 = _mm_load_si128(&row0[v1 - 2]); x5  = _mm_load_si128(&row1[v1 - 2]); x10 = _mm_load_si128(&row2[v1 - 2]); 
			x1 = _mm_load_si128(&row0[v1 - 1]); x6  = _mm_load_si128(&row1[v1 - 1]); x11 = _mm_load_si128(&row2[v1 - 1]); 
			x2 = _mm_load_si128(&row0[v1 + 0]); x7  = _mm_load_si128(&row1[v1 + 0]); x12 = _mm_load_si128(&row2[v1 + 0]); 
			x3 = _mm_load_si128(&row0[v1 + 1]); x8  = _mm_load_si128(&row1[v1 + 1]); x13 = _mm_load_si128(&row2[v1 + 1]);
			y0 = vector_and3(x0, x5, x10);
			y1 = vector_and3(x1, x6, x11);
			y2 = vector_and3(x2, x7, x12);
			y3 = vector_and3(x3, x8, x13);
			out_row0[v1 - 1] = vector_and3_row1shift(y0, y1, y2);
			out_row0[v1 + 0] = vector_and3_row1shift(y1, y2, y3);
			
		}
		break;
	case 1:
		for (row = nrl; row < nrh + 1; row++) {
			row0 = X[row - 1];
			row1 = X[row + 0];
			row2 = X[row + 1];

			out_row0 = Y[row + 0];
			x0 = _mm_load_si128(&row0[v1 - 1]); x5  = _mm_load_si128(&row1[v1 - 1]); x10 = _mm_load_si128(&row2[v1 - 1]); 
			x1 = _mm_load_si128(&row0[v1 + 0]); x6  = _mm_load_si128(&row1[v1 + 0]); x11 = _mm_load_si128(&row2[v1 + 0]); 
			x2 = _mm_load_si128(&row0[v1 + 1]); x7  = _mm_load_si128(&row1[v1 + 1]); x12 = _mm_load_si128(&row2[v1 + 1]); 
			y0 = vector_and3(x0, x5, x10);
			y1 = vector_and3(x1, x6, x11);
			y2 = vector_and3(x2, x7, x12);
			out_row0[v1 + 0] = vector_and3_row1shift(y0, y1, y2);;
			
		}
		break;
	default:
		break;
	}
}
void ui8matrix_erosion_SIMD_col_pipeline(vuint8** X, int nrl, int nrh, long ncl, long nch, int v0, int v1, vuint8 **vTempBuffer , vuint8 **Y) 
{
	long row = nrl, col = v0, x, y;
	// uint8 **temp = ui8matrix(nrl + (-1), nrh + 1, v0 + (-1), v1 + 1);
	vuint8 *in_row0, *in_row1, *in_row2, *temp_row, *out_row0;
	vuint8 x0, x1, x2;
	// Prologue	
	for (row = nrl; row < nrh + 1; row++) {
		in_row0 = X[row - 1];
		in_row1 = X[row + 0];
		in_row2 = X[row + 1];
		temp_row = vTempBuffer[row];
    	temp_row[v0 - 1] = vector_and3(in_row0[v0 - 1], in_row1[v0 - 1], in_row2[v0 - 1]);
    	temp_row[v0 + 0] = vector_and3(in_row0[v0 + 0], in_row1[v0 + 0], in_row2[v0 + 0]);
	}	
	

	for (row = nrl; row < nrh + 1; row++){
		temp_row = vTempBuffer[row];
		
		in_row0 = X[row - 1];
		in_row1 = X[row + 0];
		in_row2 = X[row + 1];

		out_row0 = Y[row]; 
		
		for (col = v0; col < v1 + 1; col ++) {
			x0 = temp_row[col - 1];		// LOAD => A
			x1 = temp_row[col - 0];		// LOAD => B
			x2 = vector_and3(in_row0[col + 1], in_row1[col + 1], in_row2[col + 1]);
			
			out_row0[col + 0] = vector_and3_row1shift( x0, x1, x2); // A | B | C
			temp_row[col + 1] = x2;		  
		}
	}
}
void ui8matrix_erosion_SIMD_col_pipeline_RR(vuint8** X, int nrl, int nrh, long ncl, long nch, int v0, int v1, vuint8 **vTempBuffer , vuint8 **Y) 
{
	long row = nrl, col = v0, x, y;
	// uint8 **temp = ui8matrix(nrl + (-1), nrh + 1, v0 + (-1), v1 + 1);
	vuint8 *in_row0, *in_row1, *in_row2, *temp_row, *out_row0;
	vuint8 x0, x1, x2;
	// Prologue	
	for (row = nrl; row < nrh + 1; row++) {
		in_row0 = X[row - 1];
		in_row1 = X[row + 0];
		in_row2 = X[row + 1];
		temp_row = vTempBuffer[row];
    	temp_row[v0 - 1] = vector_and3(in_row0[v0 - 1], in_row1[v0 - 1], in_row2[v0 - 1]);
    	temp_row[v0 + 0] = vector_and3(in_row0[v0 + 0], in_row1[v0 + 0], in_row2[v0 + 0]);
	}	
	

	in_row0 = X[nrl - 1];
	in_row1 = X[nrl + 0];
	for (row = nrl; row < nrh + 1; row++){
		temp_row = vTempBuffer[row];
		
		in_row2 = X[row + 1];

		out_row0 = Y[row]; 
		
		x0 = temp_row[v0 - 1];		// LOAD => A
		x1 = temp_row[v0 - 0];		// LOAD => B
		for (col = v0; col < v1 + 1; col ++) {
			x2 = vector_and3(in_row0[col + 1], in_row1[col + 1], in_row2[col + 1]);
			out_row0[col + 0] = vector_and3_row1shift( x0, x1, x2); // A | B | C
			temp_row[col + 1] = x2;
			_mm_store_si128(&x0, x1);
			_mm_store_si128(&x1, x2);
		}
		in_row0 = in_row1;
		in_row1 = in_row2;

	}
}

/*-------------------------------------------------------------------------------------------------*/
void ui8matrix_dilation_SIMD_RR_row (vuint8** X, int nrl, int nrh, long ncl, long nch, int v0, int v1, vuint8 **vTempBuffer , vuint8 **Y) {
/*-------------------------------------------------------------------------------------------------*/
	
	long row, col, x, y, i;
	long snrl = SE_NRL, snrh = SE_NRH, sncl = SE_NCL, snch = SE_NCH;
	
	vuint8 v0_0, v0_1, v0_2; //X[row-1][col-1], X[row-1][col  ], X[row-1][col+1]
	vuint8 v1_0, v1_1, v1_2; //X[row  ][col-1], X[row  ][col  ], X[row  ][col+1]
	vuint8 v2_0, v2_1, v2_2; //X[row+1][col-1], X[row+1][col  ], X[row+1][col+1]

	vuint8 vY0, vY1, vY2;
	vuint8 vTMP;

	// Erosion
	for (row = nrl; row <= nrh; row++) {

		col = 0;
		v0_0 = _mm_load_si128((vuint8*) &X[row-1][col-1]);
		v1_0 = _mm_load_si128((vuint8*) &X[row  ][col-1]);
		v2_0 = _mm_load_si128((vuint8*) &X[row+1][col-1]);

		v0_1 = _mm_load_si128((vuint8*) &X[row-1][col  ]);
		v1_1 = _mm_load_si128((vuint8*) &X[row  ][col  ]);
		v2_1 = _mm_load_si128((vuint8*) &X[row+1][col  ]);

		//Column -1 "or" operator
		vY0 = vector_or3(v0_0, v1_0, v2_0);
		//Column  0 "or" operator 
		vY1  = vector_or3(v0_1, v1_1, v2_1);

		for (col = v0; col <= v1; col++) {	

			v0_2 = _mm_load_si128((vuint8*) &X[row-1][col+1]);
			v1_2 = _mm_load_si128((vuint8*) &X[row  ][col+1]);
			v2_2 = _mm_load_si128((vuint8*) &X[row+1][col+1]);

			//Column  1 "or" operator
			vY2  = vector_or3(v0_2, v1_2, v2_2);

			//Row operator
			vTMP = vector_or3_row1shift(vY0, vY1, vY2);

			_mm_store_si128((vuint8*) &Y[row][col], vTMP);

			vY0 = vY1;
			vY1 = vY2;

		}
	}

}

/*------------------------------------------------------------------------------------------*/
void ui8matrix_erosion_erosion_SIMD_FO(vuint8** X, int nrl, int nrh, long ncl, long nch, int v0, int v1, vuint8 **vTempBuffer , vuint8 **Y) {
/*------------------------------------------------------------------------------------------*/

	long row, col, x, y, i;
	long snrl = SE_NRL, snrh = SE_NRH, sncl = SE_NCL, snch = SE_NCH;

	vuint8  v0_0, v0_1, v0_2;
	vuint8  v1_0, v1_1, v1_2;
	vuint8  v2_0, v2_1, v2_2;
	vuint8  v3_0, v3_1, v3_2;
	vuint8  v4_0, v4_1, v4_2;

	vuint8 vY0, vY1, vY2;
	vuint8 vSHIFT1, vSHIFT2;


	// Erode
	for (row = nrl; row <= nrh; row++) {
		for (col = v0; col <= v1; col++) {

			v0_0 = _mm_load_si128((vuint8*) &X[row-2][col-1]);
			v1_0 = _mm_load_si128((vuint8*) &X[row-1][col-1]);
			v2_0 = _mm_load_si128((vuint8*) &X[row  ][col-1]);
			v3_0 = _mm_load_si128((vuint8*) &X[row+1][col-1]);
			v4_0 = _mm_load_si128((vuint8*) &X[row+2][col-1]);

			v0_1 = _mm_load_si128((vuint8*) &X[row-2][col  ]);
			v1_1 = _mm_load_si128((vuint8*) &X[row-1][col  ]);
			v2_1 = _mm_load_si128((vuint8*) &X[row  ][col  ]);
			v3_1 = _mm_load_si128((vuint8*) &X[row+1][col  ]);
			v4_1 = _mm_load_si128((vuint8*) &X[row+2][col  ]);

			v0_2 = _mm_load_si128((vuint8*) &X[row-2][col+1]);
			v1_2 = _mm_load_si128((vuint8*) &X[row-1][col+1]);
			v2_2 = _mm_load_si128((vuint8*) &X[row  ][col+1]);
			v3_2 = _mm_load_si128((vuint8*) &X[row+1][col+1]);
			v4_2 = _mm_load_si128((vuint8*) &X[row+2][col+1]);



			//Column -1 "and" operator
			vY0 = vector_and5(v0_0, v1_0, v2_0, v3_0, v4_0);
			//Column  0 "and" operator 
			vY1  = vector_and5(v0_1, v1_1, v2_1, v3_1, v4_1);
			//Column  1 "and" operator
			vY2  = vector_and5(v0_2, v1_2, v2_2, v3_2, v4_2);



			//Row operator
			vSHIFT1 = vector_and3_row1shift(vY0, vY1, vY2);
			vSHIFT2 = vector_and3_row2shift(vY0, vY1, vY2);
			vSHIFT1 = _mm_and_si128(vSHIFT1, vSHIFT2);
			_mm_store_si128((vuint8*) &Y[row][col], vSHIFT1);	
		}
	}

} 

/*------------------------------------------------------------------------------------------*/
void ui8matrix_dilation_dilation_SIMD_FO(vuint8** X, int nrl, int nrh, long ncl, long nch, int v0, int v1, vuint8 **vTempBuffer , vuint8 **Y) {
/*------------------------------------------------------------------------------------------*/

	long row, col, x, y, i;
	long snrl = SE_NRL, snrh = SE_NRH, sncl = SE_NCL, snch = SE_NCH;

	vuint8  v0_0, v0_1, v0_2;
	vuint8  v1_0, v1_1, v1_2;
	vuint8  v2_0, v2_1, v2_2;
	vuint8  v3_0, v3_1, v3_2;
	vuint8  v4_0, v4_1, v4_2;

	vuint8 vY0, vY1, vY2;
	vuint8 vSHIFT1, vSHIFT2;
	vuint8 vV0[5] = {_mm_setzero_si128(), _mm_setzero_si128(), _mm_setzero_si128(), _mm_setzero_si128(), _mm_setzero_si128()};
	vuint8 vV1[5] = {_mm_setzero_si128(), _mm_setzero_si128(), _mm_setzero_si128(), _mm_setzero_si128(), _mm_setzero_si128()};
	vuint8 vV2[5] = {_mm_setzero_si128(), _mm_setzero_si128(), _mm_setzero_si128(), _mm_setzero_si128(), _mm_setzero_si128()};
	vuint8 vV3[5] = {_mm_setzero_si128(), _mm_setzero_si128(), _mm_setzero_si128(), _mm_setzero_si128(), _mm_setzero_si128()};
	vuint8 vV4[5] = {_mm_setzero_si128(), _mm_setzero_si128(), _mm_setzero_si128(), _mm_setzero_si128(), _mm_setzero_si128()};
	vuint8 vX0;// = &vV0[2]; //X[row-1][col-1], X[row-1][col  ], X[row-1][col+1]
 	vuint8 vX1;// = &vV1[2]; //X[row  ][col-1], X[row  ][col  ], X[row  ][col+1]
 	vuint8 vX2;// = &vV2[2]; //X[row+1][col-1], X[row+1][col  ], X[row+1][col+1]
 	vuint8 vX3;// = &vV2[2]; //X[row+1][col-1], X[row+1][col  ], X[row+1][col+1]
 	vuint8 vX4;// = &vV2[2]; //X[row+1][col-1], X[row+1][col  ], X[row+1][col+1]


	// Dilation
	for (row = nrl; row <= nrh; row++) {
		for (col = v0; col <= v1; col++) {

			v0_0 = _mm_load_si128((vuint8*) &X[row-2][col-1]);
			v1_0 = _mm_load_si128((vuint8*) &X[row-1][col-1]);
			v2_0 = _mm_load_si128((vuint8*) &X[row  ][col-1]);
			v3_0 = _mm_load_si128((vuint8*) &X[row+1][col-1]);
			v4_0 = _mm_load_si128((vuint8*) &X[row+2][col-1]);

			v0_1 = _mm_load_si128((vuint8*) &X[row-2][col  ]);
			v1_1 = _mm_load_si128((vuint8*) &X[row-1][col  ]);
			v2_1 = _mm_load_si128((vuint8*) &X[row  ][col  ]);
			v3_1 = _mm_load_si128((vuint8*) &X[row+1][col  ]);
			v4_1 = _mm_load_si128((vuint8*) &X[row+2][col  ]);

			v0_2 = _mm_load_si128((vuint8*) &X[row-2][col+1]);
			v1_2 = _mm_load_si128((vuint8*) &X[row-1][col+1]);
			v2_2 = _mm_load_si128((vuint8*) &X[row  ][col+1]);
			v3_2 = _mm_load_si128((vuint8*) &X[row+1][col+1]);
			v4_2 = _mm_load_si128((vuint8*) &X[row+2][col+1]);
			// printf("\nSIMD Dilation FO\n");
			// //Column -1 "or" operator
			vY0 = vector_or5(v0_0, v1_0, v2_0, v3_0, v4_0);
			//Column  0 "or" operator 
			vY1 = vector_or5(v0_1, v1_1, v2_1, v3_1, v4_1);
			//Column  1 "or" operator
			vY2 = vector_or5(v0_2, v1_2, v2_2, v3_2, v4_2);
			// vV0[2] = vY0;
			// vV1[2] = vY1;
			// vV2[2] = vY2;
			
			// if (row == 9){
			// 	printf("Test %ld\n", col);
			// 	print_vui8vector(&v0_0, nrl, nrh, 0, 15, "%4u", NULL); printf(" "); print_vui8vector(&v0_1, nrl, nrh, 0, 15, "%4u", NULL); printf(" "); print_vui8vector(&v0_2, nrl, nrh, 0, 15, "%4u", NULL);printf("\n");
			// 	print_vui8vector(&v1_0, nrl, nrh, 0, 15, "%4u", NULL); printf(" "); print_vui8vector(&v1_1, nrl, nrh, 0, 15, "%4u", NULL); printf(" "); print_vui8vector(&v1_2, nrl, nrh, 0, 15, "%4u", NULL);printf("\n");
			// 	print_vui8vector(&v2_0, nrl, nrh, 0, 15, "%4u", NULL); printf(" "); print_vui8vector(&v2_1, nrl, nrh, 0, 15, "%4u", NULL); printf(" "); print_vui8vector(&v2_2, nrl, nrh, 0, 15, "%4u", NULL);printf("\n");
			// 	print_vui8vector(&v3_0, nrl, nrh, 0, 15, "%4u", NULL); printf(" "); print_vui8vector(&v3_1, nrl, nrh, 0, 15, "%4u", NULL); printf(" "); print_vui8vector(&v3_2, nrl, nrh, 0, 15, "%4u", NULL);printf("\n");
			// 	print_vui8vector(&v4_0, nrl, nrh, 0, 15, "%4u", NULL); printf(" "); print_vui8vector(&v4_1, nrl, nrh, 0, 15, "%4u", NULL); printf(" "); print_vui8vector(&v4_2, nrl, nrh, 0, 15, "%4u", NULL);printf("\n");

			// 	print_vui8vector(&vY0, nrl, nrh, 0, 15, "%4u", NULL);
			// 	print_vui8vector(&vY1, nrl, nrh, 0, 15, "%4u", NULL);
			// 	print_vui8vector(&vY2, nrl, nrh, 0, 15, "%4u", NULL);printf("\n");
			// }

			// printf("\n");
			
			//Row operator
			vSHIFT1 = vector_or3_row1shift(vY0, vY1, vY2);
			vSHIFT2 = vector_or3_row2shift(vY0, vY1, vY2);
			vSHIFT1 = _mm_or_si128(vSHIFT1, vSHIFT2);

			_mm_store_si128((vuint8*) &Y[row][col], vSHIFT1);	
		}
	}

} 

/*---------------------------------------------------------------------------------------------------------*/
void ui8matrix_erosion_erosion_SIMD_FO_RR_row (vuint8** X, int nrl, int nrh, long ncl, long nch, int v0, int v1, vuint8 **vTempBuffer , vuint8 **Y) {
/*---------------------------------------------------------------------------------------------------------*/
	
	long row, col, x, y, i;
	long snrl = SE_NRL, snrh = SE_NRH, sncl = SE_NCL, snch = SE_NCH;
	
	vuint8  v0_0, v0_1, v0_2;
	vuint8  v1_0, v1_1, v1_2;
	vuint8  v2_0, v2_1, v2_2;
	vuint8  v3_0, v3_1, v3_2;
	vuint8  v4_0, v4_1, v4_2;

	vuint8 vY0, vY1, vY2;
	vuint8 vSHIFT1, vSHIFT2;

	// Erosion
	for (row = nrl; row <= nrh; row++) {

		col = 0;
		v0_0 = _mm_load_si128((vuint8*) &X[row-2][col-1]);
		v1_0 = _mm_load_si128((vuint8*) &X[row-1][col-1]);
		v2_0 = _mm_load_si128((vuint8*) &X[row  ][col-1]);
		v3_0 = _mm_load_si128((vuint8*) &X[row+1][col-1]);
		v4_0 = _mm_load_si128((vuint8*) &X[row+2][col-1]);

		v0_1 = _mm_load_si128((vuint8*) &X[row-2][col  ]);
		v1_1 = _mm_load_si128((vuint8*) &X[row-1][col  ]);
		v2_1 = _mm_load_si128((vuint8*) &X[row  ][col  ]);
		v3_1 = _mm_load_si128((vuint8*) &X[row+1][col  ]);
		v4_1 = _mm_load_si128((vuint8*) &X[row+2][col  ]);

		//Column -1 "and" operator
		vY0 = vector_and5(v0_0, v1_0, v2_0, v3_0, v4_0);
		//Column  0 "and" operator 
		vY1  = vector_and5(v0_1, v1_1, v2_1, v3_1, v4_1);

		for (col = v0; col <= v1; col++) {	

			v0_2 = _mm_load_si128((vuint8*) &X[row-2][col+1]);
			v1_2 = _mm_load_si128((vuint8*) &X[row-1][col+1]);
			v2_2 = _mm_load_si128((vuint8*) &X[row  ][col+1]);
			v3_2 = _mm_load_si128((vuint8*) &X[row+1][col+1]);
			v4_2 = _mm_load_si128((vuint8*) &X[row+2][col+1]);

			//Column  1 "and" operator
			vY2  = vector_and5(v0_2, v1_2, v2_2, v3_2, v4_2);

			//Row operator
			vSHIFT1 = vector_and3_row1shift(vY0, vY1,   vY2);
			vSHIFT2 = vector_and3_row2shift(vY0, vY1,   vY2);
			vSHIFT1 = _mm_and_si128(vSHIFT1, vSHIFT2);
			_mm_store_si128((vuint8*) &Y[row][col], vSHIFT1);

			vY0 = vY1;
			vY1 = vY2;

		}
	}

}

/*-------------------------------------------------------------------------------------------------*/
void ui8matrix_dilation_dilation_SIMD_FO_RR_row (vuint8** X, int nrl, int nrh, long ncl, long nch, int v0, int v1, vuint8 **vTempBuffer , vuint8 **Y) {
/*-------------------------------------------------------------------------------------------------*/
	
	long row, col, x, y, i;
	long snrl = SE_NRL, snrh = SE_NRH, sncl = SE_NCL, snch = SE_NCH;
	
	vuint8  v0_0, v0_1, v0_2;
	vuint8  v1_0, v1_1, v1_2;
	vuint8  v2_0, v2_1, v2_2;
	vuint8  v3_0, v3_1, v3_2;
	vuint8  v4_0, v4_1, v4_2;

	vuint8 vY0, vY1, vY2;
	vuint8 vSHIFT1, vSHIFT2;

	// Dilation
	for (row = nrl; row <= nrh; row++) {

		col = 0;
		v0_0 = _mm_load_si128((vuint8*) &X[row-2][col-1]);
		v1_0 = _mm_load_si128((vuint8*) &X[row-1][col-1]);
		v2_0 = _mm_load_si128((vuint8*) &X[row  ][col-1]);
		v3_0 = _mm_load_si128((vuint8*) &X[row+1][col-1]);
		v4_0 = _mm_load_si128((vuint8*) &X[row+2][col-1]);

		v0_1 = _mm_load_si128((vuint8*) &X[row-2][col  ]);
		v1_1 = _mm_load_si128((vuint8*) &X[row-1][col  ]);
		v2_1 = _mm_load_si128((vuint8*) &X[row  ][col  ]);
		v3_1 = _mm_load_si128((vuint8*) &X[row+1][col  ]);
		v4_1 = _mm_load_si128((vuint8*) &X[row+2][col  ]);

		//Column -1 "or" operator
		vY0 = vector_or5(v0_0, v1_0, v2_0, v3_0, v4_0);
		//Column  0 "or" operator 
		vY1  = vector_or5(v0_1, v1_1, v2_1, v3_1, v4_1);

		for (col = v0; col <= v1; col++) {	

			v0_2 = _mm_load_si128((vuint8*) &X[row-2][col+1]);
			v1_2 = _mm_load_si128((vuint8*) &X[row-1][col+1]);
			v2_2 = _mm_load_si128((vuint8*) &X[row  ][col+1]);
			v3_2 = _mm_load_si128((vuint8*) &X[row+1][col+1]);
			v4_2 = _mm_load_si128((vuint8*) &X[row+2][col+1]);

			//Column  1 "or" operator
			vY2  = vector_or5(v0_2, v1_2, v2_2, v3_2, v4_2);

			//Row operator
			vSHIFT1 = vector_or3_row1shift(vY0, vY1,   vY2);
			vSHIFT2 = vector_or3_row2shift(vY0, vY1,   vY2);
			vSHIFT1 = _mm_or_si128(vSHIFT1, vSHIFT2);

			_mm_store_si128((vuint8*) &Y[row][col], vSHIFT1);

			vY0 = vY1;
			vY1 = vY2;

		}
	}

}

void test_functions_morpho_SIMD() {

	// p_image   t0 = create_image ("../car3/car_3000.pgm");
	// p_vimage vt0 = create_vimage("../car3/car_3000.pgm");
	// p_vimage vt1 = create_vimage("../car3/car_3000.pgm");
	// p_vimage vt2 = create_vimage("../car3/car_3000.pgm");
	// p_vimage vt3 = create_vimage("../car3/car_3000.pgm");

	// uint8** tmp;

	// ui8matrix_erosion_naive(t0->I, t0->nrl+BORD, t0->nrh-BORD, t0->ncl+BORD, t0->nch-BORD, tmp, t0->O);
	// ui8matrix_erosion_naive(t0->O, t0->nrl+BORD, t0->nrh-BORD, t0->ncl+BORD, t0->nch-BORD, tmp, t0->E);
	// ui8matrix_erosion_SIMD_naive(vt0->I, vt0->nrl+BORD, vt0->nrh-BORD, vt0->v0+vBORD, vt0->v1-vBORD, vt0->O);
	// ui8matrix_erosion_SIMD_naive(vt0->O, vt0->nrl+BORD, vt0->nrh-BORD, vt0->v0+vBORD, vt0->v1-vBORD, vt0->E);
	// ui8matrix_erosion_SIMD_RR_row(vt1->I, vt1->nrl+BORD, vt1->nrh-BORD, vt1->v0+vBORD, vt1->v1-vBORD, vt1->O);
	// ui8matrix_erosion_erosion_SIMD_FO(vt2->I, vt2->nrl+BORD, vt2->nrh-BORD, vt2->v0+vBORD, vt2->v1-vBORD, vt2->O);
	// ui8matrix_erosion_erosion_SIMD_FO_RR_row(vt3->I, vt3->nrl+BORD, vt3->nrh-BORD, vt3->v0+vBORD, vt3->v1-vBORD, vt3->O);


	// display_ui8vector((uint8*) t0->E[1], 0, 31, "%4d", "NAIVE O[1][0-31]");
	// display_ui8vector((uint8*) t0->E[2], 0, 31, "%4d", "NAIVE O[2][0-31]");
	// display_ui8vector((uint8*) t0->E[3], 0, 31, "%4d", "NAIVE O[3][0-31]");
	// puts(""); puts("");
	// display_ui8vector((uint8*) t0->I[1], 0, 31, "%4d", "NAIVE I[1][0 - 31]");
	// display_ui8vector((uint8*) t0->I[2], 0, 31, "%4d", "NAIVE I[2][0 - 31]");
	// display_ui8vector((uint8*) t0->I[3], 0, 31, "%4d", "NAIVE I[3][0 - 31]");
	
	// puts(""); puts("");
	// display_vui8vector(vt0->E[1], 0, 1, "%4d", "SIMD O[1][0 - 1]");
	// display_vui8vector(vt0->E[2], 0, 1, "%4d", "SIMD O[2][0 - 1]");
	// display_vui8vector(vt0->E[3], 0, 1, "%4d", "SIMD O[3][0 - 1]");
	// puts(""); puts("");

	// puts(""); puts("");
	// display_vui8vector(vt2->O[1], 0, 1, "%4d", "SIMD_FO O[1][0 - 1]");
	// display_vui8vector(vt2->O[2], 0, 1, "%4d", "SIMD_FO O[2][0 - 1]");
	// display_vui8vector(vt2->O[3], 0, 1, "%4d", "SIMD_FO O[3][0 - 1]");
	// puts(""); puts("");

	// puts(""); puts("");
	// display_vui8vector(vt3->O[1], 0, 1, "%4d", "SIMD_FO_RR O[1][0 - 1]");
	// display_vui8vector(vt3->O[2], 0, 1, "%4d", "SIMD_FO_RR O[2][0 - 1]");
	// display_vui8vector(vt3->O[3], 0, 1, "%4d", "SIMD_FO_RR O[3][0 - 1]");
	// puts(""); puts("");

	// display_vui8vector(vt3->I[1], 0, 1, "%4d", "SIMD_FO_RR I[1][0 - 1]");
	// display_vui8vector(vt3->I[2], 0, 1, "%4d", "SIMD_FO_RR I[2][0 - 1]");
	// display_vui8vector(vt3->I[3], 0, 1, "%4d", "SIMD_FO_RR I[3][0 - 1]");
	// puts(""); puts("");

}