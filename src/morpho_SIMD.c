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

#include "mutil.h"

#include "mymacro.h"
#include "simd_macro.h"

#include "util.h"
#include "img.h"
#include "img_SIMD.h"
#include "morpho.h"
#include "morpho_SIMD.h"

#define SE_NRL -1
#define SE_NRH  1
#define SE_NCL -1
#define SE_NCH  1
void print_vui8vector(vuint8 *vV, int nrl, int nrh, int ncl, int nch, char* format, char *name) {
	int col_cnt;

	if (name != NULL) printf("%s",name);
	
	uint8 *p;
	int i = 0;
	col_cnt = 0;
	for (int v = -1; v <= 1; v++)  {
		p = (uint8*)&vV[v];
		for(i=0; i<16; i++){
			if (16 + ncl <= col_cnt && col_cnt <= 16 + nch)
				printf(format, p[i]);
			col_cnt++;
		}
	}
	

}

/*------------------------------------------------------------------------------------------*/
void ui8matrix_erosion_SIMD_naive(vuint8** X, int nrl, int nrh, int v0, int v1, vuint8 **vTempBuffer , vuint8 **Y) {
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
void ui8matrix_erosion_SIMD_divide_row_and_conquer(vuint8** X, int nrl, int nrh, int v0, int v1, vuint8 **vTempBuffer , vuint8 **Y) {
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

void ui8matrix_erosion_SIMD_pipeline2_LU3x3_InLU_O3_RR (vuint8** X, int nrl, int nrh, int v0, int v1, vuint8 **vTempBuffer , vuint8 **Y) 
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

		_mm_store_si128(&temp_row1[v0 - 1], vector_and3(row0[v0 - 1],  row1[v0 - 1], row2[v0 - 1]));
		_mm_store_si128(&temp_row1[v0 + 0], vector_and3(row0[v0 - 0],  row1[v0 - 0], row2[v0 - 0]));
	}

	for (row = nrl; row < nrh + 1; row ++){
		temp_row1 = vTempBuffer[row + 0];
		row0 = X[row - 1];
		row1 = X[row + 0];
		row2 = X[row + 1];

		out_row1 = Y[row]; 
		
		_mm_store_si128(&x0, temp_row1[v0 - 1]);		// LOAD => A
		_mm_store_si128(&x1, temp_row1[v0 - 0]);		// LOAD => B
		for (col = v0; col < v1 + 1 - r; col += order) {
			_mm_store_si128(&x2, vector_and3(row0[col + 1], row1[col + 1], row2[col + 1]));
			_mm_store_si128(&x3, vector_and3(row0[col + 2], row1[col + 2], row2[col + 2]));
			_mm_store_si128(&x4, vector_and3(row0[col + 3], row1[col + 3], row2[col + 3]));
			
			_mm_store_si128(&out_row1[col + 0], vector_and3_row1shift(x0, x1, x2)); // A | B | C
			_mm_store_si128(&out_row1[col + 1], vector_and3_row1shift(x1, x2, x3)); // B | C | D
			_mm_store_si128(&out_row1[col + 2], vector_and3_row1shift(x2, x3, x4)); // C | D | E

			_mm_store_si128(&x0, x3);
			_mm_store_si128(&x1, x4);
		}
	}
	
	switch (r) {
		case 2: 
			row0 = X[nrl - 1];
			row1 = X[nrl + 0];
			for (row = nrl; row < nrh + 1; row++) {
				row2 = X[row + 1];
				out_row1 = Y[row]; 

				_mm_store_si128(&x0, vector_and3(row0[v1 - 2], row1[v1 - 2], row2[v1 - 2]));
				_mm_store_si128(&x1, vector_and3(row0[v1 - 1], row1[v1 - 1], row2[v1 - 1]));
				_mm_store_si128(&x2, vector_and3(row0[v1 + 0], row1[v1 + 0], row2[v1 + 0]));
				_mm_store_si128(&x3, vector_and3(row0[v1 + 1], row1[v1 + 1], row2[v1 + 2]));
				_mm_store_si128(&out_row1[v1 - 1], vector_and3_row1shift(x0, x1, x2));
				_mm_store_si128(&out_row1[v1 + 0], vector_and3_row1shift(x1, x2, x3));
				row0 = row1;
				row1 = row2;
			}
			break;
		case 1: 
			row0 = X[nrl - 1];
			row1 = X[nrl + 0];
			for (row = nrl; row < nrh + 1; row++) {
				row2 = X[row + 1];

				_mm_store_si128(&x0, vector_and3(row0[v1 - 1], row1[v1 - 1], row2[v1 - 1]));
				_mm_store_si128(&x1, vector_and3(row0[v1 - 0], row1[v1 + 0], row2[v1 - 0]));
				_mm_store_si128(&x2, vector_and3(row0[v1 + 1], row1[v1 + 1], row2[v1 + 1]));
				_mm_store_si128(&Y[row][v1], vector_and3_row1shift(x0, x1, x2));
				row0 = row1;
				row1 = row2;
			}
			break;
		default:
			break;
	}

}
/*------------------------------------------------------------------------------------------*/
void ui8matrix_erosion_SIMD_divide_col_and_conquer(vuint8** X, int nrl, int nrh, int v0, int v1, vuint8 **vTempBuffer , vuint8 **Y) {
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
void ui8matrix_dilation_SIMD_naive(vuint8** X, int nrl, int nrh, int v0, int v1, vuint8 **vTempBuffer , vuint8 **Y) {
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
void ui8matrix_erosion_SIMD_RR_row (vuint8** X, int nrl, int nrh, int v0, int v1, vuint8 **vTempBuffer , vuint8 **Y) {
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

/*---------------------------------------------------------------------------------------------------------*/
void ui8matrix_erosion_SIMD_InLU_RR_row (vuint8** X, int nrl, int nrh, int v0, int v1, vuint8 **vTempBuffer , vuint8 **Y) {
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

/*-------------------------------------------------------------------------------------------------*/
void ui8matrix_dilation_SIMD_RR_row (vuint8** X, int nrl, int nrh, int v0, int v1, vuint8 **vTempBuffer , vuint8 **Y) {
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
void ui8matrix_erosion_erosion_SIMD_FO(vuint8** X, int nrl, int nrh, int v0, int v1, vuint8 **vTempBuffer , vuint8 **Y) {
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
void ui8matrix_dilation_dilation_SIMD_FO(vuint8** X, int nrl, int nrh, int v0, int v1, vuint8 **vTempBuffer , vuint8 **Y) {
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
void ui8matrix_erosion_erosion_SIMD_FO_RR_row (vuint8** X, int nrl, int nrh, int v0, int v1, vuint8 **vTempBuffer , vuint8 **Y) {
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
void ui8matrix_dilation_dilation_SIMD_FO_RR_row (vuint8** X, int nrl, int nrh, int v0, int v1, vuint8 **vTempBuffer , vuint8 **Y) {
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