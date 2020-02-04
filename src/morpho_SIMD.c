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


/*------------------------------------------------------------------------------------------*/
void ui8matrix_erosion_SIMD_naive(vuint8** X, int nrl, int nrh, int v0, int v1, vuint8 **Y) {
/*------------------------------------------------------------------------------------------*/

	long row, col, x, y, i;
	long snrl = SE_NRL, snrh = SE_NRH, sncl = SE_NCL, snch = SE_NCH;


	vuint8 v_1X_1, v_1X, v_1X1; //X[row-1][col-1], X[row-1][col  ], X[row-1][col+1]
	vuint8 vX_1  , vX  , vX1  ; //X[row  ][col-1], X[row  ][col  ], X[row  ][col+1]
	vuint8 v1X_1 , v1X , v1X1 ; //X[row+1][col-1], X[row+1][col  ], X[row+1][col+1]

	vuint8 vY_1, vY, vY1;
	vuint8 vTMP;

	// Erode
	for (row = nrl; row <= nrh; row++) {
		for (col = v0; col <= v1; col++) {

			v_1X_1 	= _mm_load_si128((vuint8*) &X[row-1][col-1]);
			vX_1	= _mm_load_si128((vuint8*) &X[row  ][col-1]);
			v1X_1 	= _mm_load_si128((vuint8*) &X[row+1][col-1]);
			v_1X   	= _mm_load_si128((vuint8*) &X[row-1][col  ]);
			vX   	= _mm_load_si128((vuint8*) &X[row  ][col  ]);
			v1X   	= _mm_load_si128((vuint8*) &X[row+1][col  ]);
			v_1X1  	= _mm_load_si128((vuint8*) &X[row-1][col+1]);
			vX1  	= _mm_load_si128((vuint8*) &X[row  ][col+1]);
			v1X1  	= _mm_load_si128((vuint8*) &X[row+1][col+1]);


			//Column -1 "and" operator
			vY_1 = vector_and3(v_1X_1, vX_1, v1X_1);
			//Column  0 "and" operator 
			vY   = vector_and3(v_1X,   vX,   v1X);
			//Column  1 "and" operator
			vY1  = vector_and3(v_1X1,  vX1,  v1X1);

			//Row operator
			vTMP = vector_and3_row1shift(vY_1, vY, vY1);
			/*if(row == 3 && col == 1) {
				printf("SIMD\n");
    			display_vuint8(v_1X_1, "%4d", NULL); printf("\n");
    			display_vuint8(vX_1, "%4d", NULL); printf("\n");
    			display_vuint8(v1X_1, "%4d", NULL); printf("\n");
    			display_vuint8(vY_1, "%4d", NULL); printf("\n\n");
    			display_vuint8(v_1X, "%4d", NULL); printf("\n");
    			display_vuint8(vX, "%4d", NULL); printf("\n");
    			display_vuint8(v1X, "%4d", NULL); printf("\n");
    			display_vuint8(vY, "%4d", NULL); printf("\n\n");
    			display_vuint8(v_1X1, "%4d", NULL); printf("\n");
    			display_vuint8(vX1, "%4d", NULL); printf("\n");
    			display_vuint8(v1X1, "%4d", NULL); printf("\n");
    			display_vuint8(vY1, "%4d", NULL); printf("\n\n");

    			display_vuint8(vec_right1(vY_1, vY), "%4d", NULL); printf("\n");
    			display_vuint8(vY, "%4d", NULL); printf("\n");
    			display_vuint8(vec_left1(vY, vY1), "%4d", NULL); printf("\n");
    			display_vuint8(vTMP, "%4d", NULL); printf("\n\n");
    		}*/

			_mm_store_si128((vuint8*) &Y[row][col], vTMP);	
		}
	}

} 

/*------------------------------------------------------------------------------------*/
void ui8matrix_dilation_SIMD_naive(vuint8** X, int nrl, int nrh, int v0, int v1, vuint8 **Y) {
/*------------------------------------------------------------------------------------*/

	long row, col, x, y, i;
	long snrl = SE_NRL, snrh = SE_NRH, sncl = SE_NCL, snch = SE_NCH;

	vuint8 v_1X_1, v_1X, v_1X1;
	vuint8 vX_1  , vX  , vX1  ;
	vuint8 v1X_1 , v1X , v1X1 ;

	vuint8 vY_1, vY, vY1;
	vuint8 vTMP;

	// Dilation
	for (row = nrl; row <= nrh; row++) {
		for (col = v0; col <= v1; col++) {

			v_1X_1 	= _mm_load_si128((vuint8*) &X[row-1][col-1]);
			vX_1	= _mm_load_si128((vuint8*) &X[row  ][col-1]);
			v1X_1 	= _mm_load_si128((vuint8*) &X[row+1][col-1]);
			v_1X   	= _mm_load_si128((vuint8*) &X[row-1][col  ]);
			vX   	= _mm_load_si128((vuint8*) &X[row  ][col  ]);
			v1X   	= _mm_load_si128((vuint8*) &X[row+1][col  ]);
			v_1X1  	= _mm_load_si128((vuint8*) &X[row-1][col+1]);
			vX1  	= _mm_load_si128((vuint8*) &X[row  ][col+1]);
			v1X1  	= _mm_load_si128((vuint8*) &X[row+1][col+1]);


			//Column -1 "and" operator
			vY_1 = vector_or3(v_1X_1, vX_1, v1X_1);
			//Column  0 "and" operator 
			vY   = vector_or3(v_1X,   vX,   v1X);
			//Column  1 "and" operator
			vY1  = vector_or3(v_1X1,  vX1,  v1X1);

			//Row operator
			vTMP = vector_or3_row1shift(vY_1, vY, vY1);

			_mm_store_si128((vuint8*) &Y[row][col], vTMP);		
		}
	}
}

/*---------------------------------------------------------------------------------------------------------*/
void ui8matrix_erosion_SIMD_RR_row (vuint8** X, int nrl, int nrh, int v0, int v1, vuint8 **Y) {
/*---------------------------------------------------------------------------------------------------------*/
	
	long row, col, x, y, i;
	long snrl = SE_NRL, snrh = SE_NRH, sncl = SE_NCL, snch = SE_NCH;
	
	vuint8 v_1X_1, v_1X, v_1X1;
	vuint8 vX_1  , vX  , vX1  ;
	vuint8 v1X_1 , v1X , v1X1 ;

	vuint8 vY_1, vY, vY1;
	vuint8 vTMP;

	// Erosion
	for (row = nrl; row <= nrh; row++) {

		col = 0;
		v_1X_1 	= _mm_load_si128((vuint8*) &X[row-1][col-1]);
		vX_1	= _mm_load_si128((vuint8*) &X[row  ][col-1]);
		v1X_1 	= _mm_load_si128((vuint8*) &X[row+1][col-1]);
		v_1X   	= _mm_load_si128((vuint8*) &X[row-1][col  ]);
		vX   	= _mm_load_si128((vuint8*) &X[row  ][col  ]);
		v1X   	= _mm_load_si128((vuint8*) &X[row+1][col  ]);

		//Column -1 "and" operator
		vY_1 = vector_and3(v_1X_1, vX_1, v1X_1);
		//Column  0 "and" operator 
		vY   = vector_and3(v_1X,   vX,   v1X);	

		for (col = v0; col <= v1; col++) {	

			v_1X1  	= _mm_load_si128((vuint8*) &X[row-1][col+1]);
			vX1  	= _mm_load_si128((vuint8*) &X[row  ][col+1]);
			v1X1  	= _mm_load_si128((vuint8*) &X[row+1][col+1]);

			//Column  1 "and" operator
			vY1  = vector_and3(v_1X1,  vX1,  v1X1);

			//Row operator
			vTMP = vector_and3_row1shift(vY_1, vY, vY1);

			_mm_store_si128((vuint8*) &Y[row][col], vTMP);

			vY_1 = vY;
			vY = vY1;

		}
	}

}

/*-------------------------------------------------------------------------------------------------*/
void ui8matrix_dilation_SIMD_RR_row (vuint8** X, int nrl, int nrh, int v0, int v1, vuint8 **Y) {
/*-------------------------------------------------------------------------------------------------*/
	
	long row, col, x, y, i;
	long snrl = SE_NRL, snrh = SE_NRH, sncl = SE_NCL, snch = SE_NCH;
	
	vuint8 v_1X_1, v_1X, v_1X1;
	vuint8 vX_1  , vX  , vX1  ;
	vuint8 v1X_1 , v1X , v1X1 ;

	vuint8 vY_1, vY, vY1;
	vuint8 vTMP;

	// Dilation
	for (row = nrl; row <= nrh; row++) {

		col = 0;
		v_1X_1 	= _mm_load_si128((vuint8*) &X[row-1][col-1]);
		vX_1	= _mm_load_si128((vuint8*) &X[row  ][col-1]);
		v1X_1 	= _mm_load_si128((vuint8*) &X[row+1][col-1]);
		v_1X   	= _mm_load_si128((vuint8*) &X[row-1][col  ]);
		vX   	= _mm_load_si128((vuint8*) &X[row  ][col  ]);
		v1X   	= _mm_load_si128((vuint8*) &X[row+1][col  ]);

		//Column -1 "and" operator
		vY_1 = vector_or3(v_1X_1, vX_1, v1X_1);
		//Column  0 "and" operator 
		vY   = vector_or3(v_1X,   vX,   v1X);

		for (col = v0; col <= v1; col++) {	

			v_1X1  	= _mm_load_si128((vuint8*) &X[row-1][col+1]);
			vX1  	= _mm_load_si128((vuint8*) &X[row  ][col+1]);
			v1X1  	= _mm_load_si128((vuint8*) &X[row+1][col+1]);

			//Column  1 "and" operator
			vY1  = vector_or3(v_1X1,  vX1,  v1X1);

			//Row operator
			vTMP = vector_or3_row1shift(vY_1, vY, vY1);

			_mm_store_si128((vuint8*) &Y[row][col], vTMP);

			vY_1 = vY;
			vY = vY1;

		}
	}

}

/*------------------------------------------------------------------------------------------*/
void ui8matrix_erosion_erosion_SIMD_FO(vuint8** X, long nrl, long nrh, long v0, long v1, vuint8 **Y) {
/*------------------------------------------------------------------------------------------*/

	long row, col, x, y, i;
	long snrl = SE_NRL, snrh = SE_NRH, sncl = SE_NCL, snch = SE_NCH;

	vuint8 v_2X_1, v_2X, v_2X1;
	vuint8 v_1X_1, v_1X, v_1X1;
	vuint8   vX_1,   vX,   vX1;
	vuint8  v1X_1,  v1X,  v1X1;
	vuint8  v2X_1,  v2X,  v2X1;

	vuint8 vY_1, vY, vY1;
	vuint8 vSHIFT1, vSHIFT2;


	// Erode
	for (row = nrl; row <= nrh; row++) {
		for (col = v0; col <= v1; col++) {

			v_2X_1  = _mm_load_si128((vuint8*) &X[row-2][col-1]);
			v_1X_1 	= _mm_load_si128((vuint8*) &X[row-1][col-1]);
			vX_1	= _mm_load_si128((vuint8*) &X[row  ][col-1]);
			v1X_1 	= _mm_load_si128((vuint8*) &X[row+1][col-1]);
			v2X_1	= _mm_load_si128((vuint8*) &X[row+2][col-1]);
			v_2X    = _mm_load_si128((vuint8*) &X[row-2][col  ]);
			v_1X   	= _mm_load_si128((vuint8*) &X[row-1][col  ]);
			vX   	= _mm_load_si128((vuint8*) &X[row  ][col  ]);
			v1X   	= _mm_load_si128((vuint8*) &X[row+1][col  ]);
			v2X 	= _mm_load_si128((vuint8*) &X[row+2][col  ]);
			v_2X1   = _mm_load_si128((vuint8*) &X[row-2][col+1]);
			v_1X1  	= _mm_load_si128((vuint8*) &X[row-1][col+1]);
			vX1  	= _mm_load_si128((vuint8*) &X[row  ][col+1]);
			v1X1  	= _mm_load_si128((vuint8*) &X[row+1][col+1]);
			v2X1    = _mm_load_si128((vuint8*) &X[row+2][col+1]);



			//Column -1 "and" operator
			vY_1 = vector_and5(v_2X_1, v_1X_1, vX_1, v1X_1, v2X_1);
			//Column  0 "and" operator 
			vY   = vector_and5(v_2X,   v_1X,   vX,   v1X,   v2X);
			//Column  1 "and" operator
			vY1  = vector_and5(v_2X1,  v_1X1,  vX1,  v1X1,  v2X1);



			//Row operator
			vSHIFT1 = vector_and3_row1shift(vY_1, vY,   vY1);
			/*if(row == 3 && col == 1) {
				printf("SIMD_FO\n");
    			display_vuint8(v_2X_1, "%4d", NULL); printf("\n");
    			display_vuint8(v_1X_1, "%4d", NULL); printf("\n");
    			display_vuint8(vX_1, "%4d", NULL); printf("\n");
    			display_vuint8(v1X_1, "%4d", NULL); printf("\n");
    			display_vuint8(v2X_1, "%4d", NULL); printf("\n\n");
    			display_vuint8(vY_1, "%4d", NULL); printf("\n\n");
    			display_vuint8(v_2X, "%4d", NULL); printf("\n");
    			display_vuint8(v_1X, "%4d", NULL); printf("\n");
    			display_vuint8(vX, "%4d", NULL); printf("\n");
    			display_vuint8(v1X, "%4d", NULL); printf("\n");
    			display_vuint8(v2X, "%4d", NULL); printf("\n\n");
    			display_vuint8(vY, "%4d", NULL); printf("\n\n");
    			display_vuint8(v_2X1, "%4d", NULL); printf("\n");
    			display_vuint8(v_1X1, "%4d", NULL); printf("\n");
    			display_vuint8(vX1, "%4d", NULL); printf("\n");
    			display_vuint8(v1X1, "%4d", NULL); printf("\n");
    			display_vuint8(v2X1, "%4d", NULL); printf("\n\n");
    			display_vuint8(vY1, "%4d", NULL); printf("\n\n");

    			display_vuint8(vec_right1(vY_1, vY), "%4d", NULL); printf("\n");
    			display_vuint8(vY, "%4d", NULL); printf("\n");
    			display_vuint8(vec_left1(vY, vY1), "%4d", NULL); printf("\n");
    			display_vuint8(vSHIFT1, "%4d", NULL); printf("\n\n");
    		}*/
			vSHIFT2 = vector_and3_row2shift(vY_1, vY, vY1);
			/*if(row == 3 && col == 1) {

	   			display_vuint8(vec_right2(vY_1, vY_1), "%4d", NULL); printf("\n");
    			display_vuint8(vY, "%4d", NULL); printf("\n");
    			display_vuint8(vec_left2(vY, vY1), "%4d", NULL); printf("\n");
    			display_vuint8(vSHIFT2, "%4d", NULL); printf("\n\n");
    		}*/
			vSHIFT1 = _mm_and_si128(vSHIFT1, vSHIFT2);
			/*if(row == 3 && col == 1) {
    			display_vuint8(vSHIFT1, "%4d", NULL); printf("\n\n");
    		}*/

			_mm_store_si128((vuint8*) &Y[row][col], vSHIFT1);	
		}
	}

} 

/*------------------------------------------------------------------------------------------*/
void ui8matrix_dilation_dilation_SIMD_FO(vuint8** X, long nrl, long nrh, long v0, long v1, vuint8 **Y) {
/*------------------------------------------------------------------------------------------*/

	long row, col, x, y, i;
	long snrl = SE_NRL, snrh = SE_NRH, sncl = SE_NCL, snch = SE_NCH;

	vuint8 v_2X_1, v_2X, v_2X1;
	vuint8 v_1X_1, v_1X, v_1X1;
	vuint8   vX_1,   vX,   vX1;
	vuint8  v1X_1,  v1X,  v1X1;
	vuint8  v2X_1,  v2X,  v2X1;

	vuint8 vY_1, vY, vY1;
	vuint8 vSHIFT1, vSHIFT2;


	// Dilation
	for (row = nrl; row <= nrh; row++) {
		for (col = v0; col <= v1; col++) {

			v_2X_1  = _mm_load_si128((vuint8*) &X[row-2][col-1]);
			v_1X_1 	= _mm_load_si128((vuint8*) &X[row-1][col-1]);
			vX_1	= _mm_load_si128((vuint8*) &X[row  ][col-1]);
			v1X_1 	= _mm_load_si128((vuint8*) &X[row+1][col-1]);
			v_2X_1	= _mm_load_si128((vuint8*) &X[row+2][col-1]);
			v_2X    = _mm_load_si128((vuint8*) &X[row-2][col  ]);
			v_1X   	= _mm_load_si128((vuint8*) &X[row-1][col  ]);
			vX   	= _mm_load_si128((vuint8*) &X[row  ][col  ]);
			v1X   	= _mm_load_si128((vuint8*) &X[row+1][col  ]);
			v2X 	= _mm_load_si128((vuint8*) &X[row+2][col  ]);
			v_2X1   = _mm_load_si128((vuint8*) &X[row-2][col+1]);
			v_1X1  	= _mm_load_si128((vuint8*) &X[row-1][col+1]);
			vX1  	= _mm_load_si128((vuint8*) &X[row  ][col+1]);
			v1X1  	= _mm_load_si128((vuint8*) &X[row+1][col+1]);
			v2X1    = _mm_load_si128((vuint8*) &X[row+2][col+1]);


			//Column -1 "or" operator
			vY_1 = vector_or5(v_2X_1, v_1X_1, vX_1, v1X_1, v2X_1);
			//Column  0 "or" operator 
			vY   = vector_or5(v_2X,   v_1X,   vX,   v1X,   v2X);
			//Column  1 "or" operator
			vY1  = vector_or5(v_2X1,  v_1X1,  vX1,  v1X1,  v2X1);

			//Row operator
			vSHIFT1 = vector_or3_row1shift(vY_1, vY,   vY1);
			vSHIFT2 = vector_or3_row2shift(vY_1, vY,   vY1);
			vSHIFT1 = _mm_and_si128(vSHIFT1, vSHIFT2);

			_mm_store_si128((vuint8*) &Y[row][col], vSHIFT1);	
		}
	}

} 

/*---------------------------------------------------------------------------------------------------------*/
void ui8matrix_erosion_erosion_SIMD_FO_RR_row (vuint8** X, long nrl, long nrh, long v0, long v1, vuint8 **Y) {
/*---------------------------------------------------------------------------------------------------------*/
	
	long row, col, x, y, i;
	long snrl = SE_NRL, snrh = SE_NRH, sncl = SE_NCL, snch = SE_NCH;
	
	vuint8 v_2X_1, v_2X, v_2X1;
	vuint8 v_1X_1, v_1X, v_1X1;
	vuint8   vX_1,   vX,   vX1;
	vuint8  v1X_1,  v1X,  v1X1;
	vuint8  v2X_1,  v2X,  v2X1;

	vuint8 vY_1, vY, vY1;
	vuint8 vSHIFT1, vSHIFT2;

	// Erosion
	for (row = nrl; row <= nrh; row++) {

		col = 0;
		v_2X_1  = _mm_load_si128((vuint8*) &X[row-2][col-1]);
		v_1X_1 	= _mm_load_si128((vuint8*) &X[row-1][col-1]);
		vX_1	= _mm_load_si128((vuint8*) &X[row  ][col-1]);
		v1X_1 	= _mm_load_si128((vuint8*) &X[row+1][col-1]);
		v2X_1	= _mm_load_si128((vuint8*) &X[row+2][col-1]);
		v_2X    = _mm_load_si128((vuint8*) &X[row-2][col  ]);
		v_1X   	= _mm_load_si128((vuint8*) &X[row-1][col  ]);
		vX   	= _mm_load_si128((vuint8*) &X[row  ][col  ]);
		v1X   	= _mm_load_si128((vuint8*) &X[row+1][col  ]);
		v2X 	= _mm_load_si128((vuint8*) &X[row+2][col  ]);

		//Column -1 "and" operator
		vY_1 = vector_and5(v_2X_1, v_1X_1, vX_1, v1X_1, v2X_1);
		//Column  0 "and" operator 
		vY   = vector_and5(v_2X,   v_1X,   vX,   v1X,   v2X);	

		for (col = v0; col <= v1; col++) {	

			v_2X1   = _mm_load_si128((vuint8*) &X[row-2][col+1]);
			v_1X1  	= _mm_load_si128((vuint8*) &X[row-1][col+1]);
			vX1  	= _mm_load_si128((vuint8*) &X[row  ][col+1]);
			v1X1  	= _mm_load_si128((vuint8*) &X[row+1][col+1]);
			v2X1    = _mm_load_si128((vuint8*) &X[row+2][col+1]);

			//Column  1 "and" operator
			vY1  = vector_and5(v_2X1,  v_1X1,  vX1,  v1X1,  v2X1);

			//Row operator
			vSHIFT1 = vector_and3_row1shift(vY_1, vY,   vY1);
			/*if(row == 3 && col == 1) {
				printf("SIMD_FO_RR\n");
    			display_vuint8(vY_1, "%4d", NULL); printf("\n\n");
    			display_vuint8(vY, "%4d", NULL); printf("\n\n");
    			display_vuint8(v_2X1, "%4d", NULL); printf("\n");
    			display_vuint8(v_1X1, "%4d", NULL); printf("\n");
    			display_vuint8(vX1, "%4d", NULL); printf("\n");
    			display_vuint8(v1X1, "%4d", NULL); printf("\n");
    			display_vuint8(v2X1, "%4d", NULL); printf("\n\n");
    			display_vuint8(vY1, "%4d", NULL); printf("\n\n");

    			display_vuint8(vec_right1(vY_1, vY), "%4d", NULL); printf("\n");
    			display_vuint8(vY, "%4d", NULL); printf("\n");
    			display_vuint8(vec_left1(vY, vY1), "%4d", NULL); printf("\n");
    			display_vuint8(vSHIFT1, "%4d", NULL); printf("\n\n");
    		}*/
			vSHIFT2 = vector_and3_row2shift(vY_1, vY,   vY1);
			/*if(row == 3 && col == 1) {

	   			display_vuint8(vec_right2(vY_1, vY_1), "%4d", NULL); printf("\n");
    			display_vuint8(vY, "%4d", NULL); printf("\n");
    			display_vuint8(vec_left2(vY, vY1), "%4d", NULL); printf("\n");
    			display_vuint8(vSHIFT2, "%4d", NULL); printf("\n\n");
    		}*/
			vSHIFT1 = _mm_and_si128(vSHIFT1, vSHIFT2);
			/*if(row == 3 && col == 1) {
    			display_vuint8(vSHIFT1, "%4d", NULL); printf("\n\n");
    		}*/

			_mm_store_si128((vuint8*) &Y[row][col], vSHIFT1);

			vY_1 = vY;
			vY = vY1;

		}
	}

}

/*-------------------------------------------------------------------------------------------------*/
void ui8matrix_dilation_dilation_SIMD_FO_RR_row (vuint8** X, long nrl, long nrh, long v0, long v1, vuint8 **Y) {
/*-------------------------------------------------------------------------------------------------*/
	
	long row, col, x, y, i;
	long snrl = SE_NRL, snrh = SE_NRH, sncl = SE_NCL, snch = SE_NCH;
	
	vuint8 v_2X_1, v_2X, v_2X1;
	vuint8 v_1X_1, v_1X, v_1X1;
	vuint8   vX_1,   vX,   vX1;
	vuint8  v1X_1,  v1X,  v1X1;
	vuint8  v2X_1,  v2X,  v2X1;

	vuint8 vY_1, vY, vY1;
	vuint8 vSHIFT1, vSHIFT2;

	// Dilation
	for (row = nrl; row <= nrh; row++) {

		col = 0;
		v_2X_1  = _mm_load_si128((vuint8*) &X[row-2][col-1]);
		v_1X_1 	= _mm_load_si128((vuint8*) &X[row-1][col-1]);
		vX_1	= _mm_load_si128((vuint8*) &X[row  ][col-1]);
		v1X_1 	= _mm_load_si128((vuint8*) &X[row+1][col-1]);
		v2X_1	= _mm_load_si128((vuint8*) &X[row+2][col-1]);
		v_2X    = _mm_load_si128((vuint8*) &X[row-2][col  ]);
		v_1X   	= _mm_load_si128((vuint8*) &X[row-1][col  ]);
		vX   	= _mm_load_si128((vuint8*) &X[row  ][col  ]);
		v1X   	= _mm_load_si128((vuint8*) &X[row+1][col  ]);
		v2X 	= _mm_load_si128((vuint8*) &X[row+2][col  ]);

		//Column -1 "or" operator
		vY_1 = vector_or5(v_2X_1, v_1X_1, vX_1, v1X_1, v2X_1);
		//Column  0 "or" operator 
		vY   = vector_or5(v_2X,   v_1X,   vX,   v1X,   v2X);	

		for (col = v0; col <= v1; col++) {	

			v_2X1   = _mm_load_si128((vuint8*) &X[row-2][col+1]);
			v_1X1  	= _mm_load_si128((vuint8*) &X[row-1][col+1]);
			vX1  	= _mm_load_si128((vuint8*) &X[row  ][col+1]);
			v1X1  	= _mm_load_si128((vuint8*) &X[row+1][col+1]);
			v2X1    = _mm_load_si128((vuint8*) &X[row+2][col+1]);

			//Column  1 "or" operator
			vY1  = vector_or5(v_2X1,  v_1X1,  vX1,  v1X1,  v2X1);

			//Row operator
			vSHIFT1 = vector_or3_row1shift(vY_1, vY,   vY1);
			vSHIFT2 = vector_or3_row2shift(vY_1, vY,   vY1);
			vSHIFT1 = _mm_or_si128(vSHIFT1, vSHIFT2);

			_mm_store_si128((vuint8*) &Y[row][col], vSHIFT1);

			vY_1 = vY;
			vY = vY1;

		}
	}

}

void test_functions_morpho_SIMD() {

	p_image   t0 = create_image ("../car3/car_3000.pgm");
	p_vimage vt0 = create_vimage("../car3/car_3000.pgm");
	p_vimage vt1 = create_vimage("../car3/car_3000.pgm");
	p_vimage vt2 = create_vimage("../car3/car_3000.pgm");
	p_vimage vt3 = create_vimage("../car3/car_3000.pgm");

	uint8** tmp;

	ui8matrix_erosion_naive(t0->I, t0->nrl+BORD, t0->nrh-BORD, t0->ncl+BORD, t0->nch-BORD, tmp, t0->O);
	ui8matrix_erosion_naive(t0->O, t0->nrl+BORD, t0->nrh-BORD, t0->ncl+BORD, t0->nch-BORD, tmp, t0->E);
	ui8matrix_erosion_SIMD_naive(vt0->I, vt0->nrl+BORD, vt0->nrh-BORD, vt0->v0+vBORD, vt0->v1-vBORD, vt0->O);
	ui8matrix_erosion_SIMD_naive(vt0->O, vt0->nrl+BORD, vt0->nrh-BORD, vt0->v0+vBORD, vt0->v1-vBORD, vt0->E);
	ui8matrix_erosion_SIMD_RR_row(vt1->I, vt1->nrl+BORD, vt1->nrh-BORD, vt1->v0+vBORD, vt1->v1-vBORD, vt1->O);
	ui8matrix_erosion_erosion_SIMD_FO(vt2->I, vt2->nrl+BORD, vt2->nrh-BORD, vt2->v0+vBORD, vt2->v1-vBORD, vt2->O);
	ui8matrix_erosion_erosion_SIMD_FO_RR_row(vt3->I, vt3->nrl+BORD, vt3->nrh-BORD, vt3->v0+vBORD, vt3->v1-vBORD, vt3->O);


	display_ui8vector((uint8*) t0->E[1], 0, 31, "%4d", "NAIVE O[1][0-31]");
	display_ui8vector((uint8*) t0->E[2], 0, 31, "%4d", "NAIVE O[2][0-31]");
	display_ui8vector((uint8*) t0->E[3], 0, 31, "%4d", "NAIVE O[3][0-31]");
	puts(""); puts("");
	display_ui8vector((uint8*) t0->I[1], 0, 31, "%4d", "NAIVE I[1][0 - 31]");
	display_ui8vector((uint8*) t0->I[2], 0, 31, "%4d", "NAIVE I[2][0 - 31]");
	display_ui8vector((uint8*) t0->I[3], 0, 31, "%4d", "NAIVE I[3][0 - 31]");
	
	puts(""); puts("");
	display_vui8vector(vt0->E[1], 0, 1, "%4d", "SIMD O[1][0 - 1]");
	display_vui8vector(vt0->E[2], 0, 1, "%4d", "SIMD O[2][0 - 1]");
	display_vui8vector(vt0->E[3], 0, 1, "%4d", "SIMD O[3][0 - 1]");
	puts(""); puts("");

	puts(""); puts("");
	display_vui8vector(vt2->O[1], 0, 1, "%4d", "SIMD_FO O[1][0 - 1]");
	display_vui8vector(vt2->O[2], 0, 1, "%4d", "SIMD_FO O[2][0 - 1]");
	display_vui8vector(vt2->O[3], 0, 1, "%4d", "SIMD_FO O[3][0 - 1]");
	puts(""); puts("");

	puts(""); puts("");
	display_vui8vector(vt3->O[1], 0, 1, "%4d", "SIMD_FO_RR O[1][0 - 1]");
	display_vui8vector(vt3->O[2], 0, 1, "%4d", "SIMD_FO_RR O[2][0 - 1]");
	display_vui8vector(vt3->O[3], 0, 1, "%4d", "SIMD_FO_RR O[3][0 - 1]");
	puts(""); puts("");

	display_vui8vector(vt3->I[1], 0, 1, "%4d", "SIMD_FO_RR I[1][0 - 1]");
	display_vui8vector(vt3->I[2], 0, 1, "%4d", "SIMD_FO_RR I[2][0 - 1]");
	display_vui8vector(vt3->I[3], 0, 1, "%4d", "SIMD_FO_RR I[3][0 - 1]");
	puts(""); puts("");

}