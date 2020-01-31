/* ------------------------ */
/* ----- morpho_SIMD.c ---- */
/* ------------------------ */

#include <stdio.h>
#include <stdlib.h>
#include <emmintrin.h>
#include <tmmintrin.h>

#include "nrdef.h"
#include "nrutil.h"

#include "vnrdef.h"
#include "vnrutil.h"

#include "mutil.h"

#include "mymacro.h"
#include "simd_macro.h"

#include "util.h"
#include "img_SIMD.h"
#include "morpho_SIMD.h"

#define SE_NRL -1
#define SE_NRH  1
#define SE_NCL -1
#define SE_NCH  1

#define vector_and3(input, col)		(input[col - 1] & input[col + 0] & input[col + 1])
#define vector_or3(input, col)		(input[col - 1] | input[col + 0] | input[col + 1])

/*------------------------------------------------------------------------------------------*/
void ui8matrix_erosion_SIMD_naive(vuint8** X, long nrl, long nrh, long v0, long v1, vuint8 **Y) {
/*------------------------------------------------------------------------------------------*/

	long row, col, x, y, i;
	long snrl = SE_NRL, snrh = SE_NRH, sncl = SE_NCL, snch = SE_NCH;

	vuint8 v_1X_1, v_1X, v_1X1;
	vuint8 vX_1  , vX  , vX1  ;
	vuint8 v1X_1 , v1X , v1X1 ;

	vuint8 vY_1, vY, vY1;
	vuint8 vTMP;


	// Erode
	for (row = nrl+BORD; row <= nrh-BORD; row++) {
		for (col = v0+vBORD; col <= v1-vBORD; col++) {

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
			vY_1 = _mm_and_si128(v_1X_1, vX_1);
			vY_1 = _mm_and_si128(v1X_1, vY_1);
			//Column  0 "and" operator 
			vY   = _mm_and_si128(v_1X, vX);
			vY	 = _mm_and_si128(vY, v1X);
			//Column  1 "and" operator
			vY1  = _mm_and_si128(v_1X1, vX1);
			vY1  = _mm_and_si128(v1X1, vY1);

			//Row operator
			vTMP = _mm_and_si128(vec_left1(vY, vY1), vY);
			vTMP = _mm_and_si128(vTMP, vec_right1(vY_1, vY));

			_mm_store_si128((vuint8*) &Y[row][col], vTMP);	
		}
	}

} 

/*------------------------------------------------------------------------------------*/
void ui8matrix_dilation_SIMD_naive(vuint8** X, long nrl, long nrh, long v0, long v1, vuint8 **Y) {
/*------------------------------------------------------------------------------------*/

	long row, col, x, y, i;
	long snrl = SE_NRL, snrh = SE_NRH, sncl = SE_NCL, snch = SE_NCH;

	vuint8 v_1X_1, v_1X, v_1X1;
	vuint8 vX_1  , vX  , vX1  ;
	vuint8 v1X_1 , v1X , v1X1 ;

	vuint8 vY_1, vY, vY1;
	vuint8 vTMP;

	// Dilatation
	for (row = nrl+BORD; row <= nrh-BORD; row++) {
		for (col = v0+vBORD; col <= v1-vBORD; col++) {

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
			vY_1 = _mm_or_si128(v_1X_1, vX_1);
			vY_1 = _mm_or_si128(v1X_1, vY_1);
			//Column  0 "and" operator 
			vY   = _mm_or_si128(v_1X, vX);
			vY	 = _mm_or_si128(vY, v1X);
			//Column  1 "and" operator
			vY1  = _mm_or_si128(v_1X1, vX1);
			vY1  = _mm_or_si128(v1X1, vY1);

			//Row operator
			vTMP = _mm_or_si128(vec_left1(vY, vY1), vY);
			vTMP = _mm_or_si128(vTMP, vec_right1(vY_1, vY));

			_mm_store_si128((vuint8*) &Y[row][col], vTMP);		
		}
	}
}

/*---------------------------------------------------------------------------------------------------------*/
void ui8matrix_erosion_SIMD_RR_line (vuint8** X, long nrl, long nrh, long v0, long v1, vuint8 **Y) {
/*---------------------------------------------------------------------------------------------------------*/
	
	long row, col, x, y, i;
	long snrl = SE_NRL, snrh = SE_NRH, sncl = SE_NCL, snch = SE_NCH;
	
	vuint8 v_1X_1, v_1X, v_1X1;
	vuint8 vX_1  , vX  , vX1  ;
	vuint8 v1X_1 , v1X , v1X1 ;

	vuint8 vY_1, vY, vY1;
	vuint8 vTMP;

	// Dilatation
	for (row = nrl+BORD; row <= nrh-BORD; row++) {

		col = 0;
		v_1X_1 	= _mm_load_si128((vuint8*) &X[row-1][col-1]);
		vX_1	= _mm_load_si128((vuint8*) &X[row  ][col-1]);
		v1X_1 	= _mm_load_si128((vuint8*) &X[row+1][col-1]);
		v_1X   	= _mm_load_si128((vuint8*) &X[row-1][col  ]);
		vX   	= _mm_load_si128((vuint8*) &X[row  ][col  ]);
		v1X   	= _mm_load_si128((vuint8*) &X[row+1][col  ]);

		//Column -1 "and" operator
		vY_1 = _mm_and_si128(v_1X_1, vX_1);
		vY_1 = _mm_and_si128(v1X_1, vY_1);
		//Column  0 "and" operator 
		vY   = _mm_and_si128(v_1X, vX);
		vY	 = _mm_and_si128(vY, v1X);

		for (col = v0+vBORD; col <= v1-vBORD; col++) {	

			v_1X1  	= _mm_load_si128((vuint8*) &X[row-1][col+1]);
			vX1  	= _mm_load_si128((vuint8*) &X[row  ][col+1]);
			v1X1  	= _mm_load_si128((vuint8*) &X[row+1][col+1]);

			//Column  1 "and" operator
			vY1  = _mm_and_si128(v_1X1, vX1);
			vY1  = _mm_and_si128(v1X1, vY1);

			//Row operator
			vTMP = _mm_and_si128(vec_left1(vY, vY1), vY);
			vTMP = _mm_and_si128(vTMP, vec_right1(vY_1, vY));

			_mm_store_si128((vuint8*) &Y[row][col], vTMP);

			vY_1 = vY;
			vY = vY1;

		}
	}

}

/*---------------------------------------------------------------------------------------------------------*/
void ui8matrix_dilatation_SIMD_RR_line (vuint8** X, long nrl, long nrh, long v0, long v1, vuint8 **Y) {
/*---------------------------------------------------------------------------------------------------------*/
	
	long row, col, x, y, i;
	long snrl = SE_NRL, snrh = SE_NRH, sncl = SE_NCL, snch = SE_NCH;
	
	vuint8 v_1X_1, v_1X, v_1X1;
	vuint8 vX_1  , vX  , vX1  ;
	vuint8 v1X_1 , v1X , v1X1 ;

	vuint8 vY_1, vY, vY1;
	vuint8 vTMP;

	// Dilatation
	for (row = nrl+BORD; row <= nrh-BORD; row++) {

		col = 0;
		v_1X_1 	= _mm_load_si128((vuint8*) &X[row-1][col-1]);
		vX_1	= _mm_load_si128((vuint8*) &X[row  ][col-1]);
		v1X_1 	= _mm_load_si128((vuint8*) &X[row+1][col-1]);
		v_1X   	= _mm_load_si128((vuint8*) &X[row-1][col  ]);
		vX   	= _mm_load_si128((vuint8*) &X[row  ][col  ]);
		v1X   	= _mm_load_si128((vuint8*) &X[row+1][col  ]);

		//Column -1 "and" operator
		vY_1 = _mm_or_si128(v_1X_1, vX_1);
		vY_1 = _mm_or_si128(v1X_1, vY_1);
		//Column  0 "and" operator 
		vY   = _mm_or_si128(v_1X, vX);
		vY	 = _mm_or_si128(vY, v1X);

		for (col = v0+vBORD; col <= v1-vBORD; col++) {	

			v_1X1  	= _mm_load_si128((vuint8*) &X[row-1][col+1]);
			vX1  	= _mm_load_si128((vuint8*) &X[row  ][col+1]);
			v1X1  	= _mm_load_si128((vuint8*) &X[row+1][col+1]);

			//Column  1 "and" operator
			vY1  = _mm_or_si128(v_1X1, vX1);
			vY1  = _mm_or_si128(v1X1, vY1);

			//Row operator
			vTMP = _mm_or_si128(vec_left1(vY, vY1), vY);
			vTMP = _mm_or_si128(vTMP, vec_right1(vY_1, vY));

			_mm_store_si128((vuint8*) &Y[row][col], vTMP);

			vY_1 = vY;
			vY = vY1;

		}
	}

}


void test_functions_morpho_SIMD() {

	p_vimage t0 = create_vimage("../car3/car_3000.pgm");
	p_vimage t1	= create_vimage("../car3/car_3000.pgm");

	ui8matrix_erosion_SIMD_naive(t0->I, t0->nrl, t0->nrh, t0->v0, t0->v1, t0->O);
	ui8matrix_erosion_SIMD_RR_line(t1->I, t1->nrl, t1->nrh, t1->v0, t1->v1, t1->O);

	display_vui8vector(t0->O[0], t0->v0+vBORD, t0->v1+vBORD, "%4d", "SSE NAIVE O[0]");
	puts(""); puts("");
	display_vui8vector(t1->O[0], t1->v0+vBORD, t1->v1+vBORD, "%4d", "SSE NAIVE O[0]");
	puts(""); puts("");
}
