/* ------------------------ */
/* ----- morpho_SIMD.c ---- */
/* ------------------------ */

#include <stdio.h>
#include <stdlib.h>
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
			vTMP = vector_and3(vec_right1(vY_1, vY), vY, vec_left1(vY, vY1));

			if(row == 1 && col == 1) {
    			display_vuint8(vY_1, "%4d", NULL); printf(" ");
    			display_vuint8(vY, "%4d", NULL); printf("\n");
    			display_vuint8(vY1, "%4d", NULL); printf("\n\n");
    			display_vuint8(vec_right1(vY_1, vY), "%4d", NULL); printf("\n");
    			display_vuint8(vY, "%4d", NULL); printf("\n");
    			display_vuint8(vec_left1(vY, vY1), "%4d", NULL); printf("\n\n");
    		}
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
			vTMP = vector_or3(vec_right1(vY_1, vY), vY, vec_left1(vY, vY1));

			_mm_store_si128((vuint8*) &Y[row][col], vTMP);		
		}
	}
}

/*---------------------------------------------------------------------------------------------------------*/
void ui8matrix_erosion_SIMD_RR_row (vuint8** X, long nrl, long nrh, long v0, long v1, vuint8 **Y) {
/*---------------------------------------------------------------------------------------------------------*/
	
	long row, col, x, y, i;
	long snrl = SE_NRL, snrh = SE_NRH, sncl = SE_NCL, snch = SE_NCH;
	
	vuint8 v_1X_1, v_1X, v_1X1;
	vuint8 vX_1  , vX  , vX1  ;
	vuint8 v1X_1 , v1X , v1X1 ;

	vuint8 vY_1, vY, vY1;
	vuint8 vTMP;

	// Dilatation
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
			vTMP = vector_and3(vec_right1(vY_1, vY), vY, vec_left1(vY, vY1));

			_mm_store_si128((vuint8*) &Y[row][col], vTMP);

			vY_1 = vY;
			vY = vY1;

		}
	}

}

/*-------------------------------------------------------------------------------------------------*/
void ui8matrix_dilation_SIMD_RR_row (vuint8** X, long nrl, long nrh, long v0, long v1, vuint8 **Y) {
/*-------------------------------------------------------------------------------------------------*/
	
	long row, col, x, y, i;
	long snrl = SE_NRL, snrh = SE_NRH, sncl = SE_NCL, snch = SE_NCH;
	
	vuint8 v_1X_1, v_1X, v_1X1;
	vuint8 vX_1  , vX  , vX1  ;
	vuint8 v1X_1 , v1X , v1X1 ;

	vuint8 vY_1, vY, vY1;
	vuint8 vTMP;

	// Dilatation
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
			vTMP = vector_or3(vec_right1(vY_1, vY), vY, vec_left1(vY, vY1));

			_mm_store_si128((vuint8*) &Y[row][col], vTMP);

			vY_1 = vY;
			vY = vY1;

		}
	}

}


void test_functions_morpho_SIMD() {

	p_vimage vt0 = create_vimage("../car3/car_3000.pgm");
	p_vimage vt1 = create_vimage("../car3/car_3000.pgm");
	p_image   t0 = create_image ("../car3/car_3000.pgm");

	uint8** tmp;

	ui8matrix_erosion_SIMD_naive(vt0->I, vt0->nrl+BORD, vt0->nrh-BORD, vt0->v0+vBORD, vt0->v1-vBORD, vt0->O);
	ui8matrix_erosion_SIMD_RR_row(vt1->I, vt1->nrl+BORD, vt1->nrh-BORD, vt1->v0+vBORD, vt1->v1-vBORD, vt1->O);
	//printf("1;%ld %ld %ld %ld\n", t0->nrl, t0->nrh, t0->ncl, t0->nch);
	//printf("2;%ld %ld %ld %ld\n", t0->nrl+BORD, t0->nrh-BORD, t0->ncl+BORD, t0->nch-BORD);
	ui8matrix_erosion_naive(t0->I, t0->nrl+BORD, t0->nrh-BORD, t0->ncl+BORD, t0->nch-BORD, tmp, t0->O);



	display_ui8vector((uint8*) t0->O[0], 0, 31, "%4d", "NAIVE O[0][0-31]");
	display_ui8vector((uint8*) t0->O[1], 0, 31, "%4d", "NAIVE O[1][0-31]");
	display_ui8vector((uint8*) t0->O[2], 0, 31, "%4d", "NAIVE O[2][0-31]");
	puts(""); puts("");
	display_ui8vector((uint8*) t0->I[0], 0, 31, "%4d", "NAIVE I[0][0 - 31]");
	display_ui8vector((uint8*) t0->I[1], 0, 31, "%4d", "NAIVE I[1][0 - 31]");
	display_ui8vector((uint8*) t0->I[2], 0, 31, "%4d", "NAIVE I[2][0 - 31]");
	
	puts(""); puts("");
	display_vui8vector(vt0->O[0], 0, 1, "%4d", "SSE RR O[0][0 - 1]");
	display_vui8vector(vt0->O[1], 0, 1, "%4d", "SSE RR O[1][0 - 1]");
	display_vui8vector(vt0->O[2], 0, 1, "%4d", "SSE RR O[2][0 - 1]");
	puts(""); puts("");

	display_vui8vector(vt0->I[0], 0, 1, "%4d", "SSE RR I[0][0 - 1]");
	display_vui8vector(vt0->I[1], 0, 1, "%4d", "SSE RR I[1][0 - 1]");
	display_vui8vector(vt0->I[2], 0, 1, "%4d", "SSE RR I[2][0 - 1]");
	puts(""); puts("");

}