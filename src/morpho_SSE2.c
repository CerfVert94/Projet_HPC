/* ------------------------ */
/* ----- morpho_SSE2.c ---- */
/* ------------------------ */

#include <stdio.h>
#include <stdlib.h>

#include "nrdef.h"
#include "nrutil.h"

#include "vnrdef.h"
#include "vnrutil.h"

#include "mutil.h"

#include "mymacro.h"
#include "simd_macro.h"

#include "img_SSE2.h"
#include "morpho_SSE2.h"
#include "util.h"

#define SE_NRL -1
#define SE_NRH  1
#define SE_NCL -1
#define SE_NCH  1


/*------------------------------------------------------------------------------------*/
void ui8matrix_erosion_SSE(vuint8** X, long nrl, long nrh, long v0, long v1, vuint8 **Y) {
/*------------------------------------------------------------------------------------*/

	long row, col, x, y, i;
	long snrl = SE_NRL, snrh = SE_NRH, sncl = SE_NCL, snch = SE_NCH;

	vuint8 vX, vX_1, vX1, vY;

	// Erode
	for (row = nrl; row <= nrh; row++) {
		for (col = v0; col <= v1; col++) {

			// CURENT ROW
			vX   = _mm_load_si128((vuint8*) &vX[row+0][col+0]);
			vX_1 = _mm_load_si128((vuint8*) &vX[row+0][col-1]);
			vX1  = _mm_load_si128((vuint8*) &vX[row+0][col+1]);
			//LEFT SHIFT AND CENTER
			vY = _mm_and_si128(vec_left1( vX,  vX1), vX);
			//RIGHT SHIFT
			vY = _mm_and_si128(vec_right1(vX, VX_1), vY);

			// PREVIOUS ROW
			vX   = _mm_load_si128((vuint8*) &vX[row-1][col+0]);
			vX_1 = _mm_load_si128((vuint8*) &vX[row-1][col-1]);
			vX1  = _mm_load_si128((vuint8*) &vX[row-1][col+1]);
			//LEFT SHIFT
			vY = _mm_and_si128(vec_left1( vX,  vX1), vY);
			//CENTER
			vY = _mm_and_si128(vX, vY);
			//RIGHT SHIFT
			vY = _mm_and_si128(vec_right1(vX, VX_1), vY);

			// NEXT ROW
			vX   = _mm_load_si128((vuint8*) &vX[row+1][col+0]);
			vX_1 = _mm_load_si128((vuint8*) &vX[row+1][col-1]);
			vX1  = _mm_load_si128((vuint8*) &vX[row+1][col+1]);
			//LEFT SHIFT
			vY = _mm_and_si128(vec_left1( vX,  vX1), vY);
			//CENTER
			vY = _mm_and_si128(vX, vY);
			//RIGHT SHIFT
			vY = _mm_and_si128(vec_right1(vX, VX_1), vY);

			_mm_store_si128((vuint8*) &Y[row][col], vY);
		}
	}

} 

/*--------------------------------------------------------------------------------------*/
void ui8matrix_dilation_naive(uint8** X, long nrl, long nrh, long v0, long v1,  uint8 **Y) {
/*--------------------------------------------------------------------------------------*/

	long row, col, x, y, i;
	long snrl = SE_NRL, snrh = SE_NRH, sncl = SE_NCL, snch = SE_NCH;

	vuint8 vX, vX_1, vX1, vY;

	// Dilatation
	for (row = nrl; row <= nrh; row++) {
		for (col = v0; col <= v1; col++) {

			// CURRENT ROW
			vX   = _mm_load_si128((vuint8*) &vX[row+0][col+0]);
			vX_1 = _mm_load_si128((vuint8*) &vX[row+0][col-1]);
			vX1  = _mm_load_si128((vuint8*) &vX[row+0][col+1]);
			//LEFT SHIFT AND CENTER
			vY = _mm_or_si128(vec_left1( vX,  vX1), vX);
			//RIGHT SHIFT
			vY = _mm_or_si128(vec_right1(vX, VX_1), vY);

			// PREVIOUS ROW
			vX   = _mm_load_si128((vuint8*) &vX[row-1][col+0]);
			vX_1 = _mm_load_si128((vuint8*) &vX[row-1][col-1]);
			vX1  = _mm_load_si128((vuint8*) &vX[row-1][col+1]);
			//LEFT SHIFT
			vY = _mm_or_si128(vec_left1( vX,  vX1), vY);
			//CENTER
			vY = _mm_or_si128(vX, vY);
			//RIGHT SHIFT
			vY = _mm_or_si128(vec_right1(vX, VX_1), vY);

			// NEXT ROW
			vX   = _mm_load_si128((vuint8*) &vX[row+1][col+0]);
			vX_1 = _mm_load_si128((vuint8*) &vX[row+1][col-1]);
			vX1  = _mm_load_si128((vuint8*) &vX[row+1][col+1]);
			//LEFT SHIFT
			vY = _mm_or_si128(vec_left1( vX,  vX1), vY);
			//CENTER
			vY = _mm_or_si128(vX, vY);
			//RIGHT SHIFT
			vY = _mm_or_si128(vec_right1(vX, VX_1), vY);

			_mm_store_si128((vuint8*) &Y[row][col], vY);		
		}
	}
}
