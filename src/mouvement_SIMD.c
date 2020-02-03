/* ------------------------ */
/* --- mouvement_SIMD.c --- */
/* ------------------------ */

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#include "nrdef.h"
#include "nrutil.h"

#include "vnrdef.h"
#include "vnrutil.h"

#include "mutil.h"

#include "mymacro.h"
#include "simd_macro.h"

#include "img_SIMD.h"
#include "mouvement_SIMD.h"
#include "util.h"


/*---------------------------------------------------------------------------------------------*/
void copy_vui8matrix_vui8matrix(vuint8** X, long nrl, long nrh, long vmin, long vmax, vuint8** Y) {
/*---------------------------------------------------------------------------------------------*/
	long i;
	int j;

	vuint8 vX;

	for (i = nrl; i <= nrh; i++)
		for (j = vmin; j <= vmax; j++ ) {
			vX = _mm_load_si128(&X[i][j]);
			_mm_store_si128((vuint8*)&Y[i][j], vX);
		}
}

/*---------------------------------------------------------------------------------------------*/
void SigmaDelta_step0_SIMD(vuint8** M, vuint8** I, vuint8** V, long nrl, long nrh, int v0, int v1 , uint8 n_coeff, uint8 v_min, uint8 v_max) {
/*---------------------------------------------------------------------------------------------*/
	//copy_ui8matrix_ui8matrix(t0->I, t0->nrl, t0->nrh, t0->ncl, t0->nch, t0->M);
	copy_vui8matrix_vui8matrix(I, nrl, nrh, v0, v1, M);
	long i, j;
	vuint8 vec1 = init_vuint8(1);
	for (i = nrl; i <= nrh; i++)
		for (j = v0; j <= v1; j++)
			_mm_store_si128((vuint8**) &V[i][j], vec1);
}

/*-----------------------------------------------------------------------------------------------*/
void SigmaDelta_step1_SIMD(vuint8** M_1, vuint8** I, vuint8** M, long nrl, long nrh, int v0, int v1 , uint8 n_coeff, uint8 v_min, uint8 v_max) {
/*-----------------------------------------------------------------------------------------------*/

	vuint8 vI, vM_1, vM, TMP, A, B;
	vuint8 C1, C2, C3;

	vuint8 CMP = init_vuint8(128);
	vuint8 ONE = init_vuint8(1);
	vuint8 ZERO = init_vuint8(0);
	vsint8 M_ONE = init_vsint8(-1);
	vsint8 UI8MAX = init_vsint8(255);

	for(long i = nrl; i <= nrh; i++) {
		for(long j = v0; j <= v1; j++) {
			vI   = _mm_load_si128((vuint8*) &I[i][j]);
			vM_1 = _mm_load_si128((vuint8*) &M_1[i][j]);

			vec_cmplt(vM_1, vI, C1, CMP);
			// vM_1 = _mm_sub_epi8(vM_1, CMP);
			// vI = _mm_sub_epi8(vI, CMP);
			// C1 = _mm_cmplt_epi8(vM_1, vI);
			// vM_1 = _mm_add_epi8(vM_1, CMP);
			// vI = _mm_add_epi8(vI, CMP);
			C2  = _mm_cmpeq_epi8(vI, vM_1);
		    TMP = _mm_or_si128(_mm_and_si128(C1, ONE), _mm_andnot_si128(C1, M_ONE));
		    TMP = _mm_or_si128(_mm_and_si128(C2, ZERO), _mm_andnot_si128(C2, TMP));
		    vM  = _mm_add_epi8(vM_1, TMP);	
		    _mm_store_si128((vuint8*) &M[i][j], vM);

		}
	}
	
}
	

/*---------------------------------------------------------------------------------------------*/
void SigmaDelta_step2_SIMD(vuint8** M, vuint8** I, vuint8** O, long nrl, long nrh, int v0, int v1 , uint8 n_coeff, uint8 v_min, uint8 v_max) {
/*---------------------------------------------------------------------------------------------*/

	vuint8 vI, vM, vO;


	for(long i = nrl; i <= nrh; i++) {
		for(long j = v0; j <= v1; j++) {
			vI = _mm_load_si128((vuint8*) &I[i][j]);
			vM = _mm_load_si128((vuint8*) &M[i][j]);
			
			vO = vec_subabs(vI, vM);
			
			_mm_store_si128((vuint8*) &O[i][j], vO);
		}
	}
	
}

/*-----------------------------------------------------------------------------------------------*/
void SigmaDelta_step3_SIMD(vuint8** V_1, vuint8** O, vuint8** V, long nrl, long nrh, int v0, int v1 , uint8 n_coeff, uint8 v_min, uint8 v_max) {
/*-----------------------------------------------------------------------------------------------*/

	
	vuint16 vV0, Vt_1_0, vO0, NvO0, TMP0;
	vuint16 vV1, Vt_1_1, vO1, NvO1, TMP1;
	vuint16 C1_0, C2_0;
	vuint16 C1_1, C2_1;
 	vuint8  vV;



	vuint16 CMP = init_vuint16(32768);
	vuint16 ONE = init_vuint16(1);
	vuint16 ZERO = init_vuint16(0);
	vsint16 M_ONE = init_vsint16(-1);

	vuint16 vVmax = init_vuint16(v_max);
	vuint16 vVmin = init_vuint16(v_min);
	uint8 *p;
	uint16 *q, *r;
	for(long i = nrl; i <= nrh; i++) {
		for(long j = v0; j <= v1; j++) {
			
			p = (uint8*)&V_1[i][j];
			
			Vt_1_0 = _mm_set_epi16(p[7],p[6],p[5],p[4],
								  p[3],p[2],p[1],p[0]);
								
			Vt_1_1 = _mm_set_epi16(p[15], p[14], p[13], p[12],
								  p[11], p[10], p[ 9], p[ 8]);
			
			p = (uint8*)&O[i][j];
			vO0 = _mm_set_epi16(p[7],p[6],p[5],p[4],
				     		    p[3],p[2],p[1],p[0]);
								
			vO1 = _mm_set_epi16(p[15], p[14], p[13], p[12],
			 				    p[11], p[10], p[ 9], p[ 8]);
			NvO0 = vO0;		
			NvO1 = vO1;

			for (int k = 0; k < n_coeff - 1; k++) {
				NvO0 = _mm_add_epi16(NvO0, vO0);
				NvO1 = _mm_add_epi16(NvO1, vO1);
			}

			vec16_cmplt(Vt_1_0, NvO0, C1_0, CMP);
			vec16_cmplt(Vt_1_1, NvO1, C1_1, CMP);

			// NvO = _mm_sub_epi8(NvO, CMP);
			// vM_1 = _mm_sub_epi8(vM_1, CMP);
			// C1  = _mm_cmplt_epi8(vM_1, NvO);
			// NvO = _mm_add_epi8(NvO, CMP);
			// vM_1 = _mm_add_epi8(vM_1, CMP);

			C2_0  = _mm_cmpeq_epi16(NvO0, Vt_1_0);
			C2_1  = _mm_cmpeq_epi16(NvO1, Vt_1_1);
				
			
		    TMP0 = _mm_or_si128(_mm_and_si128(C1_0, ONE) , _mm_andnot_si128(C1_0, M_ONE));
			TMP0 = _mm_or_si128(_mm_and_si128(C2_0, ZERO), _mm_andnot_si128(C2_0, TMP0));
		    
			TMP1 = _mm_or_si128(_mm_and_si128(C1_1, ONE) , _mm_andnot_si128(C1_1, M_ONE));
			TMP1 = _mm_or_si128(_mm_and_si128(C2_1, ZERO), _mm_andnot_si128(C2_1, TMP1));
			// display_vui16vector(&C1_1, 0, 0, "%4u ", "C1(HI)");
			// display_vui16vector(&C2_1, 0, 0, "%4u ", "C2(HI)");
			// display_vui16vector(&NvO0, 0, 0, "%4u ", "NvO0");
			// display_vui16vector(&NvO1, 0, 0, "%4u ", "NvO1");
			// display_vui16vector(&TMP0, 0, 0, "%4u ", "TMP0");
			// display_vui16vector(&TMP1, 0, 0, "%4u ", "TMP1");
			// display_vui16vector(&Vt_1_0, 0, 0, "%4u ", "V0(Max)");
			

		    vV0  = _mm_add_epi16(Vt_1_0, TMP0);
			vV0  = _mm_max_epi16(_mm_min_epi16(vV0, vVmax), vVmin);


			vV1  = _mm_add_epi16(Vt_1_1, TMP1);
		    vV1  = _mm_max_epi16(_mm_min_epi16(vV1, vVmax), vVmin);
			q = (int16 *) &vV1;
			r = (int16 *) &vV0;
			V[i][j] = _mm_set_epi8(q[7],q[6],q[5],q[4],q[3],q[2],q[1],q[0],
					               r[7],r[6],r[5],r[4],r[3],r[2],r[1],r[0]);
		    // _mm_store_si128((vuint8*) &, vV);

		}
	}
}

/*---------------------------------------------------------------------------------------------*/
void SigmaDelta_step4_SIMD(vuint8** O, vuint8** V, vuint8** E, long nrl, long nrh, int v0, int v1 , uint8 n_coeff, uint8 v_min, uint8 v_max) {
/*---------------------------------------------------------------------------------------------*/

	vuint8 vO, vV, tE;

	vuint8 ZERO = init_vuint8(0);
	vuint8 ONE = init_vuint8(1);

	vuint8 C, CMP = init_vuint8(128);

	for(long i = nrl; i <= nrh; i++) {
		for(long j = v0; j <= v1; j++) {
			vO = _mm_load_si128((vuint8*) &O[i][j]);
			vV = _mm_load_si128((vuint8*) &V[i][j]);
			
			vec_cmplt(vO, vV, C, CMP);
			// vO = _mm_sub_epi8(vO, CMP);
			// vV = _mm_sub_epi8(vV, CMP);
			// C  = _mm_cmplt_epi8(vO, vV);
			// vO = _mm_add_epi8(vO, CMP);
			// vV = _mm_add_epi8(vV, CMP);

    		tE = _mm_or_si128(_mm_and_si128(C, ZERO), _mm_andnot_si128(C, ONE));
            _mm_store_si128(&E[i][j], tE);
		}
	}
	
}
void SigmaDelta_SIMD(p_vimage t0, p_vimage t1, uint8 n_coeff, uint8 v_min, uint8 v_max)
{
	SigmaDelta_step1_SIMD(t0->M, t1->I, t1->M, t1->nrl, t1->nrh, t1->v0, t1->v1, n_coeff, v_min, v_max);
	SigmaDelta_step2_SIMD(t1->M, t1->I, t1->O, t1->nrl, t1->nrh, t1->v0, t1->v1, n_coeff, v_min, v_max);
	SigmaDelta_step3_SIMD(t0->V, t1->O, t1->V, t1->nrl, t1->nrh, t1->v0, t1->v1, n_coeff, v_min, v_max);
	SigmaDelta_step4_SIMD(t1->O, t1->V, t1->E, t1->nrl, t1->nrh, t1->v0, t1->v1, n_coeff, v_min, v_max);
}

void test_step0_SIMD() {

	// p_vimage vT_1 = create_vimage("../car3/car_3000.pgm");
	// p_vimage vT   = create_vimage("../car3/car_3001.pgm");

	// SigmaDelta_step0_SIMD(vT_1->I, vT_1->M, vT_1->V, vT_1->nrl, vT_1->nrh, vT_1->v0, vT_1->v1);
	// SigmaDelta_step1_SIMD(vT->I, vT_1->M, vT->M, vT->nrl, vT->nrh, vT->v0, vT->v1);
}


