/* ------------------------ */
/* --- mouvement_SIMD.c --- */
/* ------------------------ */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
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
void SigmaDelta_step0_SIMD(vuint8** M, vuint8** I, vuint8** V, long nrl, long nrh, int v0, int v1 , uint8 n_coeff, uint8 v_min, uint8 v_max) {
/*---------------------------------------------------------------------------------------------*/
	//copy_ui8matrix_ui8matrix(t0->I, t0->nrl, t0->nrh, t0->ncl, t0->nch, t0->M);
	// copy_vui8matrix_vui8matrix(I, nrl, nrh, v0, v1, M);
	memcpy(&M[nrl][v0], &I[nrl][v0], (nrh - nrl + 1) * (v1 - v0 + 1) *sizeof(vuint8));
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
	vuint8 vM_1sub, vM_1add;

	for(long i = nrl; i <= nrh; i++) {
		for(long j = v0; j <= v1; j++) {
			vI   = _mm_load_si128((vuint8*) &I[i][j]);
			vM_1 = _mm_load_si128((vuint8*) &M_1[i][j]);
			vM_1sub = _mm_subs_epu8(vM_1, ONE);
			vM_1add = _mm_adds_epu8(vM_1, ONE);

			vec_cmplt(vM_1, vI, C1, CMP);
			// vM_1 = _mm_sub_epi8(vM_1, CMP);
			// vI = _mm_sub_epi8(vI, CMP);
			// C1 = _mm_cmplt_epi8(vM_1, vI);
			// vM_1 = _mm_add_epi8(vM_1, CMP);
			// vI = _mm_add_epi8(vI, CMP);
			C2  = _mm_cmpeq_epi8(vI, vM_1);
		    TMP = _mm_or_si128(_mm_and_si128(C1, vM_1add), _mm_andnot_si128(C1, vM_1sub));
		    vM = _mm_or_si128(_mm_and_si128(C2, vM_1), _mm_andnot_si128(C2, TMP));
		    _mm_store_si128((vuint8*) &M[i][j], vM);

		}
	}
	
}

/*---------------------------------------------------------------------------------------------*/
void SigmaDelta_step2_SIMD(vuint8** M, vuint8** I, vuint8** O, long nrl, long nrh, int v0, int v1 , uint8 n_coeff, uint8 v_min, uint8 v_max) {
/*---------------------------------------------------------------------------------------------*/

	vuint8 vI, vM, vO;
	vuint8 vMI_sub, vShift;

	for(long i = nrl; i <= nrh; i++) {
		for(long j = v0; j <= v1; j++) {
			vI = _mm_load_si128((vuint8*) &I[i][j]);
			vM = _mm_load_si128((vuint8*) &M[i][j]);
			
			vO = vec_subabs(vI, vM);
			
			_mm_store_si128((vuint8*) &O[i][j], vO);
		}
	}
	
}
/*---------------------------------------------------------------------------------------------*/
void SigmaDelta_step2_InLU_O3_SIMD(vuint8** M, vuint8** I, vuint8** O, long nrl, long nrh, int v0, int v1 , uint8 n_coeff, uint8 v_min, uint8 v_max) {
/*---------------------------------------------------------------------------------------------*/

	vuint8 vI, vM, vO;
	vuint8  *vI_row, *vM_row, *vO_row;
	vuint8 vMI_sub, vShift;
	long order = 3;
	long r = (v1 - v0 + 1) % order;
	long i, j;

	for(long i = nrl; i <= nrh; i++) {
		vI_row = I[i];
		vM_row = M[i];
		vO_row = O[i];
		for(long j = v0; j <= v1 - r; j+= order) {
			vI = _mm_load_si128((vuint8*) &vI_row[j + 0]);
			vM = _mm_load_si128((vuint8*) &vM_row[j + 0]);
			vO = vec_subabs(vI, vM);
			_mm_store_si128((vuint8*) &vO_row[j+0], vO);

			vI = _mm_load_si128((vuint8*) &vI_row[j + 1]);
			vM = _mm_load_si128((vuint8*) &vM_row[j + 1]);
			vO = vec_subabs(vI, vM);
			_mm_store_si128((vuint8*) &vO_row[j+1], vO);
			
			vI = _mm_load_si128((vuint8*) &vI_row[j + 2]);
			vM = _mm_load_si128((vuint8*) &vM_row[j + 2]);
			vO = vec_subabs(vI, vM);
			_mm_store_si128((vuint8*) &vO_row[j]+2, vO );
		}
	}
	switch(r) {
		case 2:

		for(long i = nrl; i <= nrh; i++) {
			
			vI = _mm_load_si128((vuint8*) &I[i][v1 - 1]);
			vM = _mm_load_si128((vuint8*) &M[i][v1 - 1]);
			vO = vec_subabs(vI, vM);
			_mm_store_si128((vuint8*) &O[i][v1 -1], vO);

			vI = _mm_load_si128((vuint8*) &I[i][v1 + 0]);
			vM = _mm_load_si128((vuint8*) &M[i][v1 + 0]);
			vO = vec_subabs(vI, vM);
			_mm_store_si128((vuint8*) &O[i][v1 - 0], vO);
		}
		break;
		case 1:

		for(long i = nrl; i <= nrh; i++) {
			vI = _mm_load_si128((vuint8*) &I[i][v1 + 0]);
			vM = _mm_load_si128((vuint8*) &M[i][v1 + 0]);
			vO = vec_subabs(vI, vM);
			_mm_store_si128((vuint8*) &O[i][v1 - 0], vO);
		}
		break;
		case 0:
		break;
	}
	
}


/*-----------------------------------------------------------------------------------------------*/
void SigmaDelta_step3_SIMD(vuint8** V_1, vuint8** O, vuint8** V, long nrl, long nrh, int v0, int v1 , uint8 n_coeff, uint8 v_min, uint8 v_max) {
/*-----------------------------------------------------------------------------------------------*/

	vuint8 vV_1, vO, NvO, vV, TMP;
	vuint8 C1, C2, C3;

	vuint8 CMP = init_vuint8(128);
	vuint8 ONE = init_vuint8(1);
	vuint8 vV_1add, vV_1sub;


	vuint8 vuint8_max = init_vuint8(255);
	vuint8 vVmax = init_vuint8(v_max);
	vuint8 vVmin = init_vuint8(v_min);
	// vuint8 vVmin = init_vuint8(v_min);
	// printf("vmin:%u\n",v_min);
	// printf("vmax:%u\n",v_max);

	for(long i = nrl; i <= nrh; i++) {
		for(long j = v0; j <= v1; j++) {
			vV_1 = _mm_load_si128((vuint8*) &V_1[i][j]); 
			vO 	 = _mm_load_si128((vuint8*) &O[i][j]);   

			vV_1sub = _mm_max_epu8(_mm_subs_epu8(vV_1, ONE), vVmin);
			vV_1add = _mm_min_epu8(_mm_adds_epu8(vV_1, ONE), vVmax);
			// display_vuint8(vV_1, " %03u", "vV_1\t\t"); 			printf("\n");
			// display_vuint8(vV_1sub, " %03u", "vV_1sub\t\t");	printf("\n");
			// display_vuint8(vV_1add, " %03u", "vV_1add\t\t");	printf("\n");
			NvO = vO;
			for (int k = 1; k < n_coeff; k++) NvO =_mm_adds_epu8(NvO, vO);
			// display_vuint8(NvO, " %03u", "NvO\t\t");	printf("\n");

			// V_t0 < n * O_t1
			vec_cmplt(vV_1, NvO, C1, CMP);

			// V_t0 == n * O_t1
			C2  = _mm_cmpeq_epi8(NvO, vV_1);

			// V_t0 > n * O_t1 <=> !(V_t0 == n * O_t1) && !(V_t0 < n * O_t1) && _mm_set_epi8(0xFF, ..., 0xFF)
			C3 = _mm_andnot_si128(C1, _mm_andnot_si128(C2, vuint8_max));
			// display_vuint8(C1, "%4u", "less than\t");	printf("\n");
			// display_vuint8(C2, "%4u", "equal to\t");	printf("\n");
			// display_vuint8(C3, "%4u", "greater than\t");printf("\n");
			// NvO = _mm_sub_epi8(NvO, CMP);
			// vM_1 = _mm_sub_epi8(vM_1, CMP);
			// C1  = _mm_cmplt_epi8(vM_1, NvO);
			// NvO = _mm_add_epi8(NvO, CMP);
			// vM_1 = _mm_add_epi8(vM_1, CMP);


			// If V_t0 < n * O_t1 then vV := vV_1add 
			// OR
			// If V_t0 > n * O_t1 then vV := vV_1sub
		    TMP = _mm_or_si128(_mm_and_si128(C1, vV_1add), _mm_and_si128(C3, vV_1sub));
			// display_vuint8(TMP, "%4u", "TMP\t\t");printf("\n");
		    TMP = _mm_or_si128(TMP, _mm_and_si128(C2, vV_1));
			// display_vuint8(TMP, "%4u", "TMP\t\t");printf("\n");
		    vV  = _mm_max_epu8(_mm_min_epu8(TMP, vVmax), vVmin);

		    _mm_store_si128((vuint8*) &V[i][j], vV);

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
	SigmaDelta_step2_InLU_O3_SIMD(t1->M, t1->I, t1->O, t1->nrl, t1->nrh, t1->v0, t1->v1, n_coeff, v_min, v_max);
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


