/* ------------------------ */
/* ------ img_SIMD.c ------ */
/* ------------------------ */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <assert.h>
#include <stdint.h>
#include <emmintrin.h>
#include <tmmintrin.h>

#include "nrdef.h"
#include "nrutil.h"
#include "util.h"

#include "vnrdef.h"
#include "vnrutil.h"
#include "myvnrutil.h"

#include "img.h"
#include "img_SIMD.h"

/*----------------------------------*/
p_vimage create_vimage(char* filename) {
/*----------------------------------*/

	int i, j, k, l, z, r;
	long rl, rh, cl, ch;
	long nrl, nrh, ncl, nch;

	uint8** img;
	uint8** imgb;

	int card;
	int v0, v1;
	int m0, m1;
	p_vimage tmp;
	vuint8 x;
	uint8 *p = (uint8*) &x;

	/*alloc in memory for image*/
	tmp = (p_vimage)malloc(sizeof(vimage));
	if (!tmp) {error("Malloc error of vimage in ");}

	img = LoadPGM_ui8matrix(filename, &rl, &rh, &cl, &ch);
	nrl = rl - BORD;
	nrh = rh + BORD;
	ncl = cl - BORD;
	nch = ch + BORD;

	tmp->nrl = nrl;
	tmp->nrh = nrh;
	tmp->ncl = ncl;
	tmp->nch = nch;

	card = card_vuint8();
    s2v1D(ncl, nch, card, &v0, &v1);
    v2m1D(v0, v1, card, &m0, &m1);
    // printf("rl:%3ld rh:%3ld\n", rl, rh);
    // printf("cl:%3ld ch:%3ld\n", cl, ch);
    // printf("nrl:%3ld nrh:%3ld\n", nrl, nrh);
    // printf("ncl:%3ld nch:%3ld\n", ncl, nch);
    // printf("v0:%3d v1:%3d\n", v0, v1);
    // printf("m0:%3d m1:%3d\n", m0, m1);
    tmp->v0 = v0;
    tmp->v1 = v1;
    tmp->m0 = m0;
    tmp->m1 = m1;


    z = (v1-vBORD)*card;
    r = ch-z;

    tmp->I = vui8matrix(nrl, nrh, v0, v1);
    for (i = rl; i <= rh; i++) {
        for (j = v0+vBORD; j < v1-vBORD; j++) {
            l = j*card;
            x = init_vuint8_all(img[i][l], img[i][l+1], img[i][l+2], img[i][l+3], img[i][l+4], img[i][l+5], img[i][l+6], img[i][l+7], img[i][l+8], img[i][l+9], img[i][l+10], img[i][l+11], img[i][l+12], img[i][l+13], img[i][l+14], img[i][l+15]);
            _mm_store_si128(&tmp->I[i][j], x);
        }
        x = init_vuint8(0);
        for (k = 0; k <= r; k++)
        	p[k] = img[i][z+k];
        _mm_store_si128(&tmp->I[i][v1-vBORD], x);
    }

	tmp->M = vui8matrix(nrl, nrh, v0, v1); 
    tmp->O = vui8matrix(nrl, nrh, v0, v1);
    tmp->V = vui8matrix(nrl, nrh, v0, v1);
    tmp->E = vui8matrix(nrl, nrh, v0, v1);
	tmp->Omega = vui8matrix(nrl, nrh, v0, v1);
	free_ui8matrix(img,rl,rh,cl,ch);
	return tmp;
}
p_vimage create_vimage_from_ui8matrix(uint8 **X, long nrl, long nrh, long ncl, long nch) {
	p_vimage tmp;
	int card, v0, v1, m0, m1;

	/*alloc in memory for image*/
	tmp = (p_vimage)malloc(sizeof(vimage));
	if (!tmp) {error("Malloc error of image in");}
	nrl -= BORD;
	nrh += BORD;
	ncl -= BORD;
	nch += BORD;

	
	tmp->nrl = nrl;
	tmp->nrh = nrh;
	tmp->ncl = ncl;
	tmp->nch = nch;
	card = card_vuint8();
    s2v1D(ncl, nch, card, &v0, &v1);
    v2m1D(v0, v1, card, &m0, &m1);
    tmp->v0 = v0;
    tmp->v1 = v1;
    tmp->m0 = m0;
    tmp->m1 = m1;

	// ui8matrix(nrl, nrh, ncl, nch);
	// memset_ui8matrix(tmp->I, 0, nrl, nrh, ncl, nch);
	// tmp->I = vui8matrix(nrl - BORD, nrh + BORD, v0, v1);

	// printf("%d %d %d %d", nrl, nrh, v0, v1);getchar();
    int z = (v1-vBORD)*card;
    int r = (nch - BORD) -z;
	card = card_vuint8();
	uint8 *p;
	vuint8 x;
	int l;
	tmp->I = vui8matrix(nrl, nrh, v0, v1); 

	for (int i = nrl; i < nrh + 1; i++) {
		for (int j = v0; j < v1 + 1; j++){
			tmp->I[i][j] = _mm_setzero_si128();
		}
	}
	for (int i = nrl + BORD; i < nrh - BORD + 1; i++) {
		int l;
		p = (uint8*)&x;
		for (int j = v0 + vBORD; j < v1; j++){
			// printf("%d\n", j);
			l = j * card;
			// printf("%d ~ %d\n", l, l + 15);
			for (int k = 0; k < card; k++)
				p[k] = X[i][l + k];
			_mm_store_si128((vuint8*)&tmp->I[i][j], x);
		}
		l = v1 * card;
		for (int j = 0; j < card && l < nch - BORD + 1; j++, l++) {
			p[j] = X[i][l + j];
		}
		_mm_store_si128((vuint8*)&tmp->I[i][v1 + 1], x);
	}
	tmp->M = vui8matrix(nrl, nrh, v0, v1); 
    tmp->O = vui8matrix(nrl, nrh, v0, v1);
    tmp->V = vui8matrix(nrl, nrh, v0, v1);
    tmp->E = vui8matrix(nrl, nrh, v0, v1);
	tmp->Omega = vui8matrix(nrl, nrh, v0, v1);
	return tmp;
}
/*-----------------------------*/
void free_vimage(p_vimage vimage) {
/*-----------------------------*/

	free_vui8matrix(vimage->I, vimage->nrl, vimage->nrh, vimage->v0, vimage->v1);
	free_vui8matrix(vimage->M, vimage->nrl, vimage->nrh, vimage->v0, vimage->v1);
	free_vui8matrix(vimage->O, vimage->nrl, vimage->nrh, vimage->v0, vimage->v1);
	free_vui8matrix(vimage->V, vimage->nrl, vimage->nrh, vimage->v0, vimage->v1);
	free_vui8matrix(vimage->E, vimage->nrl, vimage->nrh, vimage->v0, vimage->v1);
	free_vui8matrix(vimage->Omega, vimage->nrl, vimage->nrh, vimage->v0, vimage->v1);
	free(vimage);

}

/*---------------------*/
void test_SIMD_img() {
/*---------------------*/
	
	p_image t = create_image("../car3/car_3000.pgm");
	p_vimage vt = create_vimage("../car3/car_3000.pgm");

    display_ui8vector((uint8*) t->I[2], t->ncl, t->nch, "%4d", "v"); puts("");
    display_vui8vector(vt->I[2], vt->v0, vt->v1, "%4d", "vT"); puts("");


}
/*---------------------*/
void test2_SIMD_img() {
/*---------------------*/
	
	//p_image t = create_image("../car3/car_3000.pgm");
	p_vimage vt = create_vimage("../car3/car_3000.pgm");
	uint8** test; //ui8matrix(vt->nrl, vt->nrh, vt-> ncl, vt->nch);
	long nrl, nrh, ncl, nch;
	uint8** u8img = LoadPGM_ui8matrix("../car3/car_3000.pgm",&nrl, &nrh, &ncl, &nrh);
	long nrl_prime= 0, nrh_prime = 0, ncl_prime = 0, nch_prime = 0;


	printf("v0 ~ v1 : %d %d\n",  vt->v0, vt->v1);
	display_ui8matrix(u8img, 0, 5, 0, 31, "%4u", "image"); puts("");
	display_vui8matrix(vt->I, 0, 1, 0, 1, "%4u", "vimage");
    // display_vui8vector(vt->I[0], vt->v0, vt->v1, "%4u", "vT"); puts("");
    test = vui8matrix_to_ui8matrix(vt->I, 0, 5, 0, 1, &nrl_prime, &nrh_prime, &ncl_prime, &nch_prime);
	
    printf("%ld %ld %ld %ld\n", nrl_prime, nrh_prime, ncl_prime, nch_prime);
	display_ui8matrix(test, nrl_prime, nrh_prime, ncl_prime, nch_prime, "%4u", "image"); puts("");
	
	free_ui8matrix(u8img, nrl, nrh, ncl, nch);
	// display_vui8vector(vt->I[0], vt->v0, vt->v1, "%4u", "vT"); puts("");
    // display_ui8vector((uint8*) test[2], vt->ncl, vt->nch, "%4d", "v"); puts("");
}
