/* ------------------------ */
/* ------ img_SIMD.c ------ */
/* ------------------------ */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
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
    printf("rl:%3ld rh:%3ld\n", rl, rh);
    printf("cl:%3ld ch:%3ld\n", cl, ch);
    printf("nrl:%3ld nrh:%3ld\n", nrl, nrh);
    printf("ncl:%3ld nch:%3ld\n", ncl, nch);
    printf("v0:%3d v1:%3d\n", v0, v1);
    printf("m0:%3d m1:%3d\n", m0, m1);
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
	uint8** test = ui8matrix(vt->nrl, vt->nrh, vt-> ncl, vt->nch);


    display_vui8vector(vt->I[0], vt->v0, vt->v1, "%4d", "vT"); puts("");
    vui8matrix2ui8matrix(vt->I, test, vt->nrl, vt->nrh, vt->v0, vt->v1);
    display_ui8vector((uint8*) test[2], vt->ncl, vt->nch, "%4d", "v"); puts("");
}
