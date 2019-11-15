#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <stdint.h>
#include "nrdef.h"
#include "nrutil.h"
#include "util.h"

#include "img.h"

p_vimage create_vimage(char* filename) {

	int i, j, l;
	long rl, rh, cl, ch;
	uint8** img;

	int card;
	int v0, v1;
	int m0, v1;
	p_vimage tmp;

	/*alloc in memory for image*/
	tmp = (p_vimage)malloc(sizeof(vimage));
	if (!tmp) {error("Malloc error of vimage in ");}

	img = LoadPGM_ui8matrix(filename, &rl, &rh, &cl, &ch);

	tmp->nrl = rl;
	tmp->nrh = rh;
	tmp->ncl = cl;
	tmp->nch = ch;

	card = card_vuint8();
    s2v1D(rl, rh, card, v0, v1);
    v2M1D(v0, v1, card, m0, m1);
    tmp->v0 = v0;
    tmp->v1 = v1;
    tmp->m0 = m0;
    tmp->m1 = m1;

    for (i = rl; i < rh; i++) {
        for (j = v1; j <= v2; j++) {
            //x = init_vuint8(11);
            l = j*16;
            x = init_vuint8_all(img[i][j], img[i][l+1], img[i][l+2], img[i][l+3], img[i][l+4], img[i][l+5], img[i][l+6], img[i][l+7], img[i][l+8], img[i][l+9], img[i][l+10], img[i][l+11], img[i][l+12], img[i][l+13], img[i][l+14], img[i][l+15]);
            _mm_store_si128(&tmp->I[i][j], x);
        }
    }

	tmp->M = vui8matrix(rl, rh, v0, v1); 
    tmp->O = vui8matrix(rl, rh, v0, v1);
    tmp->V = vui8matrix(rl, rh, v0, v1);
    tmp->E = vui8matrix(rl, rh, v0, v1);
	tmp->Omega = vui8matrix(rl, rh, v0, v1);

	return tmp;
}

void free_vimage(p_image image) {

	free_vui8matrix(image->I, image->nrl, image->nrh, image->v0, image->v1);
	free_vui8matrix(image->M, image->nrl, image->nrh, image->v0, image->v1);
	free_vui8matrix(image->O, image->nrl, image->nrh, image->v0, image->v1);
	free_vui8matrix(image->V, image->nrl, image->nrh, image->v0, image->v1);
	free_vui8matrix(image->E, image->nrl, image->nrh, image->v0, image->v1);
	free_ui8matrix(image->Omega, image->nrl, image->nrh, image->v0, image->v1);
	free(image);

}
