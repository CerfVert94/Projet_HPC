/* ------------ */
/* ---img.c --- */
/* ------------ */

#include <stdbool.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <stdint.h>

#include "nrdef.h"
#include "nrutil.h"
#include "mynrutil.h"
#include "util.h"

#include "img.h"

p_image create_image(char* filename) {

	long rl, rh, cl, ch;
	long nrl, nrh, ncl, nch;
	p_image tmp;
	uint8** img;

	/*alloc in memory for image*/
	tmp = (p_image)malloc(sizeof(image));
	if (!tmp) {error("Malloc error of image in");}
	img = LoadPGM_ui8matrix(filename, &rl, &rh, &cl, &ch);
	nrl = rl - BORD;
	nrh = rh + BORD;
	ncl = cl - BORD;
	nch = ch + BORD;

	tmp->nrl = nrl;
	tmp->nrh = nrh;
	tmp->ncl = ncl;
	tmp->nch = nch;

	tmp->I = ui8matrix(nrl, nrh, ncl, nch);
	// printf("%ld %ld\n", nrh, nrl);
	memset_ui8matrix(tmp->I, 0, nrl, nrh, ncl, nch);
	
	copy_ui8matrix_ui8matrix(img, rl, rh, cl, ch, tmp->I);

	tmp->M = ui8matrix(nrl, nrh, ncl, nch); 
    tmp->O = ui8matrix(nrl, nrh, ncl, nch);
    tmp->V = ui8matrix(nrl, nrh, ncl, nch);
    tmp->E = ui8matrix(nrl, nrh, ncl, nch);
	tmp->Omega = ui8matrix(nrl, nrh, ncl, nch);

	free_ui8matrix(img, rl, rh, cl, ch);

	return tmp;
}


p_image create_image_from_ui8matrix(uint8 **X, long nrl, long nrh, long ncl, long nch) {
	p_image tmp;
	uint8** img;

	/*alloc in memory for image*/
	tmp = (p_image)malloc(sizeof(image));
	if (!tmp) {error("Malloc error of image in");}
	nrl -= BORD;
	nrh += BORD;
	ncl -= BORD;
	nch += BORD;

	tmp->nrl = nrl;
	tmp->nrh = nrh;
	tmp->ncl = ncl;
	tmp->nch = nch;

	tmp->I = ui8matrix(nrl, nrh, ncl, nch);
	memset_ui8matrix(tmp->I, 0, nrl, nrh, ncl, nch);
	copy_ui8matrix_ui8matrix(X, nrl + BORD, nrh - BORD, ncl + BORD, nch - BORD, tmp->I);

	tmp->M = ui8matrix(nrl, nrh, ncl, nch); 
    tmp->O = ui8matrix(nrl, nrh, ncl, nch);
    tmp->V = ui8matrix(nrl, nrh, ncl, nch);
    tmp->E = ui8matrix(nrl, nrh, ncl, nch);
	tmp->Omega = ui8matrix(nrl, nrh, ncl, nch);
	return tmp;
}

void free_image(p_image image) {

	free_ui8matrix(image->I, image->nrl, image->nrh, image->ncl, image->nch);
	free_ui8matrix(image->M, image->nrl, image->nrh, image->ncl, image->nch);
	free_ui8matrix(image->O, image->nrl, image->nrh, image->ncl, image->nch);
	free_ui8matrix(image->V, image->nrl, image->nrh, image->ncl, image->nch);
	free_ui8matrix(image->E, image->nrl, image->nrh, image->ncl, image->nch);
	free_ui8matrix(image->Omega, image->nrl, image->nrh, image->ncl, image->nch);	
	free(image);

}
