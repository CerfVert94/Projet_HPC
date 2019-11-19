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

p_image create_image(char* filename) {

	long rl, rh, cl, ch;
	long nrl, nrh, ncl, nch;
	p_image tmp;

	/*alloc in memory for image*/
	tmp = (p_image)malloc(sizeof(image));
	if (!tmp) {error("Malloc error of image in");}
	tmp->I = LoadPGM_ui8matrix(filename, &rl, &rh, &cl, &ch);
	nrl = rl - 2;
	nrh = rh + 2;
	ncl = rl - 2;
	nch = rh + 2;
	
	tmp->nrl = rl;
	tmp->nrh = rh;
	tmp->ncl = cl;
	tmp->nch = ch;

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
