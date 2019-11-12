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

	int i, j;
	long rl, rh, cl, ch;
	p_image tmp;

	/*alloc in memory for image*/
	tmp = (p_image)malloc(sizeof(image));
	if (!tmp)
		printf("Malloc error of image in %s\n", __func__);
	tmp->I = LoadPGM_ui8matrix(filename, &rl, &rh, &cl, &ch);
	tmp->nrl = rl;
	tmp->nrh = rh;
	tmp->ncl = cl;
	tmp->nch = ch;

	tmp->M = ui8matrix(rl, rh, cl, ch); 
    tmp->O = ui8matrix(rl, rh, cl, ch);
    tmp->V = ui8matrix(rl, rh, cl, ch);
    tmp->E = ui8matrix(rl, rh, cl, ch);
	tmp->Omega = ui8matrix(rl, rh, cl, ch);

	return tmp;
}

void free_image(p_image image) {

	free_ui8matrix(image->I, image->nrl, image->nrh, image->ncl, image->nch);
	free_ui8matrix(image->M, image->nrl, image->nrh, image->ncl, image->nch);
	free_ui8matrix(image->O, image->nrl, image->nrh, image->ncl, image->nch);
	free_ui8matrix(image->V, image->nrl, image->nrh, image->ncl, image->nch);
	free_ui8matrix(image->E, image->nrl, image->nrh, image->ncl, image->nch);
	free(image);

}
