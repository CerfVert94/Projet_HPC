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

	uint8** var = ui8matrix(rl, rh, cl, ch);
    tmp->M = var;
    var = ui8matrix(rl, rh, cl, ch);
    tmp->O = var;
    var = ui8matrix(rl, rh, cl, ch);
    tmp->V = var;
    var = ui8matrix(rl, rh, cl, ch);
    tmp->E = var;

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
