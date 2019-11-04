#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <stdint.h>
#include "nrdef.h"
#include "nrutil.h"

#include "img.h"

p_image create_image(char* filename) {

	int i, j;
	long rl, rh, cl, ch;
	p_image tmp;

	/*alloc in memory for image*/
	tmp = (p_image)malloc(sizeof(image));
	if (!tmp)
		printf("Malloc error of image in create_image\n");
	tmp->img = LoadPGM_ui8matrix(filename, &rl, &rh, &cl, &ch);
	tmp->nrl = rl;
	tmp->nrh = rh;
	tmp->ncl = cl;
	tmp->nch = ch;

	uint8** var = ui8matrix(rl, rh, cl, ch);
	tmp->source = var;
	var = ui8matrix(rl, rh, cl, ch);
    tmp->fond = var;
    var = ui8matrix(rl, rh, cl, ch);
    tmp->diff = var;
    var = ui8matrix(rl, rh, cl, ch);
    tmp->var = var;
    var = ui8matrix(rl, rh, cl, ch);
    tmp->move = var;

	return tmp;
}

void free_image(p_image image) {
	free(image);
	free_ui8matrix(image->img, image->nrl, image->nrh, image->ncl, image->nch);
	free_ui8matrix(image->source, image->nrl, image->nrh, image->ncl, image->nch);
	free_ui8matrix(image->fond, image->nrl, image->nrh, image->ncl, image->nch);
	free_ui8matrix(image->diff, image->nrl, image->nrh, image->ncl, image->nch);
	free_ui8matrix(image->var, image->nrl, image->nrh, image->ncl, image->nch);
	free_ui8matrix(image->move, image->nrl, image->nrh, image->ncl, image->nch);
}
