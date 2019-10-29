#include "img.h"

void init_pixel(p_pixel pixel, uint8 pixval) {
	pixel->pix = pixval;
	pixel->move = 0;
}

p_image create_image(char* filename) {

	int i, j;
	long t1, h, t2, w;
	uint8** mat;
	p_image tmp;

	mat = LoadPGM_ui8matrix(filename, &t1, &h, &t2, &w);
	/*alloc in memory for image*/
	tmp = (p_image)malloc(sizeof(image));
	if (!tmp)
		printf("Malloc error of image in create_image\n")
	/*alloc in memory for pixels*/
	tmp->img = (p_pixel*)malloc((int)h*sizeof(p_pixel));
	for (i = 0; i < (int)h; i++)
		(tmp->img)[i] = (p_pixel)malloc((int)w*sizeof(pixel));
	/*init pixels and variables for image*/
	for (i = 0; i < (int)h; i++) {
		for (j = 0; j < (int)w; j++)
			init_pixel((tmp->img)[i][j], tmp[i][j]); 
	}
	tmp->height = (int)h;
	tmp->width = (int)w;

	return tmp;
}

void free_image(p_image image) {

	int i;
	for (i = 0; i < image->height; i++)
		free((image->img[i]))
	free(image->img);
	free(image);
}
