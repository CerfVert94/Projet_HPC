#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "nrutil.h"

#define MOVE 1
#define STATIC 0

typedef struct pixel {
	uint8 pix;
	int move;
} pixel, *p_pixel;

typedef struct image {
	int width;
	int height;
	pixel** img;
} image, *p_image;


void init_pixel(p_pixel pixel, uint8 pixval);

p_image create_image(char* filename);

void free_image(p_image image);