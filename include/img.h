#ifndef __IMG_H__
#define __IMG_H__

#define MOVE 1
#define STATIC 0


typedef struct image {
	long nrl;
	long nrh;
	long ncl;
	long nch;
	uint8** I;
	uint8** M;
	uint8** O;
	uint8** V;
	uint8** E;
} image, *p_image;

p_image create_image(char* filename);

void free_image(p_image image);

#endif // __IMG_H__