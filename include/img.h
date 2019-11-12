#ifndef __IMG_H__
#define __IMG_H__

#define MOVE 1
#define STATIC 0


typedef struct image {
	long nrl;
	long nrh;
	long ncl;
	long nch;
	uint8** I;	// Image de source
	uint8** M;	// Image de fond
	uint8** O;	// Image de difference
	uint8** V;	// Image de variance
	uint8** E;	// Image d'etiquette binaires
	uint8** Omega; // Image de morphologie
} image, *p_image;

p_image create_image(char* filename);

void free_image(p_image image);

#endif // __IMG_H__