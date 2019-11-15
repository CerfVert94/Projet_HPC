#ifndef __IMG_H__
#define __IMG_H__

#define MOVE 1
#define STATIC 0


typedef struct vimage {
	int  v0, v1;
	int  m0, m1; 
	long nrl, nrh;
	long ncl, nch;
	vuint8** I;	// Image de source
	vuint8** M;	// Image de fond
	vuint8** O;	// Image de difference
	vuint8** V;	// Image de variance
	vuint8** E;	// Image d'etiquette binaires
	vuint8** Omega; // Image de morphologie
} vimage, *p_vimage;

p_vimage create_vimage(char* filename);

void free_vimage(p_vimage image);

#endif // __IMG_H__