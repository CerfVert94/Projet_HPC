
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

#include <xmmintrin.h>
#include <emmintrin.h>
#include <mmintrin.h>

#include "nrdef.h"
#include "nrutil.h"
#include "vnrdef.h"
#include "vnrutil.h"
#include "myvnrutil.h"
#include "mynrutil.h"
#include "util.h"
uint8** vui8matrix_to_ui8matrix(vuint8** vX, long i0, long i1, int j0, int j1, long *nrl, long *nrh, long *ncl, long *nch)
/*--------------------------------------------------------------------------------*/
{
    int i, j, card;
    vuint8 x;

    card = card_vuint8();
    *nrl = i0;        *nrh = i1;
    *ncl = j0 * card; *nch = (j0 + (j1 - j0 + 1)) * card - 1;
    
    uint8 **Y = filled_ui8matrix(*nrl, *nrh, *ncl, *nch, 0);

    
    for(i=i0; i<=i1; i++) {
        for(j=j0; j<=j1; j++) {
            vuint8 T[1];
            uint8 *U;
            x = _mm_load_si128(&vX[i][j]);
            _mm_store_si128(T, x);                  
            memmove(&Y[i][j * card], (uint8*)T, sizeof(uint8) * card);

        }
    }
    return Y;
}
vuint8 **LoadPGM_vui8matrix(char *filename, long *nrl, long *nrh, int *v0, int *v1)
{
	int i, j, k, l, z, r;
	long rl, rh, cl, ch;
	long ncl, nch;

	uint8** img;
	uint8** imgb;

	int card;
	int m0, m1;
	vuint8 x;
	uint8 *p = (uint8*) &x;

	img = LoadPGM_ui8matrix(filename, &rl, &rh, &cl, &ch);
	*nrl = rl - BORD;
	*nrh = rh + BORD;
	ncl = cl - BORD;
	nch = ch + BORD;
	card = card_vuint8();
    s2v1D(ncl, nch, card, v0, v1);
    v2m1D(*v0, *v1, card, &m0, &m1);
	// printf("Img:%d %d\n", *v0, *v1);
    // printf("rl:%3ld rh:%3ld\n", rl, rh);
    // printf("cl:%3ld ch:%3ld\n", cl, ch);
    // printf("nrl:%3ld nrh:%3ld\n", *nrl, *nrh);
    // printf("ncl:%3ld nch:%3ld\n", ncl, nch);
    // printf("v0:%3d v1:%3d\n", *v0, *v1);
    // printf("m0:%3d m1:%3d\n", m0, m1);
    
    z = (*v1-vBORD)*card;
    r = ch-z;

    vuint8** I = vui8matrix(*nrl, *nrh, *v0, *v1);
    for (i = rl; i <= rh; i++) {
        for (j = *v0+vBORD; j < *v1-vBORD; j++) {
            l = j*card;
            x = init_vuint8_all(img[i][l], img[i][l+1], img[i][l+2], img[i][l+3], img[i][l+4], img[i][l+5], img[i][l+6], img[i][l+7], img[i][l+8], img[i][l+9], img[i][l+10], img[i][l+11], img[i][l+12], img[i][l+13], img[i][l+14], img[i][l+15]);
            _mm_store_si128(&I[i][j], x);
        }
        x = init_vuint8(0);
        for (k = 0; k <= r; k++)
        	p[k] = img[i][z+k];
        _mm_store_si128(&I[i][*v1-vBORD], x);
    }
	return I;
}