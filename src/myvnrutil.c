
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


void print_vui8vector(vuint8 *vV, int ncl, int nch, char* format, char *name) {
	int col_cnt;

	if (name != NULL) printf("%s",name);
	
	uint8 *p;
	int i = 0;
	col_cnt = 0;
	for (int v = -1; v <= 1; v++)  {
		p = (uint8*)&vV[v];
		for(i=0; i<16; i++){
			if (16 + ncl <= col_cnt && col_cnt <= 16 + nch)
				printf(format, p[i]);
			col_cnt++;
		}
	}

}

void print_vui8matrix(vuint8 **vM, int i0, int i1, int ncl, int nch, char* format, char *name) {
	int col_cnt;
	uint8 *p;

	if (name != NULL) printf("%s",name);
	
	for (int i = i0; i < i1 + 1; i++) {
        for (int j = ncl; j < nch + 1; j++) {
            p = (uint8 *)&vM[i][j];
            for(int k = 0; k < 16; k++) {
				printf(format, p[k]);
                
            }
        }
        printf("\n");
    }
}
/*---------------------------------------------------------------------------------------------*/
void copy_vui8matrix_vui8matrix(vuint8** X, int nrl, int nrh, int v0, int v1, vuint8** Y) {
/*---------------------------------------------------------------------------------------------*/
	long i;
	int j;

	vuint8 vX;

	for (i = nrl; i <= nrh; i++)
		for (j = v0; j <= v1; j++ ) {
			vX = _mm_load_si128(&X[i][j]);
			_mm_store_si128((vuint8*)&Y[i][j], vX);
		}
}

vuint8 **ui8matrix_to_vui8matrix(uint8 **X , long  nrl, long  nrh, long ncl, long nch, int *i0, int *i1, int *j0, int *j1)
{
    uint8   **Y;
    vuint8 **vY;
    vuint8 x;
	uint8 *p = (uint8*) &x;

    int row, v, card, col, bord, sj1;
    *i0 = nrl;
    *i1 = nrh;
    card = card_vuint8();
    s2v1D(ncl, nch, card, j0, j1);

    vY =       vui8matrix(nrl, nrh, *j0       , *j1);
    Y  = filled_ui8matrix(nrl, nrh, *j0 * card, (*j1 + 1) * card, 0);
    copy_ui8matrix_ui8matrix(X, nrl, nrh, ncl, nch, Y);
   
    sj1 = nch - (nch % card);
    for (row = nrl; row < nrh + 1; row++) {
        for (v = *j0; v < *j1 + 1 ; v++){
            col = v * card;
            
            x = _mm_set_epi8(X[row][col + 15], X[row][col + 14], X[row][col + 13], X[row][col + 12],
                             X[row][col + 11], X[row][col + 10], X[row][col +  9], X[row][col +  8],
                             X[row][col +  7], X[row][col +  6], X[row][col +  5], X[row][col +  4],
                             X[row][col +  3], X[row][col +  2], X[row][col +  1], X[row][col +  0]);            
            _mm_store_si128(&vY[row][v], x);
        }
        x = _mm_setzero_si128();
        for (v = sj1; v < nch + 1; v++) 
        	p[v - sj1] = X[row][v];
        _mm_store_si128(&vY[row][*j1], x);
    }
    free_ui8matrix(Y, nrl, nrh, *j0 * card, (*j1 + 1) * card);
    return vY;
}
uint8** vui8matrix_to_ui8matrix(vuint8** vX, long i0, long i1, int j0, int j1, long *nrl, long *nrh, long *ncl, long *nch)
/*--------------------------------------------------------------------------------*/
{
    int i, j, card;
    int ncl_, nch_;
    vuint8 x;

    card = card_vuint8();
    *nrl = i0;        *nrh = i1;
    
    v2m1D(j0, j1, card, &ncl_, &nch_);
    *ncl = ncl_;
    *nch = nch_;
    
    uint8 **Y = filled_ui8matrix(i0, i1, ncl_, nch_, 0);

    
    for(i=i0; i<=i1; i++) {
        for(j=j0; j<=j1; j++) {
            vuint8 T;
            x = _mm_load_si128(&vX[i][j]);
            _mm_store_si128((vuint*)&T, x);                  
            memcpy(&Y[i][j * card], (uint8*)&T, sizeof(uint8) * card);

        }
    }
    return Y;
}
vuint8 **LoadPGM_vui8matrix(char *filename, int *nrl, int *nrh, int *v0, int *v1)
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